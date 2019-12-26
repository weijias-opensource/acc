# -*- coding: utf-8 -*-
# Copyright 2019- Weijia Sun, MIT license
"""Preprocessing and correlation"""
import glob
# from functools import partial
# import itertools
import logging
# import multiprocessing

import numpy as np

from scipy.signal import freqz, iirfilter
from scipy.fftpack import fft, ifft, fftshift, ifftshift, next_fast_len

from acc.util import smooth as smooth_func
from acc.util import IterMultipleComponents

from obspy import Stream
from obspy import read_inventory
from obspy.signal.rotate import rotate2zne

log = logging.getLogger('acc.core.log')


def spectral_whitening(tr, smooth=None, filter=None,
                       waterlevel=1e-8, corners=2, zerophase=True):
    """
    Apply spectral whitening to data

    Data is divided by its smoothed (Default: None) amplitude spectrum.

    :param tr: trace to manipulate
    :param smooth: length of smoothing window in Hz
        (default None -> no smoothing)
    :param filter: filter spectrum with bandpass after whitening
        (tuple with min and max frequency)
        (default None -> no filter)
    :param waterlevel: waterlevel relative to mean of spectrum
    :param mask_again: weather to mask array after this operation again and
        set the corresponding data to 0
    :param corners: parameters parsing to filter,
    :param zerophase: parameters parsing to filter

    :return: whitened data
    """

    sr = tr.stats.sampling_rate
    data = np.copy(tr.data)
    # data = _fill_array(data, fill_value=0)
    # mask = np.ma.getmask(data)

    # transform to frequency domain
    nfft = next_fast_len(len(data))
    spec = fft(data, nfft)

    # amplitude spectrum
    spec_ampl = np.abs(spec)

    # normalization
    spec_ampl /= np.max(spec_ampl)
    spec_ampl_raw = np.copy(spec_ampl)

    # smooth
    if smooth:
        smooth = int(smooth * nfft / sr)
        spec_ampl = ifftshift(smooth_func(fftshift(spec_ampl), smooth))
        spec_ampl_smth = np.copy(spec_ampl)

    # save guard against division by 0
    spec_ampl[spec_ampl < waterlevel] = waterlevel

    # make the spectrum have the equivalent amplitude before/after smooth
    if smooth:
        scale = np.max(spec_ampl_raw) / np.max(spec_ampl_smth)
        spec /= spec_ampl * scale
    else:
        spec /= spec_ampl

    # FFT back to time domain
    ret = np.real(ifft(spec, nfft)[:len(data)])
    tr.data = ret

    # filter
    if filter is not None:
        tr.filter(type="bandpass", freqmin=filter[0], freqmax=filter[1],
                  corners=corners, zerophase=zerophase)

    return tr


def _fill_array(data, mask=None, fill_value=None):
    """
    Mask numpy array and/or fill array value without demasking.

    Additionally set fill_value to value.
    If data is not a MaskedArray and mask is None returns silently data.

    :param mask: apply mask to array
    :param fill_value: fill value
    """
    if mask is not None and mask is not False:
        data = np.ma.MaskedArray(data, mask=mask, copy=False)
    if np.ma.is_masked(data) and fill_value is not None:
        data._data[data.mask] = fill_value
        np.ma.set_fill_value(data, fill_value)
    #    elif not np.ma.is_masked(data):
    #        data = np.ma.filled(data)
    return data


def time_norm(tr, method, time_length=5, filter=(0.01, 1),
              corners=2, zerophase=True, waterlevel=1.0e-8):
    """
    Calculate normalized data, see e.g. Bensen et al. (2007, GJI)

    :param tr: Trace to manipulate
    :param str method:
        1bit: reduce data to +1 if >0 and -1 if <0\n
        run_abs_mean: running absolute mean normalization
    :param float time_length: time length to be smoothed, default 5 sec
    :param tuple filter: bandpass frequency band, default (0.01, 1) Hz
    :param int corners: corners default 2
    :param bool zerophase: default True
    :param float waterlevel: waterlevel to safe guard division

    :return: normalized data
    """

    data = tr.data
    data = _fill_array(data, fill_value=0)
    mask = np.ma.getmask(data)

    if method == '1bit':
        np.sign(data, out=data)

    elif method == 'run_abs_mean':
        data = _run_abs_mean(tr, time_length=time_length, filter=filter,
                             corners=corners, zerophase=zerophase,
                             waterlevel=waterlevel)
    else:
        msg = 'The method passed to time_norm is not known: %s.' % method
        raise ValueError(msg)

    tr.data = _fill_array(data, mask=mask, fill_value=0)
    return tr


def _run_abs_mean(tr, time_length=5, filter=(0.01, 1),
                  corners=4, zerophase=True, waterlevel=1e-8):
    """
    running absolute mean normalization

    :param tr: trace to be normalized
    :param float time_length: time length to be smoothed, default 5 sec
    :param tuple filter: bandpass frequency band, default (0.01, 1) Hz
    :param int corners: corners default 4
    :param bool zerophase: default True
    :param float waterlevel: waterlevel to safe guard division

    :return: data after temporal normalization
    """

    tr2 = tr.copy()
    tr2.filter(type="bandpass", freqmin=filter[0], freqmax=filter[1],
               corners=corners, zerophase=zerophase)
    data = np.abs(tr2.data)

    nl = int(time_length / tr2.stats.delta)
    smth_data = smooth_func(x=data, window_len=nl,
                            window='flat', method='zeros')

    # save guard against division by 0
    smth_data[smth_data < waterlevel] = waterlevel

    data = np.copy(tr.data)
    data /= smth_data

    return data


def rotation(stream, method="NE->RT", acc_type="event"):
    """
    Rotatation.

    :param stream: obspy stream to be rotated
    :param method: "NE->RT" or "ZNE->LQT"
    :param acc_type: "event", "noise"
    :return stream: obspy stream after rotation.

    .. note::
        The keywords `component_azimuth` and `component_inclination` must be given in the stats.
        Currently, only acc_type="event" are well debugged.
        The input 3-component data should be in the roughly same periods or time window.
    """

    rot_stream = Stream()

    # function to extract 3-component stream.
    # if the data are given event the data, then use onset as the key word to find the 3-component data
    if acc_type == "event":
        key = "onset"
    elif acc_type == "noise":
        key = "starttime"

    def iter3c(stream):
        return IterMultipleComponents(stream, key=key, number_components=(2, 3))

    for st3c in iter3c(stream):

        tmin, tmax = _find_start_end_time(stream=st3c)
        if tmin > tmax:
            continue
        st3c.trim(starttime=tmin, endtime=tmax)

        # handling the borehole components 1, 2, Z or 1, 2, 3.
        comp = []
        for tr in st3c:
            comp.append(tr.stats.channel[-1])
        comp.sort()
        cc = "".join(comp)
        if cc == "12Z":
            cc = "Z12"
        # rotate2zne
        if cc in ["123", "Z12"]:
            zne = rotate2zne(st3c[0].data, st3c[0].stats.component_azimuth, st3c[0].stats.component_inclination,
                             st3c[1].data, st3c[1].stats.component_azimuth, st3c[1].stats.component_inclination,
                             st3c[2].data, st3c[2].stats.component_azimuth, st3c[2].stats.component_inclination
                             )

            for tr, new_data, component in zip(st3c, zne, "ZNE"):
                tr.data = new_data
                tr.stats.channel = tr.stats.channel[:-1] + component
        # rotate to various coordinates
        st3c.rotate(method=method)
        rot_stream += st3c

    return rot_stream


def _find_start_end_time(stream):
    """
    find the start and end time of a specific stream.

    :param stream:
    :return starttime, endtime:
    """

    starttime = []
    endtime = []

    for tr in stream:
        starttime.append(tr.stats.starttime)
        endtime.append(tr.stats.endtime)

    return max(starttime), min(endtime)


def remove_response(trace, resp_path, pre_filt=(0.01, 0.02, 8, 9), output="VEL", format="XML"):
    """
    Romove instrumental response

    :param trace:
    :param resp_path:
    :param pre_filt:
    :param output:
    :param format:
    :return:

    .. Note::
        currently only support XML DATALESS RESP format.
    """

    station_id = ".".join([trace.stats.network, trace.stats.station])
    trace_id = trace.id
    tr = trace.copy()

    if format.upper()in ["XML", "DATALESS", "RESP"]:
        if format.upper() == "RESP":
            resp_file = glob.glob(resp_path + "/*%s*" % trace_id)
        else:
            resp_file = glob.glob(resp_path + "/*%s*" % station_id)
        if len(resp_file) < 1:
            logging.error("Cannot find Response file: ", station_id)
        inv = read_inventory(resp_file[0])
        tr.remove_response(inventory=inv, pre_filt=pre_filt, output=output)
    elif method == "SACPZ":
        # from obspy.io.sac.sacpz import attach_paz
        # attach_paz(tr, paz_file="SACPZ.AU.WRKA.--.BHZ")
        # https://docs.obspy.org/packages/obspy.io.sac.sacpz.attach_paz.html#obspy.io.sac.sacpz.attach_paz
        logging.warn("The SACPZ format not supported currently.")
        pass
    else:
        logging.warn("other RESP format not supported: ", method)
        pass

    return tr
