# -*- coding: utf-8 -*-
import os
import glob
import commentjson
from acc.io import _load_json, _get_event_id_tr, _get_station_id, _get_noise_id
from acc.core import spectral_whitening, time_norm, rotation, remove_response
from acc.util import iter_time, IterMultipleComponents
from acc.core import _find_start_end_time

from obspy.signal.filter import envelope
from obspy import read, Stream, Trace

from scipy.signal import correlate
import numpy as np

from tqdm import tqdm
from functools import partial
from multiprocessing import Pool
import multiprocessing
import logging

logging.basicConfig(level=logging.NOTSET, format='%(asctime)s %(filename)s[line:%(lineno)d] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filename="acc.processing.log", filemode="a"
                    )


def processing_event_Z(jsonfile):
    """
    processing vertical component of event data

    :param jsonfile:
    :return:
    """

    # load imported data from files
    kwargs = _load_json(jsonfile)
    njobs = kwargs["njobs"]
    datapath = kwargs["io"]["outpath"] + "/0_raw"

    logging.info("This may take a lot of time if dataset is huge.")

    # get filenames in Unix style for further processing
    files = glob.glob(datapath + "/*/*.pkl")
    logging.info("%d files found.", len(files))

    # processing includes downsampling, detrend, demean,
    # instrument response removal, spectral whitening, temporal normalization
    # autocorrelation and filter, then output results.
    do_work = partial(_proc_event_Z, **kwargs)

    numbers = []
    if njobs == 1:
        logging.info('do work sequential (%d cores)', njobs)
        for file in tqdm(files, total=len(files)):
            num = do_work(file)
            numbers.append(num)
    else:
        logging.info('do work parallel (%d cores)', njobs)
        pool = multiprocessing.Pool(njobs)
        for num in tqdm(pool.imap_unordered(do_work, files), total=len(files)):
            numbers.append(num)
        pool.close()
        pool.join()

    logging.info("%d/%d files processed.", sum(numbers), len(files))


def processing_event_full(jsonfile):
    # load imported data from files
    kwargs = _load_json(jsonfile)
    njobs = kwargs["njobs"]
    datapath = kwargs["io"]["outpath"] + "/0_raw"

    logging.info("This may take a lot of time and consume internal memory if dataset is huge.")

    # get filenames in Unix style for further processing
    files = glob.glob(datapath + "/*/*.pkl")
    logging.info("%d files found.", len(files))

    # read files
    st = Stream()
    for file in files:
        st += read(file)
    sampling_rate = kwargs["preprocessing"]["sampling_rate"]
    # processing includes downsampling, detrend, demean, a large dataset
    st = _simple_proc(st, sampling_rate=sampling_rate, njobs=njobs)
    # instrument response removal, spectral whitening, temporal normalization
    # autocorrelation and filter, then output results.
    _proc_event_full(st, **kwargs)


def _proc_event_full(st, **kwargs):
    """
    processings including

    :param st:
    :param kwargs:
    :return:
    """
    # instrument response removal, spectral whitening, temporal normalization
    # autocorrelation and filter, then output results.
    def iter3c(stream):
        # for an event, there is always "onset"
        return IterMultipleComponents(stream, key="onset", number_components=(2, 3))

    # resp removal, rotation, spectral whitening, temporal normalization
    tasks = iter3c(st)
    # loop over streams in each stream containing the 3-component traces
    do_work = partial(_proc_event_rst, **kwargs)

    njobs = kwargs["njobs"]

    numbers = []
    logging.info("deep processing for full event correlogram.")
    print("deep processing for full event correlogram.")
    if njobs == 1:
        logging.info('do work sequential (%d cores)', njobs)
        for task in tqdm(tasks, total=len(tasks)):
            num = do_work(task)
            numbers.append(num)
    else:
        logging.info('do work parallel (%d cores)', njobs)
        pool = multiprocessing.Pool(njobs)
        for num in tqdm(pool.imap_unordered(do_work, tasks), total=len(tasks)):
            numbers.append(num)
        pool.close()
        pool.join()

    logging.info("%d/%d files processed.", sum(numbers), len(tasks))


def _proc_event_rst(st, **kwargs):
    """
    processing including rotation, spectral whitening and temporal normalization.
    and cross-correlation.

    :param st:
    :param kwargs:
    :return: 1 if processing successfully.

    .. Note::
        rst is short as Rotation Spectral whitening and Temporal normalization.
    """
    # calculation SNR
    # trace = Trace()
    for tr in st:
        if tr.stats.channel[-1] in ["Z", "L"]:
            trace = tr.copy()
    snr = _cal_event_snr(tr=trace)
    snr_threshold = kwargs["data_selection_event"]["snr_threshold"]
    if snr < snr_threshold:
        return 0
    for tr in st:
        tr.stats.update({"snr": snr})

    # remove instrumental response
    if kwargs["rm_resp_on"]:
        option_resp = kwargs["preprocessing"]["resp"]
        st2 = Stream()
        for tr in st:
            tr = remove_response(trace=tr, **option_resp)
            st2.append(tr)
        st = st2.copy()

    # st: stream containing 3-comp traces of an event recorded by a station
    # Rotation, Spectral whitening, and Temporal normalization
    rotation_method = kwargs["preprocessing"]["rotation_method"]
    st = _rot(st, method=rotation_method)

    # spectral whitening
    if kwargs["whiten_on"]:
        option_whiten = kwargs["preprocessing"]["whiten"]
        st2 = Stream()
        for tr in st:
            tr = spectral_whitening(tr=tr, **option_whiten)
            st2.append(tr)
        st = st2.copy()

    # temporal normalization
    if kwargs["time_norm_on"]:
        option_time_norm = kwargs["preprocessing"]["time_norm"]
        st2 = Stream()
        for tr in st:
            tr = time_norm(tr, **option_time_norm)
            st2.append(tr)
        st = st2.copy()

    # correlation
    # trace_l = Trace()
    # trace_q = Trace()
    # trace_t = Trace()
    for tr in st:
        if tr.stats.channel[-1] in ["Z", "L"]:
            trace_l = tr.copy()
        elif tr.stats.channel[-1] in ["R", "Q"]:
            trace_q = tr.copy()
        elif tr.stats.channel[-1] in ["T"]:
            trace_t = tr.copy()
        else:
            logging.warn("unknown trace channel.", tr.stats.channel)

    # auto and cross correlation
    # auto
    options = kwargs["correlate"]
    cc_ll = _crosscorrelation(tr1=trace_l, tr2=trace_l, **options)
    # cross
    options = kwargs["cross_correlate"]
    cc_ql = _crosscorrelation(tr1=trace_q, tr2=trace_l, **options)
    cc_tl = _crosscorrelation(tr1=trace_t, tr2=trace_l, **options)

    trace_l = cc_ll.copy()
    trace_q = cc_ql.copy()
    trace_t = cc_tl.copy()

    # rename channel name
    if rotation_method == "NE->RT":
        chn = ["ZZ", "RZ", "TZ"]
    elif rotation_method == "ZNE->LQT":
        chn = ["LL", "QL", 'TL']
    trace_l.stats.channel = chn[0]
    trace_q.stats.channel = chn[1]
    trace_t.stats.channel = chn[2]

    # write into disk
    for tr in [trace_l, trace_q, trace_t]:
        trace_id = tr.id
        station_id = _get_station_id(tr)
        event_id = _get_event_id_tr(tr)
        filename = trace_id + "_" + event_id + ".pkl"
        outpath = kwargs["io"]["outpath"] + "/1_results"
        filepath = "/".join([outpath, station_id])
        if not os.path.exists(filepath):
            os.makedirs(filepath)

        # extracting data, to have the same length of output data
        tlen = kwargs["correlate"]["window"]
        tlen = abs(tlen[1] - tlen[0])
        tr.trim(starttime=tr.stats.starttime, endtime=tr.stats.starttime + tlen, fill_value=0, pad=True)

        filen = filepath + "/" + filename
        force = kwargs["io"]["force"]
        if not force and os.path.exists(filen):
            pass
        else:
            tr.write(filen, format="PICKLE")
    return 1


def _rot(st, method="NE->RT"):
    """
    rotation.

    :param st:
    :param method: "NE->RT", "ZNE->LQT".
    :return:

    .. Note::
        the code can handle 'Z12', '123', 'ZNE'. It requires trace header of component_azimuth and component_inclination.
        if method == "NE->RT", then the two components rotation will be applied.
        But while the components are "123" or "Z12", they are rotated to "ZNE" first.
    """

    st3c = st.copy()

    # the rotating traces should have the same starttime and endtime
    tmin, tmax = _find_start_end_time(stream=st3c)
    if tmin > tmax:
        return 0
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
    # import matplotlib.pyplot as plt
    # print(st3c)
    # tr = st3c.select(channel="BHZ")[0]
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # ax[0].plot(tr.times(), tr.data)
    st3c.rotate(method=method)
    # tr = st3c.select(channel="BHL")[0]
    # ax[0].plot(tr.times(), tr.data)
    # ax[0].set_xlim([670, 700])
    # plt.tight_layout()
    # plt.show()
    # os._exit(0)
    return st3c


def _proc(tr, sampling_rate=10):
    """
    Basic processing including downsampling, detrend, and demean.

    :param tr: raw trace
    :param sampling_rate:
    :return tr: trace after processing
    """
    # deep copy
    tr2 = tr.copy()
    tr2.interpolate(sampling_rate)
    tr2.detrend(type="linear")
    tr2.detrend(type="demean")
    return tr2


def _simple_proc(st, sampling_rate=10, njobs=1):
    """
    A parallel version of `_proc`, i.e., Basic processing including downsampling, detrend, and demean.

    :param st: an obspy stream
    :param sampling_rate: expected sampling rate
    :param njobs: number of jobs or CPU to use
    :return st: stream after processing
    """

    # downsampling, detrend, demean
    do_work = partial(_proc, sampling_rate=sampling_rate)

    # trace_list = []
    # for tr in st:
    #     trace_list.append(tr)
    #
    st2 = Stream()
    logging.info("simple processing for full event correlogram.")
    print("simple processing for full event correlogram.")
    if njobs == 1:
        logging.info('do work sequential (%d cores)', njobs)
        for tr in tqdm(st, total=len(st)):
            tr2 = do_work(tr)
            st2.append(tr2)
    else:
        logging.info('do work parallel (%d cores)', njobs)
        pool = multiprocessing.Pool(njobs)
        for tr2 in tqdm(pool.imap_unordered(do_work, st), total=len(st)):
            st2.append(tr2)
        pool.close()
        pool.join()

    return st2


def _proc_event_Z(file, **kwargs):
    """
    Processing a single component Z including downsampling, detrend, demean,
    remove instrumental response, selection of data by SNR>threshold,
    spectral whitening, temporal normalization and auto-correlation.

    For event type data.

    :param file:
    :param kwargs:
    :return:
    """

    # read file and return an obspy trace
    tr = read(file)[0]
    # process z-component only.
    if tr.stats.channel[-1] is not "Z":
        logging.warning("component is not Z: %s", file)
        return 0

    # first step is downsamping data to reduce computational burden
    sampling_rate = kwargs["preprocessing"]["sampling_rate"]
    tr.interpolate(sampling_rate)
    tr.detrend(type="linear")
    tr.detrend(type="demean")

    options = kwargs["data_selection_event"]
    snr_threshold = options["snr_threshold"]
    snr = _select_data_event(tr, **options)
    if snr < snr_threshold:
        return 0
    tr.stats.update({"snr": snr})
    logging.info("%s SNR is %f", file, snr)

    # remove instrumental response
    if kwargs["rm_resp_on"]:
        option_resp = kwargs["preprocessing"]["resp"]
        tr = remove_response(trace=tr, **option_resp)

    # spectral whitening
    if kwargs["whiten_on"]:
        option_whiten = kwargs["preprocessing"]["whiten"]
        tr = spectral_whitening(tr=tr, **option_whiten)

    # temporal normalization
    if kwargs["time_norm_on"]:
        option_time_norm = kwargs["preprocessing"]["time_norm"]
        tr = time_norm(tr, **option_time_norm)

    # autocorrelation
    options = kwargs["correlate"]
    tr = _autocorrelation(tr, **options)

    # write data to disk
    tr.stats.channel = "ZZ"  # change trace channel to 'ZZ' after autocorrelation
    trace_id = tr.id
    station_id = _get_station_id(tr)
    event_id = _get_event_id_tr(tr)
    filename = trace_id + "_" + event_id + ".pkl"
    outpath = kwargs["io"]["outpath"] + "/1_results"
    filepath = "/".join([outpath, station_id])
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    # extracting data, to have the same length of output data
    tlen = kwargs["correlate"]["window"]
    tlen = abs(tlen[1] - tlen[0])
    tr.trim(starttime=tr.stats.starttime, endtime=tr.stats.starttime + tlen, fill_value=0, pad=True)

    filen = filepath + "/" + filename
    force = kwargs["io"]["force"]
    if not force and os.path.exists(filen):
        pass
    else:
        tr.write(filen, format="PICKLE")

    return 1


def _autocorrelation(tr, window=[-20, 70], filter=[0.5, 4], corners=2, zerophase=True):
    """
    Autocorrelation for event type data.

    :param tr:
    :param window:
    :param filter: A tuple or a list containing the lower and upper limit of filter.
    :param corners:
    :param zerophase:
    :return:

    .. Note::
        filter after autocorrelation.
    """
    tr2 = tr.copy()
    t1 = tr2.stats.onset + window[0]
    t2 = tr2.stats.onset + window[1]
    tr2.trim(starttime=t1, endtime=t2)
    npts = tr2.stats.npts
    data = tr2.data
    cc = correlate(in1=data, in2=data, mode="full")
    tr2.data = np.copy(cc)
    if filter is not None:
        tr2.filter(type="bandpass", freqmin=filter[0], freqmax=filter[1],
                   corners=corners, zerophase=zerophase)
    tr2.data = tr2.data[npts:]
    return tr2


def _crosscorrelation(tr1, tr2, window=[-20, 70], filter=[0.5, 4], corners=2, zerophase=True):
    """
    cross-correlation for event type data.

    :param tr1:
    :param tr2:
    :param window:
    :param filter:
    :param corners:
    :param zerophase:
    :return:
    """

    tra = tr1.copy()
    trb = tr2.copy()

    reverse = False
    if tra.stats.phase == "S":
        if tra.stats.channel == trb.stats.channel:
            # in this case, it is autocorrelation, which means SS
            window = window
        else:
            # in this case it is cross correlation over ZR/LQ or ZT/LT corresponding with Sp RF.
            # using different time window since converted P arrives earlier than the primary S.
            # if a window [-20, 70] is given, then a window [-70, 20] will be used in this case.
            w1, w2 = window
            window = [-w2, -w1]
            reverse = True

    t1 = tr2.stats.onset + window[0]
    t2 = tr2.stats.onset + window[1]
    tra.trim(starttime=t1, endtime=t2)
    trb.trim(starttime=t1, endtime=t2)
    npts = tra.stats.npts
    data = tra.data
    datb = trb.data

    cc = correlate(in1=data, in2=datb, mode="full")
    tr = tra.copy()
    tr.data = np.copy(cc)
    tr.filter(type="bandpass", freqmin=filter[0], freqmax=filter[1],
              corners=corners, zerophase=zerophase)
    if reverse:
        tr.data = np.copy(tr.data[::-1])
    tr.data = tr.data[npts:]
    return tr


def _select_data_event(tr, dist_range=[30, 90], magnitude=[5.0, 9.0],
                       snr_threshold=2.0, signal=[-10, 10], noise=[-100, -50],
                       waterlevel=1.0e-8):
    # select according to SNR
    distance = tr.stats.distance
    if distance < dist_range[0] or distance > dist_range[1]:
        return 0

    try:
        mag = tr.stats.magnitude
        if mag < magnitude[0] or mag > magnitude[1]:
            return 0
    except:
        pass

    # calculate SNR
    snr = _cal_event_snr(tr, signal=signal, noise=noise, waterlevel=waterlevel)

    return snr


def _cal_event_snr(tr, signal=[-10, 10], noise=[-100, -50], waterlevel=1.0e-8):
    """
    Calculation of SNR for event data.

    :param tr:
    :param signal:
    :param noise:
    :param waterlevel:
    :return:
    """

    tr_sig = tr.copy()
    tr_noi = tr.copy()

    t1 = tr.stats.onset + signal[0]
    t2 = tr.stats.onset + signal[1]
    if t1 < tr.stats.starttime or t2 > tr.stats.endtime:
        logging.warning("t1 < tr.stats.starttime")
        return 0
    tr_sig.trim(starttime=t1, endtime=t2)

    t1 = tr.stats.onset + noise[0]
    t2 = tr.stats.onset + noise[1]
    if t1 < tr.stats.starttime or t2 > tr.stats.endtime:
        logging.warning("t1 < tr.stats.starttime")
        return 0
    tr_noi.trim(starttime=t1, endtime=t2)

    sig = envelope(tr_sig.data)
    noi = envelope(tr_noi.data)
    snr = np.max(sig) / max(np.max(noi), waterlevel)

    return snr


def processing_noise_Z(jsonfile):
    """
    Parallel processing of noise autocorrelation.

    :param jsonfile:
    :return:
    """

    # load imported data from files
    kwargs = _load_json(jsonfile)
    njobs = kwargs["njobs"]
    datapath = kwargs["io"]["outpath"] + "/0_raw"

    logging.info("This may take a lot of time if dataset is huge.")

    # get filenames in Unix style for further processing
    files = glob.glob(datapath + "/*/*.pkl")
    logging.info("%d files found.", len(files))

    # processing includes downsampling, detrend, demean,
    # instrument response removal, spectral whitening, temporal normalization
    # autocorrelation and filter, then output results.
    do_work = partial(_proc_noise_Z, **kwargs)

    numbers = []
    if njobs == 1:
        logging.info('do work sequential (%d cores)', njobs)
        for file in tqdm(files, total=len(files)):
            num = do_work(file)
            numbers.append(num)
    else:
        logging.info('do work parallel (%d cores)', njobs)
        pool = multiprocessing.Pool(njobs)
        for num in tqdm(pool.imap_unordered(do_work, files), total=len(files)):
            numbers.append(num)
        pool.close()
        pool.join()

    logging.info("%d/%d files processed.", sum(numbers), len(files))


def _proc_noise_Z(file, **kwargs):
    """
    Internal functions of calculate z-comp autocorrelation for noise data.

    :param file:
    :param kwargs:
    :return:
    """

    # read file and return an obspy trace
    tr = read(file)[0]
    # process z-component only.
    if tr.stats.channel[-1] is not "Z":
        return 0

    # first step is downsamping data to reduce computational burden
    sampling_rate = kwargs["preprocessing"]["sampling_rate"]
    tr.interpolate(sampling_rate)
    tr.detrend(type="linear")
    tr.detrend(type="demean")

    # remove instrumental response
    if kwargs["rm_resp_on"]:
        option_resp = kwargs["preprocessing"]["resp"]
        tr = remove_response(trace=tr, **option_resp)

    # spectral whitening
    if kwargs["whiten_on"]:
        option_whiten = kwargs["preprocessing"]["whiten"]
        tr = spectral_whitening(tr=tr, **option_whiten)

    # temporal normalization
    if kwargs["time_norm_on"]:
        option_time_norm = kwargs["preprocessing"]["time_norm"]
        tr = time_norm(tr, **option_time_norm)

    # sliding autocorrelation
    options = kwargs["correlate_noise"]
    st = _sliding_autocorrelation(tr=tr, **options)

    # extracting data
    tlen = kwargs["correlate"]["window"]
    tlen = abs(tlen[1] - tlen[0])
    for tr in st:
        tr.trim(starttime=tr.stats.starttime, endtime=tr.stats.starttime + tlen, fill_value=0, pad=True)

    # write data to disk
    for tr in st:
        tr.stats.channel = "ZZ"  # change trace channel to 'ZZ' after autocorrelation
        trace_id = tr.id
        station_id = _get_station_id(tr)
        event_id = _get_noise_id(tr)
        filename = trace_id + "_" + event_id + ".pkl"
        outpath = kwargs["io"]["outpath"] + "/1_results"
        filepath = "/".join([outpath, station_id])
        if not os.path.exists(filepath):
            os.makedirs(filepath)

        filen = filepath + "/" + filename
        force = kwargs["io"]["force"]
        if not force and os.path.exists(filen):
            pass
        else:
            tr.write(filen, format="PICKLE")

    return 1


def _sliding_autocorrelation(tr, length=3600, overlap=1800,
                             filter=[0.5, 4], corners=2, zerophase=True):
    """
    Sliding autocorrelation for noise data.

    :param tr:
    :param length:
    :param overlap:
    :param filter:
    :param corners:
    :param zerophase:
    :return:
    """
    trace = tr.copy()
    time_series = iter_time(tr=trace, length=length, overlap=overlap)
    # print(time_series)
    if len(time_series) < 1:
        return 0

    st = Stream()
    for t1, t2 in time_series:
        tr2 = trace.copy()
        tr2.trim(starttime=t1, endtime=t2)

        # get gaps
        gap_st = Stream([tr2])
        gaps = gap_st.get_gaps()
        if len(gaps) > 0:
            continue

        npts = tr2.stats.npts
        data = tr2.data
        cc = correlate(in1=data, in2=data, mode="full")
        tr2.data = np.copy(cc)
        tr2.filter(type="bandpass", freqmin=filter[0], freqmax=filter[1],
                   corners=corners, zerophase=zerophase)
        tr2.data = tr2.data[npts:]
        st.append(tr2)

    if len(st) < 1:
        return 0
    else:
        return st
