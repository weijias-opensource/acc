# -*- coding: utf-8 -*-
# Copyright 2019- Weijia Sun, MIT license
"""
Utility functions and classes for Auto and Cross Correlogram calculation.
"""
import scipy.signal
import numpy as np
import collections
from obspy import UTCDateTime
from obspy import read

def smooth(x, window_len=None, window='flat', method='zeros'):
    """Smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.

    :param x: the input signal (numpy array)
    :param window_len: the dimension of the smoothing window; should be an
        odd integer
    :param window: the type of window from 'flat', 'hanning', 'hamming',
        'bartlett', 'blackman'
        flat window will produce a moving average smoothing.
    :param method: handling of border effects\n
        'zeros': zero padding on both ends (len(smooth(x)) = len(x))\n
        'reflect': pad reflected signal on both ends (same)\n
        'clip': pad signal on both ends with the last valid value (same)\n
        None: no handling of border effects
        (len(smooth(x)) = len(x) - len(window_len) + 1)
    """
    if window_len is None:
        return x
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is one of 'flat', 'hanning', 'hamming',"
                         "'bartlett', 'blackman'")
    if method == 'zeros':
        s = np.r_[np.zeros((window_len - 1) // 2), x,
                  np.zeros(window_len // 2)]
    elif method == 'reflect':
        s = np.r_[x[(window_len - 1) // 2:0:-1], x,
                  x[-1:-(window_len + 1) // 2:-1]]
    elif method == 'clip':
        s = np.r_[x[0] * np.ones((window_len - 1) // 2), x,
                  x[-1] * np.ones(window_len // 2)]
    else:
        s = x

    if window == 'flat':
        w = np.ones(window_len, 'd')
    else:
        w = getattr(np, window)(window_len)

    return scipy.signal.fftconvolve(w / w.sum(), s, mode='valid')


class IterMultipleComponents(object):
    """
    Return iterable to iterate over associated components of a stream.

    :param stream: Stream with different, possibly many traces. It is
        split into substreams with the same seed id (only last character
        i.e. component may vary)
    :type key: str or None
    :param key: Additionally, the stream is grouped by the values of
         the given stats entry to differentiate between e.g. different events
         (for example key='starttime', key='onset')
    :type number_components: int, tuple of ints or None
    :param number_components: Only iterate through substreams with
         matching number of components.
    """

    def __init__(self, stream, key=None, number_components=None):
        substreams = collections.defaultdict(stream.__class__)
        for tr in stream:
            k = (tr.id[:-1], str(tr.stats[key]) if key is not None else None)
            substreams[k].append(tr)
        n = number_components
        self.substreams = [s for _, s in sorted(substreams.items())
                           if n is None or len(s) == n or
                           (not isinstance(n, int) and len(s) in n)]

    def __len__(self):
        return len(self.substreams)

    def __iter__(self):
        for s in self.substreams:
            yield s


def iter_time(tr, length=3600, overlap=1800):

    tr2 = tr.copy()
    starttime = int(tr2.stats.starttime.timestamp / length) * length
    endtime = int(tr2.stats.endtime.timestamp / length) * length

    time_series = []
    for t in range(starttime, endtime, overlap):
        t1 = UTCDateTime(t)
        t2 = UTCDateTime(t1 + length)
        time_series.append([t1, t2])

    return time_series


def pkl2sac1(directory, suffix="pkl", fmt="SAC"):
    """
    Convert file from Pickle format with suffix of pkl to SAC format

    :param directory: the directory contains files to be converted.
    :param suffix: in this case, it should be "pkl".
    :param fmt: the target format to be converted. Support SAC, MSEED.

    Example: /the/path/hello.pkl to /the/path_SAC/hello.sac

    """
    import os
    import glob

    files = glob.glob(directory + "/*." + suffix)

    if fmt not in ["SAC", "MSEED"]:
        print("format should be 'SAC' or 'MSEED'.")
        os._exit(0)

    for file in files:
        tr = read(file)[0]
        savepath = os.path.dirname(file) + "_SAC"
        bn = os.path.basename(file)
        bn = os.path.splitext(bn)[0]
        bn = ".".join([bn, fmt.lower()])
        try:
            os.makedirs(savepath)
        except:
            pass
        fn = savepath + "/" + bn
        print(fn)
        tr.write(fn, format=fmt)
