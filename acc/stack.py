# -*- coding: utf-8 -*-
# Copyright 2019- Weijia Sun, MIT license

"""
stacking
"""

import glob
import numpy as np
from scipy.signal import hilbert
import commentjson
from obspy import read, Stream

from acc.io import _load_json

from tqdm import tqdm
from functools import partial
from multiprocessing import Pool
import multiprocessing
import logging

logging.basicConfig(level=logging.NOTSET, format='%(asctime)s %(filename)s[line:%(lineno)d] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filename="acc.processing.log", filemode="a"
                    )


def linear_stack(stream, normalize=True):
    if normalize:
        stream.normalize(global_max=False)

    stk = np.zeros_like(stream[0].data)
    for tr in stream:
        stk += tr.data

    stk /= len(stream)

    tr = stream[0]
    tr.data = np.copy(stk)
    tr.stats.starttime = 0
    tr.stats.update({"number_of_ray": len(stream), "stack_method":{"linear"}})

    return tr


def pws_stack(stream, power=2, normalize=True):
    if normalize:
        stream.normalize(global_max=False)

    # safe guard to be a positive integer number
    power = int(abs(power))

    tr_lin = linear_stack(stream, normalize)
    if power == 0:
        return tr_lin

    phase = np.zeros_like(stream[0].data)
    for tr in stream:
        # hilbert transform
        hilb = hilbert(tr.data)
        phase += np.abs(hilb / np.abs(hilb))

    phase /= len(stream)
    power = int(power)
    stk = tr_lin.data * np.power(phase, power)

    tr = stream[0]
    tr.data = np.copy(stk)
    tr.stats.starttime = 0
    tr.stats.update({"number_of_ray": len(stream), "stack_method":{"PWS", power}})

    return tr


def stack(jsonfile):
    kwargs = _load_json(jsonfile)

    path = kwargs["io"]["outpath"]
    stations = glob.glob(path + "/1_results/*")

    do_work = partial(_stack, **kwargs)


def _stack(path, **kwargs):

    files = glob.glob(path + "/*")

    st = Stream()
    for file in files:
        try:
            tr = read(file)[0]
        except:
            continue
        st.append(tr)

    # get channels
    channels = []
    for tr in st:
        chn = tr.stats.channel
        if chn not in channels:
            channels.append(chn)

    method = kwargs["stack"]["method"]
    power = kwargs["stack"]["power"]
    outpath = kwargs["io"]["outpath"] + "/1a_stack"

    stream = st.copy()
    for chn in channels:

        st = stream.select(channel=chn)

        if method == "linear":
            tr = linear_stack(st, normalize=True)
        elif method == "PWS":
            tr = pws_stack(st, power, normalize=True)
        elif method == "bootstrap_linear" or method == "bootstrap_PWS":
            # note tr here is a stream containing two traces, mean and std
            tr = _bootstrap(st, normalize=True, **kwargs)

        filen = outpath + "/" + tr.id + "_%s.pkl" % method
        tr.write(filename=filen, format="PICKLE")


def _bootstrap(st, normalize=True, **kwargs):

    if normalize:
        stream.normalize(global_max=False)

    # generate randoma number
    ntrace = len(st)
    options = kwargs["stack"]
    n_iter = options["n_iter"]
    percentage = options["percentage"]
    seed = kwargs["seed"]
    rand_list = _gen_random_number(low=0, high=ntrace,
                                   n_iter=n_iter, percentage=percentage, seed=seed)

    # stacking
    method = kwargs["stack"]["method"]
    power = kwargs["stack"]["power"]

    # method is bootstrap_linear or bootstrap_PWS
    method = method.split("_")[1]

    npts = st[0].stats.npts
    boot = np.zeros([n_iter, npts])
    m = 0
    for r in rand_list:
        # random stream
        st_rand = Stream()
        for i in r:
            tr = st[i]
            st_rand.append(tr)

        # stack
        if method == "linear":
            tr_stk = linear_stack(st_rand, normalize=True)
        elif method == "PWS":
            tr_stk = pws_stack(st_rand, power, normalize=True)

        boot[i] = np.copy(tr_stk.data)

    # mean value
    tr_mean = st[0].copy()
    tr_std = st[0].copy()
    tr_mean.data = np.mean(boot, axis=0)
    tr_mean.stats.update({"statistics":"mean"})
    tr_std.data = np.std(boot, axis=0)
    tr_st.stats.update({"statistics":"std"})

    st2 = Stream()
    st2.append(tr_mean)
    st2.append(tr_std)
    return st2


def _gen_random_number(low=0, high=100, n_iter=100, percentage=0.9, seed=59):

    np.random.seed(seed)

    n_trace = int(percentage * high)

    rand_list = []
    for i in range(n_iter):
        r = np.random.randint(low=high, high=high)
        rand_list.append(r)

        # reset random seed. use the multiplication of the first nonzero value
        # generated random number from the previous random seed as the new seed
        s = 1
        counter = 0
        for k in r:
            if k != 0:
                s *= k
            if counter > 5:
                break
            counter += 1
        np.random.seed(seed=s)

    return rand_list





