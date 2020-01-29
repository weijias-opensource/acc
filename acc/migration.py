# -*- coding: utf-8 -*-
from obspy.taup import TauPyModel
import pandas as pd
from scipy.interpolate import interp1d
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from acc.stack import pws_stack
from obspy import Stream

import os
import glob
from tqdm import tqdm
from functools import partial
from multiprocessing import Pool
import multiprocessing

from acc.io import _load_json
import logging
from obspy import read

# standalone version
# this function is supposed to be modified to parallel version.
def mig_one_station(stream, model="ak135", earth_radius=6378137.0, depth_range=(0, 300, 1)):

    try:
        # prior to the model used to calculate traveltime
        model = stream[0].stats.model
    except AttributeError:
        model = model

    taup_model = TauPyModel(model=model)
    for tr in stream:
        evdp = tr.stats.event_depth
        evla = tr.stats.event_latitude
        evlo = tr.stats.event_longitude

        stla = tr.stats.station_latitude
        stlo = tr.stats.station_longitude
        phase = tr.stats.phase

        # tr.filter(type="bandpass", freqmin=0.05, freqmax=0.5, zerophase=True)

        arrivals = taup_model.get_ray_paths_geo(source_depth_in_km=evdp,
                                                source_latitude_in_deg=evla, source_longitude_in_deg=evlo,
                                                receiver_latitude_in_deg=stla, receiver_longitude_in_deg=stlo,
                                                phase_list=(phase,), resample=True)

        arr = arrivals[0]
        # raypath coordinates. The dataframe contains six columns, ray_p, traveltime, dist in rad, depth, lat and lon.
        df = pd.DataFrame(arr.path)

        # one-way reflection time
        time = tr.times() / 2
        data = tr.data

        # ray path information: traveltime, depth, and coordinates
        # reverse time: the maximum traveltime minus traveltime at varied depths or points
        df["time_reverse"] = df["time"].iloc[-1] - df["time"]
        # get sub-dataset with traveltime (ray path) slightly longer than the one-way reflected traveltime from acc
        # for convenience of interpolation
        df2 = df[df["time_reverse"] < time[-1] * 1.2]

        # convert the raypath information from dataframe to numpy array
        ttime = df["time_reverse"].to_numpy()
        depth = df["depth"].to_numpy()
        lat = df["lat"].to_numpy()
        lon = df["lon"].to_numpy()

        # simple migration by back projection
        fdepth = interp1d(ttime, depth)
        flat = interp1d(ttime, lat)
        flon = interp1d(ttime, lon)

        # interpolation fitting in with one-way reflection time to accomplish back projection
        depths = fdepth(time)
        lats = flat(time)
        lons = flon(time)

        # calculate the distances from pierce points to station at different depths
        # or varied radius
        dists = np.zeros_like(depths)
        for i, d in enumerate(depths):
            dists[i], az, baz = gps2dist_azimuth(lat1=tr.stats.station_latitude, lon1=tr.stats.station_longitude,
                                                 lat2=lats[i], lon2=lons[i], a=earth_radius-d)
            # convert the unit from m to km. If the piere point is to the south of the station then distance is negative.
            if lats[i] <= tr.stats.station_latitude:
                dists[i] = -dists[i]
        dists /= 1000

        # the following several lines can be commented. This was firstly written to save data with
        # irregular depth sampling intervals.
        mig1sta = pd.DataFrame(columns=["time", "depth", "dist", "lat", "lon", "data"])
        mig1sta["depth"] = depths
        mig1sta["lat"] = lats
        mig1sta["lon"] = lons
        mig1sta["time"] = time
        tr.normalize()
        mig1sta["data"] = tr.data
        mig1sta["dist"] = dists

        # second interpolation to regular depth grid. after back-projection the depth sampling is irregular.
        # depth range 0-300 km with an interval of 0.5 km.
        delta = depth_range[2]
        d = np.arange(depth_range[0], depth_range[1]+delta, delta)
        tr.stats.delta = delta
        ftime = interp1d(depths, time)
        fdist = interp1d(depths, dists)
        flat = interp1d(depths, lats)
        flon = interp1d(depths, lons)
        fdata = interp1d(depths, tr.data)

        time = ftime(d)
        dists = fdist(d)
        lats = flat(d)
        lons = flon(d)
        tr.data = fdata(d)
        depths = np.copy(d)

        mig1sta = pd.DataFrame(columns=["time", "depth", "dist", "lat", "lon", "data"])
        mig1sta["depth"] = depths
        mig1sta["lat"] = lats
        mig1sta["lon"] = lons
        mig1sta["time"] = time
        tr.normalize()
        mig1sta["data"] = tr.data
        mig1sta["dist"] = dists

        header = {"path": df2, "mig":mig1sta}
        tr.stats.update(header)

    return stream


def _mig_1(path, model="ak135"):
    # read autocorrelograms of a given station saved in `path`
    st = read(path + "/*pkl")
    st_mig = mig_one_station(stream=st, model=model, earth_radius=6378137.0)

    return st_mig


def migration_1station(jsonfile):

    kwargs = _load_json(jsonfile)
    io = kwargs["io"]

    njobs = kwargs["njobs"]
    if njobs > multiprocessing.cpu_count():
        njobs = multiprocessing.cpu_count()

    # datapath containing auto-correlograms
    path = io["outpath"] + "/1_results"
    temp = glob.glob(path + "/*")
    stations = []
    for t in temp:
        if os.path.isdir(t):
            stations.append(t)

    model = kwargs["migration"]["model"]
    if model is None:
        model = kwargs["tt_model"]

    do_work = partial(_mig_1, model=model)

    st_mig_stations = []
    if njobs == 1:
        logging.info('do work sequential (%d cores)', njobs)
        for sta in tqdm(stations, total=len(stations)):
            st = do_work(sta)
            st_mig_stations.append(st)
    else:
        logging.debug('do work parallel (%d cores)', njobs)
        pool = multiprocessing.Pool(njobs)
        for st in tqdm(pool.imap_unordered(do_work, stations), total=len(stations)):
            st_mig_stations.append(st)
        pool.close()
        pool.join()

    # save data into disk
    outpath = io["outpath"] + "/migration_1station"
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    for st in st_mig_stations:
        # station id
        tr = st[0]
        if tr.stats.location == "":
            station_id = ".".join([tr.stats.network, tr.stats.station])
        else:
            station_id = ".".join([tr.stats.network, tr.stats.station, tr.stats.location])
        fn = outpath + "/" + station_id + ".pkl"
        st.write(fn, format="PICKLE")


def migration_one_station_stack(stream, method="PWS", power=2, time_range=[8, 16], coeff=0.5):
    """
    Stacking of migrated traces for one station. Currently not used.

    :param stream:
    :param method: "PWS" for phase-weighted stacking and "linear" for linear stacking
    :param power: used by "PWS" only. power=0 is linear stacking
    :param time_range: if None do simple stacking (PWS or linear). if a list contains two elements then
    it will select traces with high resemblance with the initial stacking to get final stacked trace.
    :param coeff: traces with correlation efficient higher than the value will be selected.
    :return: trace after stacking

    .. Note::
        The one more procedure is to improve signal-to-noise ratio. You may check the stats header `corrstack`.
        If the key equals zero, the stacking is the general ones, otherwise the correlated stacking are implemented.
        The value denotes the number of stacked traces.
    """
    # first initial stacking and find similar traces for final stacking.
    tr_stk = pws_stack(stream, power=power, normalize=True)

    h = {"corrstack": 0}
    tr_stk.stats.update(h)

    if time_range is None:
        return tr_stk

    delta = tr_stk.stats.delta
    i1 = int(time_range[0] / delta)
    i2 = int(time_range[1] / delta)

    data1 = np.copy(tr_stk.data[i1:i2])

    traces = []
    for tr in stream:
        data2 = np.copy(tr.data[i1:i2])
        c = np.corrcoef(data1, data2)
        if c[1][0] >= coeff:
            traces.append(tr)

    if len(traces) < 1:
        return tr_stk

    # final stacking with high coherence
    st2 = Stream(traces=traces)
    tr_stk = pws_stack(st2, power=power, normalize=True)
    h = {"corrstack": len(traces)}
    tr_stk.stats.update(h)

    return tr_stk





