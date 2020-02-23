# -*- coding: utf-8 -*-
"""
I/O modules
"""
import os
import commentjson
import glob
from obspy import read
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics import kilometers2degrees

from tqdm import tqdm
from functools import partial
from multiprocessing import Pool
import multiprocessing

import logging

logging.basicConfig(level=logging.NOTSET, format='%(asctime)s %(filename)s[line:%(lineno)d] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filename="acc.io.log", filemode="a"
                    )


def _load_json(jsonfile):
    """
    Load parameters from the json formatted file

    :param str jsonfile: json file containing parameters
    :return dict kwargs: dictionary containing parameters
    """
    with open(jsonfile, "r") as fp:
        kwargs = commentjson.load(fp)

    if kwargs["njobs"] is None:
        kwargs["njobs"] = multiprocessing.cpu_count()

    return kwargs


def _acc_read(file, outpath, acc_type, phase, force=True, tt_model="ak135", depth_unit="km"):
    """
    Read seismic data.

    :param file: file name to be read
    :param outpath: where to save the data
    :param acc_type: event or noise
    :param phase: if acc_type='event', then calculate the traveltime of the phase
    :param force: force to overwrite when file exists if force is True
    :param tt_model: 1-d theoretical model, e.g., ak135 used to calculate traveltime of a specific phase
    :param depth_unit: unit of depth, general one is km, but in some cases it could be m.
    :return:
    """

    tr = read(file)[0]

    station_id = _get_station_id(tr)
    trace_id = tr.id

    valid = 0
    if acc_type == "event":
        # tt_model = data_selection["tt_model"]
        event_id = _get_event_id(tr)
        # depth_unit = data_selection["depth_unit"]
        tr = _get_event_data(tr, tt_model, phase=phase, acc_type=acc_type, depth_unit=depth_unit)
        if tr is None:
            logging.warn("%s no valid arrival for phase %s", file, phase)
            return 0
        else:
            valid = 1
    elif acc_type == "noise":
        event_id = _get_noise_id(tr)
        tr = _get_noise_data(tr, acc_type=acc_type)
        valid = 1

    filename = trace_id + "_" + event_id + ".pkl"
    filepath = "/".join([outpath, station_id])
    # if not os.path.exists(filepath):
    try:
        os.makedirs(filepath)
    except:
        pass

    filen = filepath + "/" + filename
    if not force and os.path.exists(filen):
        pass
    else:
        tr.write(filen, format="PICKLE")
    return valid


def _get_sac_origin(tr):
    """
    Get the origin time of an event trace in sac format.

    :param tr: A event trace
    :return origin: origin time of an event

    .. Note::
        The trace should be sac formatted.

    """
    try:
        origin = tr.stats.starttime - tr.stats.sac.b + tr.stats.sac.o
    except AttributeError:
        logging.critical("no valid event origin found in header: tr.stats.sac.o")
        logging.critical("Please check acc_type='event' or 'noise'")
        raise AttributeError
    return origin


def _get_event_id(tr):
    """
    Get event id from a sac-formatted trace.

    :param tr:
    :return str event_id: event id
    """

    origin = _get_sac_origin(tr)
    event_id = origin.datetime.strftime("%Y%m%d%H%M%S")
    return event_id


def _get_event_id_tr(tr):
    """
    Get event id from a obspy trace. The obspy trace should have the event_time keyword.

    :param tr: a obspy trace
    :return str event_id: event id in "%Y%m%d%H%M%S", e.g., '20191019200909'
    """
    origin = tr.stats.event_time
    event_id = origin.datetime.strftime("%Y%m%d%H%M%S")
    return event_id


def _get_noise_id(tr):
    """
    Get trace id for noise type data. e.g., 'starttime-endtime'

    :param tr: an obspy trace
    :return str event_id: noise id
    """
    starttime = tr.stats.starttime
    endtime = tr.stats.endtime
    event_id = "-".join([starttime.strftime("%Y%m%d%H%M%S"), endtime.strftime("%Y%m%d%H%M%S")])
    return event_id


def _get_station_id(tr):
    """
    Get station id of a given trace.

    :param tr: trace
    :return station_id: station id formatted as '{newwork}.{station}.{location}'.
    """

    s = tr.id
    s = s.split(".")
    if s[2] == "":
        station_id = ".".join(s[:2])
    else:
        station_id = ".".join(s[:3])

    return station_id


def _get_noise_data(tr, acc_type):
    """
    Get update sac-trace header to obspy trace header station_latitude etc.

    :param tr:
    :param acc_type:
    :return tr:
    """

    station_longitude = tr.stats.sac.stlo
    station_latitude = tr.stats.sac.stla
    station_elevation = tr.stats.sac.stel

    header = {"type": acc_type,
              "station_latitude": station_latitude, "station_longitude": station_longitude,
              "station_elevation": station_elevation
              }

    tr.stats.update(header)

    return tr


def _get_event_data(tr, tt_model, phase, acc_type, depth_unit="km"):
    """
    Update a sac trace to a obspy trace and update trace header,
    and calculate theoretical traveltime of a specific model and phase

    :param tr:
    :param tt_model:
    :param phase:
    :param acc_type:
    :param depth_unit:
    :return:

    .. Note::
        The input trace should be read from sac-formatted files.

        depth_unit is not used. if depth>1000 then unit should be meter,
        since no events deeper than 700 km on the earth.

    """

    model = TauPyModel(model=tt_model)

    event_longitude = tr.stats.sac.evlo
    event_latitude = tr.stats.sac.evla
    event_depth = tr.stats.sac.evdp
    event_magnitude = tr.stats.sac.mag

    # if depth_unit == "m":
    #     event_depth /= 1000.0
    # in this case, the depth_unit is considered to be m.
    if event_depth > 1000:
        event_depth /= 1000

    station_longitude = tr.stats.sac.stlo
    station_latitude = tr.stats.sac.stla
    station_elevation = tr.stats.sac.stel

    try:
        component_azimuth = tr.stats.sac.cmpaz
        component_inclination = tr.stats.sac.cmpinc
    except:
        # print(tr.stats)
        if tr.stats.channel[-1] == "Z":
            component_azimuth = 0
            component_inclination = 0
        elif tr.stats.channel[-1] == "N":
            component_azimuth = 0
            component_inclination = 90
        elif tr.stats.channel[-1] == "E":
            component_azimuth = 90
            component_inclination = 90
        else:
            print("component is not ZNE. ", tr.stats.channel)
            os._exit(0)



    event_time = _get_sac_origin(tr)

    distance, azimuth, back_azimuth = gps2dist_azimuth(lat1=event_latitude, lon1=event_longitude,
                                                       lat2=station_latitude, lon2=station_longitude,
                                                       a=6378137.0, f=0.0033528106647474805)
    distance = kilometers2degrees(kilometer=distance / 1000.0)

    # travel time, slowness, inclinations
    arrivals = model.get_travel_times(source_depth_in_km=event_depth,
                                      distance_in_degree=distance,
                                      phase_list=[phase])
    if len(arrivals) < 1:
        return None

    arr = arrivals[0]

    onset = event_time + arr.time
    phase = phase
    inclination = arr.incident_angle
    slowness = arr.ray_param

    # pierce points
    # pp_latitude
    # pp_longitude
    # pp_depth

    # ray paths
    # arrivals = model.get_travel_times(source_depth_in_km=event_depth,
    #                                   distance_in_degree=distance,
    #                                   phase_list=[phase])
    header = {"model": tt_model, "type": acc_type,
              "event_latitude": event_latitude, "event_longitude": event_longitude, "event_depth": event_depth,
              "event_time": event_time, "event_magnitude": event_magnitude,
              "station_latitude": station_latitude, "station_longitude": station_longitude,
              "station_elevation": station_elevation,
              "component_azimuth": component_azimuth, "component_inclination":component_inclination,
              "onset": onset, "phase": phase, "inclination": inclination, "slowness": slowness,
              "distance": distance, "azimuth": azimuth, "back_azimuth": back_azimuth
              }

    tr.stats.update(header)

    return tr


def import_data(jsonfile):
    """
    Import data from external media

    :param jsonfile: parameter filename
    :return:
    """

    kwargs = _load_json(jsonfile)
    io = kwargs["io"]
    # data_selection = kwargs["data_selection"]
    njobs = kwargs["njobs"]

    if njobs > multiprocessing.cpu_count():
        njobs = multiprocessing.cpu_count()

    files = glob.glob(io["data"])
    outpath = io["outpath"] + "/0_raw"
    acc_type = kwargs["acc_type"]
    phase = kwargs["phase"]
    force = io["force"]
    tt_model = kwargs["tt_model"]
    depth_unit = kwargs["depth_unit"]

    if acc_type not in ["event", 'noise']:
        print("Invalid acc_type: %s. Aborted." % acc_type)
        logging.error("Invalid acc_type: %s. Aborted.", acc_type)
        exit(-1)

    do_work = partial(_acc_read, outpath=outpath, acc_type=acc_type,
                      phase=phase, force=force, tt_model=tt_model, depth_unit=depth_unit)
    numbers = []
    if njobs == 1:
        logging.info('do work sequential (%d cores)', njobs)
        for file in tqdm(files, total=len(files)):
            num = do_work(file)
            numbers.append(num)
    else:
        logging.debug('do work parallel (%d cores)', njobs)
        pool = multiprocessing.Pool(njobs)
        for num in tqdm(pool.imap_unordered(do_work, files), total=len(files)):
            numbers.append(num)
        pool.close()
        pool.join()

    logging.info("%d/%d files imported.", sum(numbers), len(files))
