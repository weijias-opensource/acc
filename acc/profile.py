# -*- coding: utf-8 -*-
# @Author: Weijia Sun
# @Date:   2020-04-10 12:09:15
# @Last Modified by:   Weijia Sun
# @Last Modified time: 2020-04-10 20:14:27

from shapely.geometry import Polygon
from geographiclib.geodesic import Geodesic
import pandas as pd
import numpy as np
import timeit


_LARGE_BOX_WIDTH = 2000


def _get_box(latlon0, azimuth, length, width=_LARGE_BOX_WIDTH, offset=0):
    """Create a single box."""

    start = direct_geodetic(latlon0, azimuth, offset)
    azis = ((azimuth - 90) % 360, azimuth,
            (azimuth + 90) % 360, (azimuth + 180) % 360)
    dists = (width/2, length, width, length)
    latlon = start
    corners = []
    for a, d in zip(azis, dists):
        latlon = direct_geodetic(latlon, a, d)
        corners.append(latlon[::-1])
    box = {'poly': Polygon(corners),
           'length': length,
           'pos': offset + length/2,
           'latlon': direct_geodetic(start, azimuth, length/2)}
    return box


def _get_box_cartesian(latlon0, azimuth, length, width=2000, offset=0):
    """Create a single box."""

    import math
    d = offset
    a = azimuth
    lat = d*math.cos(a*math.pi/180) / 111.195
    lon = d*math.sin(a*math.pi/180) / 111.195
    start = (latlon0[0]+lat, latlon0[1]+lon)
    azis = ((azimuth - 90) % 360, azimuth,
            (azimuth + 90) % 360, (azimuth + 180) % 360, (azimuth + 270) % 360)
    dists = (width/2, length/2, width, length, width)
    latlon = start
    corners = []
    for a, d in zip(azis, dists):
        # latlon = direct_geodetic(latlon, a, d)
        # in km
        lat = d*math.cos(a*math.pi/180) / 111.195
        lon = d*math.sin(a*math.pi/180) / 111.195
        latlon = (latlon[0]+lat, latlon[1]+lon)
        corners.append(latlon[::-1])

    lat = length/2*math.cos(azimuth*math.pi/180) / 111.195
    lon = length/2*math.sin(azimuth*math.pi/180) / 111.195
    box = {'poly': Polygon(corners),
           'length': length,
           'pos': offset + length/2,
           'latlon': (start[0], start[1]+lon)}
#            'latlon': direct_geodetic(start, azimuth, length/2)}
    return box


def direct_geodetic(latlon, azi, dist):
    """
    Solve direct geodetic problem with geographiclib.

    :param tuple latlon: coordinates of first point
    :param azi: azimuth of direction
    :param dist: distance in km

    :return: coordinates (lat, lon) of second point on a WGS84 globe
    """
    coords = Geodesic.WGS84.Direct(latlon[0], latlon[1], azi, dist * 1000)
    return coords['lat2'], coords['lon2']


def get_profile_boxes(latlon0, azimuth, bins, width=_LARGE_BOX_WIDTH):
    """
    Create 2D boxes for usage in `profile()` function.

    :param tuple latlon0: coordinates of starting point of profile
    :param azimuth: azimuth of profile direction
    :param tuple bins: Edges of the distance bins in km (e.g. (0, 10, 20, 30))
    :param width: width of the boxes in km (default: large)
    :return: List of box dicts. Each box has the entries
        'poly' (shapely polygon with lonlat corners),
        'length' (length in km),
        'pos' (midpoint of box in km from starting coordinates),
        'latlon' (midpoint of box as coordinates)
    """
    boxes = []
    for i in range(len(bins)-1):
        length = bins[i+1] - bins[i]
        # box = _get_box(latlon0, azimuth, length, width, offset=bins[i])
        box = _get_box_cartesian(latlon0, azimuth, length, width, offset=bins[i])
        if i == 0:
            box['profile'] = {}
            box['profile']['latlon'] = latlon0
            box['profile']['azimuth'] = azimuth
            box['profile']['length'] = bins[-1] - bins[0]
            box['profile']['width'] = width
        boxes.append(box)
    return boxes


def _find_box(latlon, boxes, crs=None):
    """Return the box which encloses the coordinates."""
    import cartopy.crs as ccrs
    from shapely.geometry import Point
    if crs is None:
        latlons = [boxes[len(boxes)//2]['latlon']]
        latlon0 = np.median(latlons, axis=0)
        crs = ccrs.AzimuthalEquidistant(*latlon0[::-1])
    pc = ccrs.PlateCarree()
    p = crs.project_geometry(Point(*latlon[::-1]), pc)
    for box in boxes:
        poly = crs.project_geometry(box['poly'], pc)
        if p.within(poly):
            return box


def _find_box_cartesian(latlon, boxes):
    """Return the box which encloses the coordinates."""
    from shapely.geometry import Point
    p = Point(latlon[::-1])
    for box in boxes:
        poly = box['poly']
        if p.within(poly):
            return box


def internal_profile(tr):
    return tr.stats.mig


def profile(stream, boxes, crs=None):
    """
    Stack traces in stream by piercing point coordinates in defined boxes.

    :param stream: stream with pre-calculated piercing point coordinates
    :param boxes: boxes created with `get_profile_boxes()`
    :param crs: cartopy projection (default: AzimuthalEquidistant)
    :return: profile stream
    """

    # read all migrated data saved in trace headers into Pandas DataFrame
    # dfmig = pd.DataFrame()
    # for k, tr in enumerate(stream):
    #     df = tr.stats.mig
    #     dfmig = dfmig.append(df)
    #---------------------------------------------------------------------
    # parallel data, since for a profile there are thousands of raypath, 
    # the process would be very time-consuming
    from tqdm import tqdm
    import multiprocessing
    from functools import partial

    traces = []
    for tr in stream:
        traces.append(tr)
    
    df_lst = []
    pool = multiprocessing.Pool()
    for df in tqdm(pool.imap_unordered(internal_profile, traces), total=len(traces)):
        df_lst.append(df)
    pool.close()
    pool.join()
    dfmig = pd.concat(df_lst, axis=0, join="outer", ignore_index=True)
    #---------------------------------------------------------------------

    # get depth samplings
    depths = dfmig["depth"].unique()
    depths = depths.tolist()

    # get lateral samplings
    pos_dist = []
    for b in boxes:
        pos_dist.append(b["pos"])
    # print(pos_dist)

    # sampling number
    npos = len(pos_dist)
    ndep = len(depths)
    # init numpy array, stack - final stacked CRP image, count - hit number in each cell
    stack = np.zeros((ndep, npos), dtype=float)
    count = np.zeros((ndep, npos), dtype=int)

    # loops
    # serial version----------------------------------------------------
    
    # for kk, depth in enumerate(depths):
    #     print("Depth loops {}/{} in stacking.".format(kk, len(depths)))
    #     df2 = dfmig[dfmig["depth"] == depth]

    #     lats = df2["lat"].to_numpy()
    #     lons = df2["lon"].to_numpy()
    #     data = df2["data"].to_numpy()

    #     for i, lat in enumerate(lats):
    #         lon = lons[i]
    #         ppoint = (lat, lon)
    #         box = _find_box_cartesian(ppoint, boxes)
    #         if box is None:
    #             continue
    #         pos = box['pos']

    #         idep = depths.index(depth)
    #         ipos = pos_dist.index(pos)

    #         stack[idep][ipos] += data[i]
    #         count[idep][ipos] += 1
    
    # end serial version----------------------------------------------------
    # parallel version----------------------------------------------------
    do_work = partial(df_depth, dfmig=dfmig, boxes=boxes, depths=depths, pos_dist=pos_dist, a=stack, b=count)

    results = []
    pool = multiprocessing.Pool()
    for res in tqdm(pool.imap_unordered(do_work, depths), total=len(depths)):
        results.append(res)
    pool.close()
    pool.join()
    for res in results:
        stack += res[0]
        count += res[1]
    
    # end parallel version----------------------------------------------------

    d = np.array(depths)
    p = np.array(pos_dist)

    return p, d, stack


def df_depth(depth, dfmig, boxes, depths, pos_dist, a, b):

    stack = np.copy(a)
    count = np.copy(b)

    df2 = dfmig[dfmig["depth"] == depth]

    lats = df2["lat"].to_numpy()
    lons = df2["lon"].to_numpy()
    data = df2["data"].to_numpy()

    for i, lat in enumerate(lats):
        lon = lons[i]
        ppoint = (lat, lon)
        box = _find_box_cartesian(ppoint, boxes)
        if box is None:
            continue
        pos = box['pos']

        idep = depths.index(depth)
        ipos = pos_dist.index(pos)
        # print(pos, ipos, ppoint)

        stack[idep][ipos] += data[i]
        count[idep][ipos] += 1
        
    return (stack, count)
    