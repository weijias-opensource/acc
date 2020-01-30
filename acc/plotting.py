# -*- coding: utf-8 -*-
"""
Plottings
"""
import matplotlib.pyplot as plt
from matplotlib.ticker import (FormatStrFormatter, MultipleLocator, AutoMinorLocator, FixedLocator, FixedFormatter, MaxNLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.markers import MarkerStyle
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from acc.stack import linear_stack
from acc.profile import profile, get_profile_boxes
from acc.io import _load_json
from acc.migration import mig_one_station

import numpy as np
import warnings
from obspy import Stream
from obspy import read
from obspy.geodetics import gps2dist_azimuth

from collections import OrderedDict

import os
import glob
import math


def plot_acc_one_station(stream, fname=None, fig_width=7., trace_height=0.5,
                         stack_height=0.5, scale=2, scale_stack=10, fillcolors=("red", "blue"),
                         # trim=None,
                         info=(('back_azimuth', u'baz (°)', 'C0'),
                               ('distance', u'dist (°)', 'C3'))):
    """
    Plot auto- or correlogram for one station. Reproduced from the `rf` package.

    :param stream: stream to plot
    :param fname: filename to save plot to. Can be None. In this case
        the figure is left open.
    :param fig_width: width of figure in inches
    :param trace_height: height of one trace in inches
    :param stack_height: height of stack axes in inches
    :param scale: scale for individual traces
    :param fillcolors: fill colors for positive and negative wiggles
    :param info: Plot one additional axes showing maximal two entries of
        the stats object. Each entry in this list is a list consisting of
        three entries: key, label and color.
        info can be None. In this case no additional axes is plotted.
    """
    # :param trim: trim stream relative to onset before plotting using
    #      `~.rfstream.RFStream.slice2()`

    if len(stream) == 0:
        return

    stream.sort(keys=["slowness"])

    # if trim:
    #     stream = stream.slice2(*trim, reftime='onset')
    N = len(stream)
    # calculate lag times
    # stats = stream[0].stats
    # print(stats)
    # times = stream[0].times() - (stats.onset - stats.starttime)
    times = stream[0].times()
    # calculate axes and figure dimensions
    # big letters: inches, small letters: figure fraction
    H = trace_height
    HS = stack_height
    FB = 0.5
    FT = 0.2
    DW = 0.2
    FH = H * (N + 2) + HS + FB + FT + DW
    h = H / FH
    hs = HS / FH
    fb = FB / FH
    ft = FT / FH
    FL = 0.5 # figure left
    FR = 0.2 # figure right
    FW = fig_width # figure width
    FW3 = 0.8
    FW2 = FW - FL - FR - (DW + FW3) * bool(info)
    fl = FL / FW
    fr = FR / FW
    fw2 = FW2 / FW
    fw3 = FW3 / FW
    # init figure and axes
    fig = plt.figure(figsize=(FW, FH))
    ax1 = fig.add_axes([fl, fb, fw2, h * (N + 2)])
    if info:
        ax3 = fig.add_axes(
            [1 - fr - fw3, fb, fw3, h * (N + 2)], sharey=ax1)
        info = list(info)
        info[0] = [ax3] + list(info[0])
        if len(info) > 1:
            ax4 = ax3.twiny()
            info[1] = [ax4] + list(info[1])
    # plot individual receiver functions

    def _plot(ax, t, d, i):
        c1, c2 = fillcolors
        if c1:
            ax.fill_between(t, d + i, i, where=d >= 0, lw=0., facecolor=c1)
        if c2:
            ax.fill_between(t, d + i, i, where=d < 0, lw=0., facecolor=c2)
        ax.plot(t, d + i, 'k')
    # max_ = max(np.max(np.abs(tr.data)) for tr in stream)
    for i, tr in enumerate(stream):
        # _plot(ax1, times, tr.data / max_ * scale, i + 1)
        tr.normalize()
        _plot(ax1, times, tr.data * scale, i + 1)
    # plot right axes with header information
    for ax, header, label, color in info:
        data = [tr.stats[header] for tr in stream]
        ax.plot(data, 1 + np.arange(len(stream)), '.' + color, mec=color)
        ax.set_xlabel(label, color=color, size='small')
        if header == 'back_azimuth':
            ax.set_xticks(np.arange(5) * 90)
            ax.set_xticklabels(['0', '', '180', '', '360'], size='small')
        else:
            ax.xaxis.set_major_locator(MaxNLocator(4))
            for l in ax.get_xticklabels():
                l.set_fontsize('small')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    # set x and y limits
    ax1.set_xlim(times[0], times[-1])
    ax1.set_xlim(times[0], times[-1]/2)
    ax1.set_ylim(-0.5, N + 1.5)
    ax1.set_yticklabels('')
    ax1.set_xlabel('time (s)')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())

    # plot stack
    # stack = stream.stack()
    trace = linear_stack(stream=stream, normalize=True)
    stack = Stream([trace])
    if len(stack) > 1:
        warnings.warn('Different stations or channels in one RF plot.')
    elif len(stack) == 1:
        stack.normalize()
        ax2 = fig.add_axes([fl, 1 - ft - hs, fw2, hs], sharex=ax1)
        _plot(ax2, times, stack[0].data * scale_stack, 0)
        for l in ax2.get_xticklabels():
            l.set_visible(False)
        ax2.yaxis.set_major_locator(MaxNLocator(4))
        for l in ax2.get_yticklabels():
            l.set_fontsize('small')
        # annotate plot with seed id
        bbox = dict(boxstyle='round', facecolor='white', alpha=0.8, lw=0)
        text = '%s traces  %s' % (len(stream), stack[0].id)
        ax2.annotate(text, (1 - 0.5 * fr, 1 - 0.5 * ft),
                     xycoords='figure fraction', va='top', ha='right',
                     bbox=bbox, clip_on=False)
        ax2.set_ylim(-1., 1.)
    # save plot
    if fname:
        fig.savefig(fname)
        plt.close(fig)
    else:
        return fig

def plot_ppoint(stream, fname="pp.pdf", depths=[30, 50, 80, 100, 150, 200]):

    colors = ["gray", "red", "blue", "orange", "green", "magenta", "cyan", "chocolate", "pink", "royalblue"]
    if len(depths) > 10:
        raise Exception("too many depths. Should be less than 10.")

    fig, ax = plt.subplots()

    for tr in stream:
        df = tr.stats.mig
        ax.scatter(tr.stats.station_longitude, tr.stats.station_latitude, c="black", marker="v", s=100)

        for i, depth in enumerate(depths):
            df2 = df[df["depth"] == depth]

            lat = df2["lat"].to_numpy()
            lon = df2["lon"].to_numpy()

            ax.scatter(lon, lat, c=colors[i], label=depth)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), title="Depth (km)")

    plt.tight_layout()
    plt.savefig(fname)



def plot_profile(jsonfile):

    kwargs = _load_json(jsonfile)
    paras = kwargs["profile"]

    latlon0 = paras["latlon0"]
    if latlon0 is None:
        # for a single station
        azimuth = 0
        dist = 400
    else:
        # for profile
        latlon1 = paras["latlon1"]
        lat1 = latlon0[0]
        lon1 = latlon0[1]
        lat2 = latlon1[0]
        lon2 = latlon1[1]
        dist, azimuth, baz = gps2dist_azimuth(lat1, lon1, lat2, lon2)
        dist /= 1000  # from m to km
        # padding 100 km on each side
        dist += 400

        lat0 = lat1 + math.cos((azimuth + 180)*math.pi/180) * 200 / 111.195
        lon0 = lon1 + math.sin((azimuth + 180)*math.pi/180) * 200 / 111.195
        latlon0 = (lat0, lon0)

    binsize = paras["binsize"]
    binwidth = paras["binwidth"]

    nbins = int(dist / binsize + 0.5) + 1
    dist = (nbins - 1) * binsize
    bins = np.linspace(0, dist, nbins)
    # print(bins)

    profile_id = paras["profile_id"]
    clip = paras["clip_percentage"] / 100.0

    path = kwargs["io"]["outpath"]

    depth_range = kwargs["plot"]["depth_range"]

    # each stations
    if latlon0 is None:
        files = glob.glob(path + "/migration_1station/*.pkl")
        files.sort()
        for file in files:
            stmig = read(file)
            lat0 = stmig[0].stats.station_latitude - 0.5 * (bins[-1] + bins[0]) / 111.195
            lon0 = stmig[0].stats.station_longitude
            latlon0 = (lat0, lon0)
            _plot_1station(stmig, latlon0, azimuth, bins, width=binwidth, clip=clip, savepath=path,
                           depth_range=depth_range)
    # profile
    else:
        wc = paras["wild_card"]
        files = glob.glob(path + "/migration_1station/*%s*.pkl" % wc)
        print("number of stations to stack: ", len(files))
        stmig = read(path + "/migration_1station/*%s*.pkl" % wc)
        _plot_stations(stmig, latlon0, azimuth, bins, width=binwidth, clip=clip, savepath=path,
                       depth_range=depth_range, profile_id=profile_id)


def _plot_1station(stmig, latlon0, azimuth, bins, width, clip, savepath, depth_range):
    # get boxes for mig-stacking
    boxes = get_profile_boxes(latlon0, azimuth=azimuth, bins=bins, width=width)

    # depth and stack array
    pos, depth, stack = profile(stream=stmig, boxes=boxes)

    # the setting makes the station location in the center of the image
    pos -= 0.5 * (bins[-1] - bins[0])

    dist_range = [-100, 100]

    # get station id
    tr = stmig[0]
    if tr.stats.location == "":
        station_id = ".".join([tr.stats.network, tr.stats.station])
    else:
        station_id = ".".join([tr.stats.network, tr.stats.station, tr.stats.location])
    print("Plotting - station id: ", station_id)

    # latitude corresponding to pos
    lats = []
    for b in boxes:
        # print(b["pos"], b["latlon"])
        lats.append(b["latlon"][0])

    extent = [np.min(pos), np.max(pos), np.min(depth), np.max(depth)]
    # extent = [np.min(lats), np.max(lats), np.min(depth), np.max(depth)]

    # normalization
    stack /= np.max(np.abs(stack))

    amp = np.sum(stack, axis=1)
    amp /= np.max(np.abs(amp))

    # print(extent)

    # return amp, stack, extent
    _plot(amp, depth, stack, extent, dist_range, depth_range, savepath, profile_id=station_id, clip=clip)


def _plot_stations(stmig, latlon0, azimuth, bins, width, clip, savepath, depth_range, profile_id):
    # get boxes for mig-stacking
    boxes = get_profile_boxes(latlon0, azimuth=azimuth, bins=bins, width=width)

    # depth and stack array
    pos, depth, stack = profile(stream=stmig, boxes=boxes)

    # the setting makes the station location in the center of the image
    # pos -= 0.5*(bins[-1]-bins[0])

    dist_range = [np.min(pos), np.max(pos)]

    # get station id
    # tr = stmig[0]
    # if tr.stats.location == "":
    #     station_id = ".".join([tr.stats.network, tr.stats.station])
    # else:
    #     station_id = ".".join([tr.stats.network, tr.stats.station, tr.stats.location])
    print("Plotting - profile id: ", profile_id)

    # latitude corresponding to pos
    lats = []
    for b in boxes:
        # print(b["pos"], b["latlon"])
        lats.append(b["latlon"][0])

    extent = [np.min(pos), np.max(pos), np.min(depth), np.max(depth)]
    # extent = [np.min(lats), np.max(lats), np.min(depth), np.max(depth)]

    # normalization
    stack /= np.max(np.abs(stack))

    amp = np.sum(stack, axis=1)
    amp /= np.max(np.abs(amp))

    # print(extent)

    # return amp, stack, extent
    _plot(amp, depth, stack, extent, dist_range, depth_range, savepath, profile_id=profile_id, clip=clip)


def _plot(amp, depth, stack, extent, dist_range, depth_range, savepath, profile_id, clip):
    # init figure and axes
    N = 256
    rdbu_r = cm.get_cmap('RdBu_r', N)
    newcolors = rdbu_r(np.linspace(0, 1, N))
    white = np.array([255 / 255, 255 / 255, 255 / 255, 1])
    a = int(N / 2 - N / 12)
    b = int(N / 2 + N / 12)
    newcolors[a:b, :] = white
    newcmp = ListedColormap(newcolors)

    fig, (ax2, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(4.5, 3), sharey=True,
                                   gridspec_kw={'width_ratios': [0.5, 2]})
    im = ax1.imshow(stack, interpolation='bilinear', aspect="auto",
                    origin="lower", extent=extent,
                    cmap=newcmp, vmin=-clip, vmax=clip)

    ax1.yaxis.set_major_locator(MultipleLocator(50))
    ax1.yaxis.set_minor_locator(MultipleLocator(10))
    ax1.xaxis.set_minor_locator(MultipleLocator(50))
    ax1.grid(which="both", axis="y", linestyle=":", linewidth=0.5)

    ax1.annotate(profile_id,
                 xy=(0.995, 0.995), xycoords='axes fraction',
                 # xytext=(0.95, 0.05), textcoords=station_id,
                 # arrowprops=dict(facecolor='black', shrink=0.05),
                 horizontalalignment='right', verticalalignment='top')

    ax1.set_ylim(depth_range)
    ax1.set_xlim(dist_range)
    ax1.invert_yaxis()
    # ax1.invert_xaxis()

    ax2.plot(amp, depth, linewidth=0.75, color="k")
    ax2.fill_betweenx(depth, 0, amp, where=amp >= 0, facecolor='red')
    ax2.fill_betweenx(depth, 0, amp, where=amp < 0, facecolor='blue')
    ax2.grid(which="both", axis="y", linestyle=":", linewidth=0.5)
    ax2.set_xlim([-0.1, 0.1])
    ax2.set_xlabel("Amplitude")
    ax2.set_ylabel("Depth [km]")

    cbar = fig.colorbar(im, ax=ax1)
    cbar.set_label("Amplitude")

    ax1.set_xlabel("CRP Distance [km]")
    # ax1.set_ylabel("Depth [km]")

    plt.tight_layout()
    path = savepath + "/figures"
    if not os.path.exists(path):
        os.makedirs(path)
    filen = path + "/mig_crp_" + profile_id
    for fmt in [".pdf", ".png"]:
        fn = filen + fmt
        plt.savefig(fn)
    plt.close()


