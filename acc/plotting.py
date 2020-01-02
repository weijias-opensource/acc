# -*- coding: utf-8 -*-
"""
Plottings
"""
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, FixedLocator, FixedFormatter, MaxNLocator)
from acc.stack import linear_stack
import numpy as np
import warnings
from obspy import Stream


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
