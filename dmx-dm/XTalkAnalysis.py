# imports
import math
import os

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gt
import readData as rddt


# TODO: Add the model / module name on the plots

def get_dump(path, perp, c_perp, config, verbose=True):

    # Searching directories
    dir_extension = "-XTALK-PERP-" + perp + "-C{:0}".format(c_perp) + config
    dirlist = [d for d in os.listdir(path) \
             if os.path.isdir(os.path.join(path, d)) \
               and d[-1 * len(dir_extension):] == dir_extension]
    print(path)
    print(dir_extension)
    print(dirlist)

    if len(dirlist) != 1:
        print("Error, wrong number of directories. Found {0:}".format(len(dirlist)))

    pathDir = os.path.join(path, dirlist[0])
    pathData = os.path.join(pathDir, cst.dataDirName)

    # Searching dump files
    files = [f for f in os.listdir(pathData) \
             if os.path.isfile(os.path.join(pathData, f)) \
             and f[-3:] == ".h5" and f[:5] == "dump_"]
    if verbose:
        print("Accumulating {0:} dump files".format(len(files)))

    # Reading and accumulating the dumps
    dump = np.zeros((cst.nColPerDemux, 2 * cst.nSamplesPerRow * cst.muxFactor))
    for file in files:
        idump, _ = rddt.read_dump_from_hdf5(os.path.join(pathData, file))
        dump += idump
    dump /= len(files)

    # Setting the baseline to 0
    for col in range(cst.nColPerDemux):
        # Defining the reference level (expected to be on pixels 1 to TMux-1)
        margin = 15
        first_sample = cst.nSamplesPerRow + margin
        last_sample = cst.nSamplesPerRow * (cst.muxFactor - 1) - margin
        ref_level = dump[col, first_sample : last_sample].mean()

        dump[col, :] = np.abs(dump[col, :] - ref_level)

    return dump


def xtalk_analysis(path, c_perp, pls_set, config):

    dump = get_dump(path, c_perp, pls_set, config)
    dump_ref = get_dump(path, c_perp, pls_set, "fbk-looped")

    col_dump_db = 20*np.log10(dump / (dump_ref[c_perp, :].max() - dump_ref[c_perp, :].min()))

    pathPlot = os.path.join(path, cst.plotDirName)
    gt.createdir(pathPlot)

    plotFileName = os.path.join(pathPlot, 'XTalk_'+config+'_pls{0:}'.format(pls_set)+'_cPerp{0:}.png'.format(c_perp))

    # Doing the plot
    title = 'Crosstalk from column {0:} in '.format(c_perp) + config + ' mode with ' + gt.pulseshapingtext(pls_set)
    xlabel = 'Time (ns)'
    ylabel = 'Error signal Dump (ADU)'
    ylabel1 = 'Ratio with max of column {0:} (dB)'.format(c_perp)
    col1 = 'k'
    col2, col2grid = 'b', 'lightblue'
    t = np.arange(2*cst.nSamplesPerRow*cst.muxFactor) * 1e9 / cst.fSamp

    fig = plt.figure(figsize=(12, 10))
    plt.suptitle(title)


    nb_plot_cols = 2 # number of columns in the plot
    for col in range(cst.nColPerDemux):

        ax1 = fig.add_subplot(int(cst.nColPerDemux/nb_plot_cols), nb_plot_cols, col+1)
        ax1.plot(t, dump[col, :], color=col1)
        tit = "Column {0:}".format(col)
        if col == c_perp:
            tit += ' (perpetrator)'
        ax1.set_title(tit)
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel, color=col1)
        nsamples = 80
        #xlims = [0, nsamples*1e9/cst.fSamp]
        #ax1.set_xlim(xlims)
        [t.set_color(col1) for t in ax1.yaxis.get_ticklabels()]

        ax2 = ax1.twinx()
        ax2.plot(t, col_dump_db[col, :], color=col2)
        maxi = col_dump_db[col, :].max()
        if math.isnan(maxi) or math.isinf(maxi):
            lbl = ''
        else:
            lbl = '{0:} dB'.format(int(maxi))
        ax2.plot([t[0], t[-1]], [maxi, maxi], ':', color=col2, label=lbl)
        ax2.set_ylabel(ylabel1, color=col2)
        [t.set_color(col2) for t in ax2.yaxis.get_ticklabels()]
        ax2.grid(color=col2grid)
        #ax2.set_xlim(xlims)

        ax2.legend(loc='best')
        #ax2.set_xlim([0, x_max])

    fig.tight_layout()

    plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
    print("results plotted in file " + plotFileName)


def fdbk_2_error_xtalk_analysis(path, c_perp, pls_set):
    dump_no_fdbk = get_dump(path, c_perp, pls_set, "not-loaded")[c_perp, :]
    dump_fdbk = get_dump(path, c_perp, pls_set, "fbk-looped")[c_perp, :]

    xtalk_db = 20 * np.log10(np.abs(dump_no_fdbk / dump_fdbk))

    pathPlot = os.path.join(path, cst.plotDirName)
    gt.createdir(pathPlot)

    plotFileName = os.path.join(pathPlot, 'XTalk_fdbk2error_pls{0:}'.format(pls_set) + '_cPerp{0:}.png'.format(c_perp))

    # Doing the plot
    t = np.arange(2 * cst.nSamplesPerRow * cst.muxFactor) * 1e9 / cst.fSamp
    title = 'Crosstalk in column {0:} between feedback and error'.format(c_perp) + ' with ' + gt.pulseshapingtext(
        pls_set)
    xlabel = 'Time (ns)'
    title1 = 'Plot 1: The feedback is looped in the Error'
    ylabel1 = 'Error (ADU)'
    title2 = 'Plot 2: The feedback is not looped in the Error'
    ylabel2 = 'Error (ADU)'
    title3 = 'Plot 3: Cross talk'
    ylabel3 = 'Plot 2 / Plot 1 (dB)'
    index_xmean = [5, 20]
    xmean = t[index_xmean]
    index_xlims = [0, 25]
    xlims = t[index_xlims]
    col1, col2, col3, col2grid = 'k', 'b', 'r', 'lightblue'

    fig = plt.figure(figsize=(10, 10))
    plt.suptitle(title)

    ax1 = fig.add_subplot(3, 1, 1)
    ax1.plot(t, dump_fdbk, color=col1)
    ax1.set_title(title1)
    ax1.set_ylabel(ylabel1, color=col1)
    ax1.set_xlim(xlims)
    ax1.grid(color=col2grid)

    ax2 = fig.add_subplot(3, 1, 2)
    ax2.plot(t, dump_no_fdbk, color=col2)
    ####
    ax2.set_title(title2)
    ax2.set_ylabel(ylabel2, color=col2)
    ax2.set_xlim(xlims)
    ax2.grid(color=col2grid)

    ax3 = fig.add_subplot(3, 1, 3)
    ax3.plot(t, xtalk_db, color=col3)
    xtalk = xtalk_db[index_xmean[0]:index_xmean[1]].mean()
    ax3.set_title(title3)
    ax3.plot(xmean, [xtalk, xtalk], '--', color='k', label='Cross talk: {0:3.0f} dB'.format(xtalk))
    ax3.set_xlabel(xlabel)
    ax3.set_ylabel(ylabel3, color=col3)
    ax3.set_xlim(xlims)
    ax3.grid(color=col2grid)
    ax3.legend(loc='best')

    fig.tight_layout()

    plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
    print("results plotted in file " + plotFileName)

#------------------------------------------------------------------

# TODO: Clarify the Xtalk test configuration and update the analysis script accordingly

# conf = ["bob-looped", "not-loaded", "harness-08"]
conf = ["fbk-looped", "not-loaded", "harness-08"]

pls_sets = [0, 1]

list_of_paths = [
    os.path.join(cst.BASE_DATA_PATH, cst.TP27_PATH, "20251117_110000_XTalk"),
    os.path.join(cst.BASE_DATA_PATH, cst.TP27_PATH, "20251118_160000_XTalk")
]

for path in list_of_paths:
    for pls_set in pls_sets:
        for cPerp in range(cst.nColPerDemux):
            fdbk_2_error_xtalk_analysis(path, cPerp, pls_set)
            # for config in conf:
            #xtalk_analysis(path, cPerp, pls_set, config)
