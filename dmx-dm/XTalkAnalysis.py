# imports
import math
import os

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gt
import readData as rddt


# TODO: Add the model / module name on the plots

def get_dump(path, c_perp, pls_set, config):

    # Searching directories
    dir_extension = "{0:}_".format(c_perp) + "pulseShape-{0:}_".format(pls_set) + config
    dirlist = [d for d in os.listdir(path) \
             if os.path.isdir(os.path.join(path, d)) \
             and d[-25:] == dir_extension]

    if len(dirlist) != 1:
        print("Error, wrong number of directories. Found {0:}".format(len(dirlist)))

    pathDir = os.path.join(path, dirlist[0])
    pathData = os.path.join(pathDir, cst.dataDirName)

    # Searching dump files
    files = [f for f in os.listdir(pathData) \
             if os.path.isfile(os.path.join(pathData, f)) \
             and f[-5:] == ".fits" and f[:5] == "dump_"]

    # Reading and accumulating the dumps
    print("Reading and accumulating the data from {0} dumps".format(len(files)))
    dump = np.zeros((cst.nColPerDemux, 2 * cst.nSamplesPerRow * cst.muxFactor))
    for file in files:
        idump, _ = rddt.read_dump_file(os.path.join(pathData, file))
        dump += idump
    dump /= len(files)

    # Setting the baseline to 0
    for col in range(cst.nColPerDemux):
        # Defining the reference level (expected to be on pixels 1 to 34)
        margin = 5
        first_sample = cst.nSamplesPerRow + margin
        last_sample = cst.nSamplesPerRow * (cst.muxFactor - 1) - margin
        ref_level = dump[col, first_sample : last_sample].mean()

        dump[col, :] = np.abs(dump[col, :] - ref_level)

    return dump


def xtalk_analysis(path, c_perp, pls_set, config):

    dump = get_dump(path, c_perp, pls_set, config)
    dump_ref = get_dump(path, c_perp, pls_set, "bob-looped")

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

#------------------------------------------------------------------

path = ["/Users/laurent/Data/TestPlan21-perfo/20250217_183000_XTalk"
        ]
conf = ["bob-looped", "not-loaded", "harness-08"]

pls_sets = [0, 1]

for p in path:
    for config in conf:
        for pls_set in pls_sets:
            for cPerp in range(cst.nColPerDemux):
                xtalk_analysis(p, cPerp, pls_set, config)
