# imports
import numpy as np
import matplotlib.pyplot as plt
import os
import readData as rddt
import constants as cst
import general_tools as gt

path = ["/Users/laurent/Data/TestPlan21-perfo/20241021_183800_samplingDelay-col3", \
    "/Users/laurent/Data/TestPlan21-perfo/20241021_183825_samplingDelay-col3", \
    "/Users/laurent/Data/TestPlan21-perfo/20241021_183950_samplingDelay-col2" ]

colors = ['#FF0000', '#0000FF', '#00FF00', '#FFA500', '#800080', '#00FFFF', '#FFFF00', '#FFC0CB', \
        '#A52A2A', '#808080', '#000000', '#FFFFFF', '#40E0D0', '#4B0082', '#FFD700', '#C0C0C0', '#F0E68C', \
        '#FF7F50', '#FF00FF', '#7FFF00']

def debug(path):

    col = int(path[-1])
    pathData = os.path.join(path, "dmx_data")
    plotDirname = "PLOTS_DEBUG"
    pathPlot = os.path.join(path, plotDirname)
    gt.createdir(pathPlot)
    plotFileName1 = os.path.join(pathPlot, 'dump_col{0:}.png'.format(col))
    plotFileName2 = os.path.join(pathPlot, 'error_col{0:}.png'.format(col))

    xlabel = 'Time (ns)'
    ylabel1 = 'Error signal Dump (ADU)'.format(col)
    title1 = 'Dump data'
    x_max = 320

    t = np.arange(2*cst.nSamplesPerRow*cst.muxFactor) * 1e9 / cst.fSamp

    print("Processing data from directory " + path)
    #---------------------------------------------
    # Processing Dump data
    print("  Processing dump files")

    files = [f for f in os.listdir(pathData) \
             if os.path.isfile(os.path.join(pathData, f)) \
             and f[-5:] == ".fits" and f[:5] == "dump_"]

    dump = np.zeros(2 * cst.nSamplesPerRow * cst.muxFactor)
    for file in files:
        colDumps, errors = rddt.readDumpFile(os.path.join(pathData, file))
        dump+=colDumps[col, :]
    dump /= len(files)

    fig1 = plt.figure(figsize=(6, 8))
    plt.suptitle("Debug (column {0:})".format(col))

    ax1 = fig1.add_subplot(1, 1, 1)
    ax1.plot(t, dump, color='k', label='Dump data')

    fig1.tight_layout()

    plt.savefig(plotFileName1, dpi=300, bbox_inches='tight')
    print("results plotted in file " + plotFileName1)


    #---------------------------------------------
    # Processing Error data

    suffix = "_C{0:}.fits".format(col)

    files = [f for f in os.listdir(pathData) \
             if os.path.isfile(os.path.join(pathData, f)) \
             and f[-8:] == suffix]
    print("  Processing error file ", files[0])

    if len(files) != cst.nSamplesPerRow:
        raise ValueError('Wrong number of files ({0:} instead of {1:})'.format(len(files), cst.nSamplesPerRow))


    colData, ctrl = rddt.readScienceFile(os.path.join(pathData, files[0]), perPix=True)

    print(np.where(colData > 22000))
    print(colData[0,:])
    print(colData[20,:])


debug(path[0])
#for p in path:
#    samplingDelay(p)