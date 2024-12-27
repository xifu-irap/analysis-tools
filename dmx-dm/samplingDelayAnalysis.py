# imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os
import readData as rddt
import constants as cst
import general_tools as gt

path = ["/Users/laurent/Data/TestPlan21-perfo/20241024_140448_samplingDelay-col2-Bxl00",
        "/Users/laurent/Data/TestPlan21-perfo/20241024_140537_samplingDelay-col2-Bxl01",
        "/Users/laurent/Data/TestPlan21-perfo/20241024_140613_samplingDelay-col2-Bxl03",
        "/Users/laurent/Data/TestPlan21-perfo/20241024_140642_samplingDelay-col2-Bxl07",
        "/Users/laurent/Data/TestPlan21-perfo/20241029_102042_samplingDelay-col2-Bxl07",
        "/Users/laurent/Data/TestPlan21-perfo/20241029_102328_samplingDelay-col0-Bxl07",
        "/Users/laurent/Data/TestPlan21-perfo/20241114_105823_samplingDelay-col2-Bxl07",
        "/Users/laurent/Data/TestPlan21-perfo/20241114_155055_samplingDelay-col3-Bxl01",
        "/Users/laurent/Data/TestPlan21-perfo/20241118_174434_samplingDelay-col2-Bxl01",
        "/Users/laurent/Data/TestPlan21-perfo/20241220_115812_samplingDelay-col3-Bxl03",
        "/Users/laurent/Data/TestPlan21-perfo/20241220_115846_samplingDelay-col3-Bxl01"
        ]


colors = ['#FF0000', '#0000FF', '#006400', '#FFA500', '#800080', '#008B8B', '#8B0000',
 '#696969', '#A52A2A', '#000000', '#4682B4', '#556B2F', '#4B0082', '#B8860B',
 '#C0C0C0', '#CD5C5C', '#8B008B', '#FFD700', '#228B22', '#00008B']


def samplingDelay(path):

    col = int(path[-7])  # column index
    bxl = int(path[-2:]) # boxcar length
    pathData = os.path.join(path, "dmx_data")
    plotDirname = "PLOTS"
    pathPlot = os.path.join(path, plotDirname)
    gt.createdir(pathPlot)
    plotFileName = os.path.join(pathPlot, 'samplingDelay_col{0:}.png'.format(col))

    xlabel = 'Time (ns)'
    ylabel1 = 'Error signal Dump (ADU)'.format(col)
    title1 = 'Column {0:}, Boxcar length = {1:}'.format(col, bxl)
    x_max = 320

    t = np.arange(2*cst.nSamplesPerRow*cst.muxFactor) * 1e9 / cst.fSamp

    fig = plt.figure(figsize=(8, 6))
    plt.suptitle("Test of Sampling delay and boxcar length settings".format(col))

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
    # Simulating the impact of the boxcar
    dump_avg = dump.copy()
    for i in range(bxl):
        dump_avg += np.roll(dump, -1*(i+1))
        print(np.roll(dump, -1*(i+1))[:4])
    dump_avg /= (bxl+1)

    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(t, dump_avg, color='k', label='Dump data')

    #---------------------------------------------
    # Processing Error data
    print("  Processing error files")

    suffix = "_C{0:}.fits".format(col)

    files = [f for f in os.listdir(pathData) \
             if os.path.isfile(os.path.join(pathData, f)) \
             and f[-8:] == suffix]

    if len(files) != cst.nSamplesPerRow:
        raise ValueError('Wrong number of files ({0:} instead of {1:})'.format(len(files), cst.nSamplesPerRow))

    nVal = 256
    ratio = 4 # Ratio between dump and error data (due to unused LSB in error data)

    for index, file in enumerate (np.sort(files)):
        print("  {0:} > ".format(index) + file)

        colData, ctrl = rddt.readScienceFile(os.path.join(pathData, file))
        print(colData[0, :4])
        pix0 = colData[0, :nVal].mean() / ratio
        print(pix0)

        #if index == 14:
            #print(colData[0, :nVal])
            #print(colData[33, :nVal])

        ax1.plot(t[index], pix0, 'o', color=colors[index], label='Samp Delay = {0:} ns'.format(int(index*1e9/cst.fSamp)))

    ax1.set_ylabel(ylabel1)
    ax1.set_title(title1)
    ax1.set_xlabel(xlabel)
    ax1.legend(loc='best')
    ax1.set_xlim([0, x_max])

    # Définition des intervalles majeurs et mineurs pour la grille
    #ax1.set_xticks(np.arange(0, x_max, 8))  # Intervalles majeurs tous les 8

    # Définir le pas de la grille majeure
    ax1.xaxis.set_major_locator(MultipleLocator(32))
    ax1.yaxis.set_major_locator(MultipleLocator(1000))

    # Définir le pas de la grille mineure
    ax1.xaxis.set_minor_locator(MultipleLocator(8))
    ax1.yaxis.set_minor_locator(MultipleLocator(500))

    # Activation de la grille majeure
    ax1.grid(True, which='major', linestyle='-', linewidth='0.8', color='black')  # Grille majeure

    # Activer la grille mineure
    ax1.minorticks_on()  # Activer les ticks mineurs
    ax1.grid(True, which='minor', linestyle=':', linewidth='0.5', color='gray')

    fig.tight_layout()

    plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
    print("results plotted in file " + plotFileName)

for p in path[-2:]:
#for p in path:
    samplingDelay(p)
