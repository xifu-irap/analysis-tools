# imports
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

import constants as cst
import general_tools as gt
import readData as rddt

colors = ['#FF0000', '#0000FF', '#006400', '#FFA500', '#800080', '#008B8B', '#8B0000',
 '#696969', '#A52A2A', '#000000', '#4682B4', '#556B2F', '#4B0082', '#B8860B',
 '#C0C0C0', '#CD5C5C', '#8B008B', '#FFD700', '#228B22', '#00008B']


# TODO: Add the model / module name on the plots

def samplingDelayAnalysis(verbose):
    # Paths definition
    dir_path = os.path.join("..", "..")
    hk_path = os.path.join(dir_path, cst.hkDirName)
    data_path = os.path.join(dir_path, cst.dataDirName)
    plot_path = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(plot_path)
    session_name = os.path.realpath(dir_path)

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(hk_path)

    # Looking for test configuration parameters
    index_bxl = session_name.find("BXL")
    bxl = int(session_name[index_bxl + 3:index_bxl + 3 + 2])  # boxcar length
    index_pls = session_name.find("PLS+3")
    pls = int(session_name[index_pls])  # pulse shaping set

    if verbose:
        print("/----------------------------------------------------------")
        print("/ Feedback delay test ")
        print("/ Test session name:   " + session_name)
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:     {0:}".format(fwVersion))
        print("/ Box car length:       {0:} samples".format(bxl))
        print("/ Pulse shaping set:    {0:}".format(pls))
        print("/----------------------------------------------------------\n")

    xlabel = 'Time (ns)'
    ylabel1 = 'Error signal Dump (ADU)'

    x_min = 0
    x_max = 320

    t = np.arange(2*cst.nSamplesPerRow*cst.muxFactor) * 1e9 / cst.fSamp

    for colId in range(cst.nColPerDemux):

        plotFileName = os.path.join(plot_path, 'samplingDelay_col{0:}_bxl{1:}_pls{2:}.png'.format(colId, bxl, pls))
        title1 = 'Column {0:}, nb of averaged samples = {1:}, pulse shaping set = {2:}\n'.format(colId, bxl + 1, pls) \
                 + '(' + session_name + ')'

        fig = plt.figure(figsize=(9, 7))
        plt.suptitle("Test of Sampling delay and boxcar length settings".format(colId))

        # #######
        # Processing the dumps
        # #######

        # Searching dump files
        files = [f for f in os.listdir(data_path) \
                 if os.path.isfile(os.path.join(data_path, f)) \
                 and f[-3:] == ".h5" and f[:5] == "dump_"]

        # Reading and accumulating the dumps
        dump = np.zeros(2 * cst.nSamplesPerRow * cst.muxFactor)
        for file in files:
            colDumps, _ = rddt.read_dump_from_hdf5(os.path.join(data_path, file))
            dump += colDumps[colId, :]
        dump /= len(files)

        # Simulating the impact of the boxcar
        dump_avg = dump.copy()
        for i in range(bxl):
            dump_avg += np.roll(dump, -1 * (i + 1))
        dump_avg /= (bxl + 1)

        # Doing the plot
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.plot(t, dump, color='k', label='Dump data')
        if bxl > 0:
            ax1.plot(t, dump_avg, ':', color='k', label='Dump with moving average')

        # #######
        # Processing the error data files
        # #######

        # Searching error files
        suffix = "_C{0:}.h5".format(colId)
        files = [f for f in os.listdir(data_path) \
                 if os.path.isfile(os.path.join(data_path, f)) \
                 and f[:13] == 'SamplingDelay' and f[-6:] == suffix]

        if len(files) != cst.nSamplesPerRow:
            raise ValueError('Wrong number of files ({0:} instead of {1:})'.format(len(files), cst.nSamplesPerRow))

        nVal = 256  # Number of values (number of frames) to average over.
        # Each pixel has a constant value ( + noise ) along this test

        for index, file in enumerate(np.sort(files)):
            if verbose:
                print("  {0:} > ".format(index) + file)

            colData, _ = rddt.get_science_from_hdf5(os.path.join(data_path, file))
            pix0 = colData[0, :nVal].mean()

            # Adding error data on the plot
            ax1.plot(t[index], pix0, 'o', color=colors[index], label='Samp Delay = {0:} ns'.format(int(t[index])))

        ax1.set_ylabel(ylabel1)
        ax1.set_title(title1)
        ax1.set_xlabel(xlabel)
        ax1.legend(loc='upper right')
        ax1.set_xlim([x_min, x_max])

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
        if verbose:
            print("results plotted in file " + plotFileName)

# -------------------------------------------------------------------------------------
