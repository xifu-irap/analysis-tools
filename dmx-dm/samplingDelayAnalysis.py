# imports
import os
from dataclasses import dataclass

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

def samplingDelay(tconf, col):
    path = tconf.file_path
    session_name = tconf.session_name

    bxl = int(path[-19:-17]) # boxcar length
    pls = int(path[-1]) # pulse shaping set
    pathData = os.path.join(path, cst.dataDirName)
    pathPlot = os.path.join(path, cst.plotDirName)
    gt.createdir(pathPlot)
    plotFileName = os.path.join(pathPlot, 'samplingDelay_col{0:}_bxl{1:}_pls{2:}.png'.format(col, bxl, pls))

    xlabel = 'Time (ns)'
    ylabel1 = 'Error signal Dump (ADU)'
    title1 = 'Column {0:}, nb of averaged samples = {1:}, pulse shaping set = {2:}\n'.format(col, bxl + 1, pls) \
             + '(' + tconf.testPlanPath + '     ' + session_name + ')'

    x_min = 0
    x_max = 320

    t = np.arange(2*cst.nSamplesPerRow*cst.muxFactor) * 1e9 / cst.fSamp

    fig = plt.figure(figsize=(8, 6))
    plt.suptitle("Test of Sampling delay and boxcar length settings".format(col))

    print("Processing data from directory " + path)
    #---------------------------------------------

    # #######
    # Processing the dumps
    # #######

    print("  Processing dump files")

    # Searching dump files
    files = [f for f in os.listdir(pathData) \
             if os.path.isfile(os.path.join(pathData, f)) \
             and f[-5:] == ".fits" and f[:5] == "dump_"]

    # Reading and accumulating the dumps
    dump = np.zeros(2 * cst.nSamplesPerRow * cst.muxFactor)
    for file in files:
        colDumps, _ = rddt.read_dump_file(os.path.join(pathData, file))
        dump+=colDumps[col, :]
    dump /= len(files)

    # Simulating the impact of the boxcar
    dump_avg = dump.copy()
    for i in range(bxl):
        dump_avg += np.roll(dump, -1*(i+1))
    dump_avg /= (bxl+1)

    # Doing the plot
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(t, dump, color='k', label='Dump data')
    if bxl > 0:
        ax1.plot(t, dump_avg, ':', color='k', label='Dump with moving average')

    # #######
    # Processing the error data files
    # #######

    print("  Processing error files")

    # Searching error files
    suffix = "_C{0:}.fits".format(col)
    files = [f for f in os.listdir(pathData) \
             if os.path.isfile(os.path.join(pathData, f)) \
             and f[:13] == 'SamplingDelay' and f[-8:] == suffix]

    if len(files) != cst.nSamplesPerRow:
        raise ValueError('Wrong number of files ({0:} instead of {1:})'.format(len(files), cst.nSamplesPerRow))

    nVal = 256 # Number of values to average over
    ratio = 4 # Ratio between dump and error data (due to unused LSB in error data)

    for index, file in enumerate (np.sort(files)):
        print("  {0:} > ".format(index) + file)

        colData, ctrl = rddt.read_science_file(os.path.join(pathData, file))
        pix0 = colData[0, :nVal].mean() / ratio

        # Adding error data on the plot
        ax1.plot(t[index], pix0, 'o', color=colors[index], label='Samp Delay = {0:} ns'.format(int(t[index])))

    ax1.set_ylabel(ylabel1)
    ax1.set_title(title1)
    ax1.set_xlabel(xlabel)
    ax1.legend(loc='upper right')
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


# -------------------------------------------------------------------------------------

@dataclass
class TestConfig:
    testPlanPath: str
    session_name: str

    @property
    def file_path(self) -> str:
        return f"{cst.BASE_DATA_PATH}/{self.testPlanPath}/{self.session_name}"


list_of_configs = [
    TestConfig(cst.TP21_PATH, "20250213_144452_samplingDelay-Bxl07-pulseShapingSet0"),
    TestConfig(cst.TP21_PATH, "20250213_144532_samplingDelay-Bxl03-pulseShapingSet0"),
    TestConfig(cst.TP21_PATH, "20250213_144611_samplingDelay-Bxl01-pulseShapingSet0"),
    TestConfig(cst.TP21_PATH, "20250213_144648_samplingDelay-Bxl00-pulseShapingSet0")
]

for test_conf in list_of_configs:
    for col in range(cst.nColPerDemux):
        samplingDelay(test_conf, col)
