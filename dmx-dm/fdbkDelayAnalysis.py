# imports
import os

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gt
import readData as rddt

xZoom = 49

colors = ['#FF0000', '#0000FF', '#006400', '#FFA500', '#800080', '#008B8B', '#8B0000',
          '#696969', '#A52A2A', '#000000', '#4682B4', '#556B2F', '#4B0082', '#B8860B',
          '#C0C0C0', '#CD5C5C', '#8B008B', '#FFD700', '#228B22', '#00008B']


def fdbkDelay(col, verbose=False):
    # Data directory
    dir_path = os.path.join("..", "..")
    # dir_path = "../.."
    path_data = os.path.join(dir_path, cst.dataDirName)

    session_name = os.path.realpath(dir_path)
    start_str = session_name.split("_")[-2]
    end_str = session_name.split("_")[-1]
    start = int(start_str[1:])
    if start_str[0] == "M":
        start *= -1
    end = int(end_str[1:])
    if end_str[0] == "M":
        end *= -1
    nb_steps = end - start + 1

    # Creation of a directory for the plot files
    path_plot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(path_plot)
    plotFileName = os.path.join(path_plot, 'fdbkDelay_col{0:}.png'.format(col))
    xlabel = 'Time (ns)'
    ylabel = 'Dump of error signal (ADU)'.format(col)

    if verbose:
        print("Processing dump files of column {0:} from directory ".format(col) + path_data)
    files = [f for f in os.listdir(path_data) \
             if os.path.isfile(os.path.join(path_data, f)) \
             and f[-5:] == ".fits" and f[:5] == "dump_"]

    if len(files) != nb_steps:
        raise ValueError('Wrong number of files ({0:} instead of {1:})'.format(len(files), nb_steps))

    fig = plt.figure(figsize=(12, 10))
    ax1 = fig.add_subplot(2, 1, 1)  # global plot
    ax2 = fig.add_subplot(2, 1, 2)  # zoom plot
    plt.suptitle("Feedback delay compensation setting (column {0:})".format(col))

    xTime = np.arange(2 * cst.nSamplesPerRow * cst.muxFactor) * 1e9 / cst.fSamp
    for index, file in enumerate(np.sort(files)):
        if verbose:
            print(file)
        iDelay = index - 3
        colDumps, errors = rddt.read_dump_from_fits(os.path.join(path_data, file))

        # Doing plot
        setting = start + index
        lw = 2
        if setting == 0:
            lw = 4
        lbl = 'Delay = {0:2d} ns'.format(int(setting * 1e9 / cst.fSamp))
        ax1.plot(xTime[:], colDumps[col, :], color=colors[index], label=lbl, linewidth=1)
        ax2.plot(xTime[:xZoom], colDumps[col, :xZoom], color=colors[index], label=lbl, linewidth=lw)

    x1_max = 2 * cst.nSamplesPerRow * cst.muxFactor * 1e9 / cst.fSamp
    x2_max = (xZoom - 1) * 1e9 / cst.fSamp

    ax1.set_ylabel(ylabel)
    ax2.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)
    ax2.set_xlabel(xlabel)
    # ax1.legend(loc='upper right')
    ax2.legend(loc='upper right')
    xlims = ax2.get_xlim()
    ax2.set_xlim([xlims[0], x2_max])

    # Définition des intervalles majeurs et mineurs pour la grille
    ax1.set_xticks(np.arange(0, x1_max, 4096))  # Intervalles majeurs tous les 64
    ax1.set_xticks(np.arange(0, x1_max, 512), minor=True)  # Intervalles mineurs tous les 8
    ax2.set_xticks(np.arange(0, x2_max, 64))  # Intervalles majeurs tous les 64
    ax2.set_xticks(np.arange(0, x2_max, 8), minor=True)  # Intervalles mineurs tous les 8

    # Activation de la grille majeure et mineure
    ax1.grid(which='major', linestyle='-', linewidth='0.6', color='black')  # Grille majeure
    ax1.grid(which='minor', linestyle='--', linewidth='0.4', color='gray')  # Grille mineure
    ax2.grid(which='major', linestyle='-', linewidth='0.6', color='black')  # Grille majeure
    ax2.grid(which='minor', linestyle='--', linewidth='0.4', color='gray')  # Grille mineure

    fig.tight_layout()

    plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
    if verbose:
        print("results plotted in file " + plotFileName)


# -------------------------------------------------------------------------------------

def fdbkDelayAnalysis(verbose=False):
    for colid in range(cst.nColPerDemux):
        fdbkDelay(colid)
