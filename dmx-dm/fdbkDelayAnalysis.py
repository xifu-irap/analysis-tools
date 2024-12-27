# imports
import numpy as np
import matplotlib.pyplot as plt
import os
import readData as rddt
import constants as cst

#path = "/Users/laurent/Data/TestPlan21-perfo/fdbkDelay-20241014-1/"
#path = "/Users/laurent/Data/TestPlan21-perfo/fdbkDelay-20241014-2/"
path = "/Users/laurent/Data/TestPlan21-perfo/fdbkDelay-20241017-1/"
col = 3
xZoom = 49

colors = ['#0000FF', '#006400', '#FFA500', '#FF0000', '#800080', '#000080', '#8B4513']

def fdbkDelay(path, col):

    print("Processing dump files from directory " + path)

    plotFileName = path + 'fdbkDelay_col{0:}.png'.format(col)
    xlabel = 'Time (ns)'
    ylabel = 'Dump of error signal (ADU)'.format(col)

    files = [f for f in os.listdir(path) \
             if os.path.isfile(os.path.join(path, f)) \
             and f[-5:] == ".fits" and f[:5] == "dump_"]

    if len(files) != 7:
        raise ValueError('Wrong number of files ({0:} instead of 7)'.format(len(files)))

    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(2, 1, 1)  # global plot
    ax2 = fig.add_subplot(2, 1, 2)  # zoom plot
    plt.suptitle("Feedback delay compensation setting (column {0:})".format(col))

    xTime = np.arange(2*cst.nSamplesPerRow*cst.muxFactor) * 1e9 / cst.fSamp
    for index, file in enumerate (np.sort(files)):
        print(file)
        iDelay = index-3
        colDumps, errors = rddt.readDumpFile(os.path.join(path, file))

        # Doing plot
        if iDelay < 0:
            sign = "- "
        else:
            sign = "+ "
        if iDelay == 0:
            lw = 2
        else:
            lw = 1
        lbl = 'Delay = Ref. ' + sign + '{0:2d} ns'.format(int(np.abs(iDelay)*1e9/cst.fSamp))
        ax1.plot(xTime[:], colDumps[col,:], color=colors[index], label=lbl, linewidth=lw)
        ax2.plot(xTime[:xZoom], colDumps[col,:xZoom], color=colors[index], label=lbl, linewidth=lw)

    x1_max = 2 * cst.nSamplesPerRow * cst.muxFactor * 1e9 / cst.fSamp
    x2_max = (xZoom - 1) * 1e9 / cst.fSamp

    ylims = ax1.get_ylim()
    y_min = ylims[0]
    y_max = ylims[1]

    ax1.set_ylabel(ylabel)
    ax2.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)
    ax2.set_xlabel(xlabel)
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    xlims = ax2.get_xlim()
    ax2.set_xlim([xlims[0], x2_max])

    # Définition des intervalles majeurs et mineurs pour la grille
    ax1.set_xticks(np.arange(0, x1_max, 4096))  # Intervalles majeurs tous les 64
    ax1.set_xticks(np.arange(0, x1_max, 512), minor=True)  # Intervalles mineurs tous les 8
    ax2.set_xticks(np.arange(0, x2_max, 64))  # Intervalles majeurs tous les 64
    ax2.set_xticks(np.arange(0, x2_max, 8), minor=True)  # Intervalles mineurs tous les 8

    # Activation de la grille majeure et mineure
    ax1.grid(which='major', linestyle='-', linewidth='0.8', color='black')  # Grille majeure
    ax1.grid(which='minor', linestyle='--', linewidth='0.5', color='gray')  # Grille mineure
    ax2.grid(which='major', linestyle='-', linewidth='0.8', color='black')  # Grille majeure
    ax2.grid(which='minor', linestyle='--', linewidth='0.5', color='gray')  # Grille mineure

    fig.tight_layout()

    plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
    print("results plotted in file " + plotFileName)

fdbkDelay(path, col)