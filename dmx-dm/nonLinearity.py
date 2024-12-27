# imports
import numpy as np
import matplotlib.pyplot as plt
import os
import general_tools as gt
import readData as rddt
import constants as cst
from matplotlib.ticker import MultipleLocator

path = "/Users/laurent/Data/TestPlan21-perfo/errorLinearity_20241213-112110/"


def set_grid(ax, major_on, minor_on, xmajor, xminor, ymajor, yminor):
    """
    This function plots minor and major grids on a sub_plot

    Parameters
    ----------
    ax: sub_plot id
    major_on: boolean, if True the major grid is plotted
    minor_on: boolean, if True the minor grid is plotted
    xmajor: integer, xspacing of major grid
    xminor: integer, xspacing of minor grid
    ymajor: integer, yspacing of major grid
    yminor: integer, yspacing of minor grid

    Returns
    -------
    Nothing
    """

    # Définir le pas de la grille majeure
    ax.xaxis.set_major_locator(MultipleLocator(xmajor))
    ax.yaxis.set_major_locator(MultipleLocator(ymajor))

    # Définir le pas de la grille mineure
    ax.xaxis.set_minor_locator(MultipleLocator(xminor))
    ax.yaxis.set_minor_locator(MultipleLocator(yminor))

    # Activation de la grille majeure
    ax.grid(major_on, which='major', linestyle='-', linewidth='0.8', color='black')

    # Activer la grille mineure
    if minor_on:
        ax.minorticks_on()  # Activer les ticks mineurs
        ax.grid(True, which='minor', linestyle=':', linewidth='0.5', color='gray')


def plot_non_linearity(dir_path, pixel_id):
    """
    This functions characterises the non-linearity of DMX signals from scan data.
    This is applicable to feedback + error signals or offset + error signals.
    The type of scan is detected from the name of the fits files.
    The result is a plot saved in the /PLOTS directory of the session.
    The function searches for scans of the different columns.

    Args:
        dir_path: string, the path of the data directory
        pixel_id: integer, the id of the pixel to be considered

    Returns:
        nothing

    """

    # Creation of a directory for the plot files
    pathPlot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(pathPlot)

    ytit1 = "Error values (ADU)"
    ytit2 = "Non Linearity (% FSR)"
    ylim1 = [-cst.fsrADCErrorADU / 2, cst.fsrADCErrorADU / 2 - 1]

    print("/----------------------------------------------------------")
    print("/ Processing feedback + error non linearity from directory:")
    print("/ ", dir_path)
    print("/----------------------------------------------------------\n")

    for col in range(cst.nColPerDemux):

        # Searching fits files
        files = [f for f in os.listdir(path) \
                 if os.path.isfile(os.path.join(path, f)) \
                 and f[-10:] == "_col{0:}.fits".format(col) and f[:5] == "scan_"]

        # Checking file length
        if len(files) == 0:
            print("No fits file found for column {0:}".format(col))
        else:
            if len(files) == 1:
                print("Found 1 fits file for column {0:}".format(col))
            else:
                print("Found {0:} fits files for column {1:}".format(len(files), col))

            # Detecting the scan type (feedback or offset)
            scan_type = np.array([file[21:27] for file in files])
            if np.all(scan_type == "Feedba"):
                xtit = "Feedback values (ADU)"
                xlim = [-cst.fsrDACFdbkADU / 2, cst.fsrDACFdbkADU / 2 - 1]
                plotFileName = 'fdbkAndErrorLinearity_col{0:}.png'.format(col)
                figsuptitle = 'Feedback + Error linearity measurement for column {0:}'.format(col)
            elif np.all(scan_type == "Offset"):
                xtit = "Offset values (ADU)"
                xlim = [-cst.fsrDACOfcoCoarseADU / 2, cst.fsrDACOfcoCoarseADU / 2 - 1]
                plotFileName = 'ofcoAndErrorLinearity_col{0:}.png'.format(col)
                figsuptitle = 'Offset + Error linearity measurement for column {0:}'.format(col)
            else:
                raise ValueError("Error, found different scan types!")

            # Reading and concatenating the data of the different fits files
            error=np.array([])  # Empty array
            scan=np.array([])  # Empty array
            for file in files:
                print("Reading data from file ", file)
                xName, ctrl, file_scan, fileError = rddt.readScan(os.path.join(dir_path, file))
                fileError = fileError[pixel_id, :]  # keeping the error value of a single pixel
                error = np.append(error, fileError)
                scan = np.append(scan, file_scan)

            # Linear fit of the data
            coeffs = np.polyfit(scan, error, 1)
            fit = coeffs[1] + coeffs[0] * scan
            deviationPercentFSR = 100 * (error - fit) / cst.fsrADCErrorADU

            # Doing the plots
            fig = plt.figure(figsize=(12, 8))
            plotFullFileName = os.path.join(pathPlot, plotFileName)
            fig.suptitle(figsuptitle, fontsize=16)

            ax1 = fig.add_subplot(2, 1, 1) # output vs input

            ax1.plot(scan, error, linewidth = 2, label='Scan')
            ax1.plot(scan, fit, ':r', linewidth = 1, label='Fit')
            ax1.set_xlim(xlim)
            ax1.set_ylim(ylim1)
            ax1.set_xlabel(xtit)
            ax1.set_ylabel(ytit1)
            ax1.legend(loc='upper left')
            set_grid(ax1, True, True, 2 ** 12, 2 ** 8, 2 ** 12, 2 ** 8)

            ax2 = fig.add_subplot(2, 1, 2) # non linearity

            ax2.scatter(scan, deviationPercentFSR, label='100 * (Scan - Fit) / FSR')
            ax2.set_xlim(xlim)
            ax2.set_xlabel(xtit)
            ax2.set_ylabel(ytit2)
            ylimits = ax2.get_ylim()
            deltay = ylimits[1] - ylimits[0]
            ax2.set_ylim([ylimits[0]-deltay/2, ylimits[1]+deltay/2])
            ax2.legend(loc='upper left')

            ymajor = 10 ** np.floor(np.log10(deltay))
            yminor = ymajor / 5
            set_grid(ax2, True, True, 2 ** 12, 2 ** 8, ymajor, yminor)

            fig.tight_layout()

            plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
            print("results plotted in file " + plotFileName)

        print("/---------------")

plot_non_linearity(path, 0)