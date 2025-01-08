# imports
import numpy as np
import matplotlib.pyplot as plt
import os
import general_tools as gt
import readData as rddt
import constants as cst
from matplotlib.ticker import MultipleLocator


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


def plot_non_linearity(dir_path):
    """
    This functions characterises the non-linearity of DMX signals from scan data.
    This is applicable to feedback + error signals or offset + error signals.
    The type of scan is detected from the name of the fits files.
    The result is a plot saved in the /PLOTS directory of the session.
    The function searches for scans of the different columns.

    Args:
        dir_path: string, the path of the data directory

    Returns:
        nothing

    """

    # Pixel dont les données doivent être utilisées
    # Utiliser un id >> 0 pour ne pas être impacté par le changement rapide de niveau
    # entre le dernier et le premier pixel surtout en fin / début de scan
    pixel_id = 30

    # Data directory
    pathData = os.path.join(dir_path, cst.scanDirName)

    # Session name
    session_name = os.path.basename(dir_path)

    # Creation of a directory for the plot files
    pathPlot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(pathPlot)

    ytit1 = "Error values (ADU)"
    ytit2 = "Non Linearity (% FSR)"
    dotsize = 1

    print("/----------------------------------------------------------")
    print("/ Processing non linearity from directory:")
    print("/ ", pathData)
    print("/----------------------------------------------------------\n")

    for col in range(cst.nColPerDemux):

        # Searching fits files
        files = [f for f in os.listdir(pathData) \
                 if os.path.isfile(os.path.join(pathData, f)) \
                 and f[-10:] == "_col{0:}.fits".format(col) and f[:5] == "scan_"]

        # Checking number of files
        if len(files) == 0:
            print("No fits file found for column {0:}".format(col))
        else:
            if len(files) == 1:
                print("Found 1 fits file for column {0:}".format(col))
            else:
                print("Found {0:} fits files for column {1:}".format(len(files), col))

            # Reading and concatenating the data of the different fits files
            error=np.array([])  # Empty array
            scan=np.array([])  # Empty array
            scan_type=np.array([])  # Empty array
            for file in files:
                print("Reading data from file ", file)
                xName, ctrl, file_scan, fileError = rddt.readScan(os.path.join(pathData, file))
                fileError = fileError[pixel_id, :]  # keeping the error value of a single pixel
                error = np.append(error, fileError)
                scan = np.append(scan, file_scan)
                scan_type = np.append(scan_type, xName)

            # Checking the scan type (feedback or offset)
            if np.all(scan_type == "Feedback"):
                xtit = "Feedback values (ADU)"
                xlim = [-cst.fsrDACFdbkADU / 2, cst.fsrDACFdbkADU / 2]
                ylim1 = [-cst.fsrADCErrorADU / 2, cst.fsrADCErrorADU / 2]
                plotFileName = 'fdbkAndErrorLinearity_col{0:}.png'.format(col)
                figsuptitle = ('Feedback + Error linearity measurement for column {0:}   ('.format(col)
                               + session_name + ')')
            elif np.all(scan_type == "Offset"):
                # For the offset scans we use 4 frames per steps because of the settling time
                # We keep only the last frame
                scan = scan[3::4]
                error = error[3::4]

                xtit = "Offset values (ADU)"
                xlim = [0, cst.fsrDACOfcoCoarseADU]
                ylim1 = [0, cst.fsrADCErrorADU / 2]
                plotFileName = 'ofcoAndErrorLinearity_col{0:}.png'.format(col)
                figsuptitle = ('Offset + Error linearity measurement for column {0:}   ('.format(col)
                               + session_name + ')')
            else:
                raise ValueError("Error, found different scan types!")


            # Checking if all DAC values are used
            expected_array = np.arange(scan.min(), scan.max()+1)
            #expected_array = np.arange(xlim[0], xlim[1]+1)
            scan_copy = np.unique(scan.copy())
            scan_copy.sort()
            if not (expected_array == scan_copy).all():
                print("   Error, values are missing in the scan!")
                print(expected_array, scan_copy)

            # Linear fit of the data
            coeffs = np.polyfit(scan, error, 1)
            fit = coeffs[1] + coeffs[0] * scan
            deviationPercentFSR = 100 * (error - fit) / cst.fsrADCErrorADU

            # Doing the plots
            fig = plt.figure(figsize=(12, 8))
            plotFullFileName = os.path.join(pathPlot, plotFileName)
            fig.suptitle(figsuptitle, fontsize=14)

            ax1 = fig.add_subplot(2, 1, 1) # output vs input

            ax1.scatter(scan, error, s = dotsize, label='Scan')
            if coeffs[1]>0:
                sign_str = ' + '
            else:
                sign_str = ' - '
            lbl = 'Linear fit (Y = {0:.4} X'.format(coeffs[0]) + sign_str + '{0:.4})'.format(abs(coeffs[1]))
            ax1.plot(scan, fit, ':r', linewidth = 1, label=lbl)
            ax1.set_xlim(xlim)
            ax1.set_ylim(ylim1)
            ax1.set_xlabel(xtit)
            ax1.set_ylabel(ytit1)
            ax1.legend(loc='upper left')
            set_grid(ax1, True, True, 2 ** 12, 2 ** 8, 2 ** 12, 2 ** 8)

            ax2 = fig.add_subplot(2, 1, 2) # non linearity

            ax2.scatter(scan, deviationPercentFSR, s = dotsize, label='100 * (Scan - Fit) / FSR')
            ax2.set_xlim(xlim)
            ax2.set_xlabel(xtit)
            ax2.set_ylabel(ytit2)
            ylimits = ax2.get_ylim()
            deltay = ylimits[1] - ylimits[0]
            ylim_extension = 0.75
            ax2.set_ylim([ylimits[0]-ylim_extension*deltay, ylimits[1]+ylim_extension*deltay])
            ax2.legend(loc='upper left')

            #ymajor = 10 ** np.floor(np.log10(ylim_extension*deltay))*4
            grid_ratio = 10
            ymajor = np.round((ylim_extension*deltay)*grid_ratio)/grid_ratio
            print(ymajor)
            yminor = ymajor/5
            set_grid(ax2, True, True, 2 ** 12, 2 ** 8, ymajor, yminor)

            fig.tight_layout()

            plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
            print("results plotted in file " + plotFileName)

        print("/---------------")

pathes = ["/Users/laurent/Data/TestPlan21-perfo/20250106_fdbkAndErrorLinearity",
        "/Users/laurent/Data/TestPlan21-perfo/20250108_150256_fdbkAndErrorLinearity",
        "/Users/laurent/Data/TestPlan21-perfo/20250108_152522_fdbkAndErrorLinearity",
        "/Users/laurent/Data/TestPlan21-perfo/20250108_163003_ofcoAndErrorLinearity"
        ]

#for path in pathes[-1:]:
for path in pathes:
    plot_non_linearity(path)