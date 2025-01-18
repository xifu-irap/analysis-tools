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
    pixel_id = 20

    # Data directory
    pathData = os.path.join(dir_path, cst.scanDirName)

    # Session name
    session_name = os.path.basename(dir_path)

    # Creation of a directory for the plot files
    pathPlot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(pathPlot)

    ytit1 = "Error values (ADU)"
    ytit2 = "Non Linearity (% FSR)"
    ytit3 = "LSB"
    ymargin = 1024 # margin for ylimits in the plot
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
                ylim1 = [-cst.fsrADCErrorADU / 2 -ymargin, cst.fsrADCErrorADU / 2 +ymargin]
                plotFileName = 'fdbkAndErrorLinearity_col{0:}.png'.format(col)
                plotFileName2 = 'fdbkAndErrorDNL_col{0:}.png'.format(col)
                plotFileName3 = 'fdbkAndErrorINL_col{0:}.png'.format(col)
                figsuptitle = ('Feedback + Error linearity measurement for column {0:}   ('.format(col)
                               + session_name + ')')
                figsuptitle2 = ('Feedback + Error DNL measurement for column {0:}   ('.format(col)
                               + session_name + ')')
                figsuptitle3 = ('Feedback + Error INL measurement for column {0:}   ('.format(col)
                                + session_name + ')')
            elif np.all(scan_type == "Offset"):
                # For the offset scans we use 4 frames per steps because of the settling time
                # We keep only the last frame
                scan = scan[3::4]
                error = error[3::4]

                xtit = "Offset values (ADU)"
                xlim = [0, cst.fsrDACOfcoCoarseADU]
                ylim1 = [0-ymargin, cst.fsrADCErrorADU / 2 +ymargin]
                plotFileName = 'ofcoAndErrorLinearity_col{0:}.png'.format(col)
                plotFileName2 = 'ofcoAndErrorDNL_col{0:}.png'.format(col)
                plotFileName3 = 'ofcoAndErrorINL_col{0:}.png'.format(col)
                figsuptitle = ('Offset + Error linearity measurement for column {0:}   ('.format(col)
                               + session_name + ')')
                figsuptitle2 = ('Offset + Error DNL measurement for column {0:}   ('.format(col)
                               + session_name + ')')
                figsuptitle3 = ('Offset + Error INL measurement for column {0:}   ('.format(col)
                                + session_name + ')')
            else:
                raise ValueError("Error, found different scan types!")

            # Sorting the data wrt DAC values
            unique_values, unique_i = np.unique(scan, return_index=True)
            scan_unique = scan[unique_i]
            error_unique = error[unique_i]
            sorted_i = np.argsort(scan_unique)
            scan = scan_unique[sorted_i]
            error = error_unique[sorted_i]

            # Checking if all DAC values are used
            expected_array = np.arange(scan.min(), scan.max()+1)
            if not (expected_array == scan).all():
                print("   Error, values are missing in the scan!")
                print(expected_array, scan)

            # Ignoring saturations
            i_ok = np.where(np.abs(error) < 2**(cst.dmxNbBitsADCError-1)-1)[0]

            # Linear fit of the data (full range)
            coeffs = np.polyfit(scan[i_ok], error[i_ok], 1)
            fit = coeffs[1] + coeffs[0] * scan[i_ok]
            deviationPercentFSR = 100 * (error[i_ok] - fit) / cst.fsrADCErrorADU

            # Linear fit of the data (reduced range)
            red_factor1 = 0.85
            i_red1 = np.where(np.abs(error) < 2**(cst.dmxNbBitsADCError-1)*red_factor1)[0]
            coeffs_red1 = np.polyfit(scan[i_red1], error[i_red1], 1)
            fit_red1 = coeffs_red1[1] + coeffs_red1[0] * scan[i_red1]
            deviationPercentFSR_red1 = 100 * (error[i_red1] - fit_red1) / cst.fsrADCErrorADU

            # Linear fit of the data (reduced range)
            red_factor2 = 0.7
            i_red2 = np.where(np.abs(error) < 2**(cst.dmxNbBitsADCError-1)*red_factor2)[0]
            coeffs_red2 = np.polyfit(scan[i_red2], error[i_red2], 1)
            fit_red2 = coeffs_red2[1] + coeffs_red2[0] * scan[i_red2]
            deviationPercentFSR_red2 = 100 * (error[i_red2] - fit_red2) / cst.fsrADCErrorADU

            # Computing DNL and INL
            lsb_ideal = (error[i_ok].max() - error[i_ok].min()) / len(error[i_ok])
            print(lsb_ideal)
            dnl = (error[i_ok][1:] - error[i_ok][:-1])/lsb_ideal - 1
            inl = (error[i_ok] - error[i_ok][0])/lsb_ideal - np.arange(len(error[i_ok]))


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
            ax1.plot(scan[i_ok], fit, ':', color='red', linewidth=1, label=lbl)
            nl_threshold_for_reduce_ranges = 0.3
            if np.abs(deviationPercentFSR).max() > nl_threshold_for_reduce_ranges:
                lbl = 'Linear fit2 (Y = {0:.4} X'.format(coeffs_red1[0]) + sign_str + '{0:.4})'.format(abs(coeffs_red1[1])) + ' (fit2 is limited to {0:}% of the full scale)'.format(int(red_factor1*100))
                ax1.plot(scan[i_red1], fit_red1, ':', color='orange', linewidth=1, label=lbl)
                lbl = 'Linear fit3 (Y = {0:.4} X'.format(coeffs_red2[0]) + sign_str + '{0:.4})'.format(abs(coeffs_red2[1])) + ' (fit3 is limited to {0:}% of the full scale)'.format(int(red_factor2*100))
                ax1.plot(scan[i_red2], fit_red2, ':', color='purple', linewidth=1, label=lbl)
            ax1.set_xlim(xlim)
            ax1.set_ylim(ylim1)
            ax1.set_xlabel(xtit)
            ax1.set_ylabel(ytit1)
            ax1.legend(loc='upper left')
            set_grid(ax1, True, True, 2 ** 12, 2 ** 8, 2 ** 12, 2 ** 8)

            ax2 = fig.add_subplot(2, 1, 2) # non linearity

            ax2.scatter(scan[i_ok], deviationPercentFSR, s = dotsize, color='red', label='100 * (Scan - Fit) / FSR')
            if np.abs(deviationPercentFSR).max() > nl_threshold_for_reduce_ranges:
                lbl = '100 * (Scan - Fit2) / FSR  (limited to {0:}% of the full scale)'.format(int(red_factor1*100))
                ax2.scatter(scan[i_red1], deviationPercentFSR_red1, s=dotsize, color='orange', label=lbl)
                lbl = '100 * (Scan - Fit3) / FSR  (limited to {0:}% of the full scale)'.format(int(red_factor2 * 100))
                ax2.scatter(scan[i_red2], deviationPercentFSR_red2, s=dotsize, color='purple', label=lbl)
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
            yminor = ymajor/5
            set_grid(ax2, True, True, 2 ** 12, 2 ** 8, ymajor, yminor)

            fig.tight_layout()

            plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
            print("results plotted in file " + plotFileName)


            # Doing the DNL plots
            fig = plt.figure(figsize=(12, 7))
            plotFullFileName = os.path.join(pathPlot, plotFileName2)
            fig.suptitle(figsuptitle2, fontsize=14)

            ax = fig.add_subplot(1, 1, 1)  # output vs input
            ax.scatter(scan[i_ok][:-1], dnl, s = dotsize, color='red', label='DNL')
            ax.plot(xlim, [0,0], '-k', linewidth=0.5)
            ax.set_xlim(xlim)
            yl = max(np.abs(dnl).max()*1.2, 5)

            ax.set_ylim([-yl, yl])
            ax.set_xlabel(xtit)
            ax.set_ylabel(ytit3)

            set_grid(ax, True, True, 2 ** 12, 2 ** 8, 5, 1)

            fig.tight_layout()

            plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
            print("DNL results plotted in file " + plotFileName)


            # Doing the INL plots
            fig = plt.figure(figsize=(12, 7))
            plotFullFileName = os.path.join(pathPlot, plotFileName3)
            fig.suptitle(figsuptitle3, fontsize=14)

            ax = fig.add_subplot(1, 1, 1)  # output vs input
            ax.scatter(scan[i_ok], inl, s = dotsize, color='red', label='INL')
            ax.plot(xlim, [0,0], '-k', linewidth=0.5)
            ax.set_xlim(xlim)
            yl = max(np.abs(inl).max()*1.2, 10)

            ax.set_ylim([-yl, yl])
            ax.set_xlabel(xtit)
            ax.set_ylabel(ytit3)

            set_grid(ax, True, True, 2 ** 12, 2 ** 8, 10, 2)

            fig.tight_layout()

            plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
            print("INL results plotted in file " + plotFileName)


        print("/---------------")



pathes = [
        "/Users/laurent/Data/TestPlan21-perfo/20250108_150256_fdbkAndErrorLinearity",
        "/Users/laurent/Data/TestPlan21-perfo/20250108_152522_fdbkAndErrorLinearity",
        "/Users/laurent/Data/TestPlan21-perfo/20250108_163003_ofcoAndErrorLinearity",
        "/Users/laurent/Data/TestPlan21-perfo/20250109_181148_fdbkAndErrorLinearity",
        "/Users/laurent/Data/TestPlan21-perfo/20250109_181826_fdbkAndErrorLinearity",
        "/Users/laurent/Data/TestPlan21-perfo/20250110_112602_fdbkAndErrorLinearity",
        "/Users/laurent/Data/TestPlan21-perfo/20250110_113807_fdbkAndErrorLinearity",
        "/Users/laurent/Data/TestPlan21-perfo/20250113_110936_ofcoAndErrorLinearity"
        ]

for path in pathes[-2:]:
    plot_non_linearity(path)