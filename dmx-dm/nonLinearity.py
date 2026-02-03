# imports
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

import constants as cst
import general_tools as gt
import readData as rddt


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


def nonlinearity():
    """
    This functions characterises the non-linearity of DMX signals from scan data.
    This is applicable to feedback + error signals or offset + error signals.
    The type of scan is detected from the name of the fits files.
    The result is a plot saved in the /PLOTS directory of the session.
    The function searches for scans of the different columns.

    Args:
        tconf : test configuration

    Returns:
        nothing

    """

    # Index of the first pixel that shall be used
    pixel_id = 10
    # Number of pixels to be averaged
    nb_pix_avg = 16

    # Data directory
    dirpath = os.path.join("..", "..")

    session_name = os.path.basename(os.path.realpath(dirpath))

    # Looking for test configuration parameters
    index_bxl = session_name.find("BXL")
    bxl = int(session_name[index_bxl + 3:index_bxl + 3 + 2])  # boxcar length
    pathData = os.path.join(dirpath, cst.scanDirName)
    pathPlot = os.path.join(dirpath, cst.plotDirName)
    gt.createdir(pathPlot)

    ytit1 = "Error values (ADU)"
    ytit2 = "Non linearity (ADC LSB)"
    ytit3 = "Non Linearity (% of ADC FSR)"
    ymargin = 1024 # margin for ylimits in the plot
    dotsize = 1
    major_grid_ratio = 4
    grid_ratio = 4

    print("/----------------------------------------------------------")
    print("/ Processing non linearity from directory:")
    print("/ ", pathData)
    print("/----------------------------------------------------------\n")

    # looking for the position of the col id in the file names
    files = [f for f in os.listdir(pathData) \
             if os.path.isfile(os.path.join(pathData, f)) \
             and f[:5] == "scan_" \
             and f[-5:] == ".fits"]
    index_col = files[0].find("col")  # to find the column index in the file names

    
    for col in range(cst.nColPerDemux):

        # Searching fits files
        files = [f for f in os.listdir(pathData) \
                 if os.path.isfile(os.path.join(pathData, f)) \
                 and f[:5] == "scan_" \
                 and f[index_col:index_col + 4] == "col{0:}".format(col) \
                 and f[-5:] == ".fits"]

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
                xName, ctrl, file_scan, fileError = rddt.read_scan(os.path.join(pathData, file))
                fileError = np.mean(fileError[pixel_id:pixel_id+nb_pix_avg, :], axis=0)  # averaging the error value of few pixels
                error = np.append(error, fileError)
                scan = np.append(scan, file_scan)
                scan_type = np.append(scan_type, xName)

            # Checking the scan type (feedback or offset)
            if np.all(scan_type == "Feedback"):
                xtit = "Feedback values (ADU)"
                xlim = [-cst.fsrDACFdbkADU / 2, cst.fsrDACFdbkADU / 2]
                ylim1 = [-cst.fsrADCErrorADU / 2 -ymargin, cst.fsrADCErrorADU / 2 +ymargin]
                plotFileName = 'fdbkAndErrorINL2_col{0:}_bxl{1:}.png'.format(col, bxl)
                plotFileName2 = 'fdbkAndErrorDNL_col{0:}_bxl{1:}.png'.format(col, bxl)
                plotFileName3 = 'fdbkAndErrorINL1_col{0:}_bxl{1:}.png'.format(col, bxl)
                figsuptitle = 'Feedback + Error INL(2) measurement for column {0:}\n'.format(col) \
                              + '(' + session_name + ')'
                figsuptitle2 = 'Feedback + Error DNL measurement for column {0:}\n'.format(col) \
                               + '(' + session_name + ')'
                figsuptitle3 = 'Feedback + Error INL(1) measurement for column {0:} \n'.format(col) \
                               + '(' + session_name + ')'
            elif np.all(scan_type == "Offset"):
                # For the offset scans we use 4 frames per steps because of the settling time
                # We keep only the data of the last frame
                scan = scan[3::4]
                error = error[3::4]

                xtit = "Offset values (ADU)"
                xlim = [0, cst.fsrDACOfcoCoarseADU]
                ylim1 = [0-ymargin, cst.fsrADCErrorADU / 2 +ymargin]
                plotFileName = 'ofcoAndErrorINL1_col{0:}_bxl{1:}.png'.format(col, bxl)
                plotFileName2 = 'ofcoAndErrorDNL_col{0:}_bxl{1:}.png'.format(col, bxl)
                plotFileName3 = 'ofcoAndErrorINL2_col{0:}_bxl{1:}.png'.format(col, bxl)
                figsuptitle = 'Offset + Error INL-1 measurement for column {0:}\n'.format(col) \
                              + '(' + session_name + ')'
                figsuptitle2 = 'Offset + Error DNL measurement for column {0:}\n'.format(col) \
                               + '(' + session_name + ')'
                figsuptitle3 = 'Offset + Error INL-2 measurement for column {0:}\n'.format(col) \
                               + '(' + session_name + ')'
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
            deviation_lsb = error[i_ok] - fit

            # Linear fit of the data (reduced range)
            red_factor1 = 0.85
            i_red1 = np.where(np.abs(error) < 2**(cst.dmxNbBitsADCError-1)*red_factor1)[0]
            coeffs_red1 = np.polyfit(scan[i_red1], error[i_red1], 1)
            fit_red1 = coeffs_red1[1] + coeffs_red1[0] * scan[i_red1]
            deviation_lsb_red1 = error[i_red1] - fit_red1

            # Linear fit of the data (reduced range)
            red_factor2 = 0.7
            i_red2 = np.where(np.abs(error) < 2**(cst.dmxNbBitsADCError-1)*red_factor2)[0]
            coeffs_red2 = np.polyfit(scan[i_red2], error[i_red2], 1)
            fit_red2 = coeffs_red2[1] + coeffs_red2[0] * scan[i_red2]
            deviation_lsb_red2 = error[i_red2] - fit_red2

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
            lbl = 'Linear fit (Y = {0:.4} X'.format(coeffs[0]) + sign_str + '{0:.4})'.format(abs(coeffs[1])) + ' (fit is done on FS)'
            ax1.plot(scan[i_ok], fit, ':', color='red', linewidth=1, label=lbl)
            nl_threshold_for_reduce_ranges = 0.3
            if np.abs(deviation_lsb).max() > nl_threshold_for_reduce_ranges:
                lbl = 'Linear fit2 (Y = {0:.4} X'.format(coeffs_red1[0]) + sign_str + '{0:.4})'.format(abs(coeffs_red1[1])) + ' (fit2 is done on {0:}% of FS)'.format(int(red_factor1*100))
                ax1.plot(scan[i_red1], fit_red1, ':', color='orange', linewidth=1, label=lbl)
                lbl = 'Linear fit3 (Y = {0:.4} X'.format(coeffs_red2[0]) + sign_str + '{0:.4})'.format(abs(coeffs_red2[1])) + ' (fit3 is done on {0:}% of FS)'.format(int(red_factor2*100))
                ax1.plot(scan[i_red2], fit_red2, ':', color='purple', linewidth=1, label=lbl)
            ax1.set_xlim(xlim)
            ax1.set_ylim(ylim1)
            ax1.set_xlabel(xtit)
            ax1.set_ylabel(ytit1)
            ax1.legend(loc='upper left')
            set_grid(ax1, True, True, 2 ** 12, 2 ** 8, 2 ** 12, 2 ** 8)

            ax2 = fig.add_subplot(2, 1, 2) # non linearity

            val1 = max(np.abs(deviation_lsb))
            val2 = val1 * 100 / cst.fsrADCErrorADU
            lbl = 'Scan - Fit  (on FS the non linearity is {0:2.1f} LSB or {1:.2} %)'.format(val1, val2)
            ax2.scatter(scan[i_ok], deviation_lsb, s = dotsize, color='red', label=lbl)
            if np.abs(deviation_lsb).max() > nl_threshold_for_reduce_ranges:
                val0, val1 = int(red_factor1 * 100),  max(np.abs(deviation_lsb_red1))
                val2 = val1*100/cst.fsrADCErrorADU
                lbl = 'Scan - Fit2  (on {0:}% of FS the non linearity is {1:2.1f} LSB or {2:.2} %)'.format(val0, val1, val2)
                ax2.scatter(scan[i_red1], deviation_lsb_red1, s=dotsize, color='orange', label=lbl)
                val0, val1 = int(red_factor2 * 100),  max(np.abs(deviation_lsb_red2))
                val2 = val1*100/cst.fsrADCErrorADU
                lbl = 'Scan - Fit3  (on {0:}% of FS the non linearity is {1:2.1f} LSB or {2:.2} %)'.format(val0, val1, val2)
                ax2.scatter(scan[i_red2], deviation_lsb_red2, s=dotsize, color='purple', label=lbl)
            ax2.set_xlim(xlim)
            ax2.set_xlabel(xtit)
            ax2.set_ylabel(ytit2)
            ylimits = ax2.get_ylim()
            deltay = ylimits[1] - ylimits[0]
            ylim_extension = 0.75
            ax2.set_ylim([ylimits[0]-ylim_extension*deltay, ylimits[1]+ylim_extension*deltay])
            ax2.legend(loc='upper left')

            # second y axis for LSB units
            ax22 = ax2.twinx()
            ylims = ax2.get_ylim()
            ylims22 = [ylims[0] * 100 / cst.fsrADCErrorADU, ylims[1] * 100 / cst.fsrADCErrorADU]
            ax22.set_ylim(ylims22)
            ax22.set_ylabel(ytit3)

            ymajor = np.round((ylims[1] - ylims[0]) / major_grid_ratio)
            yminor = ymajor / grid_ratio
            set_grid(ax2, True, True, 2 ** 12, 2 ** 8, ymajor, yminor)

            fig.tight_layout()

            plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
            print("INL-2 results plotted in file " + plotFileName)


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
            ax.set_ylabel(ytit2)

            ax3 = ax.twinx()
            ylims = ax.get_ylim()
            ylims3 = [ylims[0] * 100 / cst.fsrADCErrorADU, ylims[1] * 100 / cst.fsrADCErrorADU]
            ax3.set_ylim(ylims3)
            ax3.set_ylabel(ytit3)

            ymajor = np.round((ylims[1] - ylims[0]) / major_grid_ratio)
            yminor = ymajor / grid_ratio
            set_grid(ax, True, True, 2 ** 12, 2 ** 8, ymajor, yminor)

            fig.tight_layout()

            plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
            print("DNL results plotted in file " + plotFileName2)


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
            ax.set_ylabel(ytit2)

            ax3 = ax.twinx()
            ylims = ax.get_ylim()
            ylims3 = [ylims[0] * 100 / cst.fsrADCErrorADU, ylims[1] * 100 / cst.fsrADCErrorADU]
            ax3.set_ylim(ylims3)
            # ax3.set_ylim(ylimDNL)
            ax3.set_ylabel(ytit3)

            ymajor = np.round((ylims[1] - ylims[0]) / major_grid_ratio)
            yminor = ymajor / grid_ratio
            set_grid(ax, True, True, 2 ** 12, 2 ** 8, ymajor, yminor)

            fig.tight_layout()

            plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
            print("INL-1 results plotted in file " + plotFileName3)

        print("/---------------")


# -------------------------------------------------------------------------------------
