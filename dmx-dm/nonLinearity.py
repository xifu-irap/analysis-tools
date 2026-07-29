# ---------------------------------------------------------------------------------
# !/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright (C) 2021-2030 Laurent Ravera, IRAP Toulouse.
#  This file is part of the ATHENA X-IFU DRE data analysis tools software.
#
#  analysis-tools is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# ---------------------------------------------------------------------------------
#
#  laurent.ravera@irap.omp.eu
#  nonLinearity.py
#
# ---------------------------------------------------------------------------------

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


def nonlinearity(verbose=False):
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

    colors = ['steelblue', 'turquoise', 'blue', 'teal']

    # Data directory
    dirpath = os.path.join("..", "..")

    session_name = os.path.basename(os.path.realpath(dirpath))

    # Looking for test configuration parameters
    index_bxl = session_name.find("BXL") + len("BXL")
    bxl = int(session_name[index_bxl])  # boxcar length
    pathData = os.path.join(dirpath, cst.scanDirName)
    pathHk = os.path.join(dirpath, cst.hkDirName)
    pathPlot = os.path.join(dirpath, cst.plotDirName)
    gt.createdir(pathPlot)

    ytit = "Error values (ADU)"
    ytit_ADU = "Non linearity (ADC LSB)"
    ytit_pc = "Non Linearity (% of ADC FSR)"
    ymargin = 0  # margin for ylimits in the plot
    dotsize = 1
    major_grid_ratio = 4
    grid_ratio = 4

    ylim_INL = 150
    ylim_DNL = 20

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(pathHk)

    print(pathData)
    files = [f for f in os.listdir(pathData) \
             if os.path.isfile(os.path.join(pathData, f)) \
             and f[:5] == "scan_" \
             and f[-3:] == ".h5"]

    # Checking number of files
    if len(files) == 0:
        raise ValueError("Error, no scan files found in this session !")

    # Checking the scan type (feedback or offset)
    scan_type_ini = rddt.read_scan_type(os.path.join(pathData, files[0]))

    if scan_type_ini == "Feedback":
        scan_type_text = "FEEDBACK"

    elif scan_type_ini == "Offset":
        scan_type_text = "OFFSET"

    print("/----------------------------------------------------------")
    print("/ Non linearity test: ERROR + " + scan_type_text)
    print("/ Test session name: " + session_name)
    print("/----------------------------------------------------------")
    print("/ DEMUX model:       " + dmxModel + " {0:}".format(boardId))
    print("/ Firmware version:  {0:}".format(fwVersion))
    print("/ Box car length:    {0:} samples".format(bxl))
    print("/----------------------------------------------------------\n")

    for col in range(cst.nColPerDemux):

        # Searching scan files
        files = [f for f in os.listdir(pathData) \
                 if os.path.isfile(os.path.join(pathData, f)) \
                 and f[:5] == "scan_" \
                 and f[-4] == "{0:}".format(col) \
                 and f[-3:] == ".h5"]

        if len(files) == 0:
            raise ValueError("Error, no scan files found for column {0:}!".format(col))

        print("Found {0:} HDF5 file(s) for column {1:}".format(len(files), col))

        # Reading and concatenating the data of the different scan files
        error = np.array([])  # Empty array
        scan = np.array([])  # Empty array
        scan_type = np.array([])  # Empty array

        for file in files:
            if verbose:
                print("Reading data from file ", file)
            scan_type, ctrl, fileScan, fileError = rddt.read_scan(os.path.join(pathData, file))
            fileError = np.mean(fileError[:, :], axis=0)  # averaging the error in a frame

            # keeping data of the last frame only in a scan step
            if file == files[0]:  # initialisation of the array
                error = [fileError[0]]
                scan = [fileScan[0]]
            for i in range(len(fileScan)):
                if i == 0:
                    error = np.append(error, fileError[i])
                    scan = np.append(scan, fileScan[i])
                else:
                    if (fileScan[i] != fileScan[i - 1]):
                        error = np.append(error, fileError[i])
                        scan = np.append(scan, fileScan[i])
                    else:
                        error[-1] = fileError[i]
                        scan[-1] = fileScan[i]

        if scan_type == "Feedback":
            xtit = "Feedback values (ADU)"
            xlim = [-cst.fsrDACFdbkADU / 2, cst.fsrDACFdbkADU / 2]
            ylim = [-cst.fsrADCErrorADU / 2 - ymargin, cst.fsrADCErrorADU / 2 + ymargin]
            data_plotFileName = 'fdbkAndError_col{0:}_bxl{1:}.png'.format(col, bxl)
            dnl_plotFileName = 'fdbkAndErrorDNL_col{0:}_bxl{1:}.png'.format(col, bxl)
            inl_plotFileName = 'fdbkAndErrorINL_col{0:}_bxl{1:}.png'.format(col, bxl)
            figsuptitle0 = 'Feedback + Error non linearity measurement for column {0:}\n'.format(col) \
                           + '(' + dmxModel + " {0:}".format(boardId) + ", Firmware version: {0:}, ".format(
                fwVersion) + session_name + ')'
            figsuptitleDNL = 'Feedback + Error DNL measurement for column {0:}\n'.format(col) \
                             + '(' + dmxModel + " {0:}".format(boardId) + ", Firmware version: {0:}, ".format(
                fwVersion) + session_name + ')'
            figsuptitleINL = 'Feedback + Error INL measurement for column {0:} \n'.format(col) \
                             + '(' + dmxModel + " {0:}".format(boardId) + ", Firmware version: {0:}, ".format(
                fwVersion) + session_name + ')'

        elif scan_type == "Offset":
            xtit = "Offset values (ADU)"
            xlim = [0, cst.fsrDACOfcoCoarseADU]
            ylim = [0 - ymargin, cst.fsrADCErrorADU / 2 + ymargin]
            data_plotFileName = 'ofcoAndError_col{0:}_bxl{1:}.png'.format(col, bxl)
            dnl_plotFileName = 'ofcoAndErrorDNL_col{0:}_bxl{1:}.png'.format(col, bxl)
            inl_plotFileName = 'ofcoAndErrorINL_col{0:}_bxl{1:}.png'.format(col, bxl)
            figsuptitle0 = 'Offset + Error non linearity measurement for column {0:}\n'.format(col) \
                           + '(' + dmxModel + " {0:}".format(boardId) + ", Firmware version: {0:}, ".format(
                fwVersion) + session_name + ')'
            figsuptitleDNL = 'Offset + Error DNL measurement for column {0:}\n'.format(col) \
                             + '(' + dmxModel + " {0:}".format(boardId) + ", Firmware version: {0:}, ".format(
                fwVersion) + session_name + ')'
            figsuptitleINL = 'Offset + Error INL measurement for column {0:}\n'.format(col) \
                             + '(' + dmxModel + " {0:}".format(boardId) + ", Firmware version: {0:}, ".format(
                fwVersion) + session_name + ')'

        # Sorting the data wrt DAC values
        unique_values, unique_i = np.unique(scan, return_index=True)
        scan_unique = scan[unique_i]
        error_unique = error[unique_i]
        sorted_i = np.argsort(scan_unique)
        scan = scan_unique[sorted_i]
        error = error_unique[sorted_i]

        # Checking if all DAC values are used
        expected_array = np.arange(scan.min(), scan.max() + 1)
        print(expected_array[:32], len(expected_array))
        print(scan[:32], len(scan))
        if not (expected_array == scan).all():
            print("   Error, values are missing in the scan!")

        # Ignoring saturations
        i_ok = np.where(np.abs(error) < 2 ** (cst.dmxNbBitsADCError - 1) - 1)[0]

        # Linear fit of the data (full range)
        coeffs = np.polyfit(scan[i_ok], error[i_ok], 1)
        fit = coeffs[1] + coeffs[0] * scan[i_ok]
        deviation_lsb = error[i_ok] - fit

        # Computing DNL and INL
        lsb_ideal = (error[i_ok].max() - error[i_ok].min()) / len(error[i_ok])
        dnl = (error[i_ok][1:] - error[i_ok][:-1]) / lsb_ideal - 1
        inl = (error[i_ok] - error[i_ok][0]) / lsb_ideal - np.arange(len(error[i_ok]))

        # Doing the plots
        ## Non linearity tests data (output versus input)
        fig0 = plt.figure(figsize=(9, 6))
        plotFullFileName = os.path.join(pathPlot, data_plotFileName)
        fig0.suptitle(figsuptitle0, fontsize=14)
        ax0 = fig0.add_subplot(1, 1, 1)  # output vs input

        ax0.scatter(scan, error, s=dotsize, c='k', label='Scan')
        if coeffs[1] > 0:
            sign_str = ' + '
        else:
            sign_str = ' - '
        lbl = 'Linear fit (Y = {0:.4} X'.format(coeffs[0]) + sign_str + '{0:.4})'.format(abs(coeffs[1]))
        ax0.plot(scan[i_ok], fit, ':', color=colors[0], linewidth=1, label=lbl)
        ax0.set_xlim(xlim)
        ax0.set_ylim(ylim)
        ax0.set_xlabel(xtit)
        ax0.set_ylabel(ytit)
        ax0.legend(loc='upper left')
        set_grid(ax0, True, True, 2 ** 12, 2 ** 8, 2 ** 12, 2 ** 8)

        fig0.tight_layout()
        plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
        print("Linearity data plotted in file " + data_plotFileName)

        ## INL plot
        fig1 = plt.figure(figsize=(9, 6))
        plotFullFileName = os.path.join(pathPlot, inl_plotFileName)
        fig1.suptitle(figsuptitleINL, fontsize=14)
        ax1 = fig1.add_subplot(1, 1, 1)  # non linearity

        val1 = max(np.abs(deviation_lsb))
        val2 = val1 * 100 / cst.fsrADCErrorADU
        lbl = 'Scan - linear Fit  (the INL is {0:2.1f} LSB or {1:.2} %)'.format(val1, val2)
        ax1.scatter(scan[i_ok], deviation_lsb, s=dotsize, color=colors[0], label=lbl)
        ax1.set_xlim(xlim)
        ax1.set_xlabel(xtit)
        ax1.set_ylabel(ytit_ADU)

        ax1.set_ylim([-1 * ylim_INL, ylim_INL])
        ax1.legend(loc='upper left')

        # second y axis for LSB units
        ax11 = ax1.twinx()
        ylims = ax1.get_ylim()
        ylims11 = [ylims[0] * 100 / (0.5 * cst.fsrADCErrorADU), ylims[1] * 100 / (0.5 * cst.fsrADCErrorADU)]
        ax11.set_ylim(ylims11)
        ax11.set_ylabel(ytit_pc)

        ymajor = np.round((ylims[1] - ylims[0]) / major_grid_ratio)
        yminor = ymajor / grid_ratio
        set_grid(ax1, True, True, 2 ** 12, 2 ** 8, ymajor, yminor)

        fig1.tight_layout()
        plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
        print("INL results plotted in file " + inl_plotFileName)

        ## DNL plots
        fig2 = plt.figure(figsize=(9, 6))
        plotFullFileName = os.path.join(pathPlot, dnl_plotFileName)
        fig2.suptitle(figsuptitleDNL, fontsize=14)
        ax2 = fig2.add_subplot(1, 1, 1)

        val1 = max(np.abs(dnl))
        val2 = val1 * 100 / cst.fsrADCErrorADU
        lbl = 'Scan - linear Fit  (the DNL is {0:2.1f} LSB or {1:.2} %)'.format(val1, val2)
        ax2.scatter(scan[i_ok][:-1], dnl, s=dotsize, color=colors[0], label='DNL')
        ax2.plot(xlim, [0, 0], '-k', linewidth=0.5)
        ax2.set_xlim(xlim)

        ax2.set_ylim([-1 * ylim_DNL, ylim_DNL])
        ax2.set_xlabel(xtit)
        ax2.set_ylabel(ytit_ADU)

        ax22 = ax2.twinx()
        ylims = ax2.get_ylim()
        ylims22 = [ylims[0] * 100 / (0.5 * cst.fsrADCErrorADU), ylims[1] * 100 / (0.5 * cst.fsrADCErrorADU)]
        ax22.set_ylim(ylims22)
        ax22.set_ylabel(ytit_pc)

        ymajor = np.round((ylims[1] - ylims[0]) / major_grid_ratio)
        yminor = ymajor / grid_ratio
        set_grid(ax2, True, True, 2 ** 12, 2 ** 8, ymajor, yminor)

        fig2.tight_layout()
        plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
        print("DNL results plotted in file " + dnl_plotFileName)

        print("/---------------")


# -------------------------------------------------------------------------------------
