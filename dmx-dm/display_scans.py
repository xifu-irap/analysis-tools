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
#  display_scans.py
#
# ---------------------------------------------------------------------------------

# imports
import os

import matplotlib.pyplot as plt

import constants as cst
import general_tools as gt
import readData as rddt


def plot_scan(scanFileName, plotFileName, suptitle, title, xlabel, verbose=True):
    pix = cst.nPixPerCol - 1  # considered pixel
    ylabel = 'Demux error input (ADU)'

    # Reading the scan data
    if verbose:
        print('Reading data from scan file ', scanFileName)
    xName, _, x, pixels_data = rddt.read_scan(scanFileName)

    # Searching the length of each step
    i = 0
    while (x[i] == x[0]):
        i += 1
    # keeping only the data of the last frame
    x = x[i - 1::i]
    pixels_data = pixels_data[pix, i - 1::i]

    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(1, 1, 1)  # global plot
    ax1.plot(x, pixels_data)

    plt.suptitle(suptitle)
    ax1.set_title(title, fontsize=8)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)

    fig.tight_layout()

    plt.savefig(plotFileName, dpi=300, bbox_inches='tight')

    if verbose:
        print("results plotted in file " + plotFileName)


def scan_amp_squid(verbose):
    # Paths definition
    dir_path = os.path.join("..", "..")
    hk_path = os.path.join(dir_path, cst.hkDirName)
    data_path = os.path.join(dir_path, cst.scanDirName)
    plot_path = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(plot_path)
    session_name = os.path.basename(os.path.realpath(dir_path))
    col = session_name[-1]

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(hk_path)

    if verbose:
        print("/----------------------------------------------------------")
        print("/ AMP SQUID V(Phi) measurement")
        print("/ Test session name:   " + session_name)
        print("/ Tested column:       " + col)
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:     {0:}".format(fwVersion))
        print("/----------------------------------------------------------\n")

    suptitle = 'AMP SQUID V(Phi), column {0:}'.format(col)
    title = '(session ' + session_name + ')'
    xlabel = 'Demux ofco output (ADU)'

    plotFileName = 'amp_squid_VPhi_col{0:}.png'.format(col)
    fullPlotFileName = os.path.join(plot_path, plotFileName)

    start = "scan_amp_"
    files = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[-4:] == '{0:}.h5'.format(col) and f[:len(start)] == start]

    if len(files) != 1:
        raise ValueError(
            'Wrong number of files ({0:} instead of {1:})'.format(len(files), 1))

    plot_scan(os.path.join(data_path, files[0]), fullPlotFileName, suptitle, title, xlabel, verbose=verbose)


def scan_mux_squid(verbose):
    # Paths definition
    dir_path = os.path.join("..", "..")
    hk_path = os.path.join(dir_path, cst.hkDirName)
    data_path = os.path.join(dir_path, cst.scanDirName)
    plot_path = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(plot_path)
    session_name = os.path.basename(os.path.realpath(dir_path))
    col = session_name[-1]

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(hk_path)

    if verbose:
        print("/----------------------------------------------------------")
        print("/ MUX SQUID I(Phi) measurement")
        print("/ Test session name:   " + session_name)
        print("/ Tested column:       " + col)
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:     {0:}".format(fwVersion))
        print("/----------------------------------------------------------\n")

    title = '(session ' + session_name + ')'
    xlabel = 'Demux fdbk output (ADU)'

    plotFileName1 = 'mux_squid_IPhi_fll_off_error_col{0:}.png'.format(col)
    fullPlotFileName1 = os.path.join(plot_path, plotFileName1)
    suptitle1 = 'MUX SQUID V(Phi), AMP SQUID FLL OFF, column {0:}'.format(col)
    plotFileName2 = 'mux_squid_IPhi_fll_on_error_col{0:}.png'.format(col)
    fullPlotFileName2 = os.path.join(plot_path, plotFileName2)
    suptitle2 = 'MUX SQUID V(Phi), AMP SQUID FLL ON, ERROR data, column {0:}'.format(col)
    plotFileName3 = 'mux_squid_IPhi_fll_on_science_col{0:}.png'.format(col)
    fullPlotFileName3 = os.path.join(plot_path, plotFileName3)
    suptitle3 = 'MUX SQUID V(Phi), AMP SQUID FLL ON, SCIENCE data, column {0:}'.format(col)

    start1 = "scan_mux_IPhi_FLL-OFF_error_"
    file1 = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[-4:] == '{0:}.h5'.format(col) and f[:len(start1)] == start1]

    if len(file1) != 1:
        raise ValueError(
            'Wrong number of files ({0:} instead of {1:})'.format(len(file1), 1))

    plot_scan(os.path.join(data_path, file1[0]), fullPlotFileName1, suptitle1, title, xlabel, verbose=verbose)

    start2 = "scan_mux_IPhi_FLL-ON_error_"
    file2 = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[-4:] == '{0:}.h5'.format(col) and f[:len(start2)] == start2]

    if len(file2) != 1:
        raise ValueError(
            'Wrong number of files ({0:} instead of {1:})'.format(len(file2), 1))

    plot_scan(os.path.join(data_path, file2[0]), fullPlotFileName2, suptitle2, title, xlabel, verbose=verbose)

    start3 = "scan_mux_IPhi_FLL-ON_science_"
    file3 = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[-4:] == '{0:}.h5'.format(col) and f[:len(start3)] == start3]

    if len(file3) != 1:
        raise ValueError(
            'Wrong number of files ({0:} instead of {1:})'.format(len(file3), 1))

    plot_scan(os.path.join(data_path, file3[0]), fullPlotFileName3, suptitle3, title, xlabel, verbose=verbose)

# -------------------------------------------------------------------------------------
