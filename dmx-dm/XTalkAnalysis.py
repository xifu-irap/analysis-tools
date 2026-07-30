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
#  XTalkAnalysis.py
#
# ---------------------------------------------------------------------------------

# imports
import os

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gt
import readData as rddt


# TODO: Add the model / module name on the plots

def get_dump(pathData, verbose=True):
    # Constants
    ## Number of files
    ndumps_per_measure = 32
    n_measures = 2 * cst.nColPerDemux
    n_dumps = ndumps_per_measure * n_measures
    ## Area where no signals are expected (a pulse is done on pixel 0)
    pix_min = 3
    pix_max = cst.muxFactor - 2

    # Searching dump files
    files = [f for f in os.listdir(pathData) \
             if os.path.isfile(os.path.join(pathData, f)) \
             and f[-3:] == ".h5" and f[:5] == "dump_"]
    if (len(files) != n_dumps):
        raise ValueError("Wrong number of dump files. Found {0:}, {1:} were expected".format(len(files, n_dumps)))

    files = np.sort(files)  # Sorting the files in alphabetical order (i.e. chronological order)

    # Reading and accumulating the dumps for loop back configuration
    if verbose:
        print("   Reading dump files for the loop back configuration...")
    dumps_lb = np.zeros((cst.nColPerDemux, cst.nColPerDemux, 2 * cst.nSamplesPerRow * cst.muxFactor))
    for col_id in range(cst.nColPerDemux):
        if verbose:
            print("      Column {0:}".format(col_id))
        for file in files[col_id * ndumps_per_measure:(col_id + 1) * ndumps_per_measure]:
            dump, _ = rddt.read_dump_from_hdf5(os.path.join(pathData, file))
            dumps_lb[col_id] += dump
        dumps_lb[col_id] /= ndumps_per_measure
        # Setting the baseline to 0
        for col in range(cst.nColPerDemux):
            avg = dumps_lb[col_id, col, pix_min * cst.nSamplesPerRow:pix_max * cst.nSamplesPerRow].mean()
            dumps_lb[col_id, col, :] -= avg

    # Reading and accumulating the dumps for 100 Ohm configuration
    if verbose:
        print("   Reading dump files for the 100 Ohms loaded configuration...")
    dumps_100_ohms = np.zeros((cst.nColPerDemux, cst.nColPerDemux, 2 * cst.nSamplesPerRow * cst.muxFactor))
    for col_id in range(cst.nColPerDemux):
        if verbose:
            print("      Column {0:}".format(col_id))
        i0 = cst.nColPerDemux * ndumps_per_measure  # offset to read the following files
        for file in files[i0 + col_id * ndumps_per_measure:i0 + (col_id + 1) * ndumps_per_measure]:
            dump, _ = rddt.read_dump_from_hdf5(os.path.join(pathData, file))
            dumps_100_ohms[col_id] += dump
        dumps_100_ohms[col_id] /= ndumps_per_measure
        # Setting the baseline to 0
        for col in range(cst.nColPerDemux):
            dumps_100_ohms[col_id, col, :] -= dumps_100_ohms[
                col_id, col, pix_min * cst.nSamplesPerRow:pix_max * cst.nSamplesPerRow].mean()

    return dumps_lb, dumps_100_ohms


def xtalkAnalysis(verbose=True):
    # Paths definition
    dir_path = os.path.join("..", "..")
    hk_path = os.path.join(dir_path, cst.hkDirName)
    data_path = os.path.join(dir_path, cst.dataDirName)
    plot_path = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(plot_path)
    session_name = os.path.basename(os.path.realpath(dir_path))

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(hk_path)

    # Looking for test configuration parameters
    index_bxl = session_name.find("BXL") + len("BXL")
    bxl = int(session_name[index_bxl])  # boxcar length

    # Looking for the signal characterised by the test (fdbk or ofco)
    index_sig = session_name.find("PERP-") + len("PERP-")
    signal = session_name[index_sig:index_sig + 4]

    if verbose:
        print("/----------------------------------------------------------")
        print("/ " + signal + " crosstalk characterisation ")
        print("/ Test session name:   " + session_name)
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:     {0:}".format(fwVersion))
        print("/ Box car length:       {0:} samples".format(bxl))
        print("/----------------------------------------------------------\n")

    dumps_lb, dumps_100_ohms = get_dump(data_path, verbose)

    plot_fdbk_2_error_xtalk(dumps_lb, dumps_100_ohms, signal, plot_path)


def plot_fdbk_2_error_xtalk(dumps_lb, dumps_100_ohms, sig, pathPlot):

    xlabel = 'Time (ns)'
    ylabel = 'Error (ADU)'
    col1, col2, col3, col2grid = 'k', 'b', 'r', 'lightblue'

    for c_perp in range(cst.nColPerDemux):

        plotFileName = os.path.join(pathPlot, 'XTalk_' + sig + '2ERRO_cPerp{0:}.png'.format(c_perp))

        # Doing the plot
        t = np.arange(2 * cst.nSamplesPerRow * cst.muxFactor) * 1e9 / cst.fSamp
        suptitle = 'Crosstalk between ' + sig + ' and ERROR'
        title1 = 'Perpetrator column {0:}: '.format(c_perp) + sig + ' looped back to ERROR (reference)'

        fig = plt.figure(figsize=(9, 11))
        plt.suptitle(suptitle)

        ax1 = fig.add_subplot(5, 1, 1)
        ax1.plot(t, dumps_lb[c_perp, c_perp, :], color=col1)
        ax1.set_title(title1)
        ax1.set_ylabel(ylabel, color=col1)
        ax1.grid(color=col2grid)

        for col in range(cst.nColPerDemux):
            xtalk = np.abs(dumps_100_ohms[c_perp, col, :]).max() / dumps_lb[c_perp, c_perp, :].max()
            xtalk_db = 20 * np.log10(xtalk).min()

            title = ('Victim column {0}: 100 Ohms on error. '.format(col)
                     + sig + ' loaded by 100 Ohms (Crosstalk = {0:2.1f}dB)'.format(xtalk_db))

            ax = fig.add_subplot(5, 1, col + 2)
            ax.plot(t, dumps_100_ohms[c_perp, col, :], color=col2)
            ax.set_title(title)
            ax.set_ylabel(ylabel)
            ax.grid(color=col2grid)

        ax.set_xlabel(xlabel)

        fig.tight_layout()

        plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
        print("results plotted in file " + plotFileName)

#------------------------------------------------------------------

# TODO: Clarify the Xtalk test configuration and update the analysis script accordingly
