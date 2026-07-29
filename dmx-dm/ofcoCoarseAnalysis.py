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
#  ofcoCoarseAnalysis.py
#
# ---------------------------------------------------------------------------------

# imports
import os

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gt
import readData as rddt


def ofcoCoarseReso(verbose=False):
    # Paths definition
    dir_path = os.path.join("..", "..")
    hk_path = os.path.join(dir_path, cst.hkDirName)
    data_path = os.path.join(dir_path, cst.scanDirName)
    plot_path = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(plot_path)
    session_name = os.path.basename(os.path.realpath(dir_path))

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(hk_path)

    if verbose:
        print("/----------------------------------------------------------")
        print("/ Ofco Coarse setting resolution ")
        print("/ Test session name:   " + session_name)
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:     {0:}".format(fwVersion))
        print("/----------------------------------------------------------\n")

    xlabel = 'Frame number'
    ylabel = 'Ofco --> Error signal (V)'

    for colid in range(cst.nColPerDemux):

        plotFileName = os.path.join(plot_path, 'ofcoCoarse_col{0:}.png'.format(colid))

        # Searching scan files
        files = [f for f in os.listdir(data_path) \
                 if os.path.isfile(os.path.join(data_path, f)) \
                 and f[:5] == "scan_" \
                 and f[-4] == "{0:}".format(colid) \
                 and f[-3:] == ".h5"]

        if len(files) == 0:
            raise ValueError('No scan files found!')
        else:
            print("Found {0:} scan files".format(len(files)))

        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(1, 1, 1)  # global plot
        plt.suptitle("Characterisation of the OFCO COARSE resolution (col {0:})".format(colid)
                     + '\n(' + session_name + ')')

        file = files[0]
        if verbose:
            print("Reading data from file ", file)
        scan_type, ctrl, ofco, error = rddt.read_scan(os.path.join(data_path, file))
        error = np.mean(error[:, :], axis=0)  # averaging the error value of all pixels

        # Conversion to Volts
        error *= cst.fsrADCErrorV / cst.fsrADCErrorADU

        ax1.plot(error, '.', markersize=0.5)
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)

        fig.tight_layout()

        plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
        if verbose:
            print("results plotted in file " + plotFileName)

# -------------------------------------------------------------------------------------
