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

xZoom = 400


def ofcoFineResoAndRange(verbose=False):
    # Paths definition
    dir_path = os.path.join("..", "..")
    hk_path = os.path.join(dir_path, cst.hkDirName)
    data_path = os.path.join(dir_path, cst.dataDirName)
    plot_path = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(plot_path)
    session_name = os.path.basename(os.path.realpath(dir_path))

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(hk_path)

    if verbose:
        print("/----------------------------------------------------------")
        print("/ Ofco Fine resolution and range characterisation ")
        print("/ Test session name:   " + session_name)
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:     {0:}".format(fwVersion))
        # print("/ Box car length:       {0:} samples".format(bxl))
        print("/----------------------------------------------------------\n")

    # Looking for test configuration parameters
    test_mode = session_name[10:15]

    xlabel = 'Time (ns)'
    ylabel = 'Error signal (V)'

    for colid in range(cst.nColPerDemux):

        plotFileName = os.path.join(plot_path, 'ofcoFine' + test_mode + '_col{0:}.png'.format(colid))

        files = [f for f in os.listdir(data_path) \
                 if os.path.isfile(os.path.join(data_path, f)) \
                 and f[-3:] == ".h5" and f[:5] == "dump_"]

        if len(files) == 0:
            raise ValueError('No dump files found!')

        fig = plt.figure(figsize=(12, 10))
        ax1 = fig.add_subplot(1, 1, 1)  # global plot
        if test_mode == "HRESO":
            suptit = "Characterisation of the OFCO FINE resolution (col {0:})".format(colid)
        else:
            suptit = "Characterisation of the OFCO FINE range (col {0:})".format(colid)
        plt.suptitle(suptit + '\n(' + session_name + ')')

        xTime = np.arange(2 * cst.nSamplesPerRow * cst.muxFactor) * 1e9 / cst.fSamp

        # Reading data from dump files
        colDumpsAccu = np.zeros((cst.nColPerDemux, 2 * cst.nPixPerCol * cst.nSamplesPerRow))
        for index, file in enumerate(np.sort(files)):
            colDumps, errors = rddt.read_dump_from_hdf5(os.path.join(data_path, file))
            colDumpsAccu += colDumps
        colDumpsAccu /= len(files)

        # Conversion to Volts
        colDumpsAccu *= cst.fsrADCErrorV / cst.fsrADCErrorADU

        ax1.plot(xTime[:], colDumpsAccu[colid, :], color=blue, linewidth=1)

        t_max = (xZoom - 1) * 1e9 / cst.fSamp
        ax1.set_xlim([0, t_max])

        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)

        # Définition des intervalles majeurs et mineurs pour la grille
        ax1.set_xticks(np.arange(0, t_max, 160))  # Intervalles majeurs
        ax1.set_xticks(np.arange(0, t_max, 8), minor=True)  # Intervalles mineurs

        # Activation de la grille majeure et mineure
        ax1.grid(which='major', linestyle='-', linewidth='0.6', color='black')  # Grille majeure
        ax1.grid(which='minor', linestyle='--', linewidth='0.4', color='gray')  # Grille mineure

        fig.tight_layout()

        plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
        if verbose:
            print("results plotted in file " + plotFileName)

# -------------------------------------------------------------------------------------
