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
#  ofcoDelayAnalysis.py
#
# ---------------------------------------------------------------------------------

# imports
import os
import re

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gt
import readData as rddt

xZoom = 49

colors = ['#FF0000', '#0000FF', '#006400', '#FFA500', '#800080', '#008B8B', '#8B0000',
          '#696969', '#A52A2A', '#000000', '#4682B4', '#556B2F', '#4B0082', '#B8860B',
          '#C0C0C0', '#CD5C5C', '#8B008B', '#FFD700', '#228B22', '#00008B']


def ofcoMuxDelayAnalysis(verbose=False):
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
        print("/ Ofco MUX delay test ")
        print("/ Test session name:   " + session_name)
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:     {0:}".format(fwVersion))
        # print("/ Box car length:       {0:} samples".format(bxl))
        print("/----------------------------------------------------------\n")

    # Looking for test configuration parameters
    split = re.split(r'[_-]+', session_name)  # splitting the session name with '_' and '-'
    start_str = split[-2]
    end_str = split[-1]
    start = int(start_str[1:])
    if start_str[0] == "M":
        start *= -1
    end = int(end_str[1:])
    if end_str[0] == "M":
        end *= -1
    nb_steps = end - start + 1

    xlabel = 'Time (ns)'

    for colid in range(cst.nColPerDemux):

        plotFileName = os.path.join(plot_path, 'ofcoMuxDelay_col{0:}.png'.format(colid))
        ylabel = 'Dump of error signal (ADU)'.format(colid)

        files = [f for f in os.listdir(data_path) \
                 if os.path.isfile(os.path.join(data_path, f)) \
                 and f[-3:] == ".h5" and f[:5] == "dump_"]

        if len(files) != nb_steps:
            raise ValueError('Wrong number of files ({0:} instead of {1:})'.format(len(files), nb_steps))

        fig = plt.figure(figsize=(12, 10))
        ax1 = fig.add_subplot(2, 1, 1)  # global plot
        ax2 = fig.add_subplot(2, 1, 2)  # zoom plot
        plt.suptitle("Characterisation of the ofco MUX delay compensation (column {0:})\n".format(colid) \
                     + '(' + session_name + ')')

        xTime = np.arange(2 * cst.nSamplesPerRow * cst.muxFactor) * 1e9 / cst.fSamp
        for index, file in enumerate(np.sort(files)):
            colDumps, errors = rddt.read_dump_from_hdf5(os.path.join(data_path, file))

            # Doing plot
            setting = start + index
            lw = 2
            if setting == 0:
                lw = 4
            lbl = 'Delay = {0:2d} ns'.format(int(setting * 1e9 / cst.fSamp))
            ax1.plot(xTime[:], colDumps[colid, :], color=colors[index], label=lbl, linewidth=1)
            ax2.plot(xTime[:xZoom], colDumps[colid, :xZoom], color=colors[index], label=lbl, linewidth=lw)

        x1_max = 2 * cst.nSamplesPerRow * cst.muxFactor * 1e9 / cst.fSamp
        x2_max = (xZoom - 1) * 1e9 / cst.fSamp

        ax1.set_ylabel(ylabel)
        ax2.set_ylabel(ylabel)
        ax1.set_xlabel(xlabel)
        ax2.set_xlabel(xlabel)
        # ax1.legend(loc='upper right')
        ax2.legend(loc='upper right')
        xlims = ax2.get_xlim()
        ax2.set_xlim([xlims[0], x2_max])

        # Définition des intervalles majeurs et mineurs pour la grille
        ax1.set_xticks(np.arange(0, x1_max, 4096))  # Intervalles majeurs tous les 64
        ax1.set_xticks(np.arange(0, x1_max, 512), minor=True)  # Intervalles mineurs tous les 8
        ax2.set_xticks(np.arange(0, x2_max, 64))  # Intervalles majeurs tous les 64
        ax2.set_xticks(np.arange(0, x2_max, 8), minor=True)  # Intervalles mineurs tous les 8

        # Activation de la grille majeure et mineure
        ax1.grid(which='major', linestyle='-', linewidth='0.6', color='black')  # Grille majeure
        ax1.grid(which='minor', linestyle='--', linewidth='0.4', color='gray')  # Grille mineure
        ax2.grid(which='major', linestyle='-', linewidth='0.6', color='black')  # Grille majeure
        ax2.grid(which='minor', linestyle='--', linewidth='0.4', color='gray')  # Grille mineure

        fig.tight_layout()

        plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
        if verbose:
            print("results plotted in file " + plotFileName)

# -------------------------------------------------------------------------------------
