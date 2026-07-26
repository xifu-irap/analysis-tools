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
#  tstPatternTools.py
#
# ---------------------------------------------------------------------------------

import os

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gt
import readData as rddt


def mk_tst_pattern(params):

    if len(params[:,0])!=cst.tstPattern_nb_regions or len(params[0,:])!=cst.tstPattern_nb_params_per_region:
        raise ValueError('params must be a {0:}x{1:}'.format(cst.tstPattern_nb_regions, cst.tstPattern_nb_params_per_region))

    # Computing length of one iteration
    l = 0
    for region in range(cst.tstPattern_nb_regions):
        l += (1+params[region, 1])*params[region, 3]

    # Creating data
    data = np.zeros(l)
    i_start = 0
    for region in range(cst.tstPattern_nb_regions):
        for step in range(params[region, 3]):
            data[i_start:i_start+(1+params[region, 1])] = np.ones(1+params[region, 1])*(params[region, 0] + step*params[region, 2])
            i_start += 1+params[region, 1]

    return data


def check_tst_pattern(dir_path, col, params, l=0):

    # Data directory
    pathData = os.path.join(dir_path, cst.dataDirName)

    # Session name
    session_name = os.path.basename(dir_path)

    # Creation of a directory for the plot files
    pathPlot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(pathPlot)
    plotFullFileName = os.path.join(pathPlot, "check_tst_pattern.png")

    # Reading test data
    files = [f for f in os.listdir(pathData) \
                if os.path.isfile(os.path.join(pathData, f)) \
                and f == 'tptAcqMode_C{}.fits'.format(col)]

    if len(files) == 0:
        raise ValueError('No files found')

    fits_file_path = os.path.join(pathData, files[0])
    tst_data, ctrl = rddt.get_science_from_fits(fits_file_path)

    # Building reference data
    ref_data = mk_tst_pattern(params)
    lref = len(ref_data)

    # Checking if data are shifted
    decal = np.where(ref_data == tst_data[0,0])[0][0]
    print("Data sets are shifted by {0:} frame(s)".format(decal))
    tst_data = tst_data[:, lref-decal:]

    # Keeping same length for both data set
    ratio = int(len(tst_data[0,:])/lref)
    tst_data = tst_data[:, :ratio*lref]
    ref_data = np.tile(ref_data, ratio)

    # Comparing data
    error_pix = []
    result = 'Test passed'
    color = 'g'
    plot_pix = 0

    for pix in range(cst.muxFactor):
        if np.all(ref_data != tst_data[pix,:]):
            error_pix.append(pix)
            result = 'Test failed for pixels ' + str(error_pix)
            color = 'r'
            plot_pix = pix
    print(result)

    # Doing the plots
    fig = plt.figure(figsize=(13, 10))
    fig.suptitle("Test pattern check", fontsize=14)

    if l == 0:
        l = len(ref_data)

    ax1 = fig.add_subplot(3, 1, 1)
    ax1.plot(ref_data[:l], 'k')
    ax1.set_title("Expected data")

    ax2 = fig.add_subplot(3, 1, 2)
    ax2.plot(tst_data[plot_pix,:l], 'b')
    ax2.set_title("Test data (session " + session_name +")")

    ax3 = fig.add_subplot(3, 1, 3)
    ax3.plot(tst_data[plot_pix,:l] - ref_data[:l], color=color, label = result)
    ax3.set_title("Error")
    ax3.set_xlabel("Frame number")
    ax3.legend(loc='upper right')

    fig.tight_layout()
    plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')

