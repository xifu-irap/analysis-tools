#!/usr/bin/env python
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
#  laurent.ravera@irap.omp.eu
#  measure_ki.py
#

import numpy as np
import os
import matplotlib.pyplot as plt

import general_tools, get_data

def measure_ki(ki_data):
    r"""
        This function analyses a SQ1 feedback and SQ2 feedback dump to
        measure the ki gain. The function makes plot in order to measure
        the step size on both signals.

        Parameters
        ----------
        ki_data: data object
        Contains the SQ1 feedback and SQ2 feedback dumps.

        Returns
        -------
        Nothing

        """
    """
    Defining path and file names
    """
    plotdirname = os.path.join(os.path.normcase(ki_data.config['path']), ki_data.config['dir_plots'])
    general_tools.checkdir(plotdirname)
    plotfilename=os.path.join(plotdirname,'loop_gain_measurement.png')

    print("Measuring ki ...")
    data1, name1 = ki_data.values[:, 0], "Feedback SQ1"
    data2, name2 = ki_data.values[:, 1], "Feedback return"

    # In these dumps the 5MHz data are over sampled at 20MHz
    fs = float(ki_data.config["frow"])*4 # Approx 20MHz
    t = np.arange(len(data1))/fs

    """
    Looking for step
    """
    beginning = data1[:100].mean()
    end = data1[-100:].mean()
    threshold = 0.5*(end+beginning)
    lower = data1 < threshold
    step = (lower[:-1] & ~lower[1:]) | (lower[1:] & ~lower[:-1])
    i_step = np.where(step)[0][0]

    """
    Computing ki
    """
    width = 40
    npts_delay = int(ki_data.config['npix']*4)
    range_feedback1=i_step-width-2+np.arange(width)
    range_feedback2=i_step+2+np.arange(width)
    range_return1=i_step+npts_delay-width-2+np.arange(width)
    range_return2=i_step+npts_delay+2+np.arange(width)
    step_size_feedback = data1[range_feedback2].mean() - data1[range_feedback1].mean()
    step_size_return = data2[range_return2].mean() - data2[range_return1].mean()
    ki = -1*step_size_return/step_size_feedback
    print("ki is equal to: {0:4.3f}".format(ki))

    """
    Doing the plot
    """
    xtitle = "Time (ms)"
    mkr='.'
    i1 = i_step - width - 20
    i2 = i_step + npts_delay + width + 20

    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(1000*t[i1:i2], data1[i1:i2], 'b', label="Feedback SQ1 signal")
    ax1.plot(1000*t[i1:i2], data2[i1:i2], 'k', label="Return signal")
    ax1.set_ylabel(name2)
    ax1.set_xlabel(xtitle)
    ax1.grid(color='k', linestyle=':', linewidth=0.5)

    ax1.plot(1000*t[range_feedback1], data1[range_feedback1], 'r', marker=mkr, label="Feedback step size ={0:6.0f}".format(step_size_feedback))
    ax1.plot(1000*t[range_feedback2], data1[range_feedback2], 'r', marker=mkr)
    ax1.plot(1000*t[range_return1], data2[range_return1], 'g', marker=mkr, label="Return step size ={0:6.0f} ==> loop gain ={1:4.3f}".format(step_size_return, ki))
    ax1.plot(1000*t[range_return2], data2[range_return2], 'g', marker=mkr)
    ax1.set_ylabel(name1)
    ax1.set_xlabel(xtitle)
    ax1.grid(color='k', linestyle=':', linewidth=0.5)
    ax1.legend(loc="best")

    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename, bbox_inches='tight')

