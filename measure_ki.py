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

import general_tools

def measure_ki(ki_data, verbose=False):
    r"""
        This function analyses a SQ1 feedback and SQ2 feedback dump to
        measure the ki gain. The function makes plot in order to measure
        the step size on both signals.

        Parameters
        ----------
        ki_data: data object
        Contains the SQ1 feedback and SQ2 feedback dumps.

        verbose: boolean
        if True some debug messages are printed

        Returns
        -------
        Nothing

        """
    """
    Defining path and file names
    """
    datafilename=ki_data.config['datafilename']
    plotdirname = os.path.join(os.path.normcase(ki_data.config['path']), ki_data.config['dir_plots'])
    general_tools.checkdir(plotdirname)
    plotfilename_steps=os.path.join(plotdirname, datafilename[:-4]+'_steps.png')
    plotfilename_gains=os.path.join(plotdirname, datafilename[:-4]+'_loop-gains.png')

    print("Measuring ki ...")

    npix = int(ki_data.config["npix"])

    """
    Keeping a single value per row
    (In these dumps the 5MHz data are over sampled at 20MHz : x4)
    """
    ratio=4
    n_rows = len(ki_data.values[:,0])//ratio
    data1 = np.resize(ki_data.values[:n_rows*ratio, 0], (n_rows, ratio))
    data2 = np.resize(ki_data.values[:n_rows*ratio, 1], (n_rows, ratio))
    data1 = data1[:,0]
    data2 = data2[:,0]

    """
    demultiplexing data
    """
    n_frames = len(data1)//npix
    data1 = np.transpose(np.resize(data1[:n_frames*npix], (n_frames, npix)))
    data2 = np.transpose(np.resize(data2[:n_frames*npix], (n_frames, npix)))

    """
    Looking for step for each pixel
    """
    i_step=np.empty(npix, np.int32)
    window_size=10
    for pixel in range(npix):
        beginning = data1[pixel,:window_size].mean()
        end = data1[pixel,-window_size:].mean()
        threshold = 0.5*(end+beginning)
        lower = data1[pixel,:] < threshold
        step = (lower[:-1] & ~lower[1:]) | (lower[1:] & ~lower[:-1])
        i_step[pixel] = np.arange(len(step))[step][0]
    
    """
    Computing ki for each pixel
    """
    decal = 0  # shift of pixel index between the 2 outputs (should be 0)
    delay = 1  # round-trip delay of the step is equal to 1 frame
    ki=np.empty(npix)
    for pixel in range(npix):
        step_size_feedback=data1[pixel,i_step[pixel]+1]-data1[pixel,i_step[pixel]]
        step_size_return=data2[(pixel+decal)%npix,i_step[(pixel+decal)%npix]+delay+1]\
                        -data2[(pixel+decal)%npix,i_step[(pixel+decal)%npix]+delay]
        ki[pixel] = -1*step_size_return/step_size_feedback

        if verbose:
            print("return before: ", data2[(pixel+decal)%npix,i_step[(pixel+decal)%npix]+delay])
            print("return after: ", data2[(pixel+decal)%npix,i_step[(pixel+decal)%npix]+delay+1])
            print("Return step: ", step_size_return)

    print("Pixel's loop gains are:\n", ki)

    """
    Doing the plots
    """
    # Plotting the loop gain per pixel
    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.scatter(np.arange(npix), ki)
    ax1.set_ylabel("Loop gain")
    ax1.grid(color='k', linestyle=':', linewidth=0.5)

    ax2 = fig.add_subplot(2, 1, 2)
    ax2.scatter(np.arange(npix), ki)
    ax2.set_ylabel("Loop gain")
    ax2.set_xlabel("Pixel id")
    ax2.set_ylim([0,2])
    ax2.grid(color='k', linestyle=':', linewidth=0.5)

    for item in ([ax1.yaxis.label, ax2.yaxis.label, ax2.xaxis.label]):
        item.set_weight('bold')
        item.set_fontsize(14)
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels() + ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(12)
    fig.suptitle(datafilename)

    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename_gains, bbox_inches='tight')

    # Plotting the steps on the feedback and rturn signals
    xtitle = "Time (ms)"
    fs = float(ki_data.config["frow"]) * ratio
    i1 = (i_step[0]-1)*ratio*npix
    i2 = (i_step[0]+3)*ratio*npix
    t = np.arange(i2-i1)*1e3/fs

    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(1, 1, 1)
    #print('t ', t[i1:i2])
    ax1.plot(t, ki_data.values[i1:i2,1], 'k', label="Return signal")
    ax1.plot(t, ki_data.values[i1:i2,0], 'b', label="Feedback SQ1 signal")
    ax1.set_ylabel("DACs outputs (ADU)")
    ax1.set_xlabel(xtitle)
    ax1.grid(color='k', linestyle=':', linewidth=0.5)
    ax1.set_xlabel(xtitle)
    ax1.grid(color='k', linestyle=':', linewidth=0.5)
    ax1.legend(loc="best")

    for item in ([ax1.xaxis.label, ax1.yaxis.label]):
        item.set_weight('bold')
        item.set_fontsize(14)
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(12)
    fig.suptitle(datafilename)

    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename_steps, bbox_inches='tight')
    
# -----------------------------------------------------------------------------
