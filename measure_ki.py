#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt

import general_tools, get_data

def measure_ki(ki_data):
    """
    Defining path and file names
    """
    plotdirname = os.path.join(os.path.normcase(ki_data.config['path']), ki_data.config['dir_plots'])
    general_tools.checkdir(plotdirname)
    plotfilename=os.path.join(plotdirname,'ki_measurement.png')

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
    ax1.plot(1000*t[i1:i2], data2[i1:i2], 'k')
    ax1.set_ylabel(name2)
    ax1.set_xlabel(xtitle)
    ax1.grid(color='k', linestyle=':', linewidth=0.5)

    ax1.plot(1000*t[i1:i2], data1[i1:i2], 'b')
    ax1.plot(1000*t[range_feedback1], data1[range_feedback1], 'r', marker=mkr, label="Feedback signal: step size ={0:5.0f}".format(step_size_feedback))
    ax1.plot(1000*t[range_feedback2], data1[range_feedback2], 'r', marker=mkr)
    ax1.plot(1000*t[range_return1], data2[range_return1], 'g', marker=mkr, label="Return signal: step size ={0:5.0f}\n ==>  ki={1:4.3f}".format(step_size_return, ki))
    ax1.plot(1000*t[range_return2], data2[range_return2], 'g', marker=mkr)
    ax1.set_ylabel(name1)
    ax1.set_xlabel(xtitle)
    ax1.grid(color='k', linestyle=':', linewidth=0.5)
    ax1.legend(loc="best")


    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename, bbox_inches='tight')

