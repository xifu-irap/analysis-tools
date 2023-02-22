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
#  check_error_signal.py
#


#---------------------------------------------------------------------------------    
import numpy as np
import matplotlib.pyplot as plt

import constants as cst
import fake_data
import fits_tools
from general_tools import madate

#---------------------------------------------------------------------------------    
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{pslatex,amsmath}')

#---------------------------------------------------------------------------------    
samples_label='Sample number'


#---------------------------------------------------------------------------------    
def rising_edge_fit(t, a, b, c, d):
    """This function is used to fit the rising edge of adc dump.

    Args:
        t (float): time input array
        a (float): output scaling factor
        b (float): shift in time (depends on the rising edge position)
        c (float): time constant
        d (float): dump offset
    """
    return(a*np.exp(-(t-b)/c)+d)


#---------------------------------------------------------------------------------    
def measure_cutoff_freq(file_name, plot=True, guess_plot=False):
    """This function reads an ADC dump from a fits file and measures the cut off 
    frequency. The dump shall correspond to a measuremnt with a square signal at 
    DEMUX input.

    Args:
        file_name (string): Name of the fits file.
        plot (bool, optional): If True a plot is done. Defaults to True.
        guess_plot (bool, optional): If True the fit with first-guess parameters
                                     is over plotted. Default to False.
    """
    from scipy.optimize import curve_fit

    # Getting the ADC dump data from the fits file
    data = fits_tools.read_adc_dump_fits_file(file_name)

    depth_before = 1
    depth_after = 40
    ratio=0.1
    threshold = data.min()+ratio*(data.max()-data.min())
    i_high = np.where(data > threshold)
    
    i_first = max(0, i_high[0][0]-depth_before)
    i_last = i_high[0][0]+depth_after
    
    data_rise = data[i_first: i_last]
    samples_rise = i_first + np.arange(depth_before+depth_after)
    
    # Fitting the rising edge
    guess = [-(data.max()-data.min()), i_high[0][0], 4, data_rise[-1]]
    popt, pcov = curve_fit(rising_edge_fit, samples_rise, data_rise, p0=guess, bounds=([-2**14, 0, 0, 0], [0, 500, 20, 2**14]))
    fit = rising_edge_fit(samples_rise, popt[0], popt[1], popt[2], popt[3])
    fit_guess = rising_edge_fit(samples_rise, guess[0], guess[1], guess[2], guess[3])

    # Computation of the time constant and the 3dB cut frequency    
    tau = popt[2]/cst.fsamp
    fc = 1/(2*np.pi*tau)

    # Plot
    if plot:
        fig = plt.figure(figsize=(6, 8))
        ax1 = fig.add_subplot(2, 1, 1)
        ax1.plot(data, label="Input data")
        ax1.set_xlabel(samples_label)
        ax1.legend(loc='best')
        ax2 = fig.add_subplot(2, 1, 2)
        ax2.plot(samples_rise, data_rise, '.', label="Input data")
        ax2.plot(samples_rise, fit, label=r'Fit ($\tau$={0:5.2f}ns, $f_c$={1:5.2f}MHz)'.format(tau*1e9, fc/1e6))
        if guess_plot:
            ax2.plot(samples_rise, fit_guess, '--', color='grey', label=r'Fit with first-guess parameters')
        ax2.set_xlabel(samples_label)
        extra = 0.2*(data.max()-data.min())
        ax2.set_ylim(data.min()-extra, data.max()+extra)
        ax2.legend(loc='best')
        fig.tight_layout()
    
        plot_file_name = file_name.split('.')[0]+'_cutoff_freq.png'
        plt.savefig(plot_file_name, dpi=300, bbox_inches='tight')
        print(" Cutoff frequency measurement plotted in file: ", plot_file_name)
        

#---------------------------------------------------------------------------------    
def measure_sampling_freq(file_name, threshold_ratio=0.1, plot=True):
    """This function checks the periodicity of a DEMUX ADC dump in order to 
    estimate the ADC sampling frequency. 
    For this test the DEMUX ADC has to be fed with a square wave signal at Fframe.
    The ADC dump is shifted by 1 frame and the difference (which is supposed to 
    be null) is measured.
    
    Args:
        file_name (string): Fits file containing the the ADC dump data
        threshold_ratio (float, optional): Above this threshold the difference is 
            considered as too high. The threshold is expressed as a fraction of the 
            input signal amplitude. Defaults to 0.1.
        plot (bool, optional): If True a plot is done. Defaults to True.

    Raises:
        ValueError: If the length of the input data is not equal to mux x nsamp
            a ValueError is raised.
    """
    
    # Getting the ADC dump data from the fits file
    data = fits_tools.read_adc_dump_fits_file(file_name)

    npts = 2 * cst.mux_factor * cst.n_samples_per_row

    if len(data)!=npts:
        raise ValueError("The length of the data is incorrect. {0:d2} x {1:d2} samples where expected."\
            .format(cst.mux_factor, cst.n_samples_per_row))
    else:
        data_roll=np.roll(data, -cst.mux_factor * cst.n_samples_per_row)
        diff=data-data_roll
        amp=data.max()-data.min()
        limits_ratio=1.5
        threshold=threshold_ratio*(amp)
        if plot:
            fig = plt.figure(figsize=(6, 8))
            ax1 = fig.add_subplot(2, 1, 1)
            ax1.plot(data, label="Input data (2 frames, 2 x {0:d} x {1:d} = {2:d} samples)".format(cst.mux_factor, cst.n_samples_per_row, npts))
            ax1.plot(data_roll, label='Input data rolled by 1 frame = {0:d} x {1:d} samples'.format(cst.mux_factor, cst.n_samples_per_row))
            ax1.set_xlabel(samples_label)
            ax1.set_xlim(0,npts)
            ax1.set_ylim([data.min()-limits_ratio*amp/4,data.max()+limits_ratio*amp/4])
            ax1.legend(loc='best')
            ax2 = fig.add_subplot(2, 1, 2)
            ax2.plot(diff, label='Data - Rolled data')
            ax2.plot([0, npts-1], [threshold, threshold], '--k', label='Threshold')
            ax2.plot([0, npts-1], [-threshold, -threshold], '--k')
            ax2.set_xlabel(samples_label)
            ax2.set_xlim(0, npts)
            ax2.set_ylim(-threshold*limits_ratio,threshold*limits_ratio)
            ax2.legend(loc='upper right')
            if abs(diff).max()<threshold:
                ax2.text(npts*0.02,-1.3*threshold, 'If there are 2 square periods in the upper plot, \n the sampling frequency is equal to the input square frequency x {0:d} x {1:d}'\
                    .format(cst.mux_factor, cst.n_samples_per_row), color='green')
            else:
                ax2.text(npts*0.02,-1.3*threshold, 'The number of samples per period and the sampling frequency are probably incorrect\n', color='red')
                
            fig.tight_layout()

            plot_file_name = file_name.split('.')[0]+'_measure_fs.png'
            plt.savefig(plot_file_name, dpi=300, bbox_inches='tight')
            print(" Sampling measurement plotted in file: ", plot_file_name)

            

#---------------------------------------------------------------------------------    

file_name = "fake_dump_" + madate() + ".fits"
fake_data.make_fake_adc_dump_square(file_name, ok=True)

measure_sampling_freq(file_name, plot=True)
measure_cutoff_freq(file_name)


#---------------------------------------------------------------------------------    
