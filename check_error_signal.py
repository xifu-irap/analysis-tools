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
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{pslatex,amsmath}')

#---------------------------------------------------------------------------------    
fs=125e6
time_label_us='Time (µs)'
samples_label='Sample number'

#---------------------------------------------------------------------------------    
def rising_edge_fit(t, a, b, c, d):
    return(a*np.exp(-(t-b)/c)+d)



#---------------------------------------------------------------------------------    
def check_cutoff_freq(data, fs, plot=True):

    from scipy.optimize import curve_fit

    depth_before = 1
    depth_after = 40
    ratio=0.1
    threshold = data.min()+ratio*(data.max()+data.min())
    i_high = np.where(data > threshold)
    
    i_first = max(0, i_high[0][0]-depth_before)
    i_last = i_high[0][0]+depth_after
    
    data_rise = data[i_first: i_last]
    samples_rise = i_first + np.arange(depth_before+depth_after)
    #time_rise = samples_rise / fs
    
    # Fitting the rising edge
    guess = [-(data.max()-data.min()), i_high[0][0], 5, data_rise[-1]]
    print(guess)
    popt, pcov = curve_fit(rising_edge_fit, samples_rise, data_rise, p0=guess, bounds=([-2**12, 0, 0, 0], [0, 500, 20, 2**12]))
    fit = rising_edge_fit(samples_rise, popt[0], popt[1], popt[2], popt[3])

    # Computation of the time constant and the 3dB cut frequency    
    tau = popt[2]/fs
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
        ax2.set_xlabel(samples_label)
        ax2.set_ylim(data.min(), data.max())
        ax2.legend(loc='best')
            
        fig.tight_layout()
            

#---------------------------------------------------------------------------------    
def check_sampling_freq(data, threshold_ratio=0.1, mux=34, nsamp=20, plot=True):
    """This function checks the periodicity of a DEMUX ADC dump in order to 
    estimate the ADC sampling frequency. 
    For this test the DEMUX ADC has to be fed with a square wave signal at Fframe.
    The ADC dump is shifted by 1 frame and the difference (which is supposed to 
    be null) is measured.
    
    Args:
        data (array): DEMUX ADC dump data
        threshold_ratio (float, optional): Above this threshold the difference is 
            considered as too high. The threshold is expressed as a fraction of the 
            input signal amplitude. Defaults to 0.1.
        mux (int, optional): Mux factor. Defaults to 34.
        nsamp (int, optional): Number of samples per row. Defaults to 20.
        plot (bool, optional): If True a plot is done. Defaults to True.

    Raises:
        ValueError: If the length of the input data is not equal to mux x nsamp
            a ValueError is raised.
    """
    
    npts=mux*nsamp*2

    if len(data)!=npts:
        raise ValueError("The length of the data is incorrect. {0:d2} x {1:d2} samples where expected.".format(mux, nsamp))
    else:
        data_roll=np.roll(data, -mux*nsamp)
        diff=data-data_roll
        print(diff)
        amp=data.max()-data.min()
        limits_ratio=1.5
        threshold=threshold_ratio*(amp)
        print(threshold_ratio, amp, threshold)
        if plot:
            fig = plt.figure(figsize=(6, 8))
            ax1 = fig.add_subplot(2, 1, 1)
            ax1.plot(data, label="Input data (2 frames, 2 x {0:d} x {1:d} = {2:d} samples)".format(mux, nsamp, npts))
            ax1.plot(data_roll, label='Input data rolled by 1 frame = {0:d} x {1:d} samples'.format(mux, nsamp))
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
                ax2.text(npts*0.02,-1.3*threshold, 'If there are 2 square periods in the upper plot, \n the sampling frequency is equal to the input square frequency x {0:d} x {1:d}'.format(mux, nsamp), color='green')
            else:
                ax2.text(npts*0.02,-1.3*threshold, 'The number of samples per period and the sampling frequency are probably incorrect\n', color='red')
                
            fig.tight_layout()
            

#---------------------------------------------------------------------------------    
def butter_lowpass(cutoff, fs, order=1):
    from scipy.signal import butter
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)

    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=1):
    from scipy.signal import lfilter
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)

    return y


#---------------------------------------------------------------------------------    
def mk_fake_square_error_data(fs:125e6, amp=1, cutoff=6e6, ok=True, mux=34, nsamp=20, plot=True):
    """This function computes a fake data set to simulate the content of a 
    DEMUX ADC dump when a square wave is fed into the error input.

    Args:
        fs (Boolean): Expected sampling frequency. Defaults to 125e6.
        amp (int, optional): Square wave amplitude. Defaults to 1.
        cutoff (float, optional): cutoff frequency in Hz. Defaults 6MHz.
        ok (bool, optional): If OK the data set shall correspond to a correct 
            sampling frequency. Defaults to True.
        mux (int, optional): Multiplexing factor. Defaults to 34.
        nsamp (int, optional): Number of samples per row. Defaults to 20.
        plot (bool, optional): If True a plot is done. Defaults to True.
    """
    from scipy import signal

    npts = 2*mux*nsamp
    sigma_noise=amp*0.01

    if ok:
        nb_periods=2  # The period corresponds to fs=125MHz
    else:
        nb_periods=2.1  # The period does not correspond to 125MHz

    # Computing square signal
    shift_radians=np.pi*3/2
    x=shift_radians+np.arange(npts)*2*np.pi*nb_periods/npts    
    fake_square = amp*(1+signal.square(x))

    # Time values
    t=np.arange(npts)/fs
                
    # Adding some noise
    fake_square += sigma_noise*np.random.randn(npts)

    # Simulating analog low pass filter 
    order=1    
    fake_sig = butter_lowpass_filter(fake_square, cutoff, fs, order)
        
    if plot:
        fig = plt.figure(figsize=(6, 4))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.plot(t*1e6, fake_square)
        ax1.plot(t*1e6, fake_sig)
        ax1.set_xlabel(time_label_us)
        fig.tight_layout()
    
    return(t, fake_sig)

#---------------------------------------------------------------------------------    

_, sig = mk_fake_square_error_data(fs, ok=True, plot=False)

#check_sampling_freq(sig, plot=True)
check_cutoff_freq(sig, fs, plot=True)

plt.show()

#---------------------------------------------------------------------------------    
