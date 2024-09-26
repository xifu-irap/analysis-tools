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
#  fake_data.py
#

#---------------------------------------------------------------------------------    
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

import constants as cst
import fits_tools

#---------------------------------------------------------------------------------
time_label_us='Time (µs)'


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
def make_fake_adc_dump_square(file_name, amp=1, cutoff=6e6, ok=True, plot=False):
    """This function computes a fake data set to simulate the content of a 
    DEMUX ADC dump when a square wave is fed into the error input. The dumps
    is stored in a fits file.

    Args:
        file_name (string): Name of the fits file.
        amp (int, optional): Square wave amplitude. Defaults to 1.
        cutoff (float, optional): Cutoff frequency in Hz. Defaults 6MHz.
        ok (bool, optional): If ok the dataset shall correspond to a correct 
            sampling frequency. Defaults to True.
        plot (bool, optional): If True a plot is done. Defaults to False.
    """
    from scipy import signal

    npts = 2 * cst.mux_factor * cst.n_samples_per_row
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
    t=np.arange(npts)/cst.fsamp
                
    # Adding some noise
    fake_square_noise = fake_square + sigma_noise*np.random.randn(npts)

    # Simulating analog low pass filter 
    order=1    
    fake_sig = butter_lowpass_filter(fake_square_noise, cutoff, cst.fsamp, order)
    
    # Simulating quantization
    FSR = 2**cst.adc_nbits-1
    ratio = 0.8 # percentage of FSR coverage
    fake_sig -= fake_sig.min()
    fake_sig = (fake_sig * FSR / (fake_sig.max() - fake_sig.min()) - FSR/2) * ratio
    fake_sig = fake_sig.astype(int)
    fake_square = (fake_square * FSR / fake_square.max() - FSR/2) * ratio

    fits_tools.write_adc_dump_fits_file(file_name, fake_sig)
    print(" Fake dump stored in file: ", file_name)
            
    if plot:
        fig = plt.figure(figsize=(6, 4))
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.plot(t*1e6, fake_square)
        ax1.plot(t*1e6, fake_sig)
        ax1.set_xlabel(time_label_us)
        fig.tight_layout()

        plt.show()        
            
#---------------------------------------------------------------------------------    

