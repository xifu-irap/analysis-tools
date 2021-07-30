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
#  ep_tools.py
#

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize.minpack import curve_fit
from astropy.io import fits

import general_tools, get_data, plotting_tools

JITTER_MARGIN=10

"""
Label definitions
"""
sample_nb_label = 'Sample number'

# ############################################################
# Function to read noise records
# ############################################################
def get_noise_records(noise_file, noise_p):
    '''
    This function builds noise records from row dump files  
    
    Arguments:
        - noise_file: string
          dump file name
        - config: dictionnary
          contains path and constants definitions
        - noise_p: dictionnary including several parameters
            - noise_p['record_len']: number
            length of the records to be applied
            - noise_p['pix']: number
            index of the pixel to be processed
            - noise_p['remove_raw_files']: boolean
            indicates if the raw dta files shall be removed
    '''
    print("Getting noise data from file: ", noise_file)
    noisedat = get_data.data(noise_file)
    if noise_p['remove_raw_files']:
        dirname = os.path.join(os.path.normcase(noisedat.config['path']), noisedat.config['dir_data'])
        fulldatafilename = os.path.join(dirname, noise_file)
        os.remove(fulldatafilename)
    # selecting the data of the requested pixel
    # (offset of 1 is applyed) because the first line is the readout index
    noise = noisedat.values[1+noise_p['pix'],:]
    n_records = len(noise)//noise_p['record_len']
    # doing noise records
    print("  Found {0:4d} records".format(n_records))
    noise_records = np.resize(noise, (n_records, noise_p['record_len']))
    return(noise_records)

# ############################################################
# Function to read pulse records
# ############################################################
def get_pulse_records(pulse_file, pulse_p, check=False):
    '''
    This function builds pulse records from row dump files  
    
    Arguments:
        - pulse_file: string
          dump file name
        - pulse_p: dictionnary including several parameters
            - pulse_p['record_len']: number
            length of the records to be applied
            - pulse_p['prebuffer']: number
            length of prebuffer data available before each pulse
            - pulse_p['pix']: number
            index of the pixel to be processed. if 100 the routine will look for a pixel with pulses
            - pulse_p['remove_raw_files']: boolean
            indicates if the raw dta files shall be removed
        - check: boolean
          if true debugging plots are displayed (default is False)

    Returns:
        - the records
        - the t0 (indicates the begining of the pulses)
        - the pixel reference (if the input parameyter pix is equal to 100,
        the routine returns the reference of a pixel where pulses have been found)
    '''
    print("Getting pulses data from file: ", pulse_file)
    pulsedat = get_data.data(pulse_file)
    if pulse_p['remove_raw_files']:
        dirname = os.path.join(os.path.normcase(pulsedat.config['path']), pulsedat.config['dir_data'])
        fulldatafilename = os.path.join(dirname, pulse_file)
        os.remove(fulldatafilename)
    if check:
        parameters={\
        't0': 0, \
        'duration':0.1, \
        'pix_zoom':0 \
        }
        plotting_tools.plot_science_dump(pulsedat.values, pulsedat.config, parameters)

    # selecting the pixels data 
    # (offset of 1 is applyed) because the first line is the readout index
    pulses = pulsedat.values[1:,:]
    # Pulse detection
    return(trig(pulses, pulse_p, pulsedat.config))

# ############################################################
# Function to search for pixels containing pulses
# ############################################################
def search_pix_with_pulses_from_snr(data, snr_threshold=100):
    '''
    This function search in a data set pixels containing pulses.
    the function compares the SNR of the pixels. 
    !!! If there are pulses in all the pixels the function will 
    not work !!!
    
    Arguments:
        - data: numpy array
          contains the science data
        - snr_threshold: number
          indicates the snr value above which we consider that a 
          pixel contains pulses

    Returns:
        - pixels_with_pulses: an array of pixel indexes.
        contains indexes of pixels containing pulses.
    '''

    """
    Estimating signal to noise ratio
    """
    data_std=data.std(axis=1).min()
    data_pp=data.max(axis=1)-data.min(axis=1)
    data_snr=data_pp/data_std
    print(data_snr)

    """
    Recherche des pixels avec pulses
    """
    pixels_with_pulses=np.arange(len(data[0:]))[data_snr>snr_threshold]
    return(pixels_with_pulses)

def search_pix_with_pulses(data, min_threshold=0.95):
    '''
    This function search in a data set pixels containing pulses.
    the function compares the minimum of the pixel values with 
    the maximun (should be the baseline).
    
    Arguments:
        - data: numpy array
          contains the science data
        - min_threshold: number
          if the minimum value of the pixels data is below this ratio
          we consider that this pixel contains pulses.

    Returns:
        - pixels_with_pulses: an array of pixel indexes.
        contains indexes of pixels containing pulses.
    '''

    """
    Estimating signal to noise ratio
    """
    min_to_baseline_ratio=data.min(axis=1)/data.max(axis=1)

    """
    Recherche des pixels avec pulses
    """
    pixels_with_pulses=np.arange(len(data[0:]))[min_to_baseline_ratio<min_threshold]
    return(pixels_with_pulses)

# ############################################################
# Function to detect pulses in data
# ############################################################
def trig(data, pulse_p, config):
    '''
    This function builds pulse records from row dump files  
    
    Arguments:
        - data: numpy array
          contains the science data
        - pulse_p: dictionnary including several parameters
            - pulse_p['record_len']: number
            length of the records to be applied
            - pulse_p['prebuffer']: number
            length of prebuffer data available before each pulse
            - pulse_p['pix']: number
            index of the pixel to be processed. if 100 the routine will look for a pixel with pulses
        - config: dictionnary
          contains path and constants definitions

    Returns:
        - the records
        - the t0 (indicates the begining of the pulses)
        - the pixel reference (if the input parameter pulse_p['pix'] is equal to 100,
        the routine returns the reference of a pixel where pulses have been found)
    '''
    pix=pulse_p['pix']
    if pix==100: # looking for pixels with pulses
        pixels_with_pulses=search_pix_with_pulses(data)
        if len(pixels_with_pulses)>0:
            print("  There are probably pulses in the following pixels:")
            print(pixels_with_pulses)
            pix = pixels_with_pulses[0]
        else:
            print("there are probably no pulses in these data")

    """
    Rising edge detection (assuming there are no multiple pulses)
    """
    print("  Processing data of pixel {0:2d}".format(pix))
    data=data[pix,:]
    
    ratio = 0.1
    threshold = data.max() - ratio*(data.max()-data.min())
    t=np.arange(len(data))*(config['frow']/config['npix'])
    noise_limit = 100
    if data.std() < noise_limit:
        print("there are probably no pulses in these data")
        data_records = 0
    else:
        """
        Smoothing data to avoid false detection
        """
        data_smoothed = general_tools.smooth(data, 4)

        """
        Front detection
        """
        boole = data_smoothed < threshold
        boole[1:] = boole[1:] & ~boole[:-1]

        """
        Building pulse records
        """
        ipulses = np.arange(len(boole))[boole]

        """
        Checking pulses are complete
        """
        print("  Found {0:4d} pulses".format(len(ipulses)))

        spacing1 = np.append(ipulses[1:], len(data)) - ipulses
        spacing2 = ipulses - np.append(pulse_p['prebuffer']+2, ipulses[:-1])
        condition = (spacing1 > pulse_p['record_len']) & (spacing2 > pulse_p['record_len'])
        ipulses = ipulses[condition]

        n_records = len(ipulses)
        print("    > {0:4d} records are long enough for processing".format(n_records))
        data_records = np.zeros((n_records, pulse_p['record_len']))
        for i in range(n_records):
            data_records[i,:] = data[ipulses[i]-pulse_p['prebuffer']-2: ipulses[i]-pulse_p['prebuffer']-2+pulse_p['record_len']]
    return(data_records, t[ipulses], pix)

# ############################################################
# Function to create pulse average
# ############################################################
def pulse_average(pulse_list,remove_outlayers=True):
    """Creates pulse template from a set of pulses. Outlayers can be rejected.
    
    Arguments:
        - pulse_list:
        - remove_outlayers: if True outlayers are rejected
        
    Returns: pulse_template 
        - pulse_template: template created from the average of detected pulses
        - i_pulse_list: indexes of the pulse list cleaned from the outlayers
    """
    
    # Compute the average of the detected pulses and reject 1% worst if requested
    if remove_outlayers:
        mean_pulse = np.mean(pulse_list,0)
        diff_list = np.sum(abs(pulse_list-mean_pulse),1)
        pulse_list = pulse_list[(diff_list<np.percentile(diff_list,99))]
            
    return mean_pulse, (diff_list<np.percentile(diff_list,99))


# ############################################################
# Function to estimate a noise spectrum average
# ############################################################
def accumulate_noise_spectra(noise_records,abs_mean=False,rebin=1,normalize=False,dt=6.4e-6):
    '''Accumulates noise spectra from pulse free data streams
    
    Arguments:
        - noise_records: input pulse-free records
        - pulse_length: length of the records, fft to perform
        - abs_mean: option to average the abs values instead of the squares
        - rebin: rebinning factor to apply on the final spectrum
        - normalize: option to normalize the spectrum in proper A/rHz units
        - dt: sampling rate of the data (only relevant in case of normalize option)
        
    Returns: average noise spectra in a np vector of length pulse_length 
    '''
    nb_records=len(noise_records)
    pulse_length=len(noise_records[0])
    noise_spectrum_tot = np.zeros(int(pulse_length/rebin))
    nb_noise_estimates=0
    for i_record in range(nb_records):
        record=noise_records[i_record]
        if abs_mean:
            noise_spectrum_tot+=(abs(np.fft.fft(record))).reshape(-1,rebin).mean(1)
        else:
            noise_spectrum_tot+=(abs(np.fft.fft(record))**2).reshape(-1,rebin).mean(1)
        nb_noise_estimates+=1
    print("  Number of records used in noise spectrum calibration:",nb_noise_estimates)
    noise_spectrum_tot/=nb_noise_estimates
    if not abs_mean:
        noise_spectrum_tot = np.sqrt(noise_spectrum_tot)
    if normalize:
        noise_spectrum_tot*=np.sqrt(2*dt/pulse_length)
    return noise_spectrum_tot


# ############################################################
# Function to compute an optimal filter from a pulse template and a noise spectrum
# ############################################################
def compute_optimal_filter(pulse_template,noise_spectrum,energy):
    """Function to compute an optimal filter from a pulse template and a noise spectrum.
    Optimal filter is normalized to return energy when dot producted with the template
    
    Arguments:
        - pulse_template: pulse template to use in optimal filter generation
        - noise_spectrum: noise spectrum to use in optimal filter generation
        - energy: energy of the pulse template 
    """
    pulse_spectrum = np.fft.fft(pulse_template)
    time_optimal_filter = pulse_spectrum/noise_spectrum**2
    time_optimal_filter[0]=0
    cutted_time_filter=np.real(np.fft.ifft(time_optimal_filter))
    cutted_time_filter/=np.dot(pulse_template,cutted_time_filter)/energy
    
    return cutted_time_filter


# ############################################################
# Function to perform the "jitter parabola fit"
# ############################################################
def do_pulse_jitter(opt_filter, pulse_record):
    '''Performs the +-1 pulse jitter parabola technique and returns both the maximum and its phase and
    the potentially corrected pulse time
    
    Arguments:
        - opt_filter: optimal filter
        - pulse_record: data_stream to analyze
    '''
    phase_offset = 1 
    pulse_length = len(opt_filter)
    while True:
        if phase_offset<0 or phase_offset+pulse_length>len(pulse_record):
            print('Problem to find the phase of the pulse!')
            energy=-1
            break
        energy1 = np.dot(opt_filter,pulse_record[phase_offset:phase_offset+pulse_length])
        energy2 = np.dot(opt_filter,pulse_record[phase_offset+1:phase_offset+pulse_length+1])
        energy3 = np.dot(opt_filter,pulse_record[phase_offset+2:phase_offset+pulse_length+2])
        A = 0.5*(energy1-2*energy2+energy3)
        B = 0.5*(energy3-energy1)
        C = energy2
        energy = C-.25*B**2/A
        pulse_phase = -.5*B/A

        if energy1>energy2:
            phase_offset-=1
        elif energy3>energy2:
            phase_offset+=1
        else:
            pulse_phase = -.5*B/A
            break

    return energy,pulse_phase+phase_offset


# ############################################################
# Gaussian function
# ############################################################
def gauss(x,*p):
    '''
    Gaussian function for fits
    '''
    A,mu,sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


# ############################################################
# Generic function to fit a histogram
# ############################################################
def hist_and_fit(array_to_fit,bins,fit_function=gauss,c=(71./255,144./255,195./255)):
    '''Plots an histogram and fit it with a gaussian
    
    Arguments:
        - array_to_fit: array containing the data to study
        - bins: bins parameter of np.histogram
        - xlabel: label of x axis for plot
        - fit_function: callable for the fit function
        - c: RGB face color of the histogram
    '''
    hist,bins = np.histogram(array_to_fit,bins=bins)
    bin_centers = (bins[:-1] + bins[1:])/2
    p0=[len(array_to_fit)/(np.sqrt(2*np.pi)*array_to_fit.std()),array_to_fit.mean(),array_to_fit.std()]
    coeff = curve_fit(fit_function, bin_centers, hist, p0=p0,maxfev=10000000)[0]

    axe_fit = np.arange(bins[0],bins[-1],0.01*(bins[-1]-bins[0]))
    hist_fit = fit_function(axe_fit, *coeff)

    print("Fitted function with parameters:",coeff)
    if fit_function==gauss:
        print("  - mu: {0:5.3f}".format(coeff[1]))
        print("  - sigma: {0:5.3f} (FWHM = {1:5.3f})".format(coeff[2],2.355*coeff[2]))
        print("  - A: {0:5.3f}".format(coeff[0]))
        
    return coeff, array_to_fit, bins, axe_fit, hist_fit


# ############################################################
# Function for phase correction
# ############################################################
def apply_phase_correction(energies,phase,phase_correction):
    '''Applies a phase correction (should be merged with baseline correction)
    
    Arguments:
        - energies: input energies
        - phase: phase data
        - phase_correction: order of the phase correction polynom
        - fig: figure id to include the phase correction plots
        - i_subplot: index of the sub_plot in the figure
    '''
    phase_correc_poly_coeff = np.polyfit(phase,energies,phase_correction)
    phase_correc_poly = np.poly1d(phase_correc_poly_coeff)
    corrected_energies=energies.copy()-phase_correc_poly(phase)+phase_correc_poly(phase.mean())
   
    return corrected_energies, phase_correc_poly

# ############################################################
# Function for drift correction
# ############################################################
def apply_drift_correction(energies,time,drift_correction):
    '''Applies a phase correction
    
    Arguments:
        - energies: input energies
        - time: time data
        - drift_correction: order of the drift correction polynom
    '''
    #Fit drift correction
    drift_correc_poly_coeff = np.polyfit(time,energies,drift_correction)
    drift_correc_poly = np.poly1d(drift_correc_poly_coeff)
    energies_corrected=energies.copy()-drift_correc_poly(time)+energies.mean()
    
    return energies_corrected, drift_correc_poly

# ############################################################
# Function for baseline correction
# ############################################################
def apply_baseline_correction(energies,baseline,baseline_correction):
    '''Applies a phase correction
    
    Arguments:
        - energies: input energies
        - baseline: baseline data
        - baseline_correction: order of the baseline correction polynom
    '''
    #Fit baseline correction
    baseline_correc_poly_coeff = np.polyfit(baseline,energies,baseline_correction)
    baseline_correc_poly = np.poly1d(baseline_correc_poly_coeff)
    energies_corrected=energies.copy()-baseline_correc_poly(baseline)+baseline_correc_poly(baseline.mean())
    
    return energies_corrected, baseline_correc_poly

# ############################################################
# Function to perform energy reconstruction
# ############################################################
def energy_reconstruction(pulse_list,optimal_filter,prebuffer,prebuff_exclusion=10,
                          verbose=False):
    """Perform energy reconstruction on list of pulses, including jitter correction.
    Also compute baseline value before each pulse
    
    Arguments:
        - pulse_list: list of pulses to reconstruct
        - optimal_filter: filter to use for reconstruction
        - prebuffer: length of prebuffer data available before each pulse
        - prebuff_exclusion: number of points to remove from the buffer in baseline estimation 
        
    Returns: pulse_template 
        - 
    """
    # Iterate over available pulses
    energies = []
    phases = []
    baselines = []
    for i_pulse in range (len(pulse_list)):
        # Perform reconstruction with jitter
        #e,p = do_pulse_jitter(optimal_filter, pulse_list[i_pulse])
        e,p = do_pulse_jitter(optimal_filter, pulse_list[i_pulse])
        if e!=-1:  # an optimal phase has been found
            energies.append(e)
            phases.append(p)
            # Estimate baseline
            baselines.append(pulse_list[i_pulse][:prebuffer-prebuff_exclusion].mean())
             
    energies = np.array(energies)
    phases = np.array(phases)
    baselines = np.array(baselines)
    
    # Return energies, phases and times
    return energies, phases, baselines


# ############################################################
# Function to perform all the operations needed to compute the optimal filters
# ############################################################
def do_ep_filter(noise, pulses, file_xifusim_template, file_xifusim_tes_noise, prebuffer, config, plotdirname, verbose=False, do_plots=True):
    """Perform the operations to compute the optimal filters (with and without tes noise)
    and compares the results with xifusim.
    
    Arguments:
        - noise: array containing noise data of one pixel
        - pulses: array containing DRE pulses data
        - file_xifusim_template: file containing xifusim template
        - file_xifusim_tes_noise: file containing tes noise compute by xifusim
        - prebuffer: length of prebuffer data available before each pulse 
        - config: contains path and constants definitions
        - plotdirname: location of plotfiles
        - verbose: if True some informations are printed (Default=False)
        - do_plots: if True a plot is done with all the intermediate results (Default=True)
        
    Returns: optimal_filter (no TES noise), optimal_filter_tot (with TES noise)       
    """
 
    # ############################################################
    # Noise spectrum calibration
    # ############################################################
    print("\nPerforming noise spectrum calibration...")

    # Load pulse-free data to generate noise spectrum
    record_length=len(noise[0])
    print("  Record length = {0:4d}".format(record_length))
        
    # Compute average noise spectrum
    noise_spectrum = accumulate_noise_spectra(noise, normalize=True)
    #np.save('noise_spectrum.npy', noise_spectrum)
    frequencies = np.fft.fftfreq(record_length,6.4e-6)[1:int(record_length/2)]
    if verbose:
        print("  Noise spectrum:",noise_spectrum)

    # ############################################################
    # Pulse template calibration
    # ############################################################
    print("\nPerforming pulse template calibration...")

    # Generate pulse template as average of detected pulses (and removing outlayers)
    pulse_template, i_pulse_list = pulse_average(pulses,remove_outlayers=True)
    pulses=pulses[i_pulse_list]
    if verbose:
        print("  Pulse template: ", pulse_template)

    # ############################################################
    # Comparison with xifusim data
    # ############################################################
    print("\nComparing xifusim and dre data...")
    
    # Load xifusim template
    _ , xifusim_template = np.load(file_xifusim_template)
    # Re-sample xifusim template
    xifusim_template_decimated = np.convolve(xifusim_template,1/128.*np.ones(128),mode="valid")[79::128][:1000] # Fitted best phase by hand
    
    # Determine scaling factor between xifusim and DRE units with the baseline value
    dre_template = pulse_template
    baseline_scaling = dre_template[0]/xifusim_template[0]
    ph_scaling = (dre_template[0]-dre_template.min())/(xifusim_template[0]-xifusim_template.min())
    print("  Difference between average baseline and PH scaling: {0:5.3f}%".format((ph_scaling/baseline_scaling-1)*100))
                    
    # Computing power spectra
    xifusim_ps = abs(np.fft.rfft(xifusim_template_decimated*baseline_scaling))**2
    decal=172  # Shift of DRE template to match xifusim template 
    DRE_PS = abs(np.fft.rfft(dre_template[decal:1000+decal]))**2
    ps_freq = np.fft.fftfreq(1000,6.4e-6)[:501]

    # ############################################################
    # Addition of TES noise ontop DRE noise
    # ############################################################
    print("\nComparing TES noise and DRE noise...")
    # Load TES noise spectrum
    tes_noise = fits.getdata(file_xifusim_tes_noise,"SPEC{0:04d}".format(record_length))["CSD"]*baseline_scaling
    total_noise = np.sqrt(tes_noise**2+noise_spectrum**2)

    # ############################################################
    # Compute optimal filters
    # ############################################################
    print('\nComputing optimal filter...')
    optimal_filter = compute_optimal_filter(pulse_template,noise_spectrum,energy=7.)
    #np.save('optimal_filter.npy', optimal_filter)
    print("Check optimal filter normalization: {0:5.3f}keV".format(np.dot(optimal_filter,pulse_template)))

    # Compute optimal filter, now including TES noise
    print('\nComputing optimal filter including TES noise...')
    optimal_filter_tot = compute_optimal_filter(pulse_template,total_noise,7.)
    print("Check optimal filter normalization: {0:5.3f}keV".format(np.dot(optimal_filter_tot,pulse_template)))

    # ############################################################
    # Do plots
    # ############################################################
    if do_plots:
        fig = plt.figure(figsize=(14, 9))    

        # Show noise spectrum
        ax1=fig.add_subplot(2,3,2)
        ax1.loglog(frequencies,noise_spectrum[1:int(record_length/2)])
        ax1.set_title('Average noise spectrum')
        ax1.set_xlabel('Frequency [Hz]')
        ax1.set_ylabel("ADU/rHz")
        for item in ([ax1.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax1.xaxis.label, ax1.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            item.set_fontsize(6)
                
        """
        # Show pulse template
        ax2=fig.add_subplot(2,4,2)
        #ax2.plot((np.arange(len(pulse_template))-prebuffer)*6.4e-3,pulse_template, label='Pulse template')
        #ax2.plot((np.arange(len(pulse_template[:prebuffer]))-200)*6.4e-3,pulse_template[:prebuffer],'r',label='Pre-buffer')
        ax2.plot((np.arange(len(pulse_template))-prebuffer)*config['npix']*1e-3/config['frow'],pulse_template, label='Pulse template')
        ax2.plot((np.arange(len(pulse_template[:prebuffer]))-prebuffer)*config['npix']*1e-3/config['frow'],pulse_template[:prebuffer],'r',label='Pre-buffer')
        ax2.set_title('Pulse template')
        ax2.set_xlabel('Time [ms]')
        ax2.set_ylabel("ADU")
        ax2.set_ylim(0, 2**16)
        ax2.legend(loc='best', prop=dict(size=7))
        for item in ([ax2.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax2.xaxis.label, ax2.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(6)
        """

        # Show difference between xifusim template and DRE template
        ax3=fig.add_subplot(2,3,1)
        ax3.plot(xifusim_template_decimated*baseline_scaling,label='xifusim')
        ax3.plot(dre_template[decal:1000+decal],label='dre')
        ax3.plot(xifusim_template_decimated*baseline_scaling - dre_template[decal:1000+decal],label='difference')
        ax3.set_title('Pulse template')
        ax3.set_xlabel(sample_nb_label)
        ax3.set_ylabel('Module (ADU)')
        ax3.legend(loc="best", prop=dict(size=7))
        for item in ([ax3.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax3.xaxis.label, ax3.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax3.get_xticklabels() + ax3.get_yticklabels()):
            item.set_fontsize(6)

        # Relative difference
        ax4=fig.add_subplot(2,3,4)
        ax4.plot((xifusim_template_decimated*baseline_scaling - dre_template[decal:1000+decal])/dre_template[decal:1000+decal]*100,label='difference')
        ax4.set_title('Relative difference between both templates')
        ax4.set_xlabel(sample_nb_label)
        ax4.set_ylabel('Relative difference [%]')
        for item in ([ax4.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax4.xaxis.label, ax4.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax4.get_xticklabels() + ax4.get_yticklabels()):
            item.set_fontsize(6)

        # Show different power spectra
        ax5=fig.add_subplot(2,3,5)
        ax5.loglog(ps_freq,xifusim_ps,label="xifusim")
        ax5.loglog(ps_freq,DRE_PS,label="dre")
        ax5.set_title('Comparison of Pulses power spectra')
        ax5.set_xlabel("Frequency [Hz]")
        ax5.set_ylabel("PSD [AU]")
        ax5.legend(loc="best", prop=dict(size=7))
        for item in ([ax5.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax5.xaxis.label, ax5.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax5.get_xticklabels() + ax5.get_yticklabels()):
            item.set_fontsize(6)

        # Compare DAC noise with TES noise
        ax6=fig.add_subplot(2,3,3)
        ax6.loglog(frequencies,tes_noise[1:int(record_length/2)],label="TES noise")
        ax6.loglog(frequencies,noise_spectrum[1:int(record_length/2)],label="DAC noise")
        ax6.loglog(frequencies,total_noise[1:int(record_length/2)],label="Total noise")
        ax6.set_title("DAC noise vs TES noise")
        ax6.set_xlabel("Frequency [Hz]")
        ax6.set_ylabel("PSD [AU]")
        ax6.legend(loc="best", prop=dict(size=7))
        for item in ([ax6.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax6.xaxis.label, ax6.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax6.get_xticklabels() + ax6.get_yticklabels()):
            item.set_fontsize(6)

        """
        # Show optimal filter without TES noise
        ax7=fig.add_subplot(2,4,7)
        ax7.plot(optimal_filter_tot)
        ax7.set_title('Optimal filter')
        ax7.set_xlabel(sample_nb_label)
        for item in ([ax7.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax7.xaxis.label, ax7.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax7.get_xticklabels() + ax7.get_yticklabels()):
            item.set_fontsize(6)
        """

        # Show optimal filter with TES noise
        ax8=fig.add_subplot(2,3,6)
        ax8.plot(optimal_filter_tot)
        ax8.set_title('Optimal filter including TES noise')
        ax8.set_xlabel(sample_nb_label)
        for item in ([ax8.title]):
            item.set_weight('bold')
            item.set_fontsize(8)
        for item in ([ax8.xaxis.label, ax8.yaxis.label]):
            item.set_fontsize(7)
        for item in (ax8.get_xticklabels() + ax8.get_yticklabels()):
            item.set_fontsize(6)
        #plt.show()
        fig.tight_layout()
        plt.savefig(os.path.join(plotdirname,'PLOT_E_RESOL_TEMPLATES.png'),bbox_inches='tight')

    return(optimal_filter, optimal_filter_tot)


# ############################################################
# Plot of energy resolution measurements
# ############################################################
def plot_er(non_linear_factor,array_to_fit1,bins1,coeffs1,axe_fit1,hist_fit1,baselines,energies,bl_correct_poly1,energies_c_bl, \
        array_to_fit2,bins2,coeffs2,axe_fit2,hist_fit2,phases,ph_correct_poly1,energies_c_ph,\
        array_to_fit3,bins3,coeffs3,axe_fit3,hist_fit3,c,tes_text,plotfilename):
    fig = plt.figure(figsize=(8, 14))

    xlabel_baseline = 'Baseline'
    xlabel_e = 'Energy [eV]'
    ylabel_e7 = 'Energy - 7000 [eV]'
    ylabel_occ = "Occurences"
    loc_ul = "upper left"
    fit_er_txt = 'Fit : Er={0:5.3f} x {1:5.3f} = {2:5.3f}+-{3:5.3f}eV'

    # Show initial energy error
    ax1=fig.add_subplot(5, 1, 1)
    ax1.hist(array_to_fit1,bins=bins1,histtype='stepfilled',facecolor=c)
    ax1.plot(axe_fit1,hist_fit1,c='r',linewidth=2, label='Fit: Er={0:5.3f} x {1:5.3f} = {2:5.3f}+-{3:5.3f}eV' \
        .format(non_linear_factor, 2.355*coeffs1[2], 2.355*coeffs1[2]*non_linear_factor, 2.355*coeffs1[2]*non_linear_factor/(np.sqrt(2.*len(energies)))))
    ax1.legend(loc=loc_ul, prop=dict(size=7))
    ax1.set_title('Initial energy resolution '+tes_text)
    ax1.set_xlabel(xlabel_e)
    ax1.set_ylabel(ylabel_occ)
    for item in ([ax1.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax1.xaxis.label, ax1.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(6)

    ax2=fig.add_subplot(5,2,3)
    ax2.plot(baselines,(energies-7)*1000,marker='.',linestyle='', c=c)
    ax2.plot(np.sort(baselines),(bl_correct_poly1(np.sort(baselines))-7)*1000,c='r',linewidth=2, label='Fit')
    ax2.legend(loc=loc_ul, prop=dict(size=7))
    ax2.set_title('Before baseline correction')
    ax2.set_xlabel(xlabel_baseline)
    ax2.set_ylabel(ylabel_e7)
    for item in ([ax2.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax2.xaxis.label, ax2.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(6)
    
    #Correct for correlation
    ax3=fig.add_subplot(5,2,4)
    ax3.plot(baselines,(energies_c_bl-7)*1000,marker='.',linestyle='', c=c)
    ax3.set_title('After baseline correction')
    ax3.set_xlabel(xlabel_baseline)
    ax3.set_ylabel(ylabel_e7)
    for item in ([ax3.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax3.xaxis.label, ax3.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(6)

    # Show energy error after baseline correction
    ax4=fig.add_subplot(5, 1, 3)
    ax4.hist(array_to_fit2,bins=bins2,histtype='stepfilled',facecolor=c)
    ax4.plot(axe_fit2,hist_fit2,c='r',linewidth=2, label=fit_er_txt\
        .format(non_linear_factor, 2.355*coeffs2[2], 2.355*coeffs2[2]*non_linear_factor, 2.355*coeffs2[2]*non_linear_factor/(np.sqrt(2.*len(energies)))))
    ax4.legend(loc=loc_ul, prop=dict(size=7))
    ax4.set_title('Energy resolution after baseline correction '+tes_text)
    ax4.set_xlabel(xlabel_e)
    ax4.set_ylabel(ylabel_occ)
    for item in ([ax4.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax4.xaxis.label, ax4.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax4.get_xticklabels() + ax4.get_yticklabels()):
        item.set_fontsize(6)

    ax5=fig.add_subplot(5,2,7)
    ax5.plot(phases,(energies_c_bl-7)*1000,marker='.',linestyle='',c=c)
    ax5.plot(np.sort(phases),(ph_correct_poly1(np.sort(phases))-7)*1000,c='r',linewidth=2, label='Fit')
    ax5.legend(loc=loc_ul, prop=dict(size=7))
    ax5.set_title('Before phase correction')
    ax5.set_xlabel('Phase (samples)')
    ax5.set_ylabel(ylabel_e7)
    for item in ([ax5.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax5.xaxis.label, ax5.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax5.get_xticklabels() + ax5.get_yticklabels()):
        item.set_fontsize(6)
    
    ax6=fig.add_subplot(5,2,8)
    ax6.plot(phases,(energies_c_ph-7)*1000,marker='.',linestyle='', c=c)
    ax6.set_title('After phase correction')
    ax6.set_xlabel('Phase (samples)')
    ax6.set_ylabel(ylabel_e7)
    for item in ([ax6.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax6.xaxis.label, ax6.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax6.get_xticklabels() + ax6.get_yticklabels()):
        item.set_fontsize(6)

    # Show energy error after phase correction
    ax7=fig.add_subplot(5, 1, 5)
    ax7.hist(array_to_fit3,bins=bins3,histtype='stepfilled',facecolor=c, label="")
    ax7.plot(axe_fit3,hist_fit3,c='r',linewidth=2, label=fit_er_txt\
        .format(non_linear_factor, 2.355*coeffs3[2], 2.355*coeffs3[2]*non_linear_factor, 2.355*coeffs3[2]*non_linear_factor/(np.sqrt(2.*len(energies)))))
    ax7.legend(loc=loc_ul, prop=dict(size=7))
    ax7.set_title('Energy resolution after baseline and phase corrections '+tes_text)
    ax7.set_xlabel(xlabel_e)
    ax7.set_ylabel(ylabel_occ)
    for item in ([ax7.title]):
        item.set_weight('bold')
        item.set_fontsize(8)
    for item in ([ax7.xaxis.label, ax7.yaxis.label]):
        item.set_fontsize(7)
    for item in (ax7.get_xticklabels() + ax7.get_yticklabels()):
        item.set_fontsize(6)
    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename,bbox_inches='tight')

    fig = plt.figure(figsize=(10, 8))
    ax7=fig.add_subplot(1, 1, 1)
    # Show energy error after phase correction
    ax7.hist(array_to_fit3,bins=bins3,histtype='stepfilled',facecolor=c, label="")
    ax7.plot(axe_fit3,hist_fit3,c='r',linewidth=2, label=fit_er_txt\
        .format(non_linear_factor, 2.355*coeffs3[2], 2.355*coeffs3[2]*non_linear_factor, 2.355*coeffs3[2]*non_linear_factor/(np.sqrt(2.*len(energies)))))
    ax7.legend(loc=loc_ul, fontsize=12)
    ax7.set_title('Energy resolution after baseline and phase corrections '+tes_text)
    ax7.set_xlabel(xlabel_e)
    ax7.set_ylabel(ylabel_occ)
    for item in ([ax7.title, ax7.xaxis.label, ax7.yaxis.label]):
        item.set_weight('bold')
        item.set_fontsize(14)
    for item in (ax7.get_xticklabels() + ax7.get_yticklabels()):
        item.set_fontsize(12)
    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename[:-4]+'_zoom.png',bbox_inches='tight')


# ############################################################
# Pulse reconstruction
# ############################################################
def measure_er(pulse_list, t_list, optimal_filter, optimal_filter_tot, pixeldirname, plotdirname, index, prebuffer, verbose=False, do_plots=True):
    """Perform the operations to measure the energy resolution (with and without tes noise).
    
    Arguments:
        - pulse_lis: Demux pulse records
        - optimal_filter: optimal filter computed without the TES noise
        - optimal_filter_tot: optimal filter computed with the TES noise
        - pixeldirname: directory containing pixel's informations
        - plotdirname: location of plotfiles
        - index: to be included in the plot filename 
        - prebuffer: size of pre buffer
        - verbose: if True some informations are printed (Default=False)
        - do_plots: if True a plot is done with all the intermediate results (Default=True)
        
    Returns: Nothing       
    """
    bcorr=3 # Order of polynomial baseline correction
    pcorr=8 # Order of polynomial arrival phase correction
    dcorr=4 # Order of polynomial drift correction
    print('\nReconstructing pulses...')

    # Load pixel non-linearity factor from an information file
    #non_linear_factor=get_non_linear_factor(pixeldirname, verbose=verbose)
    non_linear_factor=1.361
    print("  Loading pixel non linearity factor at 7keV: ", non_linear_factor)

    print("  Record length = {0:4d}".format(len(pulse_list[0])))

    # Removing outlayers
    _, i_pulse_list = pulse_average(pulse_list,remove_outlayers=True)
    pulse_list=pulse_list[i_pulse_list]
    t_list=t_list[i_pulse_list]

    # Compute first energy estimates
    energies,phases,baselines = energy_reconstruction(pulse_list, optimal_filter, prebuffer, JITTER_MARGIN)

    coeffs1, array_to_fit1, bins1, axe_fit1, hist_fit1 = hist_and_fit((energies-7)*1000, 100)
    print("Resolution without TES noise (prior correction): {0:5.3f}eV".format(coeffs1[2]*2.355*non_linear_factor))

    # Apply baseline correction
    energies_c_bl, bl_correct_poly1 = apply_baseline_correction(energies, baselines, bcorr)
    coeffs2, array_to_fit2, bins2, axe_fit2, hist_fit2  = hist_and_fit((energies_c_bl-7)*1000, 100)
    print("Resolution without TES noise (after baseline correction): {0:5.3f}eV".format(coeffs2[2]*2.355*non_linear_factor))

    # Apply phase correction
    energies_c_ph, ph_correct_poly1 = apply_phase_correction(energies_c_bl, phases, pcorr)
    coeffs3, array_to_fit3, bins3, axe_fit3, hist_fit3  = hist_and_fit((energies_c_ph-7)*1000, 100)

    eres_notesnoise = coeffs3[2]*2.355*non_linear_factor
    print("Final resolution without TES noise (after phase correction): {0:5.3f}+-{1:5.3f}eV".format(eres_notesnoise,eres_notesnoise/(np.sqrt(2.*len(energies)))))

    if do_plots:
        plotfilename=os.path.join(plotdirname,'PLOT_E_RESOL_NO_TES_NOISE_{0:d}.png'.format(index))
        plot_er(non_linear_factor,array_to_fit1,bins1,coeffs1,axe_fit1,hist_fit1,baselines,energies,bl_correct_poly1,energies_c_bl, \
            array_to_fit2,bins2,coeffs2,axe_fit2,hist_fit2,phases,ph_correct_poly1,energies_c_ph,\
            array_to_fit3,bins3,coeffs3,axe_fit3,hist_fit3,'g','(no TES noise,{0:6d} counts)'.format(len(energies)),plotfilename)

    # ############################################################
    # Pulse reconstruction with OF also containing TES noise
    # ############################################################
    print("\nReconstruction with OF including TES noise...")
        
    # Compute first energy estimates
    energies,phases,baselines = energy_reconstruction(pulse_list, optimal_filter_tot, prebuffer, JITTER_MARGIN)

    coeffs1, array_to_fit1, bins1, axe_fit1, hist_fit1 = hist_and_fit((energies-7)*1000, 100)
    print("Resolution without TES noise (prior correction): {0:5.3f}eV".format(coeffs1[2]*2.355*non_linear_factor))

    # Apply baseline correction
    energies_c_bl, bl_correct_poly1 = apply_baseline_correction(energies, baselines, bcorr)
    coeffs2, array_to_fit2, bins2, axe_fit2, hist_fit2  = hist_and_fit((energies_c_bl-7)*1000, 100)
    print("Resolution without TES noise (after baseline correction): {0:5.3f}eV".format(coeffs2[2]*2.355*non_linear_factor))

    # Apply phase correction
    energies_c_ph, ph_correct_poly1 = apply_phase_correction(energies_c_bl, phases, pcorr)
    coeffs3, array_to_fit3, bins3, axe_fit3, hist_fit3  = hist_and_fit((energies_c_ph-7)*1000, 100)

    # Apply drift correction
    energies_c_dr, dr_correct_poly1 = apply_drift_correction(energies_c_ph, t_list, dcorr)
    coeffs4, _, _, _, _  = hist_and_fit((energies_c_dr-7)*1000, 100)
    print('with drift correction:', coeffs4[2]*2.355*non_linear_factor)

    eres_tesnoise = coeffs3[2]*2.355*non_linear_factor
    print("Final resolution with TES noise (after phase correction): {0:5.3f}+-{1:5.3f}eV".format(eres_tesnoise,eres_tesnoise/(np.sqrt(2.*len(energies)))))
        
    if do_plots:
        plotfilename=os.path.join(plotdirname,'PLOT_E_RESOL_WITH_TES_NOISE_{0:d}.png'.format(index))
        plot_er(non_linear_factor,array_to_fit1,bins1,coeffs1,axe_fit1,hist_fit1,baselines,energies,bl_correct_poly1,energies_c_bl, \
            array_to_fit2,bins2,coeffs2,axe_fit2,hist_fit2,phases,ph_correct_poly1,energies_c_ph,\
            array_to_fit3,bins3,coeffs3,axe_fit3,hist_fit3,'b','(with TES noise,{0:6d} counts)'.format(len(energies)),plotfilename)

        plotdriftfilename=os.path.join(plotdirname,'er_drift_{0:d}.png'.format(index))
        fig = plt.figure(figsize=(8, 5))
        ax1=fig.add_subplot(1, 2, 1)
        ax1.plot(t_list, (energies_c_ph-7)*1000,'.')
        ax1.plot(t_list, (dr_correct_poly1(t_list)-7)*1000, 'r')
        ax1.set_title('Residual drift along time')
        ax1.set_ylabel('Energy delta (eV)')
        ax1.set_xlabel('time (A.U.)')
        ax2=fig.add_subplot(1, 2, 2)
        ax2.plot(t_list, (energies_c_dr-7)*1000,'.')
        ax2.set_title('After drift correction')
        ax2.set_ylabel('Energy delta (eV)')
        ax2.set_xlabel('time (A.U.)')
        fig.tight_layout()
        plt.savefig(plotdriftfilename,bbox_inches='tight')

    # np.save('energies.npy', energies)
    return(eres_tesnoise, eres_tesnoise/(np.sqrt(2.*len(energies))))


# ############################################################
def ep(noise_rec_filename, calib_rec_filename, measu_rec_filename, config, prebuffer, verbose=False):
    """Perform the operations to measure the energy resolution (with and without tes noise).
    
    Arguments:
        noise_rec_filename: string
        name of the npy file containing the noise records

        calib_rec_filename: string
        name of the npy file containing the pulse records of calibration data

        measu_rec_filename: string
        name of the npy file containing the pulse records of measurement data

        config: dictionnary
        Contains path and constants definitions

        prebuffer: number
        size of pre buffer

        verbose: boolean
        If True some informations are printed (Default=False)

    Returns:
        eres_ok: boolean
        True if DRE energy resolution contribution is below the requirement       
    """

    """
    Defining path and file names
    """
    dirname = os.path.join(os.path.normcase(config['path']), config['dir_data'])
    full_noise_rec_filename = os.path.join(dirname, noise_rec_filename)
    full_calib_rec_filename = os.path.join(dirname, calib_rec_filename)
    full_measu_rec_filename = os.path.join(dirname, measu_rec_filename)
    plotdirname = os.path.join(os.path.normcase(config['path']), config['dir_plots'])
    general_tools.checkdir(plotdirname)
    pixeldirname = os.path.normcase("./Pixel_data_LPA75um_AR0.5/")
    file_xifusim_template = os.path.join(pixeldirname,"pulse_withBBFB.npy")
    file_xifusim_tes_noise = os.path.join(pixeldirname,"noise_spectra_bbfb_noFBDAC.fits")

    with open(full_noise_rec_filename, 'rb') as file:
        noise = np.load(file)
    with open(full_calib_rec_filename, 'rb') as file:
        _ = np.load(file)
        calib = np.load(file)
    with open(full_measu_rec_filename, 'rb') as file:
        t_measu = np.load(file)
        measu = np.load(file)

    #plotting_tools.over_plot_records(calib)

    """
    Computing EP filter
    """
    optimal_filter, optimal_filter_tot=do_ep_filter(noise, calib, file_xifusim_template, file_xifusim_tes_noise, prebuffer, config, plotdirname, verbose)
    
    """
    Measuring energies
    """
    eres, _=measure_er(measu, t_measu, optimal_filter, optimal_filter_tot, pixeldirname, plotdirname, 0, prebuffer, verbose)

    return(eres, config['eres_req_cbe_dre_7kev'])

# ############################################################

