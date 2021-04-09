#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from numpy.fft import rfft
import os
import matplotlib.pyplot as plt

import general_tools, get_data

"""
Some general text definitions
"""
plotting_message="Plotting dump data..."
time_label="Time (ms)"
frequency_label="Frequency (Hz)"
spectral_power_density_label=r"Power spectral density (dB/$\sqrt{Hz}$)"
error_t0_message="Wrong zoom parameter, t0 changed to 0s."

# -----------------------------------------------------------------------
def mosaic_labels(ax, box, n_cols, n_lines, x_lab, y_lab):
    r"""
        This function defines the xlabel and ylabel for a plot in a plot mosaic.
        The xlabel is displayed only for the plots at the bottom of the mosaic.
        The ylabel is displayed only for the plots at the left of the mosaic.

        Parameters
        ----------
        ax : matplotlib axis

        box : number
        the number of the plot in the mosaic

        n_cols, n_lines : numbers
        the number of columns and lines in the mosaic

        x_lab, y_lab : string
        the x and y labels

        Returns
        -------
        Nothing

        """
    if box//n_cols == n_lines-1:
        ax.set_xlabel(x_lab)
    else:
        plt.xticks(visible=False)
    if box%n_cols == 0:
        ax.set_ylabel(y_lab)
    else:
        plt.yticks(visible=False)


# -----------------------------------------------------------------------------
def plot_science_dump(data, plotfilename, config, t0=0, duration=0, pix_zoom=0, noise=False, check_noise_measurement=False):
    r"""
        This function checks the data of a DRE-DEMUX science data dump.

        Parameters
        ----------
        data: numpy array
        The data of the dump 

        plotfilename: string
        Name of the plot file (with the path)

        config: dictionnary
        contains different informations such as path and directory names

        t0: number, optional
        begining of zoom in seconds (default is 0)

        duration: number, optional
        zoom duration in seconds. If 0 all the data are plotted (default is 0)

        pix_zoom: number (integer)
        pixel id refering to the pixel for which we plot a zoom (default=0)

        noise: boolean
        Indicates if a noise analysis shall be done (default=False)

        check_noise_measurement: boolean
        Indicates if a noise analysis shall be done on fake data (default=False)        

        Returns
        -------
        Nothing

        """
    print(plotting_message)
    plotfilename_zoom = plotfilename[:-4]+"_zoom.png"
    plotfilename_all = plotfilename[:-4]+"_all.png"
    print("  >> " + plotfilename_zoom)
    print("  >> " + plotfilename_all)

    npix = int(config["npix"])
    fs = float(config["frow"]/npix)
    t = np.arange(len(data[0,:]))/fs

    packet_nb = data[0,:]
    reference = packet_nb[0] + np.arange(len(packet_nb))
    if np.abs((packet_nb-reference)%2**32).max()>0:
        print(" Error in the packet number !")
        plt.plot(packet_nb-reference)
        plt.show()

    data = data[1:,:]

    """
    Plotting data in time domain
    """
    imin = int(t0*fs)
    if imin > len(data[0,:]):
        print(error_t0_message)
        imin = 0
    if duration == 0 :
        imax = len(data[0,:])
    else :
        imax = min(int((t0+duration)*fs), len(data[0,:]))
    mkr=''
    if imax - imin < 200:
        mkr='.'

    fig = plt.figure(figsize=(18, 12))
    xtitle = time_label

    ncols, nlines = 9, 4
    ymin, ymax = 0.98*data.min(), 1.02*data.max()
    for pix in range(npix):
        ax = fig.add_subplot(nlines, ncols, 1+pix)
        ax.plot(1000*t[imin:imax], data[pix,imin:imax], 'b')
        ax.set_title("Pixel {0:2d}".format(1+pix))
        ax.grid(color='k', linestyle=':', linewidth=0.5)
        ax.set_xlim(t[imin]*1000, t[imax-1]*1000)
        ax.set_ylim(ymin, ymax)
        mosaic_labels(ax, pix, ncols, nlines, xtitle, r'ADU')
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
            item.set_weight('bold')
            item.set_fontsize(12)
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(10)
    ax = fig.add_subplot(nlines, ncols, ncols*nlines)
    ax.plot(1000*t[imin:imax], packet_nb[imin:imax], 'b')
    ax.set_title("Packet number")
    ax.set_xlabel(xtitle)
    ax.grid(color='k', linestyle=':', linewidth=0.5)

    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename_all, bbox_inches='tight')

    """
    doing zoom
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(1000*t[imin:imax], data[pix_zoom,imin:imax], 'b', marker=mkr)
    ax.set_title("Pixel {0:2d}".format(1+pix_zoom))
    ax.set_xlabel(xtitle)
    ax.set_ylabel(r'ADU')
    ax.set_xlim(t[imin]*1000, t[imax-1]*1000)
    ax.grid(color='k', linestyle=':', linewidth=0.5)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_weight('bold')
        item.set_fontsize(14)
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)

    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename_zoom, bbox_inches='tight')

    """
    Plotting spectra
    """
    if noise:
        plot_science_dump_noise(data, config, plotfilename[:-4]+"_noise.png", pix_zoom, 2000, 8192)
        if check_noise_measurement:
            plotfilename_test = plotfilename[:-4]+"_FAKE.png"
            data_test = make_fake_science_noise(config, 2000, 8192)
            plot_science_dump_noise(data_test, config, plotfilename_test, 0, 2000, 8192)

# -----------------------------------------------------------------------------
def plot_science_dump_noise(data, config, plotfilename, pix_zoom, n_records=2000, record_len=8192):
    r"""
        This function measures noise spectra on science data

        Parameters
        ----------
        data: 2-dimensional numpy array
        The input science data

        config: dictionnary
        contains different informations such as path and directory names

        plotfilename: string
        Name of the plot file (with the path)

        pix_zoom: number (integer)
        pixel id refering to the pixel for which we plot a zoom (default=0)

        n_records: number (integer)
        number of spectra to be averaged (default = 2000)

        record_len: number (integer)
        number of samples per records (default = 8192)

        Returns
        -------
        noise_spectra: 2-dimensional numpy array
        The noise spectra
        """
    print("Plotting noise from dump data...")
    plotfilename_zoom = plotfilename[:-4]+"_zoom.png"
    plotfilename_all = plotfilename[:-4]+"_all.png"
    print("  >> " + plotfilename_zoom)
    print("  >> " + plotfilename_all)

    n_pix = int(config["npix"])
    fs = float(config["frow"]/n_pix)
    ylabel=r"Power density (dB/$\sqrt{Hz}$)"

    if record_len > len(data[0]):
        print("  Requested record length is too high ({0:5d})...".format(record_len))
        record_len = 2**(np.log10(len(data[0]))//np.log10(2))
        print("  Record length reduced to {0:5d}".format(record_len))
    else:
        print("  Record length: {0:5d}".format(record_len))

    n_records_max = int(len(data[0])/record_len)
    if n_records > n_records_max:
        print("  Too much records requested ({0:5d})...".format(n_records))
        n_records = n_records_max
        print("  Number of records reduced to {0:5d}".format(n_records))
    else:
        print("  Number of records: {0:5d}".format(n_records))

    noise_spectra = np.zeros((n_pix, int(record_len/2)+1))
    for pix in range(n_pix):
        noise_spectra[pix,:]=general_tools.do_spectrum(data[pix,0:n_records*record_len], record_len)

    noise_spectra_db = noise_spectra*0
    inotzero = np.where(noise_spectra != 0)
    noise_spectra_db[inotzero] = 20.*np.log10(noise_spectra[inotzero])
    n = len(noise_spectra_db[0])
    f = np.arange(n)*(fs/2)/n

    # plotting zoom
    fig = plt.figure(figsize=(10, 8))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.semilogx(f[3:-1], noise_spectra_db[pix_zoom,3:-1], marker='.')
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(frequency_label)
    ax1.set_title("Pixel {0:2d}".format(1+pix_zoom))
    ax1.grid(color='k', linestyle=':', linewidth=0.5)
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label]):
        item.set_weight('bold')
        item.set_fontsize(14)
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(12)
    #plt.show()
    fig.tight_layout()
    plt.savefig(plotfilename_zoom, bbox_inches='tight')

    fig = plt.figure(figsize=(18, 12))
    ncols, nlines = 9, 4
    ymin, ymax = 0.9*noise_spectra_db[:,3:-1].min(), 1.1*noise_spectra_db[:,3:-1].max()
    for pix in range(n_pix):
        ax = fig.add_subplot(nlines, ncols, 1+pix)
        ax.semilogx(f[3:-1], noise_spectra_db[pix,3:-1], marker='.')
        ax.set_title("Pixel {0:2d}".format(1+pix))
        ax.grid(color='k', linestyle=':', linewidth=0.5)
        ax.set_ylim(ymin, ymax)
        mosaic_labels(ax, pix, ncols, nlines, frequency_label, ylabel)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
            item.set_weight('bold')
            item.set_fontsize(12)
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(10)
    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename_all, bbox_inches='tight')

    return(noise_spectra)

# -----------------------------------------------------------------------------
def make_fake_science_noise(config, n_records=2000, record_len=8192):
    r"""
        This function builds a fake science data dump with gaussian noise.
        Two sine waves are added in the data of pixel 0.
        s
        config : dictionnary
        contains different informations such as path and directory names
        """
    n_pix = int(config["npix"])
    fs = float(config["frow"]/n_pix)
    t = np.arange(n_records*record_len)/fs
    freq1, freq2 = 500, 2000

    print('\n#---------------------')
    print('# Checking noise routine')
    n_pix = int(config["npix"])
    data = np.random.normal(loc=2**15, scale=20, size=(n_pix, n_records*record_len))

    # Adding a sine wave on pix 20 data
    data[0,:]+= 20*np.sin(2*np.pi*freq1*t) + 5*np.sin(2*np.pi*freq2*t)

    return(data)

# -----------------------------------------------------------------------------
def plot_adc_dump(data, plotfilename, config, t0=0, duration=0, spectral=False):
    r"""
        This function checks the data of a DRE-DEMUX ADC data dump.

        Parameters
        ----------
        data: numpy array
        The data of the dump 

        plotfilename: string
        Name of the plot file (with the path)

        config: dictionnary
        contains different informations such as path and directory names

        t0: number, optional
        begining of zoom in seconds (default is 0)

        duration: number, optional
        zoom duration in seconds. If 0 all the data are plotted (default is 0)

        spectral: boolean
        If True a spectral nalysis shall be done (default=False)

        Returns
        -------
        Nothing

        """
    ylabel_adu = "ADC values (ADU)"
    ylabel_v = "ADC values (V)"

    print(plotting_message)
    print("  >> " + plotfilename)

    fs = float(config["fs_ADC"])
    t = np.arange(len(data))/fs

    """
    Plotting data in time domain
    """
    imin = int(t0*fs)
    if imin > len(data[:]):
        print(error_t0_message)
        imin = 0
    if duration == 0 :
        imax = len(data[:])
    else :
        imax = min(int((t0+duration)*fs), len(data[:]))

    fig = plt.figure(figsize=(6, 8))
    xtitle = time_label
    mkr=''
    if imax-imin<300:
        mkr='.'

    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(1000*t[imin:imax], data[imin:imax]/float(config["adc_1volt_in_adu"]), 'b', marker=mkr)
    ax1.set_ylabel(ylabel_v)
    ax1.set_xlabel(xtitle)
    ax1.grid(color='k', linestyle=':', linewidth=0.5)
    y1min, y1max=ax1.get_ylim()
    ax11=ax1.twinx()
    ax11.set_ylabel(ylabel_adu)
    ax11.set_ylim(y1min*float(config["adc_1volt_in_adu"]), y1max*float(config["adc_1volt_in_adu"]))

    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(1000*t, data/float(config["adc_1volt_in_adu"]), 'b')
    ax2.set_ylabel(ylabel_v)
    ax2.set_xlabel(xtitle)
    ax2.grid(color='k', linestyle=':', linewidth=0.5)
    y2min, y2max=ax2.get_ylim()
    ax22=ax2.twinx()
    ax22.set_ylabel(ylabel_adu)
    ax22.set_ylim(y2min*float(config["adc_1volt_in_adu"]), y2max*float(config["adc_1volt_in_adu"]))

    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename, bbox_inches='tight')

    if spectral:
        # Plotting data in time domain
        fig = plt.figure(figsize=(8, 6))
        spt = general_tools.do_spectrum(data, int(2**20))
        spt_db = spt*0
        inotzero = np.where(spt != 0)[0]
        spt_db[inotzero] = 20*np.log10(spt[inotzero])

        f = np.arange(len(spt_db))*(fs/2)/len(spt_db)

        ax1 = fig.add_subplot(1, 1, 1)
        ax1.semilogx(f, spt_db, 'b')
        ax1.set_ylabel(spectral_power_density_label)
        ax1.set_xlabel(frequency_label)
        ax1.grid(color='k', linestyle=':', linewidth=0.5)
        ax1.set_xlim(1e3, f[-1])

        fig.tight_layout()
        #plt.show()
        plt.savefig(plotfilename[:-4]+"_F.png", bbox_inches='tight')

# -----------------------------------------------------------------------------
def plot_5mega_dump(data, plotfilename, config, title1, title2, t0=0, duration=0):
    r"""
        This function checks the data of a DRE-DEMUX ADC data dump.

        Parameters
        ----------
        data : numpy array
        Contains the 2 sets of data

        plotfilename : string
        Name of the plot file (with the path)

        config: dictionnary
        contains different informations such as path and directory names

        title1 : string
        Name of the first data set

        title2 : string
        Name of the second data set

        t0: number, optional
        begining of zoom in seconds (default is 0)

        duration: number, optional
        zoom duration in seconds. If 0 all the data are plotted (default is 0)

        Returns
        -------
        Nothing

        """
    print(plotting_message)
    print("  >> " + plotfilename)

    data1 = data[:, 0]
    data2 = data[:, 1]

    # In these dumps the 5MHz data are over sampled at 20MHz
    fs = float(config["frow"])*4 # Approx 20MHz
    t = np.arange(len(data1))/fs

    """
    Plotting data
    """
    imin = int(t0*fs)
    if imin > len(data1):
        print(error_t0_message)
        imin = 0
    if duration == 0 :
        imax = len(data1)
    else :
        imax = min(int((t0+duration)*fs), len(data1))

    fig = plt.figure(figsize=(10, 8))
    xtitle = "Time (ms)"
    mkr=''
    if imax-imin<300:
        mkr='.'

    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(1000*t[imin:imax], data1[imin:imax], 'b', marker=mkr)
    ax1.set_ylabel(title1)
    ax1.set_xlabel(xtitle)
    ax1.grid(color='k', linestyle=':', linewidth=0.5)

    ax3 = fig.add_subplot(2, 2, 3)
    ax3.plot(1000*t, data1, 'b')
    ax3.set_ylabel(title1)
    ax3.set_xlabel(xtitle)
    ax3.grid(color='k', linestyle=':', linewidth=0.5)

    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(1000*t[imin:imax], data2[imin:imax], 'b', marker=mkr)
    ax2.set_ylabel(title2)
    ax2.set_xlabel(xtitle)
    ax2.grid(color='k', linestyle=':', linewidth=0.5)

    ax4 = fig.add_subplot(2, 2, 4)
    ax4.plot(1000*t, data2, 'b')
    ax4.set_ylabel(title2)
    ax4.set_xlabel(xtitle)
    ax4.grid(color='k', linestyle=':', linewidth=0.5)

    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename, bbox_inches='tight')

# -----------------------------------------------------------------------------
def plot_counter_dump(data, plotfilename, config):
    r"""
        This function checks the data of a DRE-DEMUX COUNTER data dump.

        Parameters
        ----------
        data : numpy array
        The data of the dump 

        plotfilename : string
        Name of the plot file (with the path)

        config: dictionnary
        contains different informations such as path and directory names

        Returns
        -------
        Nothing

        """
    print(plotting_message)
    print("  >> " + plotfilename)

    fs = float(config["frow"])*4 # Approx 20MHz
    t = np.arange(len(data))/fs

    reference = np.arange(len(data))
    if np.abs(data-reference).max()>0:
        print(" Error in the counter data.")

    # Plotting data
    fig = plt.figure(figsize=(6, 8))
    fig.suptitle("Counter32 dump file")
    xtitle = time_label

    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(1000*t, data, 'b')
    ytitle = "32-bit Counter values"
    ax1.set_ylabel(ytitle)
    ax1.set_xlabel(xtitle)
    ax1.grid(color='k', linestyle=':', linewidth=0.5)

    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(1000*t, data-reference, 'b')
    ytitle = "32-bit counter - reference"
    ax2.set_ylabel(ytitle)
    ax2.set_xlabel(xtitle)
    ax2.grid(color='k', linestyle=':', linewidth=0.5)

    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename, bbox_inches='tight')

# -----------------------------------------------------------------------------

"""
Testing the routines with test data.
"""
if __name__ == "__main__":

    print("Testing the dumps analysis...")
    config=general_tools.configuration("demux_tools_cfg")
    config.config["path"]="."
    config.config["dir_data"]="test_data"
    config.config["dir_plots"]="test_plots"
    datadirname = os.path.join(os.path.normcase(config.config['dir_data'])) 
    dumpfilenames = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]==".dat"]

    for file in dumpfilenames:
        print('\n#---------------------')
        d=get_data.data(file, config.config)
        d.print_dumptype()
        d.plot(config, t0=0, duration=0, pix_zoom=0, spectral=True, noise=True, check_noise_measurement=True)

# -----------------------------------------------------------------------------
