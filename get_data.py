import numpy as np
from numpy.fft import rfft
import os
import matplotlib.pyplot as plt
import general_tools

# Some general text definitions
plotting_message="Plotting dump data..."
time_label="Time (ms)"
frequency_label="Frequency (Hz)"
spectral_power_density_label=r"Power spectral density (dB/$\sqrt{Hz}$)"

# -----------------------------------------------------------------------------
def dumptype2string(dumptype):
    r"""
        This function convertsbthe dumptype id into a human readable string

        Parameters
        ----------
        dumptype : integer
        The dumptype id

        Returns
        -------
        dumptype_str : string
        Description of the dumptype.
        """

    if dumptype==0:
        dumptype_str = "DACFB1, SCIENCE"
    elif dumptype==1:
        dumptype_str = "ERROR5MHZ, SCIENCE"
    elif dumptype==2:
        dumptype_str = "DACFB2, SCIENCE"
    elif dumptype==4:
        dumptype_str = "DACFB1, DACFB2"        
    elif dumptype==5:
        dumptype_str = "ERROR40MHZ"        
    elif dumptype==6:
        dumptype_str = "DACFB1"        
    elif dumptype==8:
        dumptype_str = "SCIENCE PACKETIZED"        
    elif dumptype==9:
        dumptype_str = "ERROR5MHz, DACFB1"                
    elif dumptype==15 or dumptype==10:
        dumptype_str = "COUNTER32"                
    else:
        raise ValueError('Wrong dump type')

    return(dumptype_str)

# -----------------------------------------------------------------------------
def read_dumptype(dumpfilename):
    r"""
        This function reads the dumptype of a dumpfile

        Parameters
        ----------
        dumpfilename : string
        The name of the dump file (with the path and the extension)

        Returns
        -------
        dumptype : number
        the dumptype id.
        
        """
    fdat=open(dumpfilename, 'rb')
    data=np.fromfile(fdat, dtype='<h')
    fdat.close()
        
    DADA=-9510      # 0xDADA interpreted as int16
    if data[0] != DADA:
        raise ValueError('Problem with file format!')
    header2=data[1].astype('uint16')
    header24=int(header2/2**12)
    header23=int((header2-header24*2**12)/2**8)
    header22=int((header2-header24*2**12-header23*2**8)/2**4)
    dumptype=header22
    print("# Dumptype is:", dumptype2string(dumptype))
    return(dumptype)

# -----------------------------------------------------------------------------
def read_dump(dumpfilename):
    r"""
        This function reads data from a DRE adc dump file, and returns 1 array
        (format <h).
        The file header (first 32-bit) word is removed from the data

        Parameters
        ----------
        dumpfilename : string
        The name of the dump file (with the path and the extension)

        Returns
        -------
        data : array type
        contains the values of the file (format int16).
                    
        """
    
    fdat=open(dumpfilename, 'rb')
    data=np.fromfile(fdat, dtype='<h')
    fdat.close()
        
    DADA=-9510      # 0xDADA interpreted as int16
    if data[0] != DADA:
        raise ValueError('Problem with file format!')
    # removing header
    data=data[2:]
    return(data)

# -----------------------------------------------------------------------------
def check_dump(config, filename, max_duration=0.2, spectral=False, quiet=True):
    r"""
        This function checks the data of a DRE-DEMUX ADC data dump.

        Parameters
        ----------
        config : dictionnary
        contains different informations such as path and directory names

        filename : string
        The name of the dump file (without the path)

        Max_duration : number, optional
        Maximum duration in seconds to be considered for analysis (default is 0.2)

        spectral : boolean
        If True a spectral nalysis shall be done (default=False)

        quiet : boolean
        Defines if text info shall be written by the routine

        Returns
        -------
        Nothing

        """

    fullfilename = os.path.join(os.path.normcase(config['dir_data']), filename)
    plotdirname = os.path.normcase(config['dir_plots'])
    general_tools.checkdir(plotdirname)
    plotfilename = os.path.join(plotdirname, filename[:-4]+".png")

    # verifying the dumptype
    dumptype=read_dumptype(fullfilename)
    data = read_dump(fullfilename)

    if dumptype == 0:
        plot_5mega_dump(data, plotfilename, "DACFB1", "Science", max_duration)
    if dumptype == 1:
        plot_5mega_dump(data, plotfilename, "ERROR", "Science", max_duration)
    if dumptype == 2:
        plot_5mega_dump(data, plotfilename, "DACFB2", "Science", max_duration)
    if dumptype == 4:
        plot_5mega_dump(data, plotfilename, "DACFB1", "DACFB2", max_duration)
    if dumptype == 5:
        plot_adc_dump(data, plotfilename, config, max_duration, spectral)
    if dumptype == 8:
        plot_science_dump(data, plotfilename, config, max_duration)
    if dumptype == 9:
        plot_5mega_dump(data, plotfilename, "ERROR", "DACFB1", max_duration)
    if dumptype == 15:
        plot_counter_dump(data, plotfilename)

# -----------------------------------------------------------------------------
def plot_science_dump(data, plotfilename, config, max_duration=0.2):
    r"""
        This function checks the data of a DRE-DEMUX science data dump.

        Parameters
        ----------
        data : numpy array
        The data of the dump 

        plotfilename : string
        Name of the plot file (with the path)

        config : dictionnary
        contains different informations such as path and directory names

        Max_duration : number, optional
        Maximum duration in seconds to be considered for analysis (default is 0.2)

        Returns
        -------
        Nothing

        """
    print(plotting_message)

    # merging 16-bit words in 32-bit words
    data = np.resize(data, (len(data)//2, 2)).astype('uint16')
    data = data[:, 1]*2.**16 + data[:, 0]

    # the first haeder has been removed by read_dump function. I need it here. I add a fake one.
    data = np.append(np.zeros(1), data)

    # reordering data
    npix = int(config["npix"])
    data = np.transpose(np.resize(data, (len(data)//(npix+2), npix+2)))

    fs = 5e6/npix
    t = np.arange(len(data[0,:]))/fs

    # Plotting data in time domain
    fig = plt.figure(figsize=(8, 10))
    xtitle = time_label

    ncols, nlines = 4, 9
    for pix in range(npix):
        ax = fig.add_subplot(nlines, ncols, 1+pix)
        ax.plot(1000*t, data[2+pix,:], 'b')
        ax.set_ylabel("Pixel {0:2d}".format(pix))
        ax.grid(color='k', linestyle=':', linewidth=0.5)
        if pix/ncols>=(nlines-1):
            ax.set_xlabel(xtitle)
    ax = fig.add_subplot(nlines, ncols, npix+1)
    ax.plot(1000*t, data[1,:], 'b')
    ax.set_ylabel("Packet number")
    ax.grid(color='k', linestyle=':', linewidth=0.5)
    ax.set_xlabel(xtitle)


    fig.tight_layout()
    #plt.show()
    plt.savefig(plotfilename, bbox_inches='tight')


# -----------------------------------------------------------------------------
def plot_adc_dump(data, plotfilename, config, max_duration=0.2, spectral=False):
    r"""
        This function checks the data of a DRE-DEMUX ADC data dump.

        Parameters
        ----------
        data : numpy array
        The data of the dump 

        plotfilename : string
        Name of the plot file (with the path)

        config : dictionnary
        contains different informations such as path and directory names

        Max_duration : number, optional
        Maximum duration in seconds to be considered for analysis (default is 0.2)

        spectral : boolean
        If True a spectral nalysis shall be done (default=False)

        Returns
        -------
        Nothing

        """
    ylabel_adu = "ADC values (ADU)"
    ylabel_v = "ADC values (V)"

    # input values are in the wrong order
    data2 = np.resize(data, (len(data)//2, 2))
    data[0::2]=data2[:,1]
    data[1::2]=data2[:,0]

    print(plotting_message)
    fs = 40e6
    t = np.arange(len(data))/fs

    # Plotting data in time domain
    fig = plt.figure(figsize=(6, 8))
    xtitle = time_label
    x_zoom_max = int(max_duration * fs)

    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(1000*t[0:x_zoom_max], data[0:x_zoom_max]/float(config["adc_1volt_in_adu"]), 'b')
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
        spt_db = 20*np.log10(spt)

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
def plot_5mega_dump(data, plotfilename, title1, title2, max_duration=0.2):
    r"""
        This function checks the data of a DRE-DEMUX ADC data dump.

        Parameters
        ----------
        data : numpy array
        Contains the 2 sets of data

        plotfilename : string
        Name of the plot file (with the path)

        title1 : string
        Name of the first data set

        title2 : string
        Name of the second data set

        Max_duration : number, optional
        Maximum duration in seconds to be considered for analysis (default is 0.2)

        Returns
        -------
        Nothing

        """
    print(plotting_message)

    data = np.resize(data, (len(data)//2, 2))
    data1 = data[:, 0]
    data2 = data[:, 1]

    fs = 5e6
    t = np.arange(len(data1))/fs

    # Plotting data
    fig = plt.figure(figsize=(10, 8))
    xtitle = "Time (ms)"
    x_zoom_max = int(max_duration * fs)

    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(1000*t[0:x_zoom_max], data1[0:x_zoom_max], 'b')
    ax1.set_ylabel(title1)
    ax1.set_xlabel(xtitle)
    ax1.grid(color='k', linestyle=':', linewidth=0.5)

    ax3 = fig.add_subplot(2, 2, 3)
    ax3.plot(1000*t, data1, 'b')
    ax3.set_ylabel(title1)
    ax3.set_xlabel(xtitle)
    ax3.grid(color='k', linestyle=':', linewidth=0.5)

    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(1000*t[0:x_zoom_max], data2[0:x_zoom_max], 'b')
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
def plot_counter_dump(data, plotfilename):
    r"""
        This function checks the data of a DRE-DEMUX COUNTER data dump.

        Parameters
        ----------
        data : numpy array
        The data of the dump 

        plotfilename : string
        Name of the plot file (with the path)

        Returns
        -------
        Nothing

        """
    print(plotting_message)

    data = np.resize(data, (len(data)//2, 2)).astype('uint16')
    data = data[:, 1]*2.**16 + data[:, 0]

    fs = 20e6
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
