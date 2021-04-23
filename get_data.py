#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os

import general_tools, plotting_tools

# Some general text definitions
plotting_message="Plotting dump data..."
time_label="Time (ms)"
frequency_label="Frequency (Hz)"
spectral_power_density_label=r"Power spectral density (dB/$\sqrt{Hz}$)"
DADA=-9510      # 0xDADA interpreted as int16

# -----------------------------------------------------------------------------
class data:
    def __init__(self, filename, config):
        self.filename=filename
        self.config=config
        self.values, self.dumptype=read_data(filename, config)

    def print_dumptype(self):
        if self.dumptype==0:
            dumptype_str = "DACFB1, SCIENCE"
        elif self.dumptype==1:
            dumptype_str = "ERROR5MHZ, SCIENCE"
        elif self.dumptype==2:
            dumptype_str = "DACFB2, SCIENCE"
        elif self.dumptype==4:
            dumptype_str = "DACFB1, DACFB2"        
        elif self.dumptype==5:
            dumptype_str = "ERROR40MHZ"        
        elif self.dumptype==6:
            dumptype_str = "DACFB1"        
        elif self.dumptype==8:
            dumptype_str = "SCIENCE PACKETIZED"        
        elif self.dumptype==9:
            dumptype_str = "ERROR5MHz, DACFB1"                
        elif self.dumptype==15 or self.dumptype==10:
            dumptype_str = "COUNTER32"                
        else:
            raise ValueError('Wrong dump type')
        print(dumptype_str)

    def plot(self, t0=0, duration=0, pix_zoom=0, spectral=False, noise=False, sav_spectra=False):
        r"""
            This function checks the data of a DRE-DEMUX ADC data dump.

            Parameters
            ----------
            t0: number, optional
            begining of zoom in seconds (default is 0)

            duration: number, optional
            zoom duration in seconds. If 0 all the data are plotted (default is 0)

            pix_zoom: number (integer)
            pixel id refering to the pixel for which we plot a zoom (default=0)

            spectral : boolean
            If True a spectral nalysis shall be done (default=False)

            noise: boolean
            Indicates if a noise analysis shall be done (default=False)

            check_noise_measurement: boolean
            indicates if the function shall tested on fake data (default=False)

            sav_spectra: boolean
            indicates if spectra shall be saved in npy file (default=False)

            Returns
            -------
            Nothing

            """

        plotdirname = os.path.join(os.path.normcase(self.config['path']), self.config['dir_plots'])
        general_tools.checkdir(plotdirname)
        plotfilename = os.path.join(plotdirname, self.filename[:-4]+".png")

        if self.dumptype == 0:
            data1 = self.values[:, 0]
            data2 = self.values[:, 1]+2**16 # Convertion to 16-bit unsigned format 
            plotting_tools.plot_5mega_dump(data1, data2, plotfilename, self.config, "DACFB1", "Science", t0, duration)
        if self.dumptype == 1:
            data1 = self.values[:, 0]
            data2 = self.values[:, 1]+2**16 # Convertion to 16-bit unsigned format 
            plotting_tools.plot_5mega_dump(data1, data2, plotfilename, self.config, "ERROR", "Science", t0, duration)
        if self.dumptype == 2:
            data1 = self.values[:, 0]
            data2 = self.values[:, 1]+2**16 # Convertion to 16-bit unsigned format 
            plotting_tools.plot_5mega_dump(data1, data2, plotfilename, self.config, "DACFB2", "Science", t0, duration)
        if self.dumptype == 4:
            data1 = self.values[:, 0]
            data2 = self.values[:, 1]
            plotting_tools.plot_5mega_dump(data1, data2, plotfilename, self.config, "DACFB1", "DACFB2", t0, duration)
        if self.dumptype == 5:
            plotting_tools.plot_adc_dump(self.values, plotfilename, self.config, t0, duration, spectral, sav_spectra)
        if self.dumptype == 8:
            plotting_tools.plot_science_dump(self.values, plotfilename, self.config, t0, duration, pix_zoom, noise, sav_spectra)
        if self.dumptype == 9:
            data1 = self.values[:, 0]
            data2 = self.values[:, 1]
            plotting_tools.plot_5mega_dump(data1, data2, plotfilename, self.config, "ERROR", "DACFB1", t0, duration)
        if self.dumptype == 15:
            plotting_tools.plot_counter_dump(self.values, plotfilename, self.config)


# -----------------------------------------------------------------------------
def read_data(dumpfilename, config):
    r"""
        This function reads data from a DRE dump file, and returns 1 array
        (format <h).
        The file header and the dumptype (2 firsts 32-bit words) are removed from the data

        Parameters
        ----------
        dumpfilename : string
        The name of the dump file (with the path and the extension)

        config : dictionnary
        contains different informations such as path and directory names

        Returns
        -------
        data : array type
        contains the values of the file (format int16).

        dumptype : number
        Dumptype id.
        """

    dirname = os.path.join(os.path.normcase(config['path']), config['dir_data'])
    fullfilename = os.path.join(dirname, dumpfilename)

    fdat=open(fullfilename, 'rb')
    data=np.fromfile(fdat, dtype='<h')
    fdat.close()
        
    DADA=-9510      # 0xDADA interpreted as int16
    if data[0] != DADA:
        raise ValueError('Problem with file format!')

    """
    reading dumptype
    """
    header24=int(data[1].astype('uint16')/2**12)
    header23=int((data[1].astype('uint16')-header24*2**12)/2**8)
    dumptype=int((data[1].astype('uint16')-header24*2**12-header23*2**8)/2**4)

    """
    removing header
    """
    data=data[2:]

    """
    decommutation of data according to dumptype
    """
    ## 5MHz dump data
    if dumptype == 0 or dumptype == 1 or dumptype == 2 or dumptype == 4 or dumptype == 9:
        data = np.resize(data, (len(data)//2, 2))

    ## ADC dump data
    if dumptype == 5:
        # input values are in the wrong order
        data2 = np.resize(data, (len(data)//2, 2))
        data[0::2]=data2[:,1]
        data[1::2]=data2[:,0]

    ## science data
    if dumptype == 8:
        # merging 16-bit words in 32-bit words
        data = np.resize(data, (len(data)//2, 2)).astype('uint16')
        data = data[:, 1]*2.**16 + data[:, 0]

        # the first header has been removed by read_dump function. I need it here. I add a fake one.
        data = np.append(np.zeros(1), data)

        # reordering data
        npix = int(config["npix"])
        data = np.transpose(np.resize(data, (len(data)//(npix+2), npix+2)))

        # removing periodic headers
        # keeping packet_number
        data = data[1:,:]

    ## 32-bit counter data
    if dumptype == 15:
        data = np.resize(data, (len(data)//2, 2)).astype('uint16')
        data = data[:, 1]*2.**16 + data[:, 0]

    return(data, dumptype)

# -----------------------------------------------------------------------------

