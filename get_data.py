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
#  get_data.py
#

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
    def __init__(self, filename):
        self.filename=filename
        self.values, self.dumptype, self.config=read_data(filename)

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

    def plot(self, p):
        r"""
            This function checks the data of a DRE-DEMUX ADC data dump.

            Parameters
            ----------
 
            p: dictionnary including all the parameters

                p['t0']: number, optional
                begining of zoom in seconds (default is 0)

                p['duration']: duration: number, optional
                zoom duration in seconds. If 0 all the data are plotted (default is 0)

                p['pix_zoom']: number (integer)
                pixel id refering to the pixel for which we plot a zoom (default=0)

            Returns
            -------
            Nothing

            """

        if self.dumptype == 0:
            data1 = self.values[:, 0]
            data2 = self.values[:, 1].astype('float')
            data2[data2<0]+=2**16 # Convertion to 16-bit unsigned format 
            plotting_tools.plot_5mega_dump(data1, data2, self.config, "DACFB1", "Science", p)
        if self.dumptype == 1:
            data1 = self.values[:, 0]
            data2 = self.values[:, 1].astype('float')
            data2[data2<0]+=2**16 # Convertion to 16-bit unsigned format 
            plotting_tools.plot_5mega_dump(data1, data2, self.config, "ERROR", "Science", p)
        if self.dumptype == 2:
            data1 = self.values[:, 0]
            data2 = self.values[:, 1].astype('float')
            data2[data2<0]+=2**16 # Convertion to 16-bit unsigned format 
            plotting_tools.plot_5mega_dump(data1, data2, self.config, "DACFB2", "Science", p)
        if self.dumptype == 4:
            data1 = self.values[:, 0]
            data2 = self.values[:, 1]
            plotting_tools.plot_5mega_dump(data1, data2, self.config, "DACFB1", "DACFB2", p)
        if self.dumptype == 5:
            plotting_tools.plot_adc_dump(self.values, self.config, p)
        if self.dumptype == 8:
            plotting_tools.plot_science_dump(self.values, self.config, p)
        if self.dumptype == 9:
            data1 = self.values[:, 0]
            data2 = self.values[:, 1]
            plotting_tools.plot_5mega_dump(data1, data2, self.config, "ERROR", "DACFB1", p)
        if self.dumptype == 15:
            plotting_tools.plot_counter_dump(self.values, self.config)


# -----------------------------------------------------------------------------
def read_data(datafilename):
    r"""
        This function reads data from a DRE dump file, and returns 1 array
        (format <h).
        The file header and the dumptype (2 firsts 32-bit words) are removed from the data

        Parameters
        ----------
        datafilename : string
        The name of the dump file (with the path and the extension)

        Returns
        -------
        data : array type
        contains the values of the file (format int16).

        dumptype : number
        Dumptype id.

        config : dictionnary
        contains different informations such as path and directory names

        """

    config=general_tools.configuration("demux_tools_cfg")
    config=config.config

    """
    Managing filenames
    """
    dirname = os.path.join(os.path.normcase(config['path']), config['dir_data'])
    fulldatafilename = os.path.join(dirname, datafilename)
    config['fulldatafilename']=fulldatafilename
    config['datafilename']=datafilename

    plotdirname = os.path.join(os.path.normcase(config['path']), config['dir_plots'])
    general_tools.checkdir(plotdirname)
    plotfilename = datafilename[int(config['length_of_date_in_filenames']):-4]+".png"
    fullplotfilename = os.path.join(plotdirname, plotfilename)
    file_index=1
    # Managing multiple plot files with same name
    while os.path.isfile(fullplotfilename):
        print(fullplotfilename + " already exists")
        if file_index==1:
            length_to_be_removed=4
        else:
            length_to_be_removed=len(plotfilename.split('-')[-1])+1
        file_index+=1
        plotfilename = plotfilename[:-length_to_be_removed]+"-{0:1d}.png".format(file_index)
        fullplotfilename = os.path.join(plotdirname, plotfilename)

    config['fullplotfilename']=fullplotfilename
    config['plotfilename']=plotfilename

    """
    Opening file
    """
    fdat=open(fulldatafilename, 'rb')
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

    return(data, dumptype, config)

# -----------------------------------------------------------------------------
