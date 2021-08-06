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
#  process_files.py
#

import numpy as np
import os

import general_tools, get_data, ep_tools, measure_ki

config=general_tools.configuration("demux_tools_cfg")

# Creation of the plot directory
plotdirname = os.path.join(os.path.normcase(config.config['path']), config.config['dir_plots'])
general_tools.purgedir(plotdirname)

parameters={\
        't0': 0, \
        'duration':0.1, \
        'pix_zoom':0 \
        }

datadirname = os.path.join(os.path.normcase(config.config['path']), config.config['dir_data'])

drawline = '\n#--------------------------------------'

"""
Processing dumps files ###############################################
"""
dumpfilenames = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]==".dat" \
                and f[-13:]!="_er_calib.dat" \
                and f[-13:]!="_er_noise.dat" \
                and f[-13:]!="_er_measu.dat" ]

for filename in dumpfilenames:
    print(drawline)
    print('Analysing file ', filename)
    d=get_data.data(filename)
    print('Dumptype is: ', end="")
    d.print_dumptype()
    d.plot(parameters)

"""
Processing ki measurements data #####################################
"""
kifilename = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-13:]=="_ki_check.dat" ]
if len(kifilename)>0:
    print(drawline)
    print('Analysing file ', kifilename[0])
    d=get_data.data(kifilename[0])
    print('Dumptype is: ', end="")
    d.print_dumptype()
    measure_ki.measure_ki(d)

"""
Processing energy resolution data ###################################
"""
remove_raw_files=True

"""
Making noise records
"""
record_len=4096
prebuffer=180
rec_extension='_rec.npy' 

filename_er_noise = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-13:]=="_er_noise.dat"]
filename_er_noise_rec = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-17:]=="_er_noise_rec.npy"]
if len(filename_er_noise)>0 and len(filename_er_noise_rec)==0:
    print(drawline)
    noise_parameters={\
        'record_len':record_len, \
        'pix':0, \
        'remove_raw_files':remove_raw_files \
        }
    noise = ep_tools.get_noise_records(filename_er_noise[0], noise_parameters)
    # saving noise records
    noise_npy_file = os.path.join(datadirname, filename_er_noise[0][:-4]+rec_extension)
    with open(noise_npy_file, 'wb') as file:
        np.save(file, noise)
    
"""
Triggering calibration pulses
"""

filename_er_calib = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-13:]=="_er_calib.dat"]
filename_er_calib_rec = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-17:]=="_er_calib_rec.npy"]
if len(filename_er_calib)>0 and len(filename_er_calib_rec)==0:
    print(drawline)
    calib_parameters={\
        'record_len':record_len, \
        'prebuffer':prebuffer, \
        'pix':100, \
        'remove_raw_files':remove_raw_files \
        }
    # if pix==100 the routine will look for pixels with pulses

    calib, t_calib, pix_calib = ep_tools.get_pulse_records(filename_er_calib[0], calib_parameters)
    # saving pulse records
    calib_npy_file = os.path.join(datadirname, filename_er_calib[0][:-4]+rec_extension)
    with open(calib_npy_file, 'wb') as file:
        np.save(file, t_calib)
        np.save(file, calib)

"""
Triggering pulses
"""
filename_er_measu = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-13:]=="_er_measu.dat"]
filename_er_measu_rec = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-17:]=="_er_measu_rec.npy"]
if len(filename_er_measu)>0 and len(filename_er_measu_rec)==0:
    print(drawline)
    measu_parameters={\
        'record_len':record_len+4, \
        'prebuffer':prebuffer+2, \
        'pix':pix_calib, \
        'remove_raw_files':remove_raw_files \
        }

    measu, t_measu, _ = ep_tools.get_pulse_records(filename_er_measu[0], measu_parameters)
    # saving pulse records
    measu_npy_file = os.path.join(datadirname, filename_er_measu[0][:-4]+rec_extension)
    with open(measu_npy_file, 'wb') as file:
        np.save(file, t_measu)
        np.save(file, measu)

"""
Computing energy resolution
"""
filename_er_noise_rec = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-17:]=="_er_noise_rec.npy"]
filename_er_calib_rec = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-17:]=="_er_calib_rec.npy"]
filename_er_measu_rec = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-17:]=="_er_measu_rec.npy"]

if len(filename_er_noise_rec)>0 and len(filename_er_calib_rec)>0 and len(filename_er_measu_rec)>0:
    print(drawline)
    print("Processing energy resolution data...")
    verbose=False
    _, _= ep_tools.ep(filename_er_noise_rec[0], filename_er_calib_rec[0], filename_er_measu_rec[0], config.config, prebuffer, verbose)

print(drawline)

"""
#####################################################################
"""
