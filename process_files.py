#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

import general_tools, get_data, plotting_tools, ep_tools

config=general_tools.configuration("demux_tools_cfg")

datadirname = os.path.join(os.path.normcase(config.config['path']), config.config['dir_data'])

"""
Processing dumps files
"""
dumpfilenames = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]==".dat" and f[-7:-4]!="_er"]

for file in dumpfilenames:
    print('\n#---------------------')
    d=get_data.data(file, config.config)
    d.print_dumptype()
    d.plot(t0=0, duration=0, pix_zoom=0, spectral=True, noise=True, check_noise_measurement=False)

"""
Processing energy resolution data
"""
filename_noise_er = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-12:]=="noise_er.dat"]
filename_calib_er = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-12:]=="calib_er.dat"]
filename_measu_er = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-12:]=="measu_er.dat"]

if len(filename_noise_er)>0 and len(filename_calib_er)>0 and len(filename_measu_er)>0:
    pix = 0
    print('\n#---------------------')
    print("Processing energy resolution data of pixel {0:2d}...".format(pix))
    record_len = 4096
    verbose=False
    _, _= ep_tools.ep(filename_noise_er[0], filename_calib_er[0], filename_measu_er[0], pix, config.config, record_len, verbose)
