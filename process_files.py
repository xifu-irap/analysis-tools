#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

import general_tools, get_data, plotting_tools

config=general_tools.configuration("demux_tools_cfg")

datadirname = os.path.join(os.path.normcase(config.config['dir_data']))

dumpfilenames = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]==".dat"]

for file in dumpfilenames:
    print('\n#---------------------')
    d=get_data.data(file, config.config)
    d.print_dumptype()
    d.plot(config, t0=0, duration=0, pix_zoom=0, spectral=True, noise=True, check_noise_measurement=False)
