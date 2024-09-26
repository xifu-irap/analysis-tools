# ---------------------------------------------------------------------------------
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
# ---------------------------------------------------------------------------------
#
#  laurent.ravera@irap.omp.eu
#  fakes.py
#
# ---------------------------------------------------------------------------------

import h5py
import os
import numpy as np
import constants as cst
import general_tools as gentools


def fake_dump(path, fileName, shift=0, sigmaNoise=0, shape='single'):
    # Crée un tableau de 1360 valeurs avec des plateaux successifs

    gentools.createdir(cst.dataDirName)
    fullFileName = os.path.join(path, cst.dataDirName, fileName)

    nSamplesPerRow = cst.nSamplesPerRow
    muxFactor = cst.muxFactor

    low_level = 0
    high_level = 2**10
    nPts = nSamplesPerRow*(muxFactor+1)*2
    data = np.zeros(nPts)

    if shape == 'single':
        # nSamplesPerRowx2 points are added at the beginning to apply the filter and the shift. 
        # They are removed afterward.
        data = ([low_level]*(nSamplesPerRow*2) 
                + [high_level]*nSamplesPerRow 
                + [low_level]*nSamplesPerRow*(muxFactor-1) 
                + [high_level]*nSamplesPerRow 
                + [low_level]*nSamplesPerRow*(muxFactor-1))

    if shape == 'alternate':
        nTiles = muxFactor + 1 # +1 to make easier the filtering, will be removed afterward
        data = np.tile([low_level] * nSamplesPerRow + [high_level] * nSamplesPerRow, nTiles)

    # Addition of noise
    data = data + np.random.normal(loc=0, scale=sigmaNoise, size=nPts)

    # Shifting the data
    data = np.roll(data, shift)

    # Filtering the data
    data = gentools.first_order_low_pass_filter(data, 1/(2*np.pi*20e6))

    # Selection of the dump data
    data = data[-nSamplesPerRow*muxFactor*2:]

    # Writing the data in a HDF5 file
    with h5py.File(fullFileName, 'w') as file:
        file.create_dataset('dump', data=data)

    print("Dump data saved to ", fullFileName)

    return data

def fake_error(fileName, nPatterns, shape = 'square', sigmaNoise=0):
    import matplotlib.pyplot as plt

    gentools.createdir(cst.dataDirName)
    fullFileName = os.path.join(cst.dataDirName, fileName)
    addingNoise = True

    if shape == 'sine':
        f=1/(64*cst.muxFactor)*1.347982159
        t=np.arange(0, 256*cst.muxFactor*nPatterns)
        data = 2**13*np.sin(2*np.pi*f*t)

    if shape == 'saw':
        # Making a ramp
        a = -2**11  # Value at the beginning of the ramp
        b = 4       # Number of frames per step
        c = 64      # Increment between two steps
        N = 1024    # Number of steps

        ramp = a + np.repeat(np.arange(N)*c, b*cst.muxFactor)
        pattern = np.concatenate([ramp, np.flip(ramp)])

        # repeating the pattern
        data = np.tile(pattern, nPatterns)

    if shape == 'square':
        a = -2**13  # lowest values
        b = 8       # Number of frames per step
        c = 2**14-1 # Increment between two steps
        N = 2       # Number of steps

        pattern = a + np.repeat([0, c], b*cst.muxFactor)

        # repeating the pattern
        data = np.tile(pattern, nPatterns)

    nPts = data.size
    totalDuration = (1/cst.fRow) * nPts
    print(" Error data file corresponds to {0:6.2f} seconds".format(totalDuration))

    if addingNoise:
        # Addition of noise
        noise = np.random.normal(loc=0, scale=sigmaNoise, size=nPts)

        # Filtering the noise (antialiasing filter)
        # noise = gentools.first_order_low_pass_filter(noise, 1/(2*np.pi*cst.fFilterAdc))

        data += noise

    # Quantization
    data = np.round(data)

    print("Error data saved to ", fullFileName)
    print(" Data corresponds to a duration of {0:3.2f} s".format(nPts/cst.fRow))

    return data

# ---------------------------------------------------------------------------------
