# ---------------------------------------------------------------------------------
# !/usr/bin/env python
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
#  constants.py
#
# ---------------------------------------------------------------------------------

import numpy as np

# ---------------------------------------------------------------------------------
nColPerDemux = 4
nPixPerCol = 34
muxFactor = nPixPerCol
nSamplesPerRow = 20

# Frequencies
fSamp = 125e6
fRow = fSamp / nSamplesPerRow
fFrame = fRow / muxFactor

# Error signal ADC
dmxNbBitsADCError = 14
fsrADCErrorADU = 2 ** dmxNbBitsADCError
fsrADCErrorV = 2

# Feedback signal DAC
dmxNbBitsDACFdbk = 14
fsrDACFdbkADU = 2 ** dmxNbBitsDACFdbk
fsrDACFdbkV = 2

# Offset compensation signal DAC
dmxNbBitsDACOfcoCoarse = 12
fsrDACOfcoCoarseADU = 2 ** dmxNbBitsDACOfcoCoarse
fsrDACOfcoCoarseV = 1

# Test patterns
tstPattern_nb_regions = 5
tstPattern_nb_params_per_region = 4

fFilterAdc = 20e6  # Anti aliasing filter cutoff frequency in Hz

# Directories
dataDirName = 'dmx_data'
scanDirName = 'scans'
hkDirName = 'hks'
plotDirName = 'PLOTS'
errorSpectraDirname = 'error_spectra'
spectraDirname = 'spectra'

# DEMUX models
dmx_models = {0: "DM DMX", 1: "EM DMX", 2: "PFM DMX", 3: "FM DMX", 7: "DevKit"}

BASE_DATA_PATH = "."

# Noise requirements
## Wide band noise spectral density in V/sqrt(Hz) with BW = 62.5 MHz
nsd_erro = 25e-9
nsd_fdbk = 20e-9
nsd_ofco = 15e-9
nsd = {"ERRO-ONLY": nsd_erro,
       "FDBK-ONLY": nsd_fdbk,
       "OFCO-ONLY": nsd_ofco,
       "FDBK-ERRO": np.sqrt(nsd_erro ** 2 + nsd_fdbk ** 2),
       "OFCO-ERRO": np.sqrt(nsd_erro ** 2 + nsd_ofco ** 2)}
one_over_f_at_1hz_erro = 4e-6
one_over_f_at_1hz_fdbk = 0.7e-6
one_over_f_at_1hz_ofco = 2.7e-6
one_over_f_at_1hz = {"ERRO-ONLY": one_over_f_at_1hz_erro,
                     "FDBK-ONLY": one_over_f_at_1hz_fdbk,
                     "OFCO-ONLY": one_over_f_at_1hz_ofco,
                     "FDBK-ERRO": np.sqrt(one_over_f_at_1hz_erro ** 2 + one_over_f_at_1hz_fdbk ** 2),
                     "OFCO-ERRO": np.sqrt(one_over_f_at_1hz_erro ** 2 + one_over_f_at_1hz_ofco ** 2)}

# ---------------------------------------------------------------------------------
