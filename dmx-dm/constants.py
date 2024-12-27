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
#  cosim_constants.py
#
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
fSamp = 125e6
nColPerDemux = 4
muxFactor = 34
nSamplesPerRow = 20
fRow = fSamp / nSamplesPerRow
fFrame = fRow / muxFactor
# Error signal ADC
dmxNbBitsADCError = 14
fsrADCErrorADU = 2**dmxNbBitsADCError
fsrADCErrorV = 2
# Feedback signal DAC
dmxNbBitsDACFdbk = 14
fsrDACFdbkADU = 2**dmxNbBitsDACFdbk
fsrDACFdbkV = 2
# Offset compensation signal DAC
dmxNbBitsDACOfcoCoarse = 12
fsrDACOfcoCoarseADU = 2**dmxNbBitsDACOfcoCoarse
fsrDACOfcoCoarseV = 1

fFilterAdc = 20e6 # Anti aliasing filter cutoff frequency in Hz

# ---------------------------------------------------------------------------------
# Directories
dataDirName = 'DATA'
plotDirName = 'PLOTS'

# ---------------------------------------------------------------------------------
