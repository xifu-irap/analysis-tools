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
#  tstPatternAnalysis.py
#
# ---------------------------------------------------------------------------------

import numpy as np

import tstPatternTools as tpt

# TODO: Add the model / module name on the plots

path = ["/Users/laurent/Data/TestPlan21-perfo/20250321_115547_check_tptAcqMode",
    "/Users/laurent/Data/TestPlan21-perfo/20250321_154151_check_tptAcqMode",
    "/Users/laurent/Data/TestPlan21-perfo/20250321_162115_check_tptAcqMode",
    "/Users/laurent/Data/TestPlan21-perfo/20250409_140759_check_tptAcqMode",
    "/Users/laurent/Data/TestPlan16_validOpalKelly/20250423_172632_check_tptAcqMode"
        ]

# Definition of the pattern
a = -2**15
b = 0
c = 2**6
N = 2**10 - 1
pattern_params = np.array([
    [a, b, c, N], # Slope upward
    [a/2, b, c/2, N], # Slope upward
    [0, 0, 0, 0], # Nothing
    [0, 0, 0, 0], # Nothing
    [0, 0, 0, 0] # Nothing
    ]).astype(int)

col = 0
for p in path[-1:]:
    tpt.check_tst_pattern(p, col, pattern_params)