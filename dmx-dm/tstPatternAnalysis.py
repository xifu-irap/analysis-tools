import numpy as np
import tstPatternTools as tpt


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