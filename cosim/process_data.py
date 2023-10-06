
#imports
import os
import scd_tools
import matplotlib.pyplot as plt

path='/Users/laurent/Data/TestPlan_10_cosim/'


#----------------------------------------------------------------
# Processing selection  |
#------------------------

#----------------------------------------------------------------

filename_list=['UT_5023/ki_1p0/DRE_DMX_UT_5023_scd', \
            'UT_5023/ki_1p5/DRE_DMX_UT_5023_scd',
            'UT_5023/ki_2p0/DRE_DMX_UT_5023_scd',
            'UT_5023b/ki_2p0/DRE_DMX_UT_5023b_scd']
for filename in filename_list:
    fullfilename = os.path.join(path, filename)
    p = scd_tools.read_scd(fullfilename, verbose=True)
    p.plot(column=0, pixel=0)

#----------------------------------------------------------------
