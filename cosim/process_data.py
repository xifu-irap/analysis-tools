
#imports
import os

import scd_tools as stl
import constants as cst

path='/Users/laurent/Data/TestPlan_10_cosim/'
col = 0 # column_index

#----------------------------------------------------------------
# Processing selection  |
#------------------------

redo = True
do_x_talk = True

#----------------------------------------------------------------

print(cst.draw_line)
if redo:
    filename_list=['UT_5023/ki_1p0/DRE_DMX_UT_5023_scd', \
                'UT_5023/ki_1p5/DRE_DMX_UT_5023_scd',
                'UT_5023/ki_2p0/DRE_DMX_UT_5023_scd',
                'UT_5023b/ki_2p0/DRE_DMX_UT_5023b_scd',
                'UT_5023c/DRE_DMX_UT_5023c_scd',
                'UT_5023d/DRE_DMX_UT_5023d_scd'
                ]
    for filename in filename_list:
        fullfilename = os.path.join(path, filename)
        p = stl.read_tm(fullfilename, verbose=True)
        p.plot(column=col, pixel=0)
        p.plot(column=col)

#----------------------------------------------------------------
if do_x_talk:
    print(cst.draw_line)
    print("Processing cross-talk simulations")
    frame_baseline = 600 # Where to measure the baseline level
    pixel_pulse = 2
    pixel_victim = 3
    filename_list = ['UT_5030a/DRE_DMX_UT_5030a_scd',
                     'UT_5030b/DRE_DMX_UT_5030b_scd',
                     'UT_5030c/DRE_DMX_UT_5030c_scd',
                     'UT_5030d/DRE_DMX_UT_5030d_scd',
                     'UT_5030e/DRE_DMX_UT_5030e_scd'
                     ]
    for filename in filename_list:
        fullfilename = os.path.join(path, filename)
        p = stl.read_tm(fullfilename, verbose=True)
        p.plot(column=col)
        p.plot(column=col, xy_limits=[0, 0, 24500, 27500])
        main_pulse = p.data[col][frame_baseline:, pixel_pulse].min() - p.data[col][frame_baseline, pixel_pulse]
        xtalk_pulse = p.data[col][frame_baseline:, pixel_victim].max() - p.data[col][frame_baseline, pixel_victim]
        x_talk = 100 * xtalk_pulse / main_pulse
        print("  Xtalk level: {0:3.3f} %".format(x_talk))

print(cst.draw_line)

#----------------------------------------------------------------
