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
#  calibHKTemp.py
#
# ---------------------------------------------------------------------------------

import glob
import os

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gentools

#path = "/Users/laurent/Data/20.001/test_hk_temp_1/"
#path = "/Users/laurent/Data/20.001/test_hk_temp_2/"
#path = "/Users/laurent/Data/TestPlan24_DM-DMX2/24.010/HK_TEMP/4Col/"
path = "/Users/laurent/Data/TestPlan24_DM-DMX2/24.010/HK_TEMP/Paliers/"


# TODO: Update calibHkTemp to take into account the new format temperature HK (since XifuStudio 2.1.2)
# TODO: Add the model / module name on the plots
def calibHkTemp(path):

    fileName1 = glob.glob(path + "*-HK_TEMP_MAX.csv")[0]
    fileName2 = glob.glob(path + "*-PT104_channel1Value.csv")[0]
    fileName3 = glob.glob(path + "*-HK_TEMP_AVE.csv")[0]
    fileName4 = glob.glob(path + "*-PT104_channel2Value.csv")[0]
    date = fileName1.split("/")[-1][:8]

    # Creation of a directory for the plot files
    path_plot = os.path.join(path, cst.plotDirName)
    gentools.createdir(path_plot)

    plotfilename = os.path.join(path_plot, date + '-calibHkTemp.png')

    print("Reading temperature data from file " + fileName1)
    time1, hk1 = gentools.readHkFromCsv(fileName1)

    print("Reading temperature data from file " + fileName2)
    time2, hk2 = gentools.readHkFromCsv(fileName2)

    print("Reading temperature data from file " + fileName3)
    time3, hk3 = gentools.readHkFromCsv(fileName3)

    print("Reading temperature data from file " + fileName4)
    time4, hk4 = gentools.readHkFromCsv(fileName4)

    # Starting at t=0
    t0 = max(time1[0], time2[0], time3[0], time4[0])
    time1 -= t0
    time2 -= t0
    time3 -= t0
    time4 -= t0
    time_max_12 = min(time1[-1], time2[-1])
    time_max_34 = min(time3[-1], time4[-1])

    hk1 = hk1.astype(float)
    hk3 = hk3.astype(float)

    # resampling data on same times
    hk2_target = np.interp(time1, time2, hk2)
    hk4_target = np.interp(time3, time4, hk4)

    # linear fit of the temperature
    p_AVE = np.polyfit(hk1, hk2_target, 1)
    p_MAX = np.polyfit(hk3, hk4_target, 1)

    hk1_cor = p_AVE[0]*hk1 + p_AVE[1]
    hk3_cor = p_MAX[0]*hk3 + p_MAX[1]

    legend_location = 'lower right'

    fig = plt.figure(figsize=(14, 16))

    title1 = 'Max temperature'
    col1, col2 = 'b', 'r'
    ax1 = fig.add_subplot(3, 2, 1)
    ax1.plot(time1, hk1, color=col1, label='Uncalibrated HK')
    ax1.set_xlim([0, time_max_12])
    ax1.set_title(title1)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Temperature (ADU)", color=col1)
    ax1.plot([], [], color=col2, label='Reference PT104')
    ax1.legend(loc=legend_location)
    ax1.tick_params(axis='y', labelcolor=col1)
    ax2 = ax1.twinx()
    ax2.plot(time1, hk2_target, color=col2)
    ax2.set_ylabel("Temperature (°C)", color=col2)
    ax2.tick_params(axis='y', labelcolor=col2)

    title3 = 'Average temperature'
    col3, col4 = 'g', 'orange'
    ax3 = fig.add_subplot(3, 2, 2)
    ax3.plot(time3, hk3, color=col3, label='Uncalibrated HK')
    ax3.set_xlim([0, time_max_34])
    ax3.set_title(title3)
    ax3.set_xlabel("Time (s)")
    ax3.set_ylabel("Temperature (ADU)", color=col3)
    ax3.plot([], [], color=col4, label='Reference PT104')
    ax3.legend(loc=legend_location)
    ax3.tick_params(axis='y', labelcolor=col3)
    ax4 = ax3.twinx()
    ax4.plot(time3, hk4_target, color=col4)
    ax4.set_ylabel("Temperature (°C)", color=col4)
    ax4.tick_params(axis='y', labelcolor=col4)

    ax5 = fig.add_subplot(3, 2, 3)
    ax5.plot(time1, hk2_target, color=col2, label='Reference PT104')
    if p_AVE[1] > 0:
        sign = '+'
    else:
        sign = ''
    lbl5 = 'Calibrated HK (Eng = {0:6.5f}xADU '.format(p_AVE[0]) + sign + '{0:6.5f})'.format(p_AVE[1])
    ax5.plot(time1, hk1_cor, color='k', label=lbl5)
    ax5.set_xlim([0, time_max_12])
    ax5.set_xlabel("Time (s)")
    ax5.set_ylabel("Temperature (°C)")
    ax5.legend(loc=legend_location)

    ax6 = fig.add_subplot(3, 2, 4)
    ax6.plot(time3, hk4_target, color=col4, label='Reference PT104')
    if p_MAX[1] > 0:
        sign = '+'
    else:
        sign = ''
    lbl6 = 'Calibrated HK (Eng = {0:6.5f}xADU '.format(p_MAX[0]) + sign + '{0:6.5f})'.format(p_MAX[1])
    ax6.plot(time3, hk3_cor, color='k', label=lbl6)
    ax6.set_xlim([0, time_max_34])
    ax6.set_xlabel("Time (s)")
    ax6.set_ylabel("Temperature (°C)")
    ax6.legend(loc=legend_location)

    deltaymax=0.5
    deltaymin=-0.5
    ax7 = fig.add_subplot(3, 2, 5)
    ax7.plot(time1, hk1_cor-hk2_target, color='grey')
    ax7.plot([0, time1[-1]], [deltaymin, deltaymin], '--', color='k', linewidth=2)
    ax7.plot([0, time1[-1]], [deltaymax, deltaymax], '--', color='k', linewidth=2)
    ax7.set_ylabel("Delta (°C)")

    ax8 = fig.add_subplot(3, 2, 6)
    ax8.plot(time3, hk3_cor-hk4_target, color='grey')
    ax8.plot([0, time1[-1]], [deltaymin, deltaymin], '--', color='k', linewidth=2)
    ax8.plot([0, time1[-1]], [deltaymax, deltaymax], '--', color='k', linewidth=2)
    ax8.set_ylabel("Delta (°C)")

    fig.tight_layout()
    plt.savefig(plotfilename, dpi=300, bbox_inches='tight')
    print("plot saved in file ", plotfilename)


if __name__ == '__main__':
    calibHkTemp(path)