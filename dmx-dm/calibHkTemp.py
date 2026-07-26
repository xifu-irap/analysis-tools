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

import os
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gentools
import readData as rddt


@dataclass
class TestConfig:
    session_name: str
    module_name: str

    @property
    def file_path(self) -> str:
        return os.path.join(cst.BASE_DATA_PATH, self.session_name)


def plot_uncal_temp(title, col1, col2, ax, time, hk, hk_target, legend_location):
    ax.plot(time, hk, color=col1, label='Uncalibrated HK')
    ax.set_title(title)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Temperature (ADU)", color=col1)
    ax.plot([], [], color=col2, label='Reference PT104')
    ax.legend(loc=legend_location)
    ax.tick_params(axis='y', labelcolor=col1)
    ax2 = ax.twinx()
    ax2.plot(time, hk_target, color=col2)
    ax2.set_ylabel("Temperature (°C)", color=col2)
    ax2.tick_params(axis='y', labelcolor=col2)


def calib_hk_temp(tconf):
    # Reading the temperature data from the csv files
    path = tconf.file_path
    hk_path = os.path.join(path, cst.hkDirName)

    fileNameDMX_list = [f for f in os.listdir(hk_path) \
                        if os.path.isfile(os.path.join(hk_path, f)) \
                        and f[:8] == 'Hks_DMXA']
    if len(fileNameDMX_list) != 1:
        print("ERROR: wrong number of DMX hk files: {0:}".format(len(fileNameDMX_list)))
        return
    fileNameDMX = os.path.join(hk_path, fileNameDMX_list[0])
    date = fileNameDMX_list[0].split("_")[-1][:8]

    fileNamePt104_list = [f for f in os.listdir(hk_path) \
                          if os.path.isfile(os.path.join(hk_path, f)) \
                          and f[:10] == 'Hks_Pt104_']
    if len(fileNamePt104_list) != 1:
        print("ERROR: wrong number of Pt100 hk files: {0:}".format(len(fileNamePt104_list)))
        return
    fileNamePt104 = os.path.join(hk_path, fileNamePt104_list[0])

    header1 = 'TEMP_MAX(raw)'
    header2 = 'PT104_channel2Value()'
    header3 = 'TEMP_AVE(raw)'
    header4 = 'PT104_channel1Value()'
    header_Date = 'Date(EGSE)'  # The same header in both files

    # Creation of a directory for the plot files
    path_plot = os.path.join(path, cst.plotDirName)
    gentools.createdir(path_plot)

    plotfilename = os.path.join(path_plot, date + '-calibHkTemp.png')

    print("Reading time and temperatures data from DMX hk file " + fileNameDMX)
    t_DMX = np.array(rddt.read_hk_name_from_csv(fileNameDMX, header_Date))
    hk1 = np.array(rddt.read_hk_name_from_csv(fileNameDMX, header1)).astype(float)
    hk3 = np.array(rddt.read_hk_name_from_csv(fileNameDMX, header3)).astype(float)

    print("Reading time and temperatures data from Pt104 hk file " + fileNamePt104)
    t_Pt104 = np.array(rddt.read_hk_name_from_csv(fileNamePt104, header_Date))
    hk2 = np.array(rddt.read_hk_name_from_csv(fileNamePt104, header2))
    hk4 = np.array(rddt.read_hk_name_from_csv(fileNamePt104, header4))

    # Convertir en secondes depuis t0
    time_DMX = (t_DMX - t_DMX[0]) / np.timedelta64(1, 's')
    time_Pt104 = (t_Pt104 - t_DMX[0]) / np.timedelta64(1, 's')

    # resampling data on same times
    hk2_target = np.interp(time_DMX, time_Pt104, hk2)
    hk4_target = np.interp(time_DMX, time_Pt104, hk4)

    # linear fit of the temperature
    p_MAX = np.polyfit(hk1, hk2_target, 1)
    p_AVE = np.polyfit(hk3, hk4_target, 1)

    hk1_cor = p_MAX[0] * hk1 + p_MAX[1]
    hk3_cor = p_AVE[0] * hk3 + p_AVE[1]

    legend_location = 'lower right'

    fig = plt.figure(figsize=(14, 16))
    suptitle = tconf.session_name + "   " + tconf.module_name
    fig.suptitle(suptitle, fontsize=12)

    title1 = 'Max temperature'
    col1, col2 = 'b', 'r'
    ax1 = fig.add_subplot(3, 2, 1)
    plot_uncal_temp(title1, col1, col2, ax1, time_DMX, hk1, hk2_target, legend_location)

    title3 = 'Average temperature'
    col3, col4 = 'g', 'orange'
    ax3 = fig.add_subplot(3, 2, 2)
    plot_uncal_temp(title3, col3, col4, ax3, time_DMX, hk3, hk4_target, legend_location)

    ax5 = fig.add_subplot(3, 2, 3)
    ax5.plot(time_DMX, hk2_target, color=col2, label='Reference PT104')
    if p_AVE[1] > 0:
        sign = '+'
    else:
        sign = ''
    lbl5 = 'Calibrated HK (Eng = {0:6.5f}xADU '.format(p_AVE[0]) + sign + '{0:6.5f})'.format(p_AVE[1])
    ax5.plot(time_DMX, hk1_cor, color='k', label=lbl5)
    ax5.set_xlabel("Time (s)")
    ax5.set_ylabel("Temperature (°C)")
    ax5.legend(loc=legend_location)

    ax6 = fig.add_subplot(3, 2, 4)
    ax6.plot(time_DMX, hk4_target, color=col4, label='Reference PT104')
    if p_MAX[1] > 0:
        sign = '+'
    else:
        sign = ''
    lbl6 = 'Calibrated HK (Eng = {0:6.5f}xADU '.format(p_MAX[0]) + sign + '{0:6.5f})'.format(p_MAX[1])
    ax6.plot(time_DMX, hk3_cor, color='k', label=lbl6)
    ax6.set_xlabel("Time (s)")
    ax6.set_ylabel("Temperature (°C)")
    ax6.legend(loc=legend_location)

    deltaymax=0.5
    deltaymin=-0.5
    ax7 = fig.add_subplot(3, 2, 5)
    ax7.plot(time_DMX, hk1_cor - hk2_target, color='grey')
    ax7.plot([0, time_DMX[-1]], [deltaymin, deltaymin], '--', color='k', linewidth=2)
    ax7.plot([0, time_DMX[-1]], [deltaymax, deltaymax], '--', color='k', linewidth=2)
    ax7.set_ylabel("Delta (°C)")

    ax8 = fig.add_subplot(3, 2, 6)
    ax8.plot(time_DMX, hk3_cor - hk4_target, color='grey')
    ax8.plot([0, time_DMX[-1]], [deltaymin, deltaymin], '--', color='k', linewidth=2)
    ax8.plot([0, time_DMX[-1]], [deltaymax, deltaymax], '--', color='k', linewidth=2)
    ax8.set_ylabel("Delta (°C)")

    fig.tight_layout()
    plt.savefig(plotfilename, dpi=300, bbox_inches='tight')
    print("plot saved in file ", plotfilename)

if __name__ == '__main__':
    conf = TestConfig('..', 'DM-DMX3')
    calib_hk_temp(conf)
