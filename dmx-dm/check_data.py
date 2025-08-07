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
#  check_data.py
#  This python script is used to check demux data, either by exporting fits file
#  to txt files or by checking the value of each bits along a test.
# ---------------------------------------------------------------------------------

import os

import matplotlib.pyplot as plt

import constants as cst
import general_tools as gtools
import readData as rddt


# Exporting data from fits files to txt files

# rddt.export_science_data_one_col_2_txt("/Users/laurent/Data/TestPlan24_DM-DMX2_elec/debug_com_ADC/20250804_error_data", 0)
# rddt.export_dump_2_txt("/Users/laurent/Data/TestPlan24_DM-DMX2_elec/debug_com_ADC/dump_20250804-150321.fits", 3)
# rddt.export_science_data_one_col_2_txt("/Users/laurent/Data/TestPlan24_DM-DMX2_elec/debug_com_ADC/20250805", 1)


def check_bits(path_data, data_types, column):
    if data_types == "scan":
        files = [f for f in os.listdir(path_data) \
                 if os.path.isfile(os.path.join(path_data, f)) \
                 and f[-5:] == ".fits" and f[:5] == "scan_"]
    elif data_types == "dump":
        files = [f for f in os.listdir(path_data) \
                 if os.path.isfile(os.path.join(path_data, f)) \
                 and f[-5:] == ".fits" and f[:5] == "dump_"]

    path_plot = os.path.join(path_data, cst.plotDirName)
    gtools.createdir(path_plot)

    for file in files:
        print("processing file: ", file)
        full_file_name = os.path.join(path_data, file)
        plot_full_file_name = os.path.join(path_plot, file[:-5] + '.png')
        if data_types == "scan":
            _, _, _, errors = rddt.read_scan(full_file_name)
            data = errors.flatten('F').astype(int)
        elif data_types == "dump":
            data, _ = rddt.read_dump_file(full_file_name)
            data = data[column, :]

        fig = plt.figure(figsize=(10, 12))
        fig.suptitle(file, fontsize=12)

        ax0 = fig.add_subplot(5, 1, 1)
        ax0.set_title("dump col {0:}".format(column), fontsize=10)
        ax0.plot(data, 'k')
        ax0.grid(True)

        for bit_nb in range(14):
            bit_val = gtools.bit_value(data, bit_nb)

            ax = fig.add_subplot(5, 4, 5 + bit_nb)
            ax.set_title("bit {0:}".format(bit_nb), fontsize=10)
            ax.plot(bit_val)
            ax.set_ylim(-0.5, 1.5)

        fig.tight_layout()
        plt.savefig(plot_full_file_name, dpi=300, bbox_inches='tight')
        plt.close()


path = "/Users/laurent/Data/TestPlan24_DM-DMX2_elec/debug_com_ADC/"

dir = "20250806_fw11"
column = 3
check_bits(os.path.join(path, dir), "dump", column)
check_bits(os.path.join(path, dir), "scan", column)

dir = "20250806_fw_test_adc3"
column = 0
check_bits(os.path.join(path, dir), "dump", column)
check_bits(os.path.join(path, dir), "scan", column)
