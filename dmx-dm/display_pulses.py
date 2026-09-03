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
#  display_pulses.py
#
# ---------------------------------------------------------------------------------

# imports
import os

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gt
import readData as rddt


def plot_pulse(pulse, plotFileName, suptitle, title, ylabel, verbose=True):
    xlabel = 'Sample'

    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(1, 1, 1)  # global plot
    ax1.plot(pulse)

    plt.suptitle(suptitle)
    ax1.set_title(title, fontsize=8)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)

    fig.tight_layout()

    plt.savefig(plotFileName, dpi=300, bbox_inches='tight')

    if verbose:
        print("results plotted in file " + plotFileName)


# ---------------------------------------------------------------------------------

def get_pulse(filename, pix=0, trig_level=10):
    nPointsBefore = 100
    record_length = 1024

    col_data = rddt.read_science_from_file(filename)
    pix_data = col_data[pix, :]

    # searching a pulse
    avg = pix_data.mean()
    std = pix_data.std()
    trig = trig_level * std
    index = np.where(abs(pix_data - avg) > trig)[0][0]
    if index < nPointsBefore:
        index = np.where(abs(pix_data - avg) > trig_level)[0][1]

    start = max(index - nPointsBefore, 0)
    end = min(index + record_length, len(pix_data))
    return (pix_data[start: end])


# ---------------------------------------------------------------------------------

def display_pulses(verbose):
    # Paths definition
    dir_path = os.path.join("..", "..")
    hk_path = os.path.join(dir_path, cst.hkDirName)
    data_path = os.path.join(dir_path, cst.dataDirName)
    plot_path = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(plot_path)
    session_name = os.path.basename(os.path.realpath(dir_path))
    col = session_name[-1]

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(hk_path)

    if verbose:
        print("/----------------------------------------------------------")
        print("/ Pulses measurement")
        print("/ Test session name:   " + session_name)
        print("/ Tested column:       " + col)
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:     {0:}".format(fwVersion))
        print("/----------------------------------------------------------\n")
    title = '(session ' + session_name + ')'

    suptitle1 = 'Pulse measurements with FLL OFF, column {0:}'.format(col)
    plotFileName1 = 'pulses_fll_off_error_col{0:}.png'.format(col)
    fullPlotFileName1 = os.path.join(plot_path, plotFileName1)
    ylabel1 = 'Error data (ADU)'

    start1 = "pulses_FLL_OFF_error"
    files = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[-4:] == '{0:}.h5'.format(col) and f[:len(start1)] == start1]

    if len(files) != 1:
        raise ValueError(
            'Wrong number of files ({0:} instead of {1:})'.format(len(files), 1))

    pulse = get_pulse(os.path.join(data_path, files[0]))
    plot_pulse(pulse, fullPlotFileName1, suptitle1, title, ylabel1, verbose)

    suptitle2 = 'Pulse measurements with FLL ON, column {0:}'.format(col)
    plotFileName2 = 'pulses_fll_on_error_col{0:}.png'.format(col)
    fullPlotFileName2 = os.path.join(plot_path, plotFileName2)
    ylabel2 = 'Error data (ADU)'

    start2 = "pulses_FLL_ON_error"
    files = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[-4:] == '{0:}.h5'.format(col) and f[:len(start2)] == start2]

    if len(files) != 1:
        raise ValueError(
            'Wrong number of files ({0:} instead of {1:})'.format(len(files), 1))

    pulse = get_pulse(os.path.join(data_path, files[0]))
    plot_pulse(pulse, fullPlotFileName2, suptitle2, title, ylabel2, verbose)

    suptitle3 = 'Pulse measurements with FLL ON, column {0:}'.format(col)
    plotFileName3 = 'pulses_fll_on_science_col{0:}.png'.format(col)
    fullPlotFileName3 = os.path.join(plot_path, plotFileName3)
    ylabel3 = 'Science data (ADU)'

    start3 = "pulses_FLL_ON_science"
    files = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[-4:] == '{0:}.h5'.format(col) and f[:len(start3)] == start3]

    if len(files) != 1:
        raise ValueError(
            'Wrong number of files ({0:} instead of {1:})'.format(len(files), 1))

    pulse = get_pulse(os.path.join(data_path, files[0]))
    plot_pulse(pulse, fullPlotFileName3, suptitle3, title, ylabel3, verbose)

# ---------------------------------------------------------------------------------
