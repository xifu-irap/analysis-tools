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
#  dmx_tools.py
#
# ---------------------------------------------------------------------------------

import h5py
import os
import numpy as np
import matplotlib.pyplot as plt
import constants as cst
import general_tools as gentools

# -----------------------------------------------------------------------
def read_dump(fileName):
    """This function reads a DEMUX dump (2 frames of ADC data).

    Args:
        fileName (string): file name including the path.
    Returns:
        dump: array with the dump data.
    """
    fullFileName = os.path.join(cst.dataDirName, fileName)
    with h5py.File(fullFileName, 'r') as fichier:
        dump = fichier['dump'][:]
    return dump


# -----------------------------------------------------------------------
def plot_dump(fileName):
    """This function plots a DEMUX dump (2 frames of ADC data).

    Args:
        fileName (string): file name including the path.
    """
    plotFileName = os.path.join(cst.plotDirName, fileName[:-5]+'.png')

    dump = read_dump(fileName)
    # Trace le tableau
    fig = plt.figure(figsize=(6, 10))
    ax1 = fig.add_subplot(3, 1, 1)
    ax1.plot(dump)
    ax1.set_title("Dump")
    ax1.set_xlabel("Index")
    ax1.set_ylabel("Value")
    ax2 = fig.add_subplot(3, 1, 2)
    ax2.plot(dump[0: int(0.5 * dump.size)])
    ax2.set_title("First frame")
    ax2.set_xlabel("Index")
    ax2.set_ylabel("Value")
    ax3 = fig.add_subplot(3, 1, 3)
    ax3.plot(dump[0: 2 * cst.nSamplesPerRow], marker='.')
    ax3.set_title("2 first pixels")
    ax3.set_xlabel("Index")
    ax3.set_ylabel("Value")

    fig.tight_layout()
    plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
    print("Dump plot saved to " + plotFileName)


# ------------------------------------------------------------
def edge_detect(tab):
    """This function detects a pulse upward or downward in a signal, and
       it measures the position of the rising or falling edge.

    Args:
        tab (numpy array): array .
    """
    trigger = tab.mean()

    ratio_high = 100 * (trigger - tab.min()) / (tab.max() - tab.min())
    pulse_high = ratio_high < 50  # Signal is low with a short pulse upward

    if pulse_high:  # case of a pulse upward
        print("  The signal is low with a pulse upward", end="")
        x_edges = np.flatnonzero((tab[:-1] < trigger) & (tab[1:] > trigger)) + 1
    else:  # case of a pulse downward
        print("  The signal is high with a pulse downward", end="")
        x_edges = np.flatnonzero((tab[:-1] > trigger) & (tab[1:] < trigger)) + 1
    print("(the signal is ~{0:2.1f}% low and ~{1:2.1f}% high)".format(100 - ratio_high, ratio_high))

    if len(x_edges) > 2:
        print(x_edges)
        raise ValueError('  --> The number of Edges in the signal is incorrect (>2)!')

    frame_length = cst.muxFactor * cst.nSamplesPerRow
    x_edge = x_edges[0] - frame_length

    return x_edge

# -----------------------------------------------------------------------
def plotDemuxPattern(param, plotFileName):
    """This function plots a DEMUX test pattern according to the pattern parameters

    Args:
        param (numpy ndarray): Test pattern parameters
            param[region, 0]: parameter a (value of the pattern at the beginning of the region: S(16,2))
            param[region, 1]: parameter b (the number of frames per step is equal to 1+b: U(16, 0))
            param[region, 2]: parameter c (Increment between two successive steps S(16,2))
            param[region, 3]: parameter N (Number of steps in the region U(16, 0))
            with region in [0, 4]
        plotFileName: name of the file to save the plot
    """

    # Checking the size of the pattern array
    if np.size(param[:,0])!=5:
        raise ValueError('Number of regions is incorrect: ' + str(param[0,:].len()) + ' instead of 5')
    if np.size(param[0,:])!=4:
        raise ValueError('Number of parameters is incorrect: ' + str(param[:,0].len()) + ' instead of 4')

    # Computing pattern lenth
    l = 0
    for region in range(np.size(param[:,0])):
        l+=(1+param[region, 1])*param[region, 3]
    x=np.arange(l)

    # Computing pattern
    pattern = np.zeros(l)
    x0 = 0
    for region in range(np.size(param[:,0])):
        a = param[region, 0] / 2 ** 2
        b = param[region, 1]
        c = param[region, 2] / 2 ** 2
        N = param[region, 3]
        for step in range(N):
            pattern[x0+step*(1+b):x0+(step+1)*(1+b)] = (a + step*c)
        x0 += (1+b)*N

    # Doing the plot
    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111)
    ax.plot(x, pattern, linewidth=2)
    xlim = ax.get_xlim()
    ax.plot([xlim[0], xlim[1]],[2**13, 2**13], 'k--', linewidth=0.5)
    ax.plot([xlim[0], xlim[1]],[-2**13, -2**13], 'k--', linewidth=0.5)
    ax.set_title('DEMUX test pattern')
    ax.set_xlabel('Frame')
    ax.set_xlim([xlim[0], xlim[1]])
    ax.set_ylim([-2**13*1.2, 2**13*1.2])
    ax.grid(True)

    ax2 = ax.twiny()
    ax2.set_xlabel('Time (ms)')
    ax2.set_xlim([xlim[0]*1e3/cst.fFrame, xlim[1]*1e3/cst.fFrame])
    fig.tight_layout()
    plt.savefig(plotFileName, dpi=300, format='png', bbox_inches='tight')

# -----------------------------------------------------------------------
def plotRasPattern(param, plotFileName):
    """This function plots a RAS test pattern according to the pattern parameters

    Args:
        param (numpy ndarray): Test pattern parameters
            param[region, 0]: parameter a (value of the pattern at the beginning of the region: U(12,0))
            param[region, 1]: parameter b (the number of frames per step is equal to 1+b: U(16, 0))
            param[region, 2]: parameter c (Increment between two successive steps S(12,0))
            param[region, 3]: parameter N (Number of steps in the region U(16, 0))
            with region in [0, 4]
        plotFileName: name of the file to save the plot
    """

    # Checking the size of the pattern array
    if np.size(param[:,0])!=5:
        raise ValueError('Number of regions is incorrect: ' + str(param[0,:].len()) + ' instead of 5')
    if np.size(param[0,:])!=4:
        raise ValueError('Number of parameters is incorrect: ' + str(param[:,0].len()) + ' instead of 4')

    # Computing pattern lenth
    l = 0
    for region in range(np.size(param[:,0])):
        l+=(1+param[region, 1])*param[region, 3]
    x=np.arange(l)
    print(l, l/cst.fFrame)

    # Computing pattern
    pattern = np.zeros(l)
    x0 = 0
    for region in range(np.size(param[:,0])):
        a = param[region, 0]
        b = param[region, 1]
        c = param[region, 2]
        N = param[region, 3]
        for step in range(N):
            pattern[x0+step*(1+b):x0+(step+1)*(1+b)] = (a + step*c)
        x0 += (1+b)*N

    # Doing the plot
    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111)
    ax.plot(x, pattern, linewidth=2)
    xlim = ax.get_xlim()
    ax.plot([xlim[0], xlim[1]],[2**12, 2**12], 'k--', linewidth=0.5)
    ax.plot([xlim[0], xlim[1]],[0, 0], 'k--', linewidth=0.5)
    ax.set_title('RAS test pattern')
    ax.set_xlabel('Frame')
    ax.set_xlim([xlim[0], xlim[1]])
    ax.set_ylim([-2**12*0.2, 2**12*1.2])
    ax.grid(True)

    ax2 = ax.twiny()
    ax2.set_xlabel('Time (ms)')
    ax2.set_xlim([xlim[0]*1e3/cst.fFrame, xlim[1]*1e3/cst.fFrame])
    fig.tight_layout()
    plt.savefig(plotFileName, dpi=300, format='png', bbox_inches='tight')

# -----------------------------------------------------------------------

# Faking a pulse with a test pattern
pulsePattern = np.array([
    [0x7FFF, 0x0000, -1*0x0400, 0x0020],    # Very steep slope downward
    [-0x0001, 0x01FF, 0x0000, 0x0001],      # Plateau
    [-0x0001, 0x000F, 0x0040, 0x0100],      # Steep slope upward
    [0x3FFF, 0x000F, 0x0010, 0x0400],       # Slow slope upward
    [0x7FFF, 0x7FFF, 0x0000, 0x0001]])      # Plateau

plotFileName = os.path.join('.', cst.plotDirName, 'pattern_'+gentools.maDate()+'.png')
plotDemuxPattern(pulsePattern, plotFileName)
print("Pattern plotted in ", plotFileName)

# Faking a saw with a test pattern
sawPattern = np.array([
    [0x0000, 0x0007, 0x0002, 0x07FF],    # Slope upward
    [0x0FFF, 0x0007, -1*0x0002, 0x07FF],      # Slope downward
    [0x0000, 0x3FFF, 0x0000, 0x0001],      # Plateau
    [0x0000, 0x0020, 0x0002, 0x07FF],    # Slope upward
    [0x0FFF, 0x0020, -1*0x0002, 0x07FF]])      # Slope downward

plotFileName = os.path.join('.', cst.plotDirName, 'pattern_'+gentools.maDate()+'2.png')
plotRasPattern(sawPattern, plotFileName)
print("Pattern plotted in ", plotFileName)