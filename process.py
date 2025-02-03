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
#  laurent.ravera@irap.omp.eu
#  process.py
#

import numpy as np
import matplotlib.pyplot as plt
import os

import general_tools as gentools
import constants as cst
import dmx_tools as dmx
import fakes

plt.rcParams['text.usetex'] = True

path = "."
gentools.purgedir(os.path.join(path, cst.plotDirName))  # Erasing the directory if it already exists

fake_dump = False
fake_error = True

if fake_dump:
    shift = 4
    fileName = 'dump_' + gentools.ma_date() + '_fake.hdf5'
    dump = fakes.fake_dump(path, fileName, shift=shift, sigmaNoise=2, shape='single')
    print(dmx.edge_detect(dmx.read_dump(fileName)))
    dmx.plot_dump(fileName)
    gentools.compute_spectrum(dump, 1 / cst.fSamp, fileName, plot=True)

if fake_error:
    sigmaNoise = 0.5
    window = 'blackmanharris'
    nPatterns = 2 ** 10
    nAcc = 4
    zoomOnLine = True

    FsrADU = 2 ** cst.nBitsAdc

    # rms of full scale square
    #  FSR/2 is the amplitude of the full scale square
    #  4/pi is the amplitude of the first harmonic
    #  1/sqrt(2) gives the rms from the sine wave amplitude
    squareFsrEffLevel = (FsrADU / 2) * (4 / np.pi) / np.sqrt(2)
    squareFsrEffLevelTxt = '{0:d}'.format(round(squareFsrEffLevel))

    sineFsrEffLevel = (FsrADU / (2 * np.sqrt(2)))
    sineFsrEffLevelTxt = '{0:d}'.format(round(sineFsrEffLevel))

    fileName = 'error_' + gentools.ma_date() + '_fake.hdf5'
    error_data = fakes.fake_error(fileName, shape='sine', sigmaNoise=sigmaNoise, nPatterns=nPatterns)
    xf, spectrum = gentools.compute_spectrum(error_data, 1 / cst.fRow, fileName, window=window)
    nPts = len(spectrum)

    for i in range(nAcc - 1):
        fileName = 'error_' + gentools.ma_date() + '_fake.hdf5'
        error_data = fakes.fake_error(fileName, shape='sine', sigmaNoise=sigmaNoise, nPatterns=nPatterns)
        xf, spectrum_i = gentools.compute_spectrum(error_data, 1 / cst.fRow, fileName, window=window)
        spectrum += spectrum_i

    spectrum /= nAcc

    # searching a line
    threshold = 1000
    lineWidth = 4
    lineADU = spectrum.max()
    if lineADU > spectrum.mean() * threshold:
        thereIsALine = True
        iLine = np.where(spectrum == lineADU)[0][0]
        lineADURms = np.sum(spectrum[iLine - lineWidth:iLine + lineWidth + 1])
    else:
        thereIsALine = False
        iLine = nPts // 100  # where to measure the noise if no lines
        lineADURms = 0

    # measuring noise level with an average next to the line
    shift = 2000
    n = 200
    noiseArray = np.concatenate(
        (spectrum[iLine - shift - n: iLine - shift], spectrum[iLine + shift: iLine + shift + n]))
    noiseADU = noiseArray.mean()

    # noise level correction according to the resolution BW
    rbw = ((cst.fRow / 2) / nPts) * gentools.enb[window]
    # print((cst.fRow/2) / nPts, rbw)
    noiseADUPerRootHertz = noiseADU / np.sqrt(rbw)
    noiseADUPerRootHertzdB = 20 * (np.log10(noiseADUPerRootHertz / sineFsrEffLevel))
    noiseADUPerRootHertzdBTxt = '{0:.1f}'.format(noiseADUPerRootHertzdB)

    # SNR = (FsrADU / noiseADUPerRootHertz) / np.sqrt(cst.fRow/2)
    SNR = (sineFsrEffLevel / noiseADUPerRootHertz) / np.sqrt(cst.fRow / 2)
    SNRdB = 20 * np.log10(SNR)
    SNRdBTxt = '{0:.2f}'.format(SNRdB)
    ENOB = (SNRdB - 1.76) / 6.02
    ENOBTxt = '{0:.1f}'.format(ENOB)

    # Plot of accumulated DSP
    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12
    HUGE_SIZE = 16

    plotFileName = os.path.join(path, cst.plotDirName, fileName[:-5] + '.png')
    fig = plt.figure(figsize=(10, 5))
    fig.suptitle('Accumulated spectrum', fontsize=HUGE_SIZE)
    ax = fig.add_subplot(111)
    ax.loglog([1e-5, 10e6], [2 ** 14, 2 ** 14], 'k-', label='FSR ($2^{14}$ ADU)')
    ax.loglog(xf, spectrum, linewidth=1)
    ax.loglog([1e-5, 10e6], [squareFsrEffLevel, squareFsrEffLevel], 'k:', linewidth=0.5,
              label=r'FSR square rms ($\frac{2^{13}}{\sqrt{2}} \frac{4}{\pi}$ = ' + squareFsrEffLevelTxt + ' ADU)')
    ax.loglog([1e-5, 10e6], [sineFsrEffLevel, sineFsrEffLevel], 'k--', linewidth=0.5,
              label=r'FSR sine rms ($\frac{2^{13}}{\sqrt{2}}$ = ' + sineFsrEffLevelTxt + ' ADU)')
    if thereIsALine:
        ax.loglog([xf[iLine]], [lineADU], color='white',
                  label='Line rms (sum over 8 points in the line): {0:d} ADU'.format(round(lineADURms)))
    ax.loglog([1e-5, 10e6], [noiseADU, noiseADU], '-', color='orange', label='Noise level')
    # ax.loglog([1e-5, 10e6], [noiseADUPerRootHertz, noiseADUPerRootHertz], '-', color='r',
    # label=r'Noise level/$\sqrt{\textrm{Hz}} \Longrightarrow$ SNR = '+ SNRdBTxt + r' dB$_{\textrm{FSR}}$')
    ax.loglog([1e-5, 10e6], [noiseADUPerRootHertz, noiseADUPerRootHertz], '-', color='r',
              label=r'Noise level/$\sqrt{\textrm{Hz}}$ = ' + noiseADUPerRootHertzdBTxt + r' dB$_{\textrm{FSR}}$')
    ax.loglog([1e-5], [noiseADUPerRootHertz], '.', color='white',
              label=r'$\Longrightarrow$ SNR = ' + SNRdBTxt + r' dB$_{\textrm{FSR}}$')
    ax.loglog([1e-5], [noiseADUPerRootHertz], '.', color='white',
              label=r'$\Longrightarrow$ ENOB = (SNR - 1.76) / 6.02 = ' + ENOBTxt)
    ax.set_xlim([1, cst.fRow / 2])
    ax.set_title(fileName, fontsize=MEDIUM_SIZE)
    ax.set_xlabel('Frequency (Hz)', fontsize=BIGGER_SIZE)
    ax.set_ylabel('Spectrum', fontsize=BIGGER_SIZE)
    ax.grid(True)
    ax.legend(loc='best', fontsize=MEDIUM_SIZE)
    ax2 = ax.twinx()
    ylim = ax.get_ylim()
    ylim2 = 20 * np.log10(np.array(ylim) / sineFsrEffLevel)
    ax2.set_ylim(ylim2)
    ax2.set_ylabel(r'dB$_{\textrm{FSR}}$', fontsize=BIGGER_SIZE)

    if thereIsALine and zoomOnLine:
        # this is an inset axes over the main axes
        axin1 = ax.inset_axes([0.65, 0.5, 0.3, 0.3])
        axin1.loglog(xf, spectrum, linewidth=1, marker='.', markersize=3)
        axin1.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=False)  # labels along the bottom edge are off
        axin1.set_yticks([])
        ratio = 1.01
        axin1.set_xlim(xf[iLine]/ratio, xf[iLine]*ratio)
        axin1.set_ylim(ylim)

    fig.tight_layout()
    plt.savefig(plotFileName, dpi=300, format='png')

    print("Spectrum plot saved to ", plotFileName)
