# imports
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import os
import readData as rddt
import constants as cst
import general_tools as gt


def power_spectrum_from_dumps_noline(dataPath, colId, npts, win):

    files = [f for f in os.listdir(dataPath) \
             if os.path.isfile(os.path.join(dataPath, f)) \
             and f[:5] == "dump_" and f[-5:] == ".fits"]

    if len(files) < 1:
        raise ValueError('Wrong number of files')


    print("Processing {0:} DUMP files from ".format(len(files))+dataPath)

    dumpDataEng = np.zeros((len(files), npts))
    spectrum_total = np.zeros(int(npts/2)+1)

    for i in range(len(files)):
        dumpData, _ =rddt.readDumpFile(os.path.join(dataPath, files[i]))

        # Converting data in engineering values
        dumpDataEng[i,:] = dumpData[colId] * cst.fsrADCErrorV / cst.fsrADCErrorADU

        # removing DC
        dumpDataEng[i,:] = dumpDataEng[i,:] - dumpDataEng[i,:].mean()

        # Computing spectrum
        xf, spectrum = gt.do_power_spectrum(dumpDataEng[i,:], cst.fSamp, npts, window=win)
        spectrum_total += spectrum

    # normalisation wrt signal rms and conversion to dbFSR
    rmsInTime = dumpDataEng.std()
    SNR = 2 / rmsInTime
    rmsInFreq = np.sqrt(spectrum_total.sum())
    power_spectrum = spectrum_total / (SNR * rmsInFreq)**2

    return xf, power_spectrum


def power_spectrum_from_error_noline(dataPath, npts, win):

    files = [f for f in os.listdir(dataPath) \
             if os.path.isfile(os.path.join(dataPath, f)) \
             and f[-5:] == ".fits"]

    if len(files) != 1:
        raise ValueError('Wrong number of files')

    fileName = files[0]

    print("Processing ERROR data file from ", dataPath)

    colData, _ = rddt.readScienceFile(os.path.join(dataPath, fileName))
    # Flattening the array to have data at Frow
    # Dividing by 4 because Error data are in S(16,2) format
    colData = colData.flatten('F') / 4

    # removing DC
    colData -= colData.mean()

    # Converting data in engineering values
    colDataEng = colData * cst.fsrADCErrorV / cst.fsrADCErrorADU

    # Computing spectrum
    xf, power_spectrum = gt.do_power_spectrum(colDataEng, cst.fRow, npts, window=win)

    # normalisation wrt signal rms and conversion to dbFSR
    rmsInTime = colDataEng.std()
    SNR = 2 / rmsInTime
    rmsInFreq = np.sqrt(power_spectrum.sum())
    power_spectrum = power_spectrum / (SNR * rmsInFreq)**2

    return xf, power_spectrum


# Plotting noise spectral density for one column
def pltNspOneCol(dir_path, win, enob=11.5):

    colId = dir_path[-1]
    xlabel = r'Frequencies (MHz)'
    ylabel = r'Error signal (V / $\sqrt{Hz}$)'
    ref_noise_lvl_nv = 27
    ylims = [1e-9, 1e-5]

    # Data directory
    pathData = os.path.join(dir_path, cst.dataDirName)

    # Session name
    session_name = os.path.basename(dir_path)

    # Creation of a directory for the plot files
    pathPlot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(pathPlot)

    # Processing science files
    if dir_path[-16:-5] == 'dump_noline':
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum = power_spectrum_from_dumps_noline(pathData, int(colId), npts, win)
        plotFileName = 'noise_error_dump_col{0:}'.format(colId)
        fs = cst.fSamp
        xlims2 = [1e5, fs / 2]
    elif dir_path[-17:-5] == 'error_noline':
        npts = 2**22
        xf, power_spectrum = power_spectrum_from_error_noline(pathData, npts, win)
        plotFileName = 'noise_error_error_col{0:}'.format(colId)
        fs = cst.fRow
        xlims2 = [1, fs / 2]
    xlims1 = [0, fs / 2 / 1e6]
    plotFullFileName = os.path.join(pathPlot, plotFileName)

    # Resolution bandwidth
    if win == "blackman":
        ENBW = 1.727
    else:
        ENBW = 1
    rbw = xf[1] * ENBW
    # converting the spectrum to V / sqrt(Hz)
    spectrum = np.sqrt(power_spectrum/rbw)

    # Computation of the SNR equivalent to the requested ENOB
    snr_db = 6.02*enob + 1.76
    snr = 10**(snr_db/20)

    # Computation of noise floor corresponding to the ENOB
    noise_floor = cst.fsrADCErrorV / (snr * np.sqrt(fs/2))

    # Doing plot
    fig = plt.figure(figsize=(8, 10))
    title = session_name
    ax1 = fig.add_subplot(2, 1, 1)

    lbl1 = "Spectrum".format(rbw/1e3)
    ax1.semilogy(xf/1e6, spectrum, label=lbl1)
    lbl3 = "Expected noise floor for ENOB={0:}\nand bandwidth = {1:} MHz".format(enob, fs/2/1e6)
    ax1.semilogy(xlims1, [noise_floor, noise_floor], ':', color='purple', label=lbl3)
    lbl4 = r'{0:}'.format(ref_noise_lvl_nv)+r' nV/$\sqrt{Hz}$'
    ax1.semilogy([-1, 1e12], [ref_noise_lvl_nv*1e-9, ref_noise_lvl_nv*1e-9], color='orange', label=lbl4) # forcer un y en échelle log

    ax1.set_xlim(xlims1)
    ax1.set_ylim(ylims)
    ax1.set_title(title)
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)
    ax1.grid()
    ax1.legend(loc='best', framealpha=1)


    ax2 = fig.add_subplot(2, 1, 2)

    ax2.loglog(xf[1:], spectrum[1:], label=lbl1)
    ax2.loglog(xlims2, [noise_floor, noise_floor], ':', color='purple', label=lbl3)
    ax2.loglog([-1, 1e12], [ref_noise_lvl_nv*1e-9, ref_noise_lvl_nv*1e-9], color='orange', label=lbl4) # forcer un y en échelle log

    ax2.set_xlim(xlims2)
    ax2.set_ylim(ylims)
    ax2.set_title(title)
    ax2.set_ylabel(ylabel)
    ax2.set_xlabel(xlabel)
    ax2.grid()
    ax2.legend(loc='best', framealpha=1)

    fig.tight_layout()

    plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
    print("Results plotted in file ", plotFullFileName)



path = [
    "/Users/laurent/Data/TestPlan21-perfo/20250113_171800_errorEnob_dump_noline-col3",
    "/Users/laurent/Data/TestPlan21-perfo/20250113_171908_errorEnob_error_noline-col3"
]

#win = "none"
win = "blackman"
for p in path[:]:
    pltNspOneCol(p, win)