# imports
import numpy as np
import matplotlib.pyplot as plt
import os
import readData as rddt
import constants as cst
import general_tools as gt


def power_spectrum_from_dumps_noline(dataPath, npts, win):

    col = 3

    files = [f for f in os.listdir(dataPath) \
             if os.path.isfile(os.path.join(dataPath, f)) \
             and f[:5] == "dump_" and f[-5:] == ".fits"]

    if len(files) < 1:
        raise ValueError('Wrong number of files')


    print("Processing {0:} dump files".format(len(files)))

    dumpDataEng = np.zeros((len(files), npts))
    spectrum_total = np.zeros(int(npts/2)+1)

    for i in range(len(files)):
        dumpData, _ =rddt.readDumpFile(os.path.join(dataPath, files[i]))

        # Converting data in engineering values
        dumpDataEng[i,:] = dumpData[col] * cst.fsrADCErrorV / cst.fsrADCErrorADU

        # removing DC
        dumpDataEng[i,:] = dumpDataEng[i,:] - dumpDataEng[i,:].mean()

        # Computing spectrum
        xf, spectrum = gt.do_power_spectrum(dumpDataEng[i,:], cst.fSamp, npts, window=win)
        spectrum_total += spectrum

    # normalisation wrt signal rms and conversion to dbFSR
    rmsInTime = dumpDataEng.std()
    SNR = 2 / rmsInTime
    rmsInFreq = np.sqrt(spectrum_total.sum())
    spectrum_total = spectrum_total / (SNR * rmsInFreq)**2

    spectrum_dbfsr = 10*np.log10(spectrum_total)

    return xf, spectrum_dbfsr


# Plotting noise spectral density for one column
def pltNspOneCol(dir_path, win, enob=11.5, npts=2**18):

    colId = dir_path[-1]
    ylabel = r'Error signal (dB FSR)'
    ylabel2 = r'Noise (V/$\sqrt{Hz}$)'
    ymin = 1e-9
    ymax = 1e-5

    # Data directory
    pathData = os.path.join(dir_path, cst.dataDirName)

    # Session name
    session_name = os.path.basename(dir_path)

    # Creation of a directory for the plot files
    pathPlot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(pathPlot)
    plotFileName = 'errorENOB_col{0:}'.format(colId)
    plotFullFileName = os.path.join(pathPlot, plotFileName)

    # Processing science files
    xf, spectrum_dbfsr = power_spectrum_from_dumps_noline(pathData, npts, win)

    # Resolution bandwidth
    if win == "blackman":
        ENBW = 1.727
    else:
        ENBW = 1
    rbw = xf[1] * ENBW

    # Computation of the SNR equivalent to the requested ENOB
    snr_db = 6.02*enob + 1.76
    noise_floor_db = -(snr_db + 10*np.log10(npts/2))
    if win=='blackman': # Blackman window enlarges the RBW
        noise_floor_db += 2.373

    # Doing plot
    fig = plt.figure(figsize=(8, 10))
    title = session_name
    xlabel = r'Frequencies (MHz)'
    ax1 = fig.add_subplot(2, 1, 1)
    xlims = [0, cst.fSamp/2/1e6]
    ylims = [-160, 0]

    lbl1 = "Spectrum (RBW = {0:.4}kHz)".format(rbw/1e3)
    lns1 = ax1.plot(xf/1e6, spectrum_dbfsr, label=lbl1)
    lbl2 = "Expected SNR for ENOB={0:}".format(enob)
    lns2 = ax1.plot(xlims, [-snr_db, -snr_db], ':r', label=lbl2)
    lbl3 = "Expected noise floor for ENOB={0:}\n(SNR + FFT gain)".format(enob)
    lns3 = ax1.plot(xlims, [noise_floor_db, noise_floor_db], ':', color='purple', label=lbl3)

    ax1.set_xlim(xlims)
    ax1.set_ylim(ylims)
    ax1.set_title(title)
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)
    ax1.grid()

    # Echelle en nV/sqrt(Hz)
    noise_spec_nv = 27
    lbl4 = r'Required noise level: {0:}'.format(noise_spec_nv)+r'nV/$\sqrt{Hz}$'
    ax11 = ax1.twinx()
    lns4 = ax11.semilogy([-1, 1e12], [noise_spec_nv*1e-9, noise_spec_nv*1e-9], color='orange', label=lbl4) # forcer un y en échelle log
    yl_high = cst.fsrADCErrorV / np.sqrt(rbw)
    yl_low = (cst.fsrADCErrorV / np.sqrt(rbw)) / (10**(-ylims[0]/20))
    ax11.set_ylim([yl_low, yl_high])
    ax11.set_ylabel(ylabel2)
    ax11.grid(False)

    # Building a single legend
    lns = lns1 + lns2 + lns3 + lns4
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc='upper left', framealpha=1)

    ax2 = fig.add_subplot(2, 1, 2)
    xlims = [1e5, cst.fSamp/2]
    ylims = [-160, 0]
    xlabel = r'Frequencies (Hz)'

    lns1 = ax2.semilogx(xf[1:], spectrum_dbfsr[1:], label=lbl1)
    lns2 = ax2.semilogx(xlims, [-snr_db, -snr_db], ':r', label=lbl2)
    lns3 = ax2.semilogx(xlims, [noise_floor_db, noise_floor_db], ':', color='purple', label=lbl3)

    ax2.set_xlim(xlims)
    ax2.set_ylim(ylims)
    ax2.set_title(title)
    ax2.set_ylabel(ylabel)
    ax2.set_xlabel(xlabel)
    ax2.grid()

    # Echelle en nV/sqrt(Hz)
    ax22 = ax2.twinx()
    ax22.semilogy([-1, 1e12], [noise_spec_nv*1e-9, noise_spec_nv*1e-9], color='orange', label=lbl4) # forcer un y en échelle log
    ax22.set_ylim([yl_low, yl_high])
    ax22.set_ylabel(ylabel2)
    ax22.grid(False)

    # Building a single legend
    lns = lns1 + lns2 + lns3 + lns4
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc='upper left', framealpha=1)

    fig.tight_layout()

    plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
    print("Results plotted in file ", plotFullFileName)


path = [
    "/Users/laurent/Data/TestPlan21-perfo/20250113_171800_errorEnob-col3_dump_noline"
]
#win = "none"
win = "blackman"
for p in path[:]:
    pltNspOneCol(p, win, npts=2*cst.nSamplesPerRow*cst.muxFactor)