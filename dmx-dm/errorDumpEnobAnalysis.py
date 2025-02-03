# imports
import numpy as np
import matplotlib.pyplot as plt
import os
import readData as rddt
import constants as cst
import general_tools as gt


def power_spectrum_from_dumps(dataPath, npts, win, remove_dc=True):

    col = 3

    files = [f for f in os.listdir(dataPath) \
             if os.path.isfile(os.path.join(dataPath, f)) \
             and f[:5] == "dump_" and f[-5:] == ".fits"]

    if len(files) < 1:
        raise ValueError('Wrong number of files')


    print("Processing {0:} dump files".format(len(files)))

    spectrum_total = np.zeros(int(npts/2)+1)

    for dump_file in files:
        dumpData, _ =rddt.readDumpFile(os.path.join(dataPath,dump_file))

        # removing DC
        if remove_dc:
            dumpData = dumpData[col] - dumpData[col].mean()

        # Computing spectrum
        xf, spectrum = gt.do_power_spectrum(dumpData, cst.fSamp, npts, window=win)

        spectrum_total += spectrum

    spectrum_db = 10*np.log10(spectrum_total)

    return xf, spectrum_db


# Plotting noise spectral density for one column
def pltNspOneCol(dir_path, win, remove_dc = True, enob=11.5, npts=2**18):

    colId = dir_path[-1]
    ylabel = r'Error signal (dB FSR)'
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
    xf, spectrum_db = power_spectrum_from_dumps(pathData, npts, win, remove_dc=remove_dc)

    # Translation to dbFSR (considering that the sine wave is at full scale)
    # Few spectral channels are considered in case a window has been applied before FFT
    # or in case the line is not well centered on a spectral channel. The extra noise power
    # added doing this is negligible.
    iline = np.where(spectrum_db == spectrum_db[1:].max())[0][0]
    HW = 3
    sine_power = 10*np.log10(np.sum(10**(spectrum_db[iline-HW:iline+HW]/10))) # True sine power

    fsr_vs_sine_power = 20*np.log10(2*np.sqrt(2))
    spectrum_dbfsr = spectrum_db - sine_power - fsr_vs_sine_power

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
    ax1.plot(xf/1e6, spectrum_dbfsr, label=lbl1)
    lbl2 = "Expected SNR for ENOB={0:}".format(enob)
    ax1.plot(xlims, [-snr_db, -snr_db], ':r', label=lbl2)
    lbl3 = "Expected noise floor for ENOB={0:}\n(SNR + FFT gain)".format(enob)
    ax1.plot(xlims, [noise_floor_db, noise_floor_db], ':', color='purple', label=lbl3)
    lbl4 = r"ADC FS / $(2\sqrt{2})$"
    ax1.plot(xlims, [-fsr_vs_sine_power, -fsr_vs_sine_power], '--k', linewidth=0.5, label=lbl4)

    ax1.set_xlim(xlims)
    ax1.set_ylim(ylims)
    ax1.set_title(title)
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)
    ax1.grid()
    ax1.legend(loc='upper right', framealpha=1)

    ax2 = fig.add_subplot(2, 1, 2)
    xlims = [1e5, cst.fSamp/2]
    ylims = [-160, 0]
    xlabel = r'Frequencies (Hz)'

    ax2.semilogx(xf[1:], spectrum_dbfsr[1:], label=lbl1)
    ax2.semilogx(xlims, [-snr_db, -snr_db], ':r', label=lbl2)
    ax2.semilogx(xlims, [noise_floor_db, noise_floor_db], ':', color='purple', label=lbl3)
    ax2.semilogx(xlims, [-fsr_vs_sine_power, -fsr_vs_sine_power], '--k', linewidth=0.5, label=lbl4)

    ax2.set_xlim(xlims)
    ax2.set_ylim(ylims)
    ax2.set_title(title)
    ax2.set_ylabel(ylabel)
    ax2.set_xlabel(xlabel)
    ax2.grid()
    ax2.legend(loc='upper right', framealpha=1)

    fig.tight_layout()

    plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
    print("Results plotted in file ", plotFullFileName)


path = [
    "/Users/laurent/Data/TestPlan21-perfo/20250113_165326_errorEnob_dump-col3"
]
#win = "none"
win = "blackman"
for p in path[:]:
    pltNspOneCol(p, win, remove_dc=True, npts=2*cst.nSamplesPerRow*cst.muxFactor)