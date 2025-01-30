# imports
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import zipfile
import readData as rddt
import constants as cst
import general_tools as gt
from cosim.cosim_constants import nb_col


def power_spectrum_from_dumps(data_path, col_id, npts, win_name):
    """
    Processes multiple dump files, extracts power spectrum data, and performs normalization.

    This function processes dump files (.fits) located in the specified directory, computes
    the power spectrum for a specified column, and averages the spectrum across all the files.
    The computed spectrum is normalized with respect to the root mean square (RMS) values of
    the signal in both time and frequency domains.

    Args:
        data_path (str): Path to the directory containing the dump files. The files must have
            names starting with 'dump_' and ending with '.fits'.
        col_id (int): Index of the data column to be selected for power spectrum computation.
        npts (int): Number of points for computing the power spectrum. Should be consistent
            with the sampling rate and signal properties.
        win_name (str): Name of the window function to be applied during the computation.

    Raises:
        ValueError: If no valid dump files are found in the specified directory.

    Returns:
        tuple: A tuple containing the following:
            - xf (numpy.ndarray): Frequency bins corresponding to the computed power spectrum.
            - power_spectrum_total (numpy.ndarray): Averaged and normalized power spectrum.
    """

    files = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[:5] == "dump_" and f[-5:] == ".fits"]

    if len(files) < 1:
        raise ValueError('Wrong number of files')


    print("Processing {0:} DUMP files from ".format(len(files)) + data_path)

    rmsInTime_ADU = 0
    power_spectrum_total = np.zeros(int(npts/2)+1)

    for i in range(len(files)):
        dumpData, _ =rddt.readDumpFile(os.path.join(data_path, files[i]))

        # keeping a single column and removing DC
        dumpData = dumpData[col_id, :] - dumpData[col_id, :].mean()

        # Computing spectrum
        xf, power_spectrum = gt.do_power_spectrum(dumpData, cst.fSamp, npts, window=win_name)

        # Averaging
        rmsInTime_ADU += dumpData.std()
        power_spectrum_total += power_spectrum

    # normalisation wrt signal rms
    rmsInTime_ADU /= len(files)
    rmsInTime_V = rmsInTime_ADU * cst.fsrADCErrorV / cst.fsrADCErrorADU
    rmsInFreq = np.sqrt(power_spectrum_total.sum())
    power_spectrum_total *= (rmsInTime_V / rmsInFreq)**2

    return xf, power_spectrum_total


def power_spectrum_from_error(data_path, col_id, npts, win_name):
    """
    Compute the power spectrum from error data.

    This function processes error data files from the specified directory, performs
    data preprocessing steps such as flattening, DC component removal, and spectral
    analysis, and then computes the normalized power spectrum. The computation
    includes steps for windowing and normalization with respect to the signal's
    root mean square (RMS). The processed frequency and power spectrum data are
    returned as output.

    Args:
        data_path: Path to the directory containing error data files.
        col_id: Column identifier to differentiate error files.
        npts: Number of points used for power spectrum computation.
        win_name: Name of the window function to use for spectral analysis.

    Returns:
        Tuple containing:
        - xf: 1D array of frequencies corresponding to the computed power spectrum.
        - power_spectrum: 1D array of the normalized power spectrum values.

    Raises:
        ValueError: If no valid error file matching the specified parameters is found
            in the provided directory.
    """

    files = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[:6] == 'error_' and f[-6:] == '{0:}.fits'.format(col_id)]

    if len(files) == 0:
        raise ValueError('Wrong number of files')

    fileName = files[0]

    print("Processing ERROR data file from ", data_path)

    colData, _ = rddt.readScienceFile(os.path.join(data_path, fileName))
    # Flattening the array to have data at Frow
    # Dividing by 4 because Error data are in S(16,2) format
    colData = colData.flatten('F') / 4

    # removing DC
    colData -= colData.mean()

    # Computing spectrum
    xf, power_spectrum = gt.do_power_spectrum(colData, cst.fRow, npts, window=win_name)

    # normalisation wrt signal rms
    rmsInTime_ADU = colData.std()
    rmsInTime_V = rmsInTime_ADU * cst.fsrADCErrorV / cst.fsrADCErrorADU
    rmsInFreq = np.sqrt(power_spectrum.sum())
    power_spectrum *= (rmsInTime_V / rmsInFreq)**2

    return xf, power_spectrum

def one_over_f(f, num):
    """
    This function calculates the result of dividing a given value 'a' by the square
    root of 'f'. It is used for the fitting of the spectra at low frequencies.

    Args:
        f (float): The denominator parameter that will be square-rooted as part of the
            division process. Should typically be a positive number.
        num (float): The numerator parameter that will be divided by the square root of 'f'.

    Returns:
        float: The resulting value after dividing 'a' by the square root of 'f'.
    """
    return num / np.sqrt(f)

# Plotting noise spectral density for one column
def plot_col_spectrum(dir_path, col_id, win_name, acq_mode, enob=11.5):
    """
    Plots a column spectrum based on the provided input parameters.

    Processes scientific data files from a specified directory path corresponding to
    a specific column, applies windowing functions, and computes power spectrum
    data. Converts the spectrum to an appropriate scale, computes equivalent
    noise bandwidth, signal-to-noise ratio (SNR), and the noise floor. Optionally
    fits 1/f noise if the acquisition mode is 'error'. Generates and saves plots
    in the specified directory.

    Args:
        dir_path (str): The directory path containing data files for processing.
        col_id (int): The column identifier for the science data to process.
        win_name (str): The window function to apply (e.g., "blackman").
        acq_mode (str): The acquisition mode of the data, either 'dump' or 'error'.
        enob (float, optional): The Effective Number Of Bits (ENOB) for the
            computation of SNR and noise floor. Defaults to 11.5.
    """

    signal = dir_path[-9:]
    xlabel1 = r'Frequencies (MHz)'
    xlabel2 = r'Frequencies (Hz)'
    ylabel = r'Error signal (V / $\sqrt{Hz}$)'
    ylims = [1e-9, 1e-4]

    # Data directory
    pathData = os.path.join(dir_path, cst.dataDirName)

    # Session name
    session_name = os.path.basename(dir_path)

    # Creation of a directory for the plot files
    pathPlot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(pathPlot)

    # Processing science files
    if acq_mode == 'dump':
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum = power_spectrum_from_dumps(pathData, col_id, npts, win_name)
        plotFileName = 'noise_'+signal+'_dumps_c{0:}'.format(col_id)
        fs = cst.fSamp
        xlims2 = [1e5, fs / 2]

    elif acq_mode == 'error':
        npts = 2**23
        xf, power_spectrum = power_spectrum_from_error(pathData, col_id, npts, win_name)
        plotFileName = 'noise_'+signal+'_error_c{0:}'.format(col_id)
        fs = cst.fRow
        xlims2 = [1, fs / 2]
    xlims1 = [0, fs / 2 / 1e6]
    plotFullFileName = os.path.join(pathPlot, plotFileName)

    # Equivalent Noise BandWidth (ENBW)
    if win_name == "blackman":
        ENBW = 1.727
    else:
        ENBW = 1
    rbw = xf[1] * ENBW

    # converting the spectrum to V/sqrt(Hz)
    spectrum = np.sqrt(power_spectrum/rbw)

    # Computation of the SNR equivalent to the requested ENOB
    snr_db = 6.02*enob + 1.76
    snr = 10**(snr_db/20)

    # Computation of the noise floor corresponding to the ENOB
    noise_floor = cst.fsrADCErrorV / (snr * np.sqrt(fs/2))

    # Measuring the noise floor at high frequencies
    noise_floor_hf = spectrum[-100:-2].mean()

    # Fit of 1/f behaviour
    if acq_mode == 'error':
        xf_start = 3 # to avoid DC perturbations
        f_stop = 1e3 # max frequency of 1/f area
        xf_stop = np.where(xf < f_stop)[0][-1]
        params, params_covariance = curve_fit(one_over_f, xf[xf_start:xf_stop], spectrum[xf_start:xf_stop])
        a = params[0]

    # Doing plot
    fig = plt.figure(figsize=(8, 10))
    title = session_name
    ax1 = fig.add_subplot(2, 1, 1)

    lbl1 = "Spectrum".format(rbw/1e3)
    ax1.semilogy(xf/1e6, spectrum, label=lbl1)
    lbl3 = "Expected noise floor for ENOB={0:}\nand bandwidth = {1:} MHz".format(enob, fs/2/1e6)
    ax1.semilogy(xlims1, [noise_floor, noise_floor], ':', color='purple', label=lbl3)
    lbl5 = r'{0:3.1f}'.format(noise_floor_hf*1e9)+r' nV/$\sqrt{Hz}$'
    ax1.semilogy([xf[0]/1e6, xf[-1]/1e6], [noise_floor_hf, noise_floor_hf], '--', color='red', label=lbl5)

    ax1.set_xlim(xlims1)
    ax1.set_ylim(ylims)
    ax1.set_title(title)
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel1)
    ax1.grid()
    ax1.legend(loc='best', framealpha=1)


    ax2 = fig.add_subplot(2, 1, 2)

    ax2.loglog(xf[1:], spectrum[1:], label=lbl1)
    ax2.loglog(xlims2, [noise_floor, noise_floor], ':', color='purple', label=lbl3)
    #ax2.loglog([-1, 1e12], [ref_noise_lvl_nv*1e-9, ref_noise_lvl_nv*1e-9], color='orange', label=lbl4)
    ax2.loglog([xf[0], xf[-1]], [noise_floor_hf, noise_floor_hf], '--', color='red', label=lbl5)
    if acq_mode == 'error':
        x = np.arange(1e4)
        lbl6 = r'1/f noise ({0:3.1f}'.format(one_over_f(1, a) * 1e6) + r' µV/$\sqrt{Hz}$ at 1 Hz)'
        ax2.loglog(x[1:], one_over_f(x[1:], a), '-.', color='k', label=lbl6)

    ax2.set_xlim(xlims2)
    ax2.set_ylim(ylims)
    ax2.set_title(title)
    ax2.set_ylabel(ylabel)
    ax2.set_xlabel(xlabel2)
    ax2.grid()
    ax2.legend(loc='best', framealpha=1)

    fig.tight_layout()

    plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')
    print("Results plotted in file ", plotFullFileName)


def process_list_of_dir(dir_list, win_name, process_dump, process_error):
    """
    Processes a list of directories, extracts data from zip files if required, applies specific
    data processing functions, and optionally cleans up extracted files. The function handles
    two main tasks: unzipping and plotting datasets for the specified directory paths, depending
    on the supplied parameters.

    Args:
        dir_list (list[str]): List of directory paths to process.
        win_name (str): Window name or identifier used in the data processing.
        process_dump (bool): Flag indicating whether to process dump data.
        process_error (bool): Flag indicating whether to process error data.
    """

    for p in dir_list:
        print("Processing data from ", p)

        data_from_zip_file = False
        if process_error:
            # Looking for zip files if any
            dataPath = os.path.join(p, cst.dataDirName)
            zipfiles = [f for f in os.listdir(dataPath) \
                        if os.path.isfile(os.path.join(dataPath, f)) \
                        and f[-4:] == '.zip']

            # Unzipping zip files
            if len(zipfiles) == 0:
                print('There is no zip file...')
            else:
                data_from_zip_file = True
                for z in zipfiles:
                    print('Unzipping Error data from ', z)
                    with zipfile.ZipFile(os.path.join(dataPath, z), 'r') as zip_ref:
                        zip_ref.extractall(dataPath)

        for col in range(nb_col):
            if process_dump:
                plot_col_spectrum(p, col, win_name, "dump")
            if process_error:
                plot_col_spectrum(p, col, win_name, "error")

        # Removing fits files if extracted from a zip file
        if data_from_zip_file:
            os.remove(os.path.join(dataPath, "error_noise_C0.fits"))
            os.remove(os.path.join(dataPath, "error_noise_C1.fits"))
            os.remove(os.path.join(dataPath, "error_noise_C2.fits"))
            os.remove(os.path.join(dataPath, "error_noise_C3.fits"))

#-------------------------------------------------------------------------------------

list_of_dir = [
    "/Users/laurent/Data/TestPlan21-perfo/20250113_170000_noise_erro-only_col3",
    "/Users/laurent/Data/TestPlan21-perfo/20250121_110000_noise_erro-fdbk_col3",
    "/Users/laurent/Data/TestPlan21-perfo/20250124_170723_noise_erro-only",
    "/Users/laurent/Data/TestPlan21-perfo/20250124_171221_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_153842_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_175122_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_175444_noise_erro-only",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_175802_noise_erro-ofco",
    "/Users/laurent/Data/TestPlan21-perfo/20250113_165326_errorEnob_dump-col3"
]

#win = "none"
win = "blackman"

dump = True
error = False
process_list_of_dir(list_of_dir[-1:], win, dump, error)

#-------------------------------------------------------------------------------------
