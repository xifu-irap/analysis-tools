# imports
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import zipfile
import readData as rddt
import constants as cst
import general_tools as gt
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter


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
            - power_spectrum (numpy.ndarray): Averaged and normalized power spectrum.
    """

    files = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[:5] == "dump_" and f[-5:] == ".fits"]

    if len(files) < 1:
        raise ValueError('Wrong number of files')


    print("Processing {0:} DUMP files from ".format(len(files)) + data_path)

    power_spectrum = np.zeros(int(npts/2)+1)

    for i in range(len(files)):
        dumpData, _ =rddt.readDumpFile(os.path.join(data_path, files[i]))

        # keeping a single column and removing DC
        dumpData = dumpData[col_id, :] - dumpData[col_id, :].mean()

        # Computing spectrum
        xf, power_spectrum_i = gt.do_power_spectrum(dumpData, cst.fSamp, npts, window=win_name)

        # Averaging
        power_spectrum += power_spectrum_i

    power_spectrum = power_spectrum / len(files)

    # passage V => ADU
    power_spectrum *= (cst.fsrADCErrorV/cst.fsrADCErrorADU)**2

    return xf, power_spectrum


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
        raise ValueError('No fits file for column {0:}'.format((col_id)))

    file_name = files[0]

    print("Processing ERROR data file from ", data_path)

    col_data, _ = rddt.readScienceFile(os.path.join(data_path, file_name))
    # Flattening the array to have data at Frow
    # Dividing by 4 because Error data are in S(16,2) format
    col_data = col_data.flatten('F') / 4

    # removing DC
    col_data -= col_data.mean()

    # Computing spectrum
    xf, power_spectrum = gt.do_power_spectrum(col_data, cst.fRow, npts, window=win_name)

    # passage V => ADU
    power_spectrum *= (cst.fsrADCErrorV/cst.fsrADCErrorADU)**2

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


def read_two_vectors_from_file(nom_fichier):
    """
    Reads two vectors from a specified file, ignoring header rows.

    This function loads numerical data from a given file, skipping the first row that is
    presumed to contain headers. The data is expected to be structured into two columns.
    The first column is interpreted as a vector of frequencies, while the second column
    is interpreted as a corresponding vector of noise values. The two vectors are then
    returned separately.

    Args:
        nom_fichier: The path to the file containing the data to be read. The file should
            have a specific structure with two columns of numerical data. The first row
            of the file is ignored as it is assumed to be a header.

    Returns:
        A tuple containing:
            - frequency: A NumPy array representing the first vector (frequency values)
              extracted from the specified file.
            - onoise: A NumPy array representing the second vector (noise values)
              extracted from the specified file.
    """
    # Charger les données en ignorant les lignes d'en-tête
    data = np.loadtxt(nom_fichier, skiprows=1)

    # Séparer les colonnes en deux vecteurs
    frequency = data[:, 0]
    onoise = data[:, 1]

    return frequency, onoise


def undersampling_with_aliasing(frequencies, signal, new_fs, old_fs, new_n = 2048):
    """
    Computes the downsampled and aliased frequencies and amplitudes of a signal when
    undersampling is employed. This function takes a signal sampled at a higher sampling
    frequency and recalculates its frequencies and spectrum at a lower sampling
    frequency with aliasing included.

    Args:
        frequencies (ndarray): The array of original frequency components of the
            signal, corresponding to the sampling frequency `old_fs`.
        signal (ndarray): The amplitude or power values of the signal corresponding
            to the `frequencies`.
        new_fs (float): The desired lower sampling frequency to undersample the
            signal.
        old_fs (float): The original higher sampling frequency of the signal.
        new_n (int, optional): The number of equally-spaced frequency points for the
            output signal. Defaults to 2048.

    Returns:
        tuple: A tuple containing:
            - new_frequencies (ndarray): The new frequency bins after downsampling
              and aliasing, up to `new_fs/2`.
            - folded_pow (ndarray): The aliased and recalculated amplitude or power
              spectrum based on the new sampling frequency.
    """
    # Calculer le facteur entier de sous-échantillonnage
    if old_fs % new_fs != 0:
        raise ValueError("old_fs must be a multiple of new_fs")
    downsampling_factor = int(old_fs // new_fs)

    # re-sampling frequencies and signal at regularly spaced frequencies
    freq = np.arange(0, old_fs/2, old_fs/2/(downsampling_factor*new_n))
    sig = np.interp(freq, frequencies, signal)
    pow = sig**2

    # Aliasing of the signal around new_fs/2
    folded_pow = np.zeros(new_n)
    for i in range(downsampling_factor):
        if i % 2 == 0:
            folded_pow += pow[i*new_n : (i+1)*new_n]
        else:
            folded_pow += np.flip(pow[i*new_n : (i+1)*new_n])

    # New frequency array
    new_frequencies = freq[:new_n]  # Repliement spectral

    return new_frequencies, np.sqrt(folded_pow)


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
    xlabel_lin = r'Frequencies (MHz)'
    xlabel_log = r'Frequencies (Hz)'
    ylabel = r'Error signal (nV / $\sqrt{Hz}$)'
    xlims_erro = [1, cst.fRow/2]
    xlims_dump = [1e5, cst.fSamp/2]
    ylims_erro = [10, 1e4]
    ylims_dump = [10, 100]

    # Selection of noise model
    path_models = 'noise_models'
    model = True # a model exists
    if signal == 'erro-only':
        model_filename = os.path.join(path_models, "erro-only.txt")
    elif signal == 'erro-fdbk':
        model_filename = os.path.join(path_models, "erro-fdbk.txt")
    elif signal == 'erro-ofco':
        model_filename = os.path.join(path_models, "erro-ofco.txt")
    else:
        model = False # file type is unknown, no model exists

    # Data directory
    path_data = os.path.join(dir_path, cst.dataDirName)

    # Session name
    session_name = os.path.basename(dir_path)

    # Creation of a directory for the plot files
    path_plot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(path_plot)

    # Processing science files
    if acq_mode == 'dump':
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum = power_spectrum_from_dumps(path_data, col_id, npts, win_name)
        plot_file_name = 'noise_'+signal+'_dumps_c{0:}'.format(col_id)
        fs = cst.fSamp
        xlims = xlims_dump
        ylims = ylims_dump

    elif acq_mode == 'error':
        npts = 2**23
        xf, power_spectrum = power_spectrum_from_error(path_data, col_id, npts, win_name)
        plot_file_name = 'noise_'+signal+'_error_c{0:}'.format(col_id)
        fs = cst.fRow
        xlims = xlims_erro
        ylims = ylims_erro
    xlims1 = [0, fs / 2 / 1e6]
    plot_full_file_name = os.path.join(path_plot, plot_file_name)

    # converting the spectrum to V/sqrt(Hz)
    rbw = xf[1]
    spectrum = np.sqrt(power_spectrum/rbw)

    # Computation of the SNR equivalent to the requested ENOB
    snr_db = 6.02*enob + 1.76
    snr = 10**(snr_db/20)

    # Computation of the noise floor per sqrt(Hz) corresponding to the ENOB
    noise_floor = cst.fsrADCErrorV / (snr * np.sqrt(fs/2))

    # Measuring the noise floor at high frequencies
    noise_floor_hf = spectrum[-100:-2].mean()

    # Getting theoretical noise data for fs=125MHz
    if model:
        f_theo, noise_theo = read_two_vectors_from_file(model_filename)

    # Fit of 1/f behavior
    if acq_mode == 'error':
        xf_start = 3 # to avoid DC perturbations
        f_stop = 1e3 # max frequency of 1/f area
        xf_stop = np.where(xf < f_stop)[0][-1]
        params, params_covariance = curve_fit(one_over_f, xf[xf_start:xf_stop], spectrum[xf_start:xf_stop])
        a = params[0]

        # Computing theoretical noise data after aliasing at fRow
        if model:
            f_theo, noise_theo = undersampling_with_aliasing(f_theo, noise_theo, cst.fRow, cst.fSamp)

    # Doing plot
    fig = plt.figure(figsize=(8, 10))
    title = session_name
    ax1 = fig.add_subplot(2, 1, 1)

    lbl1 = "Spectrum"
    ax1.semilogy(xf/1e6, spectrum*1e9, label=lbl1)
    lbl2 = "Expected noise floor for ENOB={0:}\nand bandwidth = {1:} MHz".format(enob, fs/2/1e6)
    ax1.semilogy(xlims1, [noise_floor*1e9, noise_floor*1e9], ':', color='purple', label=lbl2)
    #lbl3 = r'{0:3.1f}'.format(noise_floor_hf*1e9)+r' nV/$\sqrt{Hz}$'
    #ax1.semilogy([xf[0]/1e6, xf[-1]/1e6], [noise_floor_hf, noise_floor_hf], '--', color='red', label=lbl3)
    if model:
        lbl4 = 'Model (from datasheets)'
        ax1.semilogy(f_theo/1e6, noise_theo*1e9, '--', color='r', label=lbl4)

    ax1.set_xlim([xlims[0]/1e6, xlims[1]/1e6])
    ax1.set_ylim(ylims)
    ax1.set_title(title)
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel_lin)

    if acq_mode == 'dump':
        # Désactiver la notation scientifique pour l'axe Y
        ax1.yaxis.set_major_formatter(ScalarFormatter())  # Utiliser ScalarFormatter
        ax1.yaxis.set_minor_formatter(ScalarFormatter())  # Idem pour les ticks mineurs
        ax1.ticklabel_format(style='plain', axis='y')  # Forcer le style "plain" (non scientifique)

        # Configurer les ticks pour afficher uniquement des entiers
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax1.grid(True, which='both', linestyle='--')
    ax1.legend(loc='best', framealpha=1)


    ax2 = fig.add_subplot(2, 1, 2)

    ax2.loglog(xf[1:], spectrum[1:]*1e9, label=lbl1)
    ax2.loglog(xlims, [noise_floor*1e9, noise_floor*1e9], ':', color='purple', label=lbl2)
    #ax2.loglog([xf[0], xf[-1]], [noise_floor_hf*1e9, noise_floor_hf*1e9], '--', color='red', linewidth=0.8, label=lbl3)
    if model:
        ax2.loglog(f_theo, noise_theo*1e9, '--', color='r', label=lbl4)

    if acq_mode == 'error':
        x = np.arange(1e4)
        lbl6 = r'1/f noise ({0:3.1f}'.format(one_over_f(1, a) * 1e6) + r' µV/$\sqrt{Hz}$ at 1 Hz)'
        ax2.loglog(x[1:], one_over_f(x[1:], a)*1e9, '-.', color='k', label=lbl6)

    ax2.set_xlim(xlims)
    ax2.set_ylim(ylims)
    ax2.set_title(title)
    ax2.set_ylabel(ylabel)
    ax2.set_xlabel(xlabel_log)

    if acq_mode == 'dump':
        # Désactiver la notation scientifique pour l'axe Y
        ax2.yaxis.set_major_formatter(ScalarFormatter())  # Utiliser ScalarFormatter
        ax2.yaxis.set_minor_formatter(ScalarFormatter())  # Idem pour les ticks mineurs
        ax2.ticklabel_format(style='plain', axis='y')  # Forcer le style "plain" (non scientifique)

        # Configurer les ticks pour afficher uniquement des entiers
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax2.grid(True, which='both', linestyle='--')
    ax2.legend(loc='best', framealpha=1)

    fig.tight_layout()

    plt.savefig(plot_full_file_name, dpi=300, bbox_inches='tight')
    print("Results plotted in file ", plot_full_file_name)


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

        for col in range(cst.nColPerDemux):
        #for col in range(1):
            if process_dump:
                plot_col_spectrum(p, col, win_name, "dump")
            if process_error:
                file_name = os.path.join(dataPath, "error_noise_C{0:}.fits".format(col))
                if os.path.isfile(file_name):
                    plot_col_spectrum(p, col, win_name, "error")
                    # Removing fits files if extracted from a zip file
                    if data_from_zip_file:
                        os.remove(file_name)

#-------------------------------------------------------------------------------------

list_of_dir = [
    "/Users/laurent/Data/TestPlan21-perfo/20250113_170000_noise_erro-only",
    "/Users/laurent/Data/TestPlan21-perfo/20250121_110000_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250124_170723_noise_erro-only",
    "/Users/laurent/Data/TestPlan21-perfo/20250124_171221_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_153842_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_175444_noise_erro-only",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_175802_noise_erro-ofco",
    "/Users/laurent/Data/TestPlan21-perfo/20250130_110005_noise_conf-FPAs"
]

win = "none"
#win = "blackman"

dump = True
error = True
process_list_of_dir(list_of_dir, win, dump, error)

#-------------------------------------------------------------------------------------
