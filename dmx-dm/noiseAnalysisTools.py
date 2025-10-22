# imports
import os

import numpy as np
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
from numpy.typing import ArrayLike
from scipy.optimize import curve_fit

import constants as cst
import general_tools as gt
import readData as rddt


def power_spectrum_from_dumps(data_path, npts, win_name):
    """
    Processes multiple dump files, extracts power spectrum data, and performs normalization.

    This function processes dump files (.fits) located in the specified directory, computes
    the power spectrum for the 4 columns, and averages the spectrum across all the files.
    The computed spectra are normalized with respect to the root-mean-square (RMS) values of
    the signal in both time and frequency domains.

    Args:
        data_path (str): Path to the directory containing the dump files. The files must have
            names starting with 'dump_' and ending with '.fits'.
        npts (int): Number of points for computing the power spectrum. It should be consistent
            with the sampling rate and signal properties.
        win_name (str): Name of the window function to be applied during the computation.

    Raises:
        ValueError: If no valid dump files are found in the specified directory.

    Returns:
        tuple: A tuple containing the following:
            - xf (numpy.ndarray): Frequency bins corresponding to the computed power spectrum.
            - power_spectrum (numpy.ndarray): Averaged and normalized power spectrum.
    """

    from joblib import Parallel, delayed

    files = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[:5] == "dump_" and f[-5:] == ".fits"]

    if len(files) < 1:
        raise ValueError('Wrong number of files')

    print("    Accumulating {0:} DUMP files... ".format(len(files)))

    power_spectrum = np.zeros((cst.nColPerDemux, int(npts / 2) + 1))

    for i in range(len(files)):
        dumpData, _ = rddt.read_dump_file(os.path.join(data_path, files[i]))

        # Computing spectrum (the 4 columns in parallel)
        # xf, power_spectrum_i = gt.do_power_spectrum(dumpData, cst.fSamp, npts, window=win_name)
        results = Parallel(n_jobs=cst.nColPerDemux)(
            delayed(gt.do_power_spectrum)(dumpData[col, :], cst.fSamp, npts, window=win_name) for col in
            range(cst.nColPerDemux))

        # Averaging
        for col_id in range(cst.nColPerDemux):
            power_spectrum[col_id, :] += results[col_id][1]
        xf = results[0][0]

    power_spectrum = power_spectrum / len(files)

    # passage V => ADU
    power_spectrum *= (cst.fsrADCErrorV/cst.fsrADCErrorADU)**2

    return xf, power_spectrum


def power_spectrum_from_1error_column(data_path, col_id, npts, win_name):
    """
    Compute the power spectrum from error data.

    This function processes error data files from the specified directory, performs
    data preprocessing steps such as flattening, DC component removal, and spectral
    analysis, and then computes the normalized power spectrum. The computation
    includes steps for windowing and normalization with respect to the signal's
    root-mean-square (RMS). The processed frequency and power spectrum data are
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

    remove_dc = True
    col_data = rddt.read_science_data_one_col(data_path, col_id, remove_dc, False)

    # Computing spectrum
    xf, power_spectrum = gt.do_power_spectrum(col_data, cst.fRow, npts, window=win_name)

    # passage V => ADU
    power_spectrum *= (cst.fsrADCErrorV/cst.fsrADCErrorADU)**2

    return xf, power_spectrum


def power_spectrum_from_error(data_path, npts, win_name):
    from joblib import Parallel, delayed

    results = Parallel(n_jobs=cst.nColPerDemux)(
        delayed(power_spectrum_from_1error_column)(data_path, col_id, npts, win_name) for col_id in
        range(cst.nColPerDemux))

    power_spectrum = np.zeros((cst.nColPerDemux, int(npts / 2) + 1))

    for col_id in range(cst.nColPerDemux):
        power_spectrum[col_id, :] = results[col_id][1]
    xf = results[0][0]

    return xf, power_spectrum


def low_pass_filter_1(f, a_dc, fc):
    """
    Modèle d'une fonction de transfert pour un filtre passe-bas de premier ordre
    Args:
        f: frequency array
        a_dc: amplitude at DC
        fc: 3dB cutoff frequency

    Returns:

    """
    return a_dc / (1 + f / fc)


def cross_power_spectrum_from_error(data_path, col_ref, npts, win_name):
    """
    Compute the cross-power spectrum from error data.

    This function processes error data files from the specified directory, performs
    data preprocessing steps such as flattening, DC component removal, and spectral
    analysis, and then computes the normalized power spectrum. The computation
    includes steps for windowing and normalization with respect to the signal's
    root-mean-square (RMS). The processed frequency and power spectrum data are
    returned as output.

    Args:
        data_path: Path to the directory containing error data files.
        col_ref: Column identifier to differentiate error files.
        npts: Number of points used for power spectrum computation.
        win_name: Name of the window function to use for spectral analysis.

    Returns:
        Tuple containing:
        - xf: 1D array of frequencies corresponding to the computed power spectrum.
        - cross_power_spectrum: 1D array of the normalized power spectrum values.

    Raises:
        ValueError: If no valid error file matching the specified parameters is found
            in the provided directory.
    """

    cross_power_spectrum = np.zeros((cst.nColPerDemux, int(npts/2+1)))

    remove_dc = True
    ref_data = rddt.read_science_data_one_col(data_path, col_ref, remove_dc)

    for col_id in range(cst.nColPerDemux):
        if col_id == col_ref:
            col_data = ref_data
        else:
            col_data = rddt.read_science_data_one_col(data_path, col_id, remove_dc)

        # Computing cross-spectrum
        xf, cross_power_spectrum[col_id,:] = gt.do_cross_power_spectrum(ref_data, col_data, cst.fRow, npts, window=win_name)

    # passage V => ADU
    cross_power_spectrum *= (cst.fsrADCErrorV/cst.fsrADCErrorADU)**2

    return xf, cross_power_spectrum


def one_over_f(f: ArrayLike, num: ArrayLike) -> np.ndarray:
    """
    Calculate num / sqrt(f), supports scalars, lists, and arrays.
    Always returns a NumPy array.
    """
    return np.asarray(num) / np.sqrt(np.asarray(f))

def fit_one_over_f(x, y):
    """
    Fit a one-over-f model to the provided data.

    This function takes in two arrays, `x` and `y`, and fits a one-over-f
    model to the finite values of `y`. It uses the `curve_fit` function
    from the scipy library to estimate the parameters of the model.

    Parameters:
    ----------
    x : array_like
        The independent variable data (input values).

    y : array_like
        The dependent variable data (output values) to be fitted.
        Only finite values will be considered for fitting.

    Returns:
    -------
    float
        The estimated parameter of the one-over-f model based on the fit
        to the finite values of `y`.

    Raises:
    ------
    ValueError
        If there are no finite values in `y` to fit the model.

    Notes:
    -----
    The function assumes that `one_over_f` is defined elsewhere in the code
    and is the model function used for fitting.
    """

    indexes = np.where(np.isfinite(y))[0]
    params, params_covariance = curve_fit(one_over_f, x[indexes], y[indexes])
    return params[0]


def combine_noises_data(nom_fichier1, fs1, nom_fichier2, fs2, nom_fichier3, npts=2048):
    """
    Combine noise data from two files by resampling and summing quadratically.

    This function loads data from two files, resamples them to the same frequency range
    based on the lower sampling frequency, and combines the noise data quadratically.
    It assumes the input data files are formatted in two columns: frequency and data,
    with one header line to be skipped.

    Args:
        nom_fichier1: Path to the first file containing frequency and data as two columns.
        fs1: Sampling frequency of the first data.
        nom_fichier2: Path to the second file containing frequency and data as two columns.
        fs2: Sampling frequency of the second data.
        nom_fichier3: Output file.
        npts: Number of points for resampling. Default is 2048.

    Returns:
        Tuple containing:
        - freq: Resampled frequencies array.
        - dt: Combined noise data as an array after quadratic summation.
    """

    # Charger les données en ignorant les lignes d'en-tête
    freq1, data1 = gt.read_two_vectors_from_file(nom_fichier1)
    freq2, data2 = gt.read_two_vectors_from_file(nom_fichier2)

    # Keeping the lowest sampling frequency
    fs = min(fs1, fs2)

    # re-sampling data at same frequencies
    freq = np.arange(0, fs/2, fs/2/npts)
    dt1 = np.interp(freq, freq1, data1)
    dt2 = np.interp(freq, freq2, data2)

    # combining noise quadratically
    dt = np.sqrt(dt1**2 + dt2**2)

    gt.save_two_vectors_to_file(freq, dt, nom_fichier3)


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
        new_n (int, optional): The number of equally spaced frequency points for the
            output signal. Default to 2048.

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

    debug = False
    if debug:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        ax.semilogx(frequencies, signal * 1e9)
        ax.semilogx(new_frequencies, np.sqrt(folded_pow) * 1e9)
        ax.set_xlim([1e5, 125e6])

    return new_frequencies, np.sqrt(folded_pow)


def plot_spectrum(ax, xf, spectrum, acq_mode, model_filename, enob=11.5):
    """
    Plot the spectrum and noise floor in both linear and logarithmic scales.

    This function generates a spectrum plot in logarithmic scale,
    displaying the spectrum of a signal and the expected noise floor based on the
    specified Effective Number of Bits (ENOB). The function also supports
    theoretical noise data from a model file and fits a 1/f behavior to the spectrum
    if the acquisition mode is set to 'error'.

    Parameters:
    ----------
    ax : list
        A list of two matplotlib Axes objects for plotting the spectrum and noise floor.

    xf : array_like
        The frequency data (in Hz) corresponding to the spectrum.

    spectrum : array_like
        The spectrum data to be plotted (in nV/√Hz).

    acq_mode : str
        The acquisition mode, which can be either 'dump' or 'error'. This affects the
        limits and the fitting behavior of the function.

    model_filename : str
        The filename of the model data to be used for theoretical noise calculations.
        If an empty string is provided, the model will not be used.

    enob : float, optional
        The Effective Number of Bits for the measurement. Default is 11.5.

    Returns:
    -------
    None
        This function does not return any value. It modifies the provided Axes
        objects directly to display the plots.

    Raises:
    ------
    FileNotFoundError
        If the specified model_filename does not exist.

    Notes:
    -----
    - The function assumes that the constant values and functions such as `cst.fRow`,
      `cst.fSamp`, and `one_over_f` are defined elsewhere in the code.
    - The noise floor is calculated based on the SNR derived from the provided ENOB.
    - The function handles both linear and logarithmic axes for better visualization of
      the spectrum and noise floor.
    """

    xlabel_log: str = r'Frequencies (Hz)'
    ylabel: str = r'Error signal (nV / $\sqrt{Hz}$)'
    xlims_erro: list[int | float] = [1, cst.fRow/2]
    xlims_dump: list[float] = [1e5, cst.fSamp/2]
    ylims_erro: list[float] = [1e1, 1e4]
    ylims_dump: list[float] = [1e1, 1e2]

    # Processing science files
    if acq_mode == 'dump':
        fs = cst.fSamp
        xlims = xlims_dump
        ylims = ylims_dump

    else:  # acq_mode == 'error'
        fs = cst.fRow
        xlims = xlims_erro
        ylims = ylims_erro

    # Computation of the SNR equivalent to the requested ENOB
    snr_db = 6.02 * enob + 10.792
    snr = 10**(snr_db/20)

    # Computation of the noise floor per sqrt(Hz) corresponding to the ENOB
    noise_floor = cst.fsrADCErrorV / (snr * np.sqrt(fs/2))

    if acq_mode == 'error':
        # Fit of 1/f behavior
        xf_start = 3 # to avoid DC perturbations
        f_stop = 1e3 # max frequency of 1/f area
        xf_stop = np.where(xf < f_stop)[0][-1]
        a = fit_one_over_f(xf[xf_start:xf_stop], spectrum[xf_start:xf_stop])

        # Fit of white noise
        f_white_start = 30e3  # to avoid 1/f contribution
        xf_white_start = np.where(xf > f_white_start)[0][0]
        white_noise = spectrum[xf_white_start:].mean()

    if model_filename != '':
        # Getting theoretical noise data from 0 to 125 MHz
        f_theo, noise_theo = gt.read_two_vectors_from_file(model_filename)
        # Aliasing the noise
        f_theo, noise_theo = undersampling_with_aliasing(f_theo, noise_theo, fs, 2 * cst.fSamp)
        print("     Noise from model at low frequencies: {0:6.3f} nV/sqrt(Hz)".format(noise_theo[0] * 1e9))

    lbl1 = "Spectrum"
    ax.loglog(xf[1:], spectrum[1:]*1e9, label=lbl1)

    if model_filename != '':
        lbl2 = 'Model (from datasheets)'
        ax.loglog(f_theo, noise_theo * 1e9, '--', color='r', label=lbl2)

    if acq_mode == 'error':
        x = np.arange(1e4)
        lbl3 = r'White noise level ({0:3.0f}'.format(white_noise * 1e9) + r' nV/$\sqrt{Hz}$)'
        ax.loglog([xf[xf_white_start], xf[-1]], [white_noise * 1e9, white_noise * 1e9], '--', color='k', label=lbl3)
        lbl4 = r'1/f noise ({0:3.1f}'.format(one_over_f(1, a) * 1e6) + r' µV/$\sqrt{Hz}$ at 1 Hz)'
        ax.loglog(x[1:], one_over_f(x[1:], a) * 1e9, '-.', color='k', label=lbl4)

    if enob != 0:
        lbl5 = "Expected noise floor for ENOB={0:}\nand bandwidth = {1:} MHz".format(enob, fs / 2 / 1e6)
        ax.loglog(xlims, [noise_floor * 1e9, noise_floor * 1e9], ':', color='purple', label=lbl5)


    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel_log)

    if acq_mode == 'dump':
        # Désactiver la notation scientifique pour l'axe Y
        ax.yaxis.set_major_formatter(ScalarFormatter())  # Utiliser ScalarFormatter
        ax.yaxis.set_minor_formatter(ScalarFormatter())  # Idem pour les ticks mineurs
        ax.ticklabel_format(style='plain', axis='y')  # Forcer le style "plain" (non scientifique)

        # Configurer les ticks pour afficher uniquement des entiers
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax.grid(True, which='both', linestyle='--')
    ax.legend(loc='best', framealpha=1)


def plot_spectrum2(ax, xf, spectrum, acq_mode, model_filename, enob=11.5):
    """
    Plot the spectrum and noise floor in both linear and logarithmic scales.

    This function generates a spectrum plot in logarithmic scale,
    displaying the spectrum of a signal and the expected noise floor based on the
    specified Effective Number of Bits (ENOB). The function also supports
    theoretical noise data from a model file and fits a 1/f behavior to the spectrum
    if the acquisition mode is set to 'error'.

    Parameters:
    ----------
    ax : list
        A list of two matplotlib Axes objects for plotting the spectrum and noise floor.

    xf : array_like
        The frequency data (in Hz) corresponding to the spectrum.

    spectrum : array_like
        The spectrum data to be plotted (in nV/√Hz).

    acq_mode : str
        The acquisition mode, which can be either 'dump' or 'error'. This affects the
        limits and the fitting behavior of the function.

    model_filename : str
        The filename of the model data to be used for theoretical noise calculations.
        If an empty string is provided, the model will not be used.

    enob : float, optional
        The Effective Number of Bits for the measurement. Default is 11.5.

    Returns:
    -------
    None
        This function does not return any value. It modifies the provided Axes
        objects directly to display the plots.

    Raises:
    ------
    FileNotFoundError
        If the specified model_filename does not exist.

    Notes:
    -----
    - The function assumes that the constant values and functions such as `cst.fRow`,
      `cst.fSamp`, and `one_over_f` are defined elsewhere in the code.
    - The noise floor is calculated based on the SNR derived from the provided ENOB.
    - The function handles both linear and logarithmic axes for better visualization of
      the spectrum and noise floor.
    """

    model_filename = os.path.join('noise_models',
                                  "fdbk-awaxe.txt")  # to subtract this model from erro + feedback measurements
    erro_model_filename = os.path.join('noise_models', "erro-only.txt")  # to overplot the moise model
    xlabel_log: str = r'Frequencies (Hz)'
    ylabel: str = r'Error signal (nV / $\sqrt{Hz}$)'
    xlims_erro: list[int | float] = [1, cst.fRow / 2]
    xlims_dump: list[float] = [1e5, cst.fSamp / 2]
    ylims_erro: list[float] = [1e1, 1e4]
    ylims_dump: list[float] = [1e1, 1e2]

    # Processing science files
    if acq_mode == 'dump':
        fs = cst.fSamp
        xlims = xlims_dump
        ylims = ylims_dump

    else:  # acq_mode == 'error'
        fs = cst.fRow
        xlims = xlims_erro
        ylims = ylims_erro

    # Computation of the SNR equivalent to the requested ENOB
    snr_db = 6.02 * enob + 1.76
    snr = 10 ** (snr_db / 20)

    # Computation of the noise floor per sqrt(Hz) corresponding to the ENOB
    noise_floor = cst.fsrADCErrorV / (snr * np.sqrt(fs / 2))

    if acq_mode == 'error':
        # Fit of 1/f behavior
        xf_start = 3  # to avoid DC perturbations
        f_stop = 1e3  # max frequency of 1/f area
        xf_stop = np.where(xf < f_stop)[0][-1]
        a = fit_one_over_f(xf[xf_start:xf_stop], spectrum[xf_start:xf_stop])

    # Getting theoretical noise data from 0 to 125 MHz
    f_theo, noise_theo = gt.read_two_vectors_from_file(model_filename)
    # Aliasing the noise
    f_theo, noise_theo = undersampling_with_aliasing(f_theo, noise_theo, fs, 2 * cst.fSamp)
    print("     Noise from model at low frequencies: {0:6.3f} nV/sqrt(Hz)".format(noise_theo[0] * 1e9))

    # Getting theoretical noise data from 0 to 125 MHz
    f_theo2, noise_erro = gt.read_two_vectors_from_file(erro_model_filename)
    # Aliasing the noise
    f_theo2, noise_erro = undersampling_with_aliasing(f_theo2, noise_erro, fs, 2 * cst.fSamp)
    print("     Noise from model at low frequencies: {0:6.3f} nV/sqrt(Hz)".format(noise_theo[0] * 1e9))

    lbl1 = "Spectrum"
    lbl2 = "Expected noise floor for ENOB={0:}\nand bandwidth = {1:} MHz".format(enob, fs / 2 / 1e6)
    lbl4 = 'Model of error chain (from datasheets)'

    noise_theo_resamp = np.interp(xf, f_theo, noise_theo)

    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    # ax.loglog(xf[1:], noise_theo_resamp[1:])
    # plt.close()

    spectrum_erro_estimated = np.sqrt(spectrum ** 2 - noise_theo_resamp ** 2)
    ax.loglog(xf[1:], spectrum_erro_estimated[1:] * 1e9, label=lbl1)
    # ax.loglog(xf[1:], spectrum[1:]*1e9, label=lbl1)
    ax.loglog(f_theo2, noise_erro * 1e9, '--', color='r', label=lbl4)

    if acq_mode == 'error':
        x = np.arange(1e4)
        lbl6 = r'1/f noise ({0:3.1f}'.format(one_over_f(1, a) * 1e6) + r' µV/$\sqrt{Hz}$ at 1 Hz)'
        ax.loglog(x[1:], one_over_f(x[1:], a) * 1e9, '-.', color='k', label=lbl6)

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel_log)

    if acq_mode == 'dump':
        # Désactiver la notation scientifique pour l'axe Y
        ax.yaxis.set_major_formatter(ScalarFormatter())  # Utiliser ScalarFormatter
        ax.yaxis.set_minor_formatter(ScalarFormatter())  # Idem pour les ticks mineurs
        ax.ticklabel_format(style='plain', axis='y')  # Forcer le style "plain" (non scientifique)

        # Configurer les ticks pour afficher uniquement des entiers
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax.grid(True, which='both', linestyle='--')
    ax.legend(loc='best', framealpha=1)


def plot_cross_spectrum(ax, xf, spectrum, col_ref):
    """
    Plot the spectrum and noise floor in both linear and logarithmic scales.

    This function generates two plots: one in linear scale and one in logarithmic scale,
    displaying the spectrum of a signal and the expected noise floor based on the
    specified Effective Number of Bits (ENOB). The function also supports
    theoretical noise data from a model file and fits a 1/f behavior to the spectrum
    if the acquisition mode is set to 'error'.

    Parameters:
    ----------
    ax : list
        A list of two matplotlib Axes objects for plotting the spectrum and noise floor.

    xf : array_like
        The frequency data (in Hz) corresponding to the spectrum.

    spectrum : array_like
        The spectrum data to be plotted (in nV/√Hz).

    col_ref : int
        Column id

    Returns:
    -------
    None
        This function does not return any value. It modifies the provided Axes
        objects directly to display the plots.

    Raises:
    ------
    FileNotFoundError
        If the specified model_filename does not exist.

    Notes:
    -----
    - The function assumes that the constant values and functions such as `cst.fRow`,
      `cst.fSamp`, and `one_over_f` are defined elsewhere in the code.
    - The noise floor is calculated based on the SNR derived from the provided ENOB.
    - The function handles both linear and logarithmic axes for better visualization of
      the spectrum and noise floor.
    """

    # normalisation en A.U.
    #spectrum /= spectrum.max()

    xlabel: str = r'Frequencies (Hz)'
    ylabel: str = r'AU'
    xlims_erro: list[int | float] = [1, cst.fRow/2]

    fs = cst.fRow
    xlims = xlims_erro

    colors = ['k', 'b', 'g', 'r']
    # Fit of 1/f behavior
    xf_start = 1 # to avoid DC perturbations
    f_stop = 1e3 # max frequency of 1/f area
    xf_stop = np.where(xf < f_stop)[0][-1]
    x = np.arange(1e4)
    # plotting columns data (col_ref below other columns)
    ax.loglog(xf[1:], spectrum[col_ref, 1:], color=colors[col_ref], label='correlation with column {0:}'.format(col_ref))
    a = fit_one_over_f(xf[xf_start:xf_stop], spectrum[col_ref, xf_start:xf_stop])
    ax.loglog(x[1:], one_over_f(x[1:], a), '-.', color='k')
    for col in range(cst.nColPerDemux):
        if col != col_ref:
            ax.loglog(xf[1:], spectrum[col, 1:], color=colors[col], label='correlation with column {0:}'.format(col))
            a = fit_one_over_f(xf[xf_start:xf_stop], spectrum[col, xf_start:xf_stop])
            ax.loglog(x[1:], one_over_f(x[1:], a), '-.', color='k')

    ax.set_xlim(xlims)
    #ax.set_ylim(ylims)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    ax.grid(True, which='both', linestyle='--')
    ax.legend(loc='best', framealpha=1)

#-------------------------------------------------------------------------------------
