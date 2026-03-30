# imports
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
from numpy.typing import ArrayLike
from scipy.optimize import curve_fit

import constants as cst
import general_tools as gt
import readData as rddt


def power_spectrum_from_dumps(data_path, npts, win_name="none"):
    """
    Processes multiple dump files, extracts power spectrum data, and performs normalization.

    This function processes dump files located in the specified directory, computes
    the power spectrum for the 4 columns, and averages the spectrum across all the files.
    The computed spectra are normalized with respect to the root-mean-square (RMS) values of
    the signal in both time and frequency domains.

    Args:
        data_path (str): Path to the directory containing the dump files. The files must have
            names starting with 'dump_' and ending with '.h5'.
        npts (int): Number of points for computing the power spectrum. It should be consistent
            with the sampling rate and signal properties.
        win_name (str): Name of the window function to be applied during the computation.
            default="none"

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
             and f[:5] == "dump_" and f[-3:] == ".h5"]

    if len(files) < 1:
        raise ValueError('Wrong number of files')

    print("    Accumulating {0:} DUMP files... ".format(len(files)))

    power_spectrum = np.zeros((cst.nColPerDemux, int(npts / 2) + 1))

    for i in range(len(files)):
        dumpData, _ = rddt.read_dump_from_hdf5(os.path.join(data_path, files[i]))

        # Computing spectrum (the 4 columns in parallel)
        # xf, power_spectrum_i = gt.do_power_spectrum(dumpData, cst.fSamp, npts, window=win_name)
        results = Parallel(n_jobs=cst.nColPerDemux)(
            delayed(gt.do_power_spectrum)(dumpData[col, :], cst.fSamp, npts, window=win_name) for col in
            range(cst.nColPerDemux))

        # Accumulating power spectra
        for col_id in range(cst.nColPerDemux):
            power_spectrum[col_id, :] += results[col_id][1]
        xf = results[0][0]

    power_spectrum = power_spectrum / len(files)

    # passage V => ADU
    power_spectrum *= (cst.fsrADCErrorV/cst.fsrADCErrorADU)**2

    return xf, power_spectrum


def power_spectrum_from_1error_column(data_path, col_id, npts, rate, win_name):
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

    if rate == 'FROW':  # Column data
        data, data_exists = rddt.read_col_science_from_dir(data_path, col_id, flatten=True, remove_dc=True,
                                                           verbose=False)

        # Computing spectrum
        if not data_exists:
            xf, power_spectrum = 0, 0
        else:
            xf, power_spectrum = gt.do_power_spectrum(data, cst.fRow, npts, window=win_name)

    elif rate == 'FFRAME':  # Pixels data
        data, data_exists = rddt.read_col_science_from_dir(data_path, col_id, flatten=False, remove_dc=True,
                                                           verbose=False)

        # Computing spectrum
        if not data_exists:
            xf, power_spectrum = 0, 0
        else:
            power_spectrum = np.zeros(int(npts / 2) + 1)
            for pix in range(cst.nPixPerCol):
                xf, power_spectrum_pix = gt.do_power_spectrum(data[pix, :], cst.fFrame, npts, window=win_name)
                power_spectrum += power_spectrum_pix / cst.nPixPerCol

    # passage V => ADU
    power_spectrum *= (cst.fsrADCErrorV / cst.fsrADCErrorADU) ** 2

    return xf, power_spectrum, data_exists


def power_spectrum_from_error(data_path, npts, rate, win_name="none"):
    from joblib import Parallel, delayed

    results = Parallel(n_jobs=cst.nColPerDemux)(
        delayed(power_spectrum_from_1error_column)(data_path, col_id, npts, rate, win_name) for col_id in
        range(cst.nColPerDemux))

    # Construction du vecteur booléen indiquant la présence de données par colonne
    data_exists = [res[2] for res in results]

    power_spectrum = np.zeros((cst.nColPerDemux, int(npts / 2) + 1))

    xf = None
    for col_id in range(cst.nColPerDemux):
        if data_exists[col_id]:
            xf_col, ps_col, _ = results[col_id]
            # Initialisation de xf à partir de la première colonne valide
            if xf is None:
                xf = xf_col
            power_spectrum[col_id, :] = ps_col

    return xf, power_spectrum, data_exists


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


def low_pass_filter_2(f, a_dc, fc, Q=1 / np.sqrt(2)):
    """
    Calcule le module (gain) d'un filtre passe-bas d'ordre 2.

    Paramètres
    ----------
    f : float ou array-like
        Fréquence ou tableau de fréquences (Hz)
    a_dc : amplitude at DC
    fc : float
        Fréquence propre (ou de coupure) en Hz
    Q : float
        Facteur de qualité

    Retour
    ------
    |H| : float ou ndarray
        Module du gain en valeur absolue (non en dB)
    """
    f = np.asarray(f)
    w = 2 * np.pi * f
    w0 = 2 * np.pi * fc

    num = a_dc * w0 ** 2
    den = np.sqrt((w0 ** 2 - w ** 2) ** 2 + (w * w0 / Q) ** 2)
    H = num / den
    return H


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

    ref_data, _ = rddt.read_col_science_from_dir(data_path, col_ref, flatten=True, remove_dc=True, verbose=False)

    for col_id in range(cst.nColPerDemux):
        if col_id == col_ref:
            col_data = ref_data
        else:
            col_data = rddt.read_col_science_from_dir(data_path, col_id, flatten=True, remove_dc=True, verbose=False)

        # Computing cross-spectrum
        xf, cross_power_spectrum[col_id,:] = gt.do_cross_power_spectrum(ref_data, col_data, cst.fRow, npts, window=win_name)

    # passage V => ADU
    cross_power_spectrum *= (cst.fsrADCErrorV/cst.fsrADCErrorADU)**2

    return xf, cross_power_spectrum


def one_over_f(freq: ArrayLike, num: ArrayLike) -> np.ndarray:
    """
    Calculate num / sqrt(f), supports scalars, lists, and arrays.
    Always returns a NumPy array.
    """
    return np.asarray(num) / np.sqrt(np.asarray(freq))


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


def undersampling_with_aliasing(frequencies, signal, new_fs, new_n=2048):
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
    old_fs = frequencies[-1] * 2
    if old_fs % new_fs != 0:
        raise ValueError("old_fs must be a multiple of new_fs. The ratio is {0:}".format(old_fs / new_fs))
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
        ax.set_xlim(1e5, 125e6)

    return new_frequencies, np.sqrt(folded_pow)


def enob_to_nsd(enob, fsr, fs, reference="picpic"):
    """
    Computes the noise spectral density in V/sqrt(Hz) from the ENOB, the FSR the BW

    Args:
        enob (float): The effective number of bits.
        fsr (float): The full scale range (in volts).
        fs (float): The sampling frequency (in Hz).
        reference (string):
            if "picpic" the fsr parameter corresponds to the picpic value of the FSR sine wave
            if "rms" the sfr parameter corresponds to the rms value of the FSR sine wave
            default is "picpic"
    """
    if reference == "picpic":
        snr_db = 6.02 * enob + 10.792
    elif reference == "rms":
        snr_db = 6.02 * enob + 1.76
    else:
        raise ValueError('Wrong reference type. Should be picpic or rms.')
    snr = 10 ** (snr_db / 20)
    nsd = fsr / (snr * np.sqrt(fs / 2))
    return nsd


def plot_spectrum(xf, spectrum, acq_mode, col_id, config, lpf=0, verbose=False):
    """
    Plot the spectrum and noise floor in both linear and logarithmic scales.

    This function generates one_over_f_at_one_hz spectrum plot in logarithmic scale,
    displaying the spectrum of one_over_f_at_one_hz signal and the expected noise floor. The function also supports
    theoretical noise data from one_over_f_at_one_hz model file and fits one_over_f_at_one_hz 1/f behavior to the spectrum
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
        The acquisition mode, which can be either 'DUMP' or 'ERRO'. This affects the
        limits and the fitting behavior of the function.

    model_filename : str
        The filename of the model data to be used for theoretical noise calculations.
        If an empty string is provided, the model will not be used.

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
    plt.rcParams['agg.path.chunksize'] = 200

    xlabel_log: str = r'Frequencies (Hz)'
    ylabel: str = r'Error signal noise (nV / $\sqrt{\mathrm{Hz}}$)'

    plot_path = os.path.join(config["dir_path"], cst.plotDirName)

    plot_file_name_base = 'noise_' + config["setup"] + '_' + config["rate"] + '_BXL' + \
                          "{0:}".format(config["bxl"]) + '_FDBK' + config["fdbk"] + '_OFCO' + config["ofco"]

    col_plot_file_name = plot_file_name_base + '_c{0}'.format(col_id)
    plot_full_file_name = os.path.join(plot_path, col_plot_file_name)

    # Doing the plot
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    suptitle = config["setup"] + ' (Fs = ' + config["rate"] + ' in column {0:})\n'.format(col_id) \
               + ' BXL = {0:}'.format(config["bxl"]) + ', FDBK = ' + config["fdbk"] + ', OFCO = ' + config[
                   "ofco"]

    fig.suptitle(suptitle, fontsize=12)
    ax.set_title(config["session_name"], fontsize=10)

    path_models = 'noise_models'
    match config["rate"]:
        case 'FREF':
            fs = cst.fSamp
            xlims: list[float] = [xf[1], fs / 2]
            ylims: list[float] = [1e1, 1e2]
            model_file_name_end = "_Fref.txt"
            af_text = ""  # Aliasing factor description
        case 'FROW':
            fs = cst.fRow
            xlims: list[int | float] = [xf[1], fs / 2]
            ylims: list[float] = [1e1, 1e4]
            model_file_name_end = "_Frow.txt"
            af_text = r' with $AF=\sqrt{\pi {F_{Ref}}/{F_{Row}}}$'  # Aliasing factor description
        case 'FFRAME':
            fs = cst.fFrame
            xlims: list[int | float] = [xf[1], fs / 2]
            ylims: list[float] = [1e2, 1e4]
            model_file_name_end = "_Fframe.txt"
            af_text = r' with $AF=\sqrt{\pi {F_{Ref}}/{F_{Frame}}}$'  # Aliasing factor description

    match config["setup"]:
        case 'ERRO-ONLY':
            model_file_name_base = "mod-erro-only"
        case 'FDBK-ERRO':
            if col_id == 0 or col_id == 3:
                model_file_name_base = "mod-erro-fdbk-awaxe"
            else:
                model_file_name_base = "mod-erro-fdbk-rhf200"
        case 'OFCO-ERRO':
            if config["ofco"] == "P0V":
                model_file_name_base = "mod-erro-ofco_0"
            else:
                model_file_name_base = "mod-erro-ofco_1023"
        case _:
            model_file_name_base = ''

    if acq_mode[:4] == 'ERRO':
        # Fit of white noise
        f_white_start = 30e3  # to avoid 1/f contribution
        xf_white_start = np.where(xf > f_white_start)[0][0]
        white_noise = spectrum[xf_white_start:].mean()

        # Fit of 1/f behavior
        xf_start = 3 # to avoid DC perturbations
        f_stop = 1e1  # max frequency of 1/f area
        xf_stop = np.where(xf < f_stop)[0][-1]
        one_over_f_at_one_hz = fit_one_over_f(xf[xf_start:xf_stop], spectrum[xf_start:xf_stop])

        # Corner frequency
        f_corner = (one_over_f_at_one_hz / white_noise) ** 2

        # 1/f plus white noise contributions
        one_over_f_plus_white_noise = np.sqrt(
            one_over_f(xf, one_over_f_at_one_hz) ** 2 + np.ones(len(xf)) * white_noise ** 2)

        # Sauvegarde des résultats
        spectra_path = os.path.join(config["dir_path"], cst.spectraDirname)
        spectrum_file_name = os.path.join(spectra_path, config["setup"] + '_' + config["rate"] + "_col{0:}".format(
            col_id) + "_spectrum.txt")
        spectrum_fit_file_name = os.path.join(spectra_path, config["setup"] + '_' + config["rate"] + "_col{0:}".format(
            col_id) + "_spectrum_fit.txt")
        param_fit_file_name = os.path.join(spectra_path, config["setup"] + '_' + config["rate"] + "_col{0:}".format(
            col_id) + "_param_fit.txt")
        gt.save_two_vectors_to_file(xf, spectrum, spectrum_file_name, label1='frequency', label2='noise_req')
        gt.save_two_vectors_to_file(xf, one_over_f_plus_white_noise, spectrum_fit_file_name, label1='frequency',
                                    label2='noise_req')
        gt.save_two_vectors_to_file(np.array([one_over_f_at_one_hz]), np.array([white_noise]), param_fit_file_name,
                                    label1='one_over_f_at_one_hz', label2='white_noise')

    # Computing requirement limits (1/f plus white noise contributions)
    ## Applying the aliasing factor on the white noise
    if fs == cst.fSamp:
        aliasing_factor = 1
    else:
        aliasing_factor = np.sqrt(np.pi * cst.fSamp / fs)

    one_over_f_at_one_hz_req = cst.one_over_f_at_1hz[config["setup"]]
    white_noise_req = cst.nsd[config["setup"]] * aliasing_factor
    f_corner_req = (one_over_f_at_one_hz_req / white_noise_req) ** 2

    # Quadratic sum of 1/f and white noise requirements
    noise_req = np.sqrt(one_over_f(xf, one_over_f_at_one_hz_req) ** 2 + np.ones(len(xf)) * (white_noise_req ** 2))

    # Plot of the measured spectrum
    lbl1 = "Measurement"
    ax.loglog(xf, spectrum * 1e9)

    # Plot of the averaged noise in the band for dump data
    if acq_mode == "DUMP":
        f_range_limits = [500e3, 10e6]
        i_range_limit0 = int(len(xf) * f_range_limits[0] / (cst.fSamp / 2))
        i_range_limit1 = int(len(xf) * f_range_limits[1] / (cst.fSamp / 2))
        noise_avg_nV = spectrum[i_range_limit0:i_range_limit1].mean() * 1e9
        lbl_noise_avg = r'Averaged noise between {0:3.0f} kHz and {1:3.0f} MHz: {2:3.1f} nV/$\sqrt{{\mathrm{{Hz}}}}$' \
            .format(f_range_limits[0] / 1e3, f_range_limits[1] / 1e6, noise_avg_nV)
        ax.loglog([xf[i_range_limit0], xf[i_range_limit1]], [noise_avg_nV, noise_avg_nV], '-', color='k',
                  label=lbl_noise_avg)

    if acq_mode[0:4] == 'ERRO':
        lbl2 = r'Measured white noise: {0:3.0f}'.format(white_noise * 1e9) + r' nV/$\sqrt{{\mathrm{{Hz}}}}$'
        ax.loglog([f_corner, xf[-1]], [white_noise * 1e9, white_noise * 1e9], '--', color='k', label=lbl2)
        lbl3 = r'Measured 1/f noise: {0:3.1f}'.format(
            one_over_f(1, one_over_f_at_one_hz) * 1e6) + r' µV/$\sqrt{Hz}$ at 1 Hz'
        ax.loglog([1, f_corner], [one_over_f_at_one_hz * 1e9, white_noise * 1e9], '-.', color='k', label=lbl3)
        lbl4 = r'Quadratic sum of 1/f and white noise contributions'
        ax.loglog(xf, one_over_f_plus_white_noise * 1e9, '-', linewidth=2, color='k', label=lbl4)

    # Plot of the requirements
    lbl5 = r'White noise req.: {0:3.0f}'.format(white_noise_req * 1e9) + r' nV/$\sqrt{{\mathrm{{Hz}}}}$' + af_text
    ax.loglog([f_corner_req, xf[-1]], [white_noise_req * 1e9, white_noise_req * 1e9], '--', color='r', label=lbl5)
    lbl6 = r'1/f noise req.: {0:3.1f}'.format(
        one_over_f(1, one_over_f_at_one_hz_req) * 1e6) + r' µV/$\sqrt{Hz}$ at 1 Hz'
    ax.loglog([1, f_corner_req], [one_over_f_at_one_hz_req * 1e9, white_noise_req * 1e9], '-.', color='r',
              label=lbl6)
    lbl7 = r'Quadratic sum of both req.'
    ax.loglog(xf, noise_req * 1e9, '-', linewidth=2, color='r', label=lbl7)

    # Plot of the model
    # Getting theoretical noise data from 0 to fs
    if model_file_name_base != '':
        model_file_name = os.path.join(path_models, model_file_name_base + model_file_name_end)
        f_mod, noise_mod = gt.read_two_vectors_from_file(model_file_name)

        lbl8 = 'Model (from datasheets)'
        ax.loglog(f_mod, noise_mod * 1e9, '--', color='purple', label=lbl8)

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel_log)

    if acq_mode == 'DUMP':
        # Désactiver la notation scientifique pour l'axe Y
        ax.yaxis.set_major_formatter(ScalarFormatter())  # Utiliser ScalarFormatter
        ax.yaxis.set_minor_formatter(ScalarFormatter())  # Idem pour les ticks mineurs
        ax.ticklabel_format(style='plain', axis='y')  # Forcer le style "plain" (non scientifique)

        # Configurer les ticks pour afficher uniquement des entiers
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    if lpf != 0:  # comparison with a low pass filter model
        a_dc_estimated = spectrum[4:50].mean()
        popt, pcov = curve_fit(low_pass_filter_1, xf, spectrum, np.array([a_dc_estimated, lpf]))
        a_dc = popt[0]
        fc = popt[1]
        if verbose:
            print("Estimated cutoff frequency: {0:4.0f} MHz".format(fc / 1e6))
            print("Estimated DC level: {0:4.1f} nV/rtHz".format(a_dc * 1e9))

        # Plotting the 1st order LPF fit
        lbl = '1st order LPF fit (fc = {0:4.0f} MHz)'.format(fc / 1e6)
        ax.loglog(xf, low_pass_filter_1(xf, a_dc, fc) * 1e9, '--', color="orange", label=lbl)

    # Afficher la légende
    ax.legend(loc='upper right', fontsize=9)

    ax.grid(True, which='both', linestyle='--')
    fig.tight_layout()
    plt.savefig(plot_full_file_name, dpi=300, bbox_inches='tight')
    plt.close()

    if verbose:
        print("Results plotted in file ", plot_full_file_name)


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
