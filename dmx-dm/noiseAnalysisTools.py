# imports
import numpy as np
import os
import readData as rddt
import constants as cst
import general_tools as gt


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


#-------------------------------------------------------------------------------------
