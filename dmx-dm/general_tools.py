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
#  laurent.ravera@cnrs.fr
#  general_tools.py
#
# ---------------------------------------------------------------------------------

import os
import shutil
from datetime import datetime
from typing import Tuple, Union

import matplotlib.pyplot as plt
import numpy as np

import constants as cst

# definition of equivalent noise bandwidth for some well known windows
enb = {'boxcar': 1.0, 'hann': 1.5, 'hamming': 1.363, 'blackman': 1.727,
       'blackmanharris': 2.004, 'barlett': 1.333, 'flattop': 3.77,
       'triang': 1.333}


# -----------------------------------------------------------------------
def checkdir(dirname: str) -> None:
    """Check if a directory exists, creating it if it does not.

    Parameters
    ----------
    dirname : str
        Name of the directory to be verified / created.

    Returns
    -------
    None
    """
    if not os.path.isdir(dirname):
        os.mkdir(dirname)


# -----------------------------------------------------------------------
def purgedir(dirname: str) -> None:
    """Erase a directory if it exists, then create a new empty directory.

    Parameters
    ----------
    dirname : str
        Name of the directory to be erased / created.

    Returns
    -------
    None
    """
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)
    os.mkdir(dirname)


# -----------------------------------------------------------------------
def createdir(dirname: str) -> None:
    """Create a directory if it does not exist.

    Parameters
    ----------
    dirname : str
        Name of the directory to be created.

    Returns
    -------
    None
    """
    if not os.path.isdir(dirname):
        os.mkdir(dirname)


# -----------------------------------------------------------------------
def compute_spectrum(signal: np.ndarray, period: float, file_name: str,
                     window: str = 'blackman',
                     plot: bool = False) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the amplitude spectrum of a signal and optionally plot it.

    Parameters
    ----------
    signal : numpy.ndarray
        Input time series.
    period : float
        Sampling period in seconds.
    file_name : str
        Base file name used to derive the plot file name.
    window : str, optional
        Windowing function to apply ('blackman', 'blackmanharris', 'hann').
        Default is 'blackman'.
    plot : bool, optional
        If True, save the spectrum plot to disk. Default is False.

    Returns
    -------
    xf : numpy.ndarray
        Frequencies of the spectrum.
    fft_norm : numpy.ndarray
        Normalized amplitude spectrum in ADU.
    """
    from scipy.signal.windows import blackman, blackmanharris, hann

    plot_file_name = os.path.join(cst.plotDirName,
                                  file_name[:-5] + '_sp.png')
    n_pts = len(signal)

    # create the selected window
    if window == 'blackman':
        w = blackman(n_pts)
    elif window == 'blackmanharris':
        w = blackmanharris(n_pts)
    elif window == 'hann':
        w = hann(n_pts)
    else:
        w = np.ones(n_pts)

    # compute the FFT
    fft_result = np.fft.rfft(signal * w)

    # normalize the result in ADU
    fft_norm = np.abs(fft_result) / n_pts * np.sqrt(2)

    # adjust normalization for the DC and Nyquist components
    fft_norm[0] /= np.sqrt(2)
    if len(signal) % 2 == 0:
        fft_norm[-1] /= np.sqrt(2)

    # create the frequency axis
    xf = np.linspace(0, 0.5 / period, len(fft_norm))

    # display the spectrum
    if plot:
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111)
        ax.loglog(xf, fft_norm)
        ax.set_title('Spectrum')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Amplitude (ADU)')
        ax.grid(True)
        fig.tight_layout()
        plt.savefig(plot_file_name, dpi=300, format='png',
                    bbox_inches='tight')

    return xf, fft_norm


# -----------------------------------------------------------------------
def is_even(x: int) -> bool:
    """Check if a number is even.

    Parameters
    ----------
    x : int
        The value to be verified.

    Returns
    -------
    bool
        True if the input is even, False if the input is odd.
    """
    return x % 2 == 0


# -----------------------------------------------------------------------
def first_order_low_pass_filter(data_vector: np.ndarray,
                                tau: float) -> np.ndarray:
    """Apply a first-order low-pass filter to a data vector.

    Parameters
    ----------
    data_vector : numpy.ndarray
        The input data vector.
    tau : float
        Time constant of the filter (in seconds).

    Returns
    -------
    numpy.ndarray
        The filtered vector.
    """
    dt = 1.0 / cst.fSamp  # time interval
    alpha = dt / (tau + dt)  # filter coefficient

    filtered_vector = np.zeros_like(data_vector)
    filtered_vector[0] = data_vector[0]  # initial condition

    for i in range(1, len(data_vector)):
        filtered_vector[i] = (alpha * data_vector[i]
                              + (1 - alpha) * filtered_vector[i - 1])

    return filtered_vector


# -----------------------------------------------------------------------
def ma_date() -> str:
    """Return a string containing the current date and time.

    Returns
    -------
    str
        The current date and time as 'YYYY-MM-DD_HHhMMmnSSs'.
    """
    n = datetime.now()
    return ("{0:04d}-{1:02d}-{2:02d}_{3:02d}h{4:02d}mn{5:02d}s"
            .format(n.year, n.month, n.day, n.hour, n.minute, n.second))


# -----------------------------------------------------------------------
def convert_time(time_in_csv: Union[list, np.ndarray]) -> np.ndarray:
    """Convert times from a csv format into python timestamps.

    Parameters
    ----------
    time_in_csv : list or numpy.ndarray
        Array of date strings or timestamps to convert.

    Returns
    -------
    numpy.ndarray
        The converted timestamps.
    """
    n = len(time_in_csv)
    time_converted = np.zeros(n)
    for i in range(n):
        dt = parse_date_auto(time_in_csv[i])
        time_converted[i] = dt.timestamp()
    return time_converted


# -----------------------------------------------------------------------
def do_power_spectrum(x: np.ndarray, fs: float, npts: int,
                      window: str = "none",
                      verbose: bool = False) -> Tuple[np.ndarray,
                                                      np.ndarray]:
    """Compute the power spectrum of the input vector.

    If the input vector is long enough, several computations are averaged.
    A window can be applied before the rfft.

    Parameters
    ----------
    x : numpy.ndarray
        Input vector.
    fs : float
        Sampling frequency.
    npts : int
        Number of values to be used in the rfft.
    window : str, optional
        Indicates if a window shall be applied ('blackman' or 'none').
        Default is 'none'.
    verbose : bool, optional
        If True, print power conservation checks. Default is False.

    Returns
    -------
    xf : numpy.ndarray
        Frequencies.
    spectrum : numpy.ndarray
        Computed power spectrum.
    """
    from numpy.fft import rfft
    from scipy.signal.windows import blackman

    if window == "blackman":
        w = blackman(npts, False)
    else:
        w = np.ones(npts)

    if len(x) < npts:
        raise ValueError("Not enough values in input vector to compute "
                         "spectra.")

    # removing DC
    x = x.astype(float) - x.mean()

    # spectral analysis with rfft
    xf = np.fft.rfftfreq(npts, 1 / fs)

    power_spectrum = np.zeros_like(xf)
    nslices = int(len(x) / npts)
    for slice in range(nslices):
        power_spectrum += (abs(rfft(x[slice * npts:(slice + 1) * npts] * w))
                           ** 2)

    power_spectrum /= nslices  # amplitude per slice
    power_spectrum /= npts ** 2  # amplitude per sample
    power_spectrum[1:-1] *= 2  # compensation for positive components

    if verbose:
        # power conservation check
        power_time_domain = np.mean(x ** 2)
        power_freq_domain = np.sum(power_spectrum)
        print(f"Puissance temporelle : {power_time_domain:.5f}")
        print(f"Puissance fréquentielle : {power_freq_domain:.5f}")

    return xf, power_spectrum


# -----------------------------------------------------------------------
def do_cross_power_spectrum(x: np.ndarray, y: np.ndarray, fs: float,
                            npts: int, window: str = "none",
                            verbose: bool = False) -> Tuple[np.ndarray,
                                                            np.ndarray]:
    """Compute the cross-power spectrum of two input vectors.

    If the input vectors are long enough, several computations are averaged.
    A window can be applied before the rfft.

    Parameters
    ----------
    x : numpy.ndarray
        First input vector.
    y : numpy.ndarray
        Second input vector.
    fs : float
        Sampling frequency.
    npts : int
        Number of values to be used in the rfft.
    window : str, optional
        Indicates if a window shall be applied ('blackman' or 'none').
        Default is 'none'.
    verbose : bool, optional
        If True, print power conservation checks. Default is False.

    Returns
    -------
    xf : numpy.ndarray
        Frequencies.
    spectrum : numpy.ndarray
        Computed cross-power spectrum.
    """
    from numpy.fft import rfft, rfftfreq
    from scipy.signal.windows import blackman
    from scipy.signal import correlate

    if window == "blackman":
        w = blackman(npts, False)
    else:
        w = np.ones(npts)

    # ensuring that x and y have the same size
    n = min(len(x), len(y))
    x = x[:n]
    y = y[:n]

    if n < npts:
        raise ValueError("Not enough values in input vector to compute "
                         "spectra.")

    # compute the correlation of the two data vectors
    print("    Computing cross-correlation......... ", end='')
    correl = correlate(x, y, mode='full')
    print("Done")

    # spectral analysis with rfft
    print("    Doing spectral analysis......... ", end='')
    xf = rfftfreq(npts, 1 / fs)
    power_spectrum = np.zeros_like(xf)
    nslices = int(len(correl) / npts)
    for slice in range(nslices):
        power_spectrum += abs(rfft(correl[slice * npts:(slice + 1) * npts]
                                   * w))

    power_spectrum /= nslices  # amplitude per slice
    power_spectrum /= npts ** 2  # amplitude per sample
    power_spectrum[1:-1] *= 2  # compensation for positive components
    print("Done")

    if verbose:
        # power conservation check
        power_time_domain = np.mean(x ** 2)
        power_freq_domain = np.sum(power_spectrum)
        print(f"Puissance temporelle : {power_time_domain:.5f}")
        print(f"Puissance fréquentielle : {power_freq_domain:.5f}")

    return xf, power_spectrum


# -----------------------------------------------------------------------
def derivative(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Compute the derivative of a function from x and y values.

    Parameters
    ----------
    x : numpy.ndarray
        x values.
    y : numpy.ndarray
        y values.

    Returns
    -------
    numpy.ndarray
        The derivative of the function.

    Raises
    ------
    ValueError
        If x and y do not have the same length.
    """
    if len(x) != len(y):
        raise ValueError("x and y shall have the same length")

    return (y[1:] - y[:-1]) / (x[1:] - x[:-1])


# -----------------------------------------------------------------------
def unzip_files_from_dir(dir: str) -> bool:
    """Unzip all zip files found in the data directory of dir.

    Parameters
    ----------
    dir : str
        Base directory containing the data subdirectory.

    Returns
    -------
    bool
        True if at least one zip file was found and unzipped, False otherwise.
    """
    import zipfile

    data_from_zip_file = False
    # looking for zip files if any
    data_path = os.path.join(dir, cst.dataDirName)
    zipfiles = [f for f in os.listdir(data_path)
                if os.path.isfile(os.path.join(data_path, f))
                and f.endswith('.zip')]

    # unzipping zip files
    if len(zipfiles) == 0:
        print('There is no zip file...')
    else:
        data_from_zip_file = True
        for z in zipfiles:
            print('Unzipping Error data from ', z)
            with zipfile.ZipFile(os.path.join(data_path, z), 'r') as zip_ref:
                zip_ref.extractall(data_path)

    return data_from_zip_file


# -----------------------------------------------------------------------
def read_two_vectors_from_file(
        nom_fichier: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read two vectors from a two-column file, ignoring the header line.

    Parameters
    ----------
    nom_fichier : str
        Path to the file containing the data. The first line is ignored
        as it is assumed to be a header.

    Returns
    -------
    x : numpy.ndarray
        The first vector extracted from the file.
    y : numpy.ndarray
        The second vector extracted from the file.

    Raises
    ------
    ValueError
        If the total number of values is not even (2 columns expected).
    """
    # read all lines of the file
    with open(nom_fichier, 'r') as fichier:
        lines = fichier.readlines()

    # ignore the first line (header)
    data_lines = lines[1:]

    # prepare a list to store the data
    data = []
    for line in data_lines:
        # clean the line: replace multiple spaces by a single one
        line = line.strip()
        line = ' '.join(line.split())
        data.extend(line.split())

    # convert to a numpy array
    data = np.array(data, dtype=float)

    # reshape into two columns
    if len(data) % 2 != 0:
        raise ValueError("Le nombre total de valeurs doit être pair "
                         "(2 colonnes).")

    data = data.reshape((-1, 2))

    # split the columns into two vectors
    x = data[:, 0]
    y = data[:, 1]

    return x, y


# -----------------------------------------------------------------------
def save_two_vectors_to_file(vecteur1: Union[list, np.ndarray],
                             vecteur2: Union[list, np.ndarray],
                             fichier: str, label1: str = 'x',
                             label2: str = 'y') -> None:
    """Save two vectors to a text file.

    Parameters
    ----------
    vecteur1 : list or numpy.ndarray
        First vector.
    vecteur2 : list or numpy.ndarray
        Second vector.
    fichier : str
        Path of the text file to create or overwrite.
    label1 : str, optional
        Header of the first column. Default is 'x'.
    label2 : str, optional
        Header of the second column. Default is 'y'.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the two vectors do not have the same length.
    """
    if len(vecteur1) != len(vecteur2):
        raise ValueError("Les deux vecteurs doivent avoir la même longueur.")

    with open(fichier, "w") as f:
        f.write(label1 + '   ' + label2 + "\n")
        for val1, val2 in zip(vecteur1, vecteur2):
            f.write(f"{val1}  {val2}\n")


# -----------------------------------------------------------------------
def pulseshapingtext(pls_set: int) -> str:
    """Return the text description of a pulse shaping setting.

    Parameters
    ----------
    pls_set : int
        The pulse shaping setting code.

    Returns
    -------
    str
        The text description of the setting.
    """
    switch = {
        0: 'no pulse shaping',
        1: 'pulse shaping at 20 MHz',
        2: 'pulse shaping at 25 MHz',
        3: 'pulse shaping at 30 MHz'
    }
    return switch.get(pls_set, "Invalid input")


# -----------------------------------------------------------------------
def parse_date_auto(value: Union[int, float, str, None]) -> datetime:
    """Convert a value into a datetime automatically.

    Accepts either a timestamp (float or string) or a date string in
    'DD/MM/YYYY HH:MM:SS', 'DD/MM/YYYY HH:MM' or 'DD/MM/YYYY' format.

    Parameters
    ----------
    value : int, float, str or None
        The value to convert.

    Returns
    -------
    datetime
        The converted date and time.

    Raises
    ------
    ValueError
        If the value format is not recognized.
    """
    if isinstance(value, (int, float)):
        # direct timestamp format
        return datetime.fromtimestamp(value)

    if value is None:
        raise ValueError("Format de date non reconnu : valeur None")

    value_str = str(value).strip()

    try:
        # try float to timestamp conversion
        ts = float(value_str)
        return datetime.fromtimestamp(ts)
    except ValueError:
        pass  # not a float, maybe a date string

    date_formats = (
        "%d/%m/%Y %H:%M:%S",
        "%d/%m/%Y %H:%M",
        "%d/%m/%Y",
    )

    for date_format in date_formats:
        try:
            return datetime.strptime(value_str, date_format)
        except ValueError:
            pass

    raise ValueError(
        f"Format de date non reconnu : {value!r}. "
        f"Formats acceptés : timestamp, DD/MM/YYYY HH:MM:SS, "
        f"DD/MM/YYYY HH:MM, DD/MM/YYYY"
    )


# -----------------------------------------------------------------------
def bit_value(x: int, bit_id: int) -> int:
    """Retrieve the value of a specific bit from an integer.

    Parameters
    ----------
    x : int
        The integer from which to extract the bit.
    bit_id : int
        The position of the bit to retrieve (0-indexed).

    Returns
    -------
    int
        The bit value.

    Raises
    ------
    TypeError
        If bit_id is not an integer.
    """
    if not isinstance(bit_id, int):
        raise TypeError("The bit number shall be an integer")
    else:
        return (x >> bit_id) % 2


# -----------------------------------------------------------------------
def peakdetect(sig: np.ndarray, half_space: int,
               ref_ratio: float = 5) -> np.ndarray:
    """Detect peaks in a signal.

    Parameters
    ----------
    sig : numpy.ndarray
        The input signal.
    half_space : int
        Half width of the detection window.
    ref_ratio : float, optional
        Detection threshold ratio. Default is 5.

    Returns
    -------
    numpy.ndarray
        The indices of the detected peaks.
    """
    # first search of a spike (could do a double detection on the rise and
    # fall edges)
    ratio = np.zeros_like(sig)
    ratio[half_space:-1 * half_space] = (sig[half_space:-1 * half_space]
                                         / (np.abs(sig[2 * half_space:]
                                                   + sig[:-2 * half_space])
                                            / 2))
    x_ini = np.where(ratio > ref_ratio)[0]

    # second search around the first detections based on a maximum detection
    for ix in range(len(x_ini)):
        new_ix = x_ini[ix] - half_space + np.where(
            sig[x_ini[ix] - half_space:x_ini[ix] + half_space]
            == sig[x_ini[ix] - half_space:x_ini[ix] + half_space].max())[0]
        x_ini[ix] = new_ix

    # removing multiple identical results
    return np.unique(x_ini)


# -----------------------------------------------------------------------
def get_test_type_from_session_name(session_name: str) -> str:
    """Extract the test type from a session directory name.

    The test type is read from a fixed position in the session name,
    according to the naming convention of the DEMUX performance tests.

    Parameters
    ----------
    session_name : str
        Name of the session directory.

    Returns
    -------
    str
        The test type, as found in the session name.
    """
    return session_name[16:31]
