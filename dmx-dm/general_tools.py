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
#  laurent.ravera@irap.omp.eu
#  general_tools.py
#
# ---------------------------------------------------------------------------------

import csv
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import constants as cst

# definition of equivalent noise bandwidth for some well known windows
enb = {'boxcar': 1.0, 'hann': 1.5, 'hamming': 1.363, 'blackman': 1.727, 'blackmanharris': 2.004,
       'barlett': 1.333, 'flattop': 3.77, 'triang': 1.333}


# -----------------------------------------------------------------------
def checkdir(dirname):
    r"""
        This function checks if a directory exists. If not it creates it.

        Parameters:
        -----------
        dirname: String
        Name of the directory to be verified / created.

        Returns
        -------
        Nothing

        """
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    return ()


# -----------------------------------------------------------------------
def purgedir(dirname):
    r"""
        This function erase a directory if it exists, and it creates a new empty directory

        Parameters:
        -----------
        dirname: String
        Name of the directory to be erased / created.

        Returns
        -------
        Nothing

        """
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)
    os.mkdir(dirname)
    return ()


# -----------------------------------------------------------------------
def createdir(dirname):
    r"""
        This function creates a directory if it doesn't exist

        Parameters:
        -----------
        dirname: String
        Name of the directory to be erased / created.

        Returns
        -------
        Nothing

        """
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    return ()


# -----------------------------------------------------------------------
def get_csv(filename):
    r"""
        This function reads a dictionnary from a csv file.

        Parameters:
        -----------
        filename: string
        The name of the csv file

        Returns
        -------
        dictionnary: dictionnary

        """

    dictionnary = {}

    if not os.path.exists(os.path.join(filename)):
        print("File " + filename + " not found.")
    else:
        with open(filename, newline='') as csvfile:
            dict_reader = csv.reader(csvfile, delimiter='=', quotechar='|')
            for row in dict_reader:
                try:  # for numbers
                    dictionnary[row[0]] = float(row[1].replace(',', '.'))
                except Exception:  # for strings
                    dictionnary[row[0]] = row[1]
    return dictionnary


# -----------------------------------------------------------------------
def compute_spectrum(signal, period, file_name, window='blackman', plot=False):
    from scipy.signal.windows import blackman, blackmanharris, hann

    plotFileName = os.path.join(cst.plotDirName, file_name[:-5] + '_sp.png')
    nPts = len(signal)

    w = np.ones(nPts)
    # create blackman window
    if window == 'blackman':
        w = blackman(nPts)
    if window == 'blackmanharris':
        w = blackmanharris(nPts)
    if window == 'hann':
        w = hann(nPts)

    # Calculer la FFT
    fft_result = np.fft.rfft(signal*w)

    # Normaliser le résultat en ADU
    fft_norm = np.abs(fft_result) / nPts * np.sqrt(2)

    # Ajuster la normalisation pour les composantes à zéro et à la fréquence d'échantillonnage / 2
    fft_norm[0] /= np.sqrt(2)
    if len(signal) % 2 == 0:
        fft_norm[-1] /= np.sqrt(2)

    # Créer les fréquences pour l'axe des x
    xf = np.linspace(0, 0.5/period, len(fft_norm))

    # Afficher le spectre
    if plot:
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111)
        ax.loglog(xf, fft_norm)
        ax.set_title('Spectrum')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('Amplitude (ADU)')
        ax.grid(True)
        fig.tight_layout()
        plt.savefig(plotFileName, dpi=300, format='png', bbox_inches='tight')

    return xf, fft_norm


# -----------------------------------------------------------------------
def is_even(x):
    r"""
        This function checks if a number is even.

        Parameters:
        -----------
        x: number
        The value to be verified

        Returns
        -------
        result: boolean
        True if the input is even, False if the input is odd.

        """
    return x % 2 == 0


# -----------------------------------------------------------------------
def first_order_low_pass_filter(data_vector, tau):
    """
    Applies a first-order low-pass filter to a data vector.

    Args:
        data_vector (numpy.ndarray): The input data vector.
        tau (float): Time constant of the filter (in seconds).

    Returns:
        numpy.ndarray: The filtered vector.
    """
    dt = 1.0 / cst.fSamp  # Time interval (can be adjusted as needed)
    alpha = dt / (tau + dt)  # Filter coefficient

    filtered_vector = np.zeros_like(data_vector)
    filtered_vector[0] = data_vector[0]  # Initial condition

    for i in range(1, len(data_vector)):
        filtered_vector[i] = alpha * data_vector[i] + (1 - alpha) * filtered_vector[i - 1]

    return filtered_vector


# -----------------------------------------------------------------------
def ma_date():
    """This function returns a string containing the date and time"""
    from datetime import datetime
    n = datetime.now()
    return ("{0:04d}-{1:02d}-{2:02d}_{3:02d}h{4:02d}mn{5:02d}s"
            .format(n.year, n.month, n.day, n.hour, n.minute, n.second))


# -----------------------------------------------------------------------
def readHkFromCsv(filename):
    """This function reads a set of HK and a set of time from a csv file"""
    df = np.array(pd.read_csv(filename, sep=";"))
    time = df[:, 0]
    hk = df[:, 1]
    return time, hk


# -----------------------------------------------------------------------
def do_power_spectrum(x, fs, npts, window="none", verbose=False):
    r"""
        This function computes the spectrum of the input vector.
        If the input vector is long enough, several computations are averaged.
        A Blackman window is applied before the rfft.

        Parameters:
        -----------
        x: numpy array
        input vector

        fs: float
        sampling frequency

        npts: number
        Number of values to be used in the rfft.

        window: string
        indicates if a window shall be applied ("blackman" or "none")
        (default is "none")

        Returns
        -------
        xf: numpy array
        frequencies

        spectrum: numpy array
        computed spectrum.
        """
    from numpy.fft import rfft
    from scipy.signal.windows import blackman

    if window=="blackman":
        w=blackman(npts, False)
    else:
        w=np.ones(npts)

    if len(x)<npts:
        raise ValueError("Not enough values in input vector to compute spectra.")

    # Analyse spectrale avec rfft
    xf = np.fft.rfftfreq(npts, 1 / fs)  # Fréquences associées

    power_spectrum = np.zeros_like(xf)  # Initialisation du spectre
    nslices=int(len(x)/npts)
    for slice in range(nslices):
        power_spectrum += abs(rfft(x[slice*npts:(slice+1)*npts]*w))**2

    power_spectrum /= nslices  # Amplitude par slice
    power_spectrum /= npts**2  # Amplitude par échantillon
    power_spectrum[1:-1] *= 2  # Compensations pour les composantes positives

    if verbose:
        # Vérification de la conservation de la puissance
        power_time_domain = np.mean(x ** 2)  # Puissance du signal temporel
        power_freq_domain = np.sum(power_spectrum)  # Puissance spectrale
        print(f"Puissance temporelle : {power_time_domain:.5f}")
        print(f"Puissance fréquentielle : {power_freq_domain:.5f}")

    return xf, power_spectrum


def derivative(x, y):
    """Computes the derivative of a function from x and y values.

    Parameters:
        x (array): xvalues.
        y (array): yvalues.

    Returns:
        derivative (array): the derivative of the function from x to y
    """
    if len(x) != len(y):
        raise ValueError("x and y shall have the same length")

    return (y[1:] - y[:-1]) / (x[1:] - x[:-1])

