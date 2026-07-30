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
from datetime import datetime

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
    n = datetime.now()
    return ("{0:04d}-{1:02d}-{2:02d}_{3:02d}h{4:02d}mn{5:02d}s"
            .format(n.year, n.month, n.day, n.hour, n.minute, n.second))


# -----------------------------------------------------------------------
def readHkFromCsv(filename):
    """This function reads a set of HK and a set of time from a csv file"""
    df = np.array(pd.read_csv(filename, sep=";"))
    time = df[:, 0]
    # conversion du format de date en timestamp
    time_converted = np.zeros(len(time))
    for i in range(len(time)):
        dt = parse_date_auto(time[i])
        time_converted[i] = dt.timestamp()
    hk = df[:, 1]
    return time_converted, hk


# -----------------------------------------------------------------------
def convert_time(time_in_csv):
    """
    This function converts times from a csv format into a python timestamp format
    """
    # conversion du format de date en timestamp
    l = len(time_in_csv)
    time_converted = np.zeros(l)
    for i in range(l):
        dt = parse_date_auto(time_in_csv[i])
        time_converted[i] = dt.timestamp()
    return time_converted


# -----------------------------------------------------------------------
def do_power_spectrum(x, fs, npts, window="none", verbose=False):
    r"""
        This function computes the spectrum of the input vector.
        If the input vector is long enough, several computations are averaged.
        A Blackman window can be applied before the rfft.

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

    # Removing DC
    x = x.astype(float) - x.mean()

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


# -----------------------------------------------------------------------
def do_cross_power_spectrum(x, y, fs, npts, window="none", verbose = False):
    r"""
        This function computes the cross-power spectrum of two input vectors.
        If the input vector is long enough, several computations are averaged.
        A Blackman window can be applied before the rfft.

        Parameters:
        -----------
        x: numpy array
        input vector

        y: numpy array
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
    from numpy.fft import rfft, rfftfreq
    from scipy.signal.windows import blackman
    from scipy.signal import correlate

    if window=="blackman":
        w=blackman(npts, False)
    else:
        w=np.ones(npts)

    # Ensuring that x and y have the same size
    l = min(len(x), len(y))
    x = x[:l]
    y = y[:l]

    if l<npts:
        raise ValueError("Not enough values in input vector to compute spectra.")

    # Calcul de la correlation des deux vecteurs de données
    print("    Computing cross-correlation......... ", end = '')
    correl = correlate(x, y, mode='full')
    #norm_factor = np.linalg.norm(x) * np.linalg.norm(y)
    #correl /= norm_factor
    print("Done")

    # Analyse spectrale avec rfft
    print("    Doing spectral analysis......... ", end = '')
    xf = rfftfreq(npts, 1 / fs)  # Fréquences associées
    power_spectrum = np.zeros_like(xf)  # Initialisation du spectre
    nslices=int(len(correl)/npts)
    for slice in range(nslices):
        power_spectrum += abs(rfft(correl[slice*npts:(slice+1)*npts]*w))

    power_spectrum /= nslices  # Amplitude par slice
    power_spectrum /= npts**2  # Amplitude par échantillon
    power_spectrum[1:-1] *= 2  # Compensations pour les composantes positives
    print("Done")

    if verbose:
        # Vérification de la conservation de la puissance
        power_time_domain = np.mean(x ** 2)  # Puissance du signal temporel
        power_freq_domain = np.sum(power_spectrum)  # Puissance spectrale
        print(f"Puissance temporelle : {power_time_domain:.5f}")
        print(f"Puissance fréquentielle : {power_freq_domain:.5f}")

    return xf, power_spectrum

# -----------------------------------------------------------------------
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


# -----------------------------------------------------------------------
def unzip_files_from_dir(dir):

    import zipfile

    data_from_zip_file = False
    # Looking for zip files if any
    dataPath = os.path.join(dir, cst.dataDirName)
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

    return data_from_zip_file


# -----------------------------------------------------------------------
def read_two_vectors_from_file_obsolete(nom_fichier):
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
            - x: A NumPy array representing the first vector
              extracted from the specified file.
            - y: A NumPy array representing the second vector
              extracted from the specified file.
    """
    # Charger les données en ignorant les lignes d'en-tête
    data = np.loadtxt(nom_fichier, skiprows=1)

    # Séparer les colonnes en deux vecteurs
    x = data[:, 0]
    y = data[:, 1]

    return x, y


def read_two_vectors_from_file(nom_fichier):
    """
    Lit deux vecteurs depuis un fichier contenant 2 colonnes et 1 ou plusieurs lignes de données.
    Ignore la première ligne (en-tête).
    Retourne deux tableaux numpy x et y.
    """
    # Lire toutes les lignes du fichier
    with open(nom_fichier, 'r') as fichier:
        lines = fichier.readlines()

    # Ignorer la première ligne (en-tête)
    data_lines = lines[1:]

    # Préparer une liste pour stocker les données
    data = []
    for line in data_lines:
        # Nettoyer la ligne : remplacer les espaces multiples par un seul
        line = line.strip()
        line = ' '.join(line.split())
        # Ajouter les valeurs à la liste data
        data.extend(line.split())

    # Convertir en tableau numpy
    data = np.array(data, dtype=float)

    # Reshaper en deux colonnes (on suppose que le nombre total de valeurs est pair)
    if len(data) % 2 != 0:
        raise ValueError("Le nombre total de valeurs doit être pair (2 colonnes).")

    data = data.reshape((-1, 2))

    # Séparer les colonnes en deux vecteurs
    x = data[:, 0]
    y = data[:, 1]

    return x, y

def save_two_vectors_to_file(vecteur1, vecteur2, fichier, label1='x', label2='y'):
    """
    Sauvegarde deux vecteurs dans un fichier texte.

    :param fichier: Chemin du fichier texte à créer ou à écraser.
    :param vecteur1: Premier vecteur (list ou numpy array).
    :param vecteur2: Deuxième vecteur (list ou numpy array).
    :param separateur: Caractère utilisé pour séparer les données (par défaut : espace).
    """
    if len(vecteur1) != len(vecteur2):
        raise ValueError("Les deux vecteurs doivent avoir la même longueur.")

    with open(fichier, "w") as f:
        f.write(label1 + '   ' + label2 + "\n")
        for val1, val2 in zip(vecteur1, vecteur2):
            f.write(f"{val1}{"  "}{val2}\n")


def pulseshapingtext(pls_set):
    switch={
      0:'no pulse shaping',
      1:'pulse shaping at 20 MHz',
      2:'pulse shaping at 25 MHz',
      3:'pulse shaping at 30 MHz'
      }
    return switch.get(pls_set, "Invalid input")


def parse_date_auto_old(value):
    """
    Convertit automatiquement une valeur de date en datetime.
    Accepte soit :
    - une chaîne de type '28/05/2025 09:42:21'
    - un timestamp (float ou chaîne de float comme '1748347100.938')
    """

    if isinstance(value, (int, float)):
        # Format timestamp direct
        return datetime.fromtimestamp(value)

    try:
        # Essai conversion float → timestamp
        ts = float(value)
        return datetime.fromtimestamp(ts)
    except ValueError:
        pass  # Ce n'était pas un float, peut-être une chaîne de date

    try:
        # Essai conversion chaîne de date classique
        return datetime.strptime(value.strip(), "%d/%m/%Y %H:%M:%S")
    except ValueError as e:
        raise ValueError(f"Format de date non reconnu : {value}")


def parse_date_auto(value):
    """
    Convertit automatiquement une valeur de date en datetime.
    Accepte soit :
    - une chaîne de type '28/05/2025 09:42:21'
    - une chaîne de type '28/05/2025 09:42'
    - une chaîne de type '28/05/2025'
    - un timestamp (float ou chaîne de float comme '1748347100.938')
    """

    if isinstance(value, (int, float)):
        # Format timestamp direct
        return datetime.fromtimestamp(value)

    if value is None:
        raise ValueError("Format de date non reconnu : valeur None")

    value_str = str(value).strip()

    try:
        # Essai conversion float → timestamp
        ts = float(value_str)
        return datetime.fromtimestamp(ts)
    except ValueError:
        pass  # Ce n'était pas un float, peut-être une chaîne de date

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
        f"Formats acceptés : timestamp, DD/MM/YYYY HH:MM:SS, DD/MM/YYYY HH:MM, DD/MM/YYYY"
    )

def bit_value(x, bit_id):
    """
    This function retrieves the value of a specific bit from an integer.
    The function raises a TypeError if 'bit_id' is not an integer.
    This function can be applied to integers or numpy arrays.

    Args:
        x: the integer from which to extract the bit
        bit_id: the position of the bit to retrieve (0-indexed)

    Returns: The bit value

    """
    if not isinstance(bit_id, int):
        raise TypeError("The bit number shall be an integer")
    else:
        return (x >> bit_id) % 2


def peakdetect(sig, half_space, ref_ratio=5):
    # First search of a spike (could do a double detection on the rise and fall edges)
    ratio = np.zeros_like(sig)
    ratio[half_space:-1 * half_space] = sig[half_space:-1 * half_space] / (
                np.abs(sig[2 * half_space:] + sig[:-2 * half_space]) / 2)
    x_ini = np.where(ratio > ref_ratio)[0]

    # Second search around the first detections based on a maximum detection
    for ix in range(len(x_ini)):
        new_ix = x_ini[ix] - half_space + np.where(sig[x_ini[ix] - half_space:x_ini[ix] + half_space] == sig[
            x_ini[ix] - half_space:x_ini[ix] + half_space].max())[0]
        x_ini[ix] = new_ix

    # Removing multiple identical results
    return np.unique(x_ini)
