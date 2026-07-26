# ---------------------------------------------------------------------------------
# !/usr/bin/env python
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
#  readData.py
#
# ---------------------------------------------------------------------------------

# Imports
import os

import h5py
import numpy as np
import pandas as pd
from astropy.io import fits

import constants as cst


def read_event_records(fits_file):
    """
    Reads DEMUX science data (error or science) from a fits file.

    Parameters:
        fits_file (string): Name of the fits file (includes the path).

    Returns:
        firstDf (int): number of the Data Frame corresponding to the beginning of the fits file.

        data (array): event records. Each element corresponds to a record. It contains:
            - the data frame number at the beginning of the record
            - the data
    """
    with fits.open(fits_file) as hdul:
        # Reading the header
        firstDf = hdul[1].header['FIRST_DF']
        data = hdul[1].data
        return firstDf, data


def get_science_from_fits(fullFileName):
    """
    Reads DEMUX science data (error or science) from a fits file.

    Parameters:
        fullFileName (string): Name of the fits file (includes the path).

    Returns:
        error (array): The error values (one value per pixel and per step)
        ctrl (array): The control words (one value per step)
    """
    with fits.open(fullFileName) as hdul:

        # Checking the extension name 'pixels_data'
        if 'pixels_data' in hdul:
            data_ext = hdul['pixels_data']
        else:
            raise ValueError("'pixels_data' extension not found in the FITS file.")

        # Extraire les données en tant que tableau NumPy
        data = data_ext.data
        dataArray = np.array(data.tolist()).T

    ctrl = dataArray[0, :]
    data = dataArray[1:, :].astype(float) / 4
    print(data.shape)

    return data, ctrl


def get_science_from_hdf5(fullFileName):
    """
    Reads DEMUX science data (error or science) from an HDF5 file.

    Parameters:
        fullFileName (str): Path to the HDF5 file.

    Returns:
        data (np.ndarray): The science data (one value per pixel and per step, divided by 4).
        ctrl (np.ndarray): The control words (one value per step).
    """
    with h5py.File(fullFileName, 'r') as f:
        # Lire les données
        ctrl = f["ctrl"][()]
        data = np.array(f["pixels"][()]).T

        data = data.astype(float) / 4  # Conversion to s(16,2) format

    return data, ctrl


def read_science_from_file(fullFileName, flatten=False, remove_dc=True, verbose=True):
    """
    Reads DEMUX science data for one column from a fits file.

    Parameters:

        fullFileName (string): filename including the path.
        flatten (boolean): If True the data are arranged at Frow, else the data are arranged at FFrame
        remove_dc (boolean): If True, the dc of the data is removed (default is True)
        verbose: (boolean): If True, some text is displayed

    Returns:
        col_data (numpy array): The data
     """

    if verbose:
        print("    Reading TM data from ", fullFileName, ".... ")

    col_data, _ = get_science_from_hdf5(fullFileName)
    # If requested, flattening the array to have data at Frow
    if flatten:
        col_data = col_data.flatten('F')

    # removing DC
    if remove_dc:
        if verbose:
            print("     Print removing DC")
        col_data -= col_data.mean()

    if verbose:
        print("Done!")

    return col_data


def read_col_science_from_dir(data_path, col_id, flatten=False, remove_dc=True, verbose=True):
    """
    Reads DEMUX science data for one column from a fits file.

    Parameters:
        data_path (string): Path to the data files.
        col_id (int): Column ID (0 to 3)
        flatten (boolean): If True the data are arranged at Frow, else the data are arranged at FFrame
        remove_dc (boolean): If True, the dc of the data is removed (default is True)
        verbose: (boolean): If True, some text is displayed

    Returns:
        col_data (numpy array): The data
    """
    files = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             # and f[:6] == 'error_' and f[-6:] == '{0:}.fits'.format(col_id)]
             and f[:4] != 'dump' \
             and f[-4:] == '{0:}.h5'.format(col_id)]

    file_exists = (len(files) != 0)

    if file_exists:
        if len(files) > 1:
            print("   Warning, {0:3d} files in the directory, processing only one file...".format(len(files)))
        file_name = files[0]
        file_name_with_path = os.path.join(data_path, file_name)

        col_data = read_science_from_file(file_name_with_path, flatten, remove_dc, verbose)

    else:
        col_data = 0

    return col_data, file_exists


def read_dump_from_fits(fits_file):
    """Reads DEMUX dump data from a fits file.

    Parameters
        fits_file (string): Name of the fits file (includes the path).

    Returns:
        dump (array): the dump data (4 x 1360 values)
        adc_error (array): the ADC error data (1360 bytes)
    """
    with fits.open(fits_file) as hdul:

        # Checking the extension name 'Dumps'
        if 'Dumps' in hdul:
            data_ext = hdul['Dumps']
        else:
            raise ValueError("'Dumps' extension not found in the FITS file.")

        # Extraire les données en tant que tableau NumPy
        data = data_ext.data
        dump = np.array(data.tolist())[0]

        # Return columns and errors
        return dump[0:-1,:], dump[-1,:]


def read_dump_from_hdf5(hdf5_file):
    """
    Lit les données DEMUX dump et les erreurs ADC depuis un fichier HDF5.

    Args:
        hdf5_file (str): Chemin vers le fichier HDF5.

    Returns:
        tuple: (dump, adc_error)
            - dump: tableau numpy de forme (4, 1360) contenant Col0, Col1, Col2, Col3
            - adc_error: tableau numpy de forme (1360,) contenant les erreurs ADC
    """
    size = 2 * cst.nSamplesPerRow * cst.nPixPerCol
    with h5py.File(hdf5_file, 'r') as f:
        # Read the columns data Col0, Col1, Col2, Col3
        col0 = f['Col0'][0, :]
        col1 = f['Col1'][0, :]
        col2 = f['Col2'][0, :]
        col3 = f['Col3'][0, :]

        # Read the conversion errors
        adc_error = f['Errors'][0, :]

        # Check data format consistency
        if col0.shape != (size,) or col1.shape != (size,) or col2.shape != (size,) or col3.shape != (
                size,) or adc_error.shape != (size,):
            print(col0.shape)
            raise ValueError("Data have not the expected size ({0:},).".format(size))

        # Retourner les données des colonnes et les erreurs
        dump = np.array([col0, col1, col2, col3])
        return dump, adc_error


def read_scan_fits(fits_file):
    """Reads DEMUX scan data from a fits file.

    Parameters
        fits_file (string): Name of the fits file (includes the path).

    Returns:
        xName (string): the name of the signal on the X axis
        ctrl (array): array with the control words
        xValues (array): array with the x values
        error (array): array with the error values (one value per pixel and per step)
    """
    with fits.open(fits_file) as hdul:

        # Checking the extension name 'pixels_data'
        if 'pixels_data' in hdul:
            data_ext = hdul['pixels_data']
        else:
            raise ValueError("'pixels_data' extension not found in the FITS file.")

        xName = data_ext.header['TTYPE2']

        # Extraire les données en tant que tableau NumPy
        data = data_ext.data
        scan = np.array(data.tolist())
        errors = scan[:,2:].T

        # Return xName, CTRL, xValues and error per pixels
        return xName, scan[:,0], scan[:,1], errors


def read_scan(hdf5_file):
    """Reads DEMUX scan data from an hdf5 file.

    Parameters
        hdf5_file (string): Name of the hdf5 file (includes the path).

    Returns:
        xName (string): the name of the signal on the X axis
        ctrl (array): array with the control words
        xValues (array): array with the x values
        error (array): array with the error values (one value per pixel and per step)
    """
    with h5py.File(hdf5_file, 'r') as f:
        # Getting the name of the x data (feedback or offset)
        xName = f.attrs["X_LABEL"].decode("utf-8")

        ctrl = np.array(f['ctrl'])
        pixels_data = np.array(f['pixels']).T
        x = np.array(f['x'])

        # Return xName, CTRL, xValues and error per pixels
        return xName, ctrl, x, pixels_data


def read_scan_type(hdf5_file):
    """Reads DEMUX scan data from an hdf5 file.

    Parameters
        hdf5_file (string): Name of the hdf5 file (includes the path).

    Returns:
        xName (string): the name of the signal on the X axis
    """
    with h5py.File(hdf5_file, 'r') as f:
        # Getting the name of the x data (feedback or offset)
        xName = f.attrs["X_LABEL"].decode("utf-8")

        return xName


def export_dump_2_txt(fits_file, col):
    """Saves the data of a column from a dump to a text file.

    Parameters
        fits_file (string): the name of the dump fits_file
        col (integer): the col id
    """
    data, _ = read_dump_from_fits(fits_file)
    np.savetxt(fits_file[0:-5] + '_col_{0:}.txt'.format(col), data[col,:])


def export_science_data_one_col_2_txt(data_path, col, remove_dc=False, verbose=False):
    """Saves the data of a column from a science fits file to a text file.

    Parameters
        data_path (string): path to the fits_file
        col (integer): the col id
    """
    print(data_path)
    data = read_col_science_from_dir(data_path, col, flatten=True, remove_dc=False, verbose=True)
    txt_file_name = os.path.join(data_path, "error_txt_file" + '_col_{0:}.txt'.format(col))
    np.savetxt(txt_file_name, data)


# def read_hk_name_from_csv(hk_file, hk_name, encoding="latin1"):
#    # Loading the csv file
#    df = pd.read_csv(hk_file, sep=';', encoding=encoding)
#
#    if hk_name[:4] == 'Date':
#        df[hk_name] = pd.to_datetime(df[hk_name], dayfirst=True, errors="coerce")
#        df[hk_name] = pd.to_datetime(df[hk_name], format="%Y:%M:%D %H:%M:%S")
#
#    return df[hk_name]


def read_hk_name_from_csv(hk_file, hk_suffix, encoding="latin1"):
    """
    Lit une colonne d'un fichier CSV en utilisant les 'suffix_length' derniers caractères de son nom.

    Args:
        hk_file (str): Chemin vers le fichier CSV.
        hk_suffix (str): Les 'suffix_length' derniers caractères du nom de la colonne.
        suffix_length (int, optionnel): Longueur du suffixe à vérifier. Par défaut 10.
        encoding (str, optionnel): Encodage du fichier. Par défaut "latin1".

    Returns:
        pd.Series: La colonne correspondante.

    Raises:
        ValueError: Si aucune ou plusieurs colonnes correspondent au suffixe.
    """
    df = pd.read_csv(hk_file, sep=';', encoding=encoding)

    # Trouver les colonnes dont le nom se termine par hk_suffix
    suffix_length = len(hk_suffix)
    matching_columns = [
        col for col in df.columns
        if col[-suffix_length:] == hk_suffix
    ]

    if not matching_columns:
        raise ValueError(f"No HK match the name '{hk_suffix}'.")
    if len(matching_columns) > 1:
        raise ValueError(f"More than one HK match the name '{hk_suffix}': {matching_columns}")

    selected_column = matching_columns[0]

    # Conversion en datetime si le nom de la colonne commence par 'Date'
    if selected_column[:4] == 'Date':
        df[selected_column] = pd.to_datetime(df[selected_column], dayfirst=True, errors="coerce")
        df[selected_column] = pd.to_datetime(df[selected_column], format="%Y:%M:%D %H:%M:%S")

    return df[selected_column]


def read_fwVersion_dmxModel(path):
    # looking for the hk files
    files = [f for f in os.listdir(path) \
             if os.path.isfile(os.path.join(path, f)) \
             and f[:7] == "Hks_DMX" \
             and f[-4:] == ".csv"]

    if len(files) == 0:
        print("No HK found")
    else:
        fwVersion = read_hk_name_from_csv(os.path.join(path, files[0]), "Firmware Version")[0]
        # fwVersion = read_hk_name_from_csv(os.path.join(path, files[0]), "firmwareVersion")[0]
        ref = read_hk_name_from_csv(os.path.join(path, files[0]), "Hardware Version")[0]
        #ref = read_hk_name_from_csv(os.path.join(path, files[0]), "hardwareVersion")[0]

        dmxModel = (ref >> 8) & 3
        boardId = ref & (2 ** 5) - 1

        return cst.dmx_models[dmxModel], boardId, fwVersion


def read_dmxConfig_fromXml(path):
    from xml.dom import minidom

    dict_conf = {'fw_version': '',
                 'hw_version': '',
                 'boxcar_length': '',
                 'c0_pulse_shaping_set': '',
                 'c0_offset_coarse': '',
                 'c0_offset_lsb': '',
                 'c0_sampling_delay': '',
                 'c0_feedback_delay': '',
                 'c0_offset_dac_delay': '',
                 'c0_offset_mux_delay': '',
                 'relock_delay': '',
                 'relock_threshold': ''}

    # looking for the hk files
    files = [f for f in os.listdir(path) \
             if os.path.isfile(os.path.join(path, f)) \
             and f[-4:] == ".xml"]

    if len(files) != 1:
        raise ValueError(f"Wrong number of xml files: expected 1, found {0:}".format(len(files)))
    else:
        # parsing the xml file
        file = minidom.parse(files[0])
        dmx = file.getElementsByTagName('dmx')
        dict_conf[key] = models[1].attributes['name'].value

    return dmx_config
