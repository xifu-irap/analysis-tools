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

"""Functions to read DEMUX data (science, dump, scan, HK) from files."""

import os
from typing import Any, Tuple
from xml.dom import minidom

import h5py
import numpy as np
import pandas as pd

import constants as cst


def get_science_from_hdf5(
    full_file_name: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """Read DEMUX science data (error or science) from an HDF5 file.

    Parameters
    ----------
    full_file_name : str
        Path to the HDF5 file.

    Returns
    -------
    data : np.ndarray
        The science data (one value per pixel and per step, divided by 4).
    ctrl : np.ndarray
        The control words (one value per step).
    """
    with h5py.File(full_file_name, "r") as f:
        ctrl = f["ctrl"][()]
        data = np.array(f["pixels"][()]).T

        # Conversion to s(16,2) format
        data = data.astype(float) / 4

    return data, ctrl


def read_science_from_file(
    full_file_name: str,
    flatten: bool = False,
    remove_dc: bool = True,
    verbose: bool = True,
) -> np.ndarray:
    """Read DEMUX science data for one column from an HDF5 file.

    Parameters
    ----------
    full_file_name : str
        Filename including the path.
    flatten : bool, optional
        If True the data are arranged at Frow, else at FFrame.
    remove_dc : bool, optional
        If True the dc of the data is removed (default is True).
    verbose : bool, optional
        If True some text is displayed.

    Returns
    -------
    col_data : np.ndarray
        The data.
    """
    if verbose:
        print("    Reading TM data from ", full_file_name, ".... ")

    col_data, _ = get_science_from_hdf5(full_file_name)

    # If requested, flattening the array to have data at Frow
    if flatten:
        col_data = col_data.flatten("F")

    # Removing DC
    if remove_dc:
        if verbose:
            print("     Print removing DC")
        col_data -= col_data.mean()

    if verbose:
        print("Done!")

    return col_data


def read_col_science_from_dir(
    data_path: str,
    col_id: int,
    flatten: bool = False,
    remove_dc: bool = True,
    verbose: bool = True,
) -> Tuple[np.ndarray, bool]:
    """Read DEMUX science data for one column from a data directory.

    Parameters
    ----------
    data_path : str
        Path to the data files.
    col_id : int
        Column ID (0 to 3).
    flatten : bool, optional
        If True the data are arranged at Frow, else at FFrame.
    remove_dc : bool, optional
        If True the dc of the data is removed (default is True).
    verbose : bool, optional
        If True some text is displayed.

    Returns
    -------
    col_data : np.ndarray
        The data (0 if no file matches).
    file_exists : bool
        True if a matching file was found.
    """
    files = [
        f
        for f in os.listdir(data_path)
        if os.path.isfile(os.path.join(data_path, f))
        and f[:4] != "dump"
        and f[-4:] == "{0:}.h5".format(col_id)
    ]

    file_exists = len(files) != 0

    if file_exists:
        if len(files) > 1:
            print(
                "   Warning, {0:3d} files in the directory, "
                "processing only one file...".format(len(files))
            )
        file_name = files[0]
        file_name_with_path = os.path.join(data_path, file_name)

        col_data = read_science_from_file(
            file_name_with_path, flatten, remove_dc, verbose
        )
    else:
        col_data = 0

    return col_data, file_exists


def read_dump_from_hdf5(hdf5_file: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read DEMUX dump data and ADC errors from an HDF5 file.

    Parameters
    ----------
    hdf5_file : str
        Path to the HDF5 file.

    Returns
    -------
    dump : np.ndarray
        Array of shape (4, 1360) containing Col0, Col1, Col2, Col3.
    adc_error : np.ndarray
        Array of shape (1360,) containing the ADC errors.
    """
    size = 2 * cst.nSamplesPerRow * cst.nPixPerCol
    with h5py.File(hdf5_file, "r") as f:
        # Read the columns data Col0, Col1, Col2, Col3
        col0 = f["Col0"][0, :]
        col1 = f["Col1"][0, :]
        col2 = f["Col2"][0, :]
        col3 = f["Col3"][0, :]

        # Read the conversion errors
        adc_error = f["Errors"][0, :]

        # Check data format consistency
        expected = (size,)
        if (
            col0.shape != expected
            or col1.shape != expected
            or col2.shape != expected
            or col3.shape != expected
            or adc_error.shape != expected
        ):
            print(col0.shape)
            raise ValueError(f"Data have not the expected size ({size},).")

        # Return the columns data and the errors
        dump = np.array([col0, col1, col2, col3])
        return dump, adc_error


def _decode_text(value: Any) -> str:
    """Decode an HDF5 attribute to a string.

    Parameters
    ----------
    value : Any
        The attribute value, either ``str`` or ``bytes``.

    Returns
    -------
    str
        The decoded text.
    """
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def read_scan(
    hdf5_file: str,
) -> Tuple[str, np.ndarray, np.ndarray, np.ndarray]:
    """Read DEMUX scan data from an HDF5 file.

    Parameters
    ----------
    hdf5_file : str
        Name of the HDF5 file (includes the path).

    Returns
    -------
    x_name : str
        The name of the signal on the X axis.
    ctrl : np.ndarray
        Array with the control words.
    x_values : np.ndarray
        Array with the x values.
    error : np.ndarray
        Array with the error values (one value per pixel and per step).
    """
    with h5py.File(hdf5_file, "r") as f:
        # Getting the name of the x data (feedback or offset)
        x_name = _decode_text(f.attrs["X_LABEL"])

        ctrl = np.array(f["ctrl"])
        pixels_data = np.array(f["pixels"]).T
        x_values = np.array(f["x"])

        # Return xName, CTRL, xValues and error per pixels
        return x_name, ctrl, x_values, pixels_data


def read_scan_type(hdf5_file: str) -> str:
    """Read the type of a DEMUX scan from an HDF5 file.

    Parameters
    ----------
    hdf5_file : str
        Name of the HDF5 file (includes the path).

    Returns
    -------
    x_name : str
        The name of the signal on the X axis.
    """
    with h5py.File(hdf5_file, "r") as f:
        # Getting the name of the x data (feedback or offset)
        x_name = _decode_text(f.attrs["X_LABEL"])

        return x_name


def read_pulses_from_hdf5(
    hdf5_file: str,
    column: int = None,
    pixel: int = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read DEMUX pulse data from an HDF5 file.

    Parameters
    ----------
    hdf5_file : str
        Name of the HDF5 file (includes the path).
    column : int, optional
        If given, keep only the pulses of this column (0 to 3).
    pixel : int, optional
        If given, keep only the pulses of this pixel (0 to 33).

    Returns
    -------
    frame_num : np.ndarray
        Array of frame numbers (one per selected pulse).
    columns : np.ndarray
        Array of column ids (one per selected pulse).
    pixels : np.ndarray
        Array of pixel ids (one per selected pulse).
    pulses : np.ndarray
        Array of shape (M, S) with the S samples of each selected pulse.
    """
    with h5py.File(hdf5_file, "r") as f:
        columns = np.array(f["Pulses/col"][:])
        pixels = np.array(f["Pulses/pixel"][:])

        selected = np.ones(len(columns), dtype=bool)
        if column is not None:
            selected &= columns == column
        if pixel is not None:
            selected &= pixels == pixel

        indices = np.flatnonzero(selected)
        frame_num = np.array(f["Pulses/FrameNum"][indices])
        columns = columns[indices]
        pixels = pixels[indices]
        pulses = np.array(f["Pulses/pulse"][indices, :])

        return frame_num, columns, pixels, pulses


def _to_python_value(value: Any) -> Any:
    """Convert an HDF5 value into a plain Python object.

    Parameters
    ----------
    value : Any
        The value read from the HDF5 file.

    Returns
    -------
    Any
        The value converted to a plain Python type.
    """
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.ndarray):
        return [_to_python_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _read_config_group(group: Any) -> dict:
    """Recursively read an HDF5 configuration group.

    Parameters
    ----------
    group : h5py.Group
        The HDF5 group to read.

    Returns
    -------
    dict
        The group content as nested dictionaries.
    """
    result: dict = {}
    if group.attrs:
        result["_attributes"] = {
            name: _to_python_value(value)
            for name, value in group.attrs.items()
        }
    for name, item in group.items():
        if isinstance(item, h5py.Group):
            result[name] = _read_config_group(item)
        else:
            result[name] = _to_python_value(item[()])
    return result


def read_instrument_configuration(
    hdf5_filename: str, required: bool = False
) -> Any:
    """Read the instrument configuration from an HDF5 file.

    The values stored in the group ``/Configuration`` are already decoded:
    fractional parameters are exposed as floats and the others keep their
    integer type.

    Parameters
    ----------
    hdf5_filename : str
        Path to the HDF5 file.
    required : bool, optional
        If True, raise a ``KeyError`` when the file has no ``/Configuration``
        group. If False (default), return ``None`` instead (e.g. for old v1
        scans and dumps).

    Returns
    -------
    dict or None
        The content of ``/Configuration`` as a Python dictionary, or ``None``
        if the group is absent and ``required`` is False.

    Raises
    ------
    KeyError
        If the file does not contain a ``/Configuration`` group and
        ``required`` is True.
    """
    with h5py.File(hdf5_filename, "r") as h5:
        if "Configuration" not in h5:
            if required:
                raise KeyError(
                    f"{hdf5_filename!r} does not contain /Configuration"
                )
            return None
        return _read_config_group(h5["Configuration"])


def detect_hdf5_type(filename: str) -> str:
    """Detect the type of an HDF5 file.

    Parameters
    ----------
    filename : str
        Path to the HDF5 file.

    Returns
    -------
    str
        One of ``"pulses"``, ``"scan"``, ``"science"``, ``"dump"`` or
        ``"inconnu"``.
    """
    with h5py.File(filename, "r") as h5:
        if "Pulses" in h5:
            return "pulses"
        if "x" in h5 and "pixels" in h5:
            return "scan"
        if "ctrl" in h5 and "pixels" in h5:
            return "science"
        if all(f"Col{i}" in h5 for i in range(4)) and "Errors" in h5:
            return "dump"
        return "inconnu"


def hks_exist(hk_file: str, encoding: str = "latin1") -> bool:
    """Return True if the HK CSV file contains data beyond the header.

    Parameters
    ----------
    hk_file : str
        Path to the CSV file.
    encoding : str, optional
        Encoding of the file (default is "latin1").

    Returns
    -------
    bool
        True if there is more than one row in the file.
    """
    df = pd.read_csv(hk_file, sep=";", encoding=encoding)

    # If len(df) == 1 there are only the HK names in the file
    return len(df) > 1


def read_hk_name_from_csv(
    hk_file: str, hk_suffix: str, encoding: str = "latin1"
) -> pd.Series:
    """Read a column of a CSV file matching the end of its name.

    Parameters
    ----------
    hk_file : str
        Path to the CSV file.
    hk_suffix : str
        The last characters of the column name to look for.
    encoding : str, optional
        Encoding of the file (default is "latin1").

    Returns
    -------
    pd.Series
        The matching column.

    Raises
    ------
    ValueError
        If no or more than one column matches the suffix.
    """
    df = pd.read_csv(hk_file, sep=";", encoding=encoding)

    # Find the columns whose name ends with hk_suffix
    matching_columns = [
        col for col in df.columns if col.endswith(hk_suffix)
    ]

    if not matching_columns:
        raise ValueError(f"No HK matches the name '{hk_suffix}'.")
    if len(matching_columns) > 1:
        names = ", ".join(matching_columns)
        raise ValueError(
            f"More than one HK matches the name '{hk_suffix}': {names}"
        )

    selected_column = matching_columns[0]

    # Convert to datetime if the column name starts with 'Date'
    if selected_column.startswith("Date"):
        df[selected_column] = pd.to_datetime(
            df[selected_column], dayfirst=True, errors="coerce"
        )

    return df[selected_column]


def read_fwVersion_dmxModel(path: str) -> Tuple[str, int, Any]:
    """Read the firmware version, DMX model and board id from the HK files.

    Parameters
    ----------
    path : str
        Directory containing the HK files.

    Returns
    -------
    dmx_model : str
        The DEMUX model.
    board_id : int
        The board id.
    fw_version : Any
        The firmware version.
    """
    # Looking for the HK files
    files = [
        f
        for f in os.listdir(path)
        if os.path.isfile(os.path.join(path, f))
        and f.startswith("Hks_DMXA")
        and f.endswith(".csv")
    ]

    if len(files) == 0:
        raise ValueError("No HK files found")

    if hks_exist(os.path.join(path, files[0])):
        fw_version = read_hk_name_from_csv(
            os.path.join(path, files[0]), "Firmware Version"
        )[0]
        ref = read_hk_name_from_csv(
            os.path.join(path, files[0]), "Hardware Version"
        )[0]

        dmx_model_id = (ref >> 8) & 3
        dmx_model = cst.dmx_models[dmx_model_id]
        board_id = ref & (2 ** 5) - 1

    else:
        print("HK file " + files[0] + " is empty")
        dmx_model = "99999"
        board_id = 99
        fw_version = 99

    return dmx_model, board_id, fw_version


def read_dmxConfig_fromXml(path: str) -> dict:
    """Read the DEMUX configuration from an XML file.

    Parameters
    ----------
    path : str
        Directory containing the XML file.

    Returns
    -------
    dict
        The DEMUX configuration parameters.
    """
    config_keys = [
        "fw_version",
        "hw_version",
        "boxcar_length",
        "c0_pulse_shaping_set",
        "c0_offset_coarse",
        "c0_offset_lsb",
        "c0_sampling_delay",
        "c0_feedback_delay",
        "c0_offset_dac_delay",
        "c0_offset_mux_delay",
        "relock_delay",
        "relock_threshold",
    ]
    dict_conf = {key: "" for key in config_keys}

    # Looking for the XML files
    files = [
        f
        for f in os.listdir(path)
        if os.path.isfile(os.path.join(path, f)) and f.endswith(".xml")
    ]

    if len(files) != 1:
        raise ValueError(
            f"Wrong number of xml files: expected 1, found {len(files)}"
        )

    # Parsing the XML file
    file = minidom.parse(os.path.join(path, files[0]))
    dmx = file.getElementsByTagName("dmx")

    if dmx:
        attributes = dmx[0].attributes
        for key in config_keys:
            if key in attributes:
                dict_conf[key] = attributes[key].value

    return dict_conf
