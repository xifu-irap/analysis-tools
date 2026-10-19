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
#  nonLinearity.py
#
# ---------------------------------------------------------------------------------

# imports
import os
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

import constants as cst
import general_tools as gt
import readData as rddt

# ----------------------------------------------------------------------
# Plotting parameters of the non-linearity test
# ----------------------------------------------------------------------
COLORS = ['steelblue', 'turquoise', 'blue', 'teal']
REQ_LEVEL = 0.5  # required level in % of the ADC FSR
DOT_SIZE = 1
YMARGIN = 0
YLIM_INL = 150  # INL plot y limit in LSB
YLIM_DNL = 20  # DNL plot y limit in LSB

# x grid steps
XMAJ = 2 ** 12
XMIN = 2 ** 7

# y grid steps for the INL plots
YMAJ_INL_PC = 0.25
YMIN_INL_PC = 0.05

# y grid steps for the DNL plots
YMAJ_DNL_PC = 0.05
YMIN_DNL_PC = 0.01

FEEDBACK = "Feedback"
OFFSET = "Offset"


# ----------------------------------------------------------------------
def set_grid(ax, major_on: bool, minor_on: bool, xmajor: float, xminor: float,
             ymajor: float, yminor: float) -> None:
    """Plot major and minor grids on a sub-plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The sub-plot axes on which to draw the grids.
    major_on : bool
        If True, the major grid is plotted.
    minor_on : bool
        If True, the minor grid is plotted.
    xmajor : float
        Spacing of the x major grid.
    xminor : float
        Spacing of the x minor grid.
    ymajor : float
        Spacing of the y major grid.
    yminor : float
        Spacing of the y minor grid.

    Returns
    -------
    None
    """
    # set the major grid spacing
    ax.xaxis.set_major_locator(MultipleLocator(xmajor))
    ax.yaxis.set_major_locator(MultipleLocator(ymajor))

    # set the minor grid spacing
    ax.xaxis.set_minor_locator(MultipleLocator(xminor))
    ax.yaxis.set_minor_locator(MultipleLocator(yminor))

    # enable the major grid
    ax.grid(major_on, which='major', linestyle='-', linewidth='0.8',
            color='black')

    # enable the minor grid
    if minor_on:
        ax.minorticks_on()
        ax.grid(True, which='minor', linestyle=':', linewidth='0.5',
                color='gray')


# ----------------------------------------------------------------------
def get_session_info() -> Tuple[str, int, str, str, str]:
    """Extract the session configuration from the directory structure.

    The function derives the session name, the boxcar length, and the paths
    to the data, housekeeping and plot directories. It also creates the plot
    directory if it does not exist.

    Returns
    -------
    session_name : str
        The session name.
    bxl : int
        The boxcar length extracted from the session name.
    path_data : str
        Path to the scan data directory.
    path_hk : str
        Path to the housekeeping data directory.
    path_plot : str
        Path to the plot directory.
    """
    # data directory
    dirpath = os.path.join("..", "..")

    session_name = os.path.basename(os.path.realpath(dirpath))

    # looking for the boxcar length in the session name
    index_bxl = session_name.find("BXL") + len("BXL")
    bxl = int(session_name[index_bxl])

    path_data = os.path.join(dirpath, cst.scanDirName)
    path_hk = os.path.join(dirpath, cst.hkDirName)
    path_plot = os.path.join(dirpath, cst.plotDirName)
    gt.createdir(path_plot)

    return session_name, bxl, path_data, path_hk, path_plot


# ----------------------------------------------------------------------
def get_scan_files(path_data: str, col: int = None) -> List[str]:
    """List the scan HDF5 files of a given column.

    Parameters
    ----------
    path_data : str
        Path to the scan data directory.
    col : int, optional
        Column number to filter on. If None, all scan files are returned.

    Returns
    -------
    list of str
        The names of the scan files.
    """
    files = [f for f in os.listdir(path_data)
             if os.path.isfile(os.path.join(path_data, f))
             and f.startswith("scan_")
             and f.endswith(".h5")]
    if col is not None:
        files = [f for f in files if f.endswith("{0:}.h5".format(col))]
    return files


# ----------------------------------------------------------------------
def last_of_each_plateau(scan_values: np.ndarray) -> np.ndarray:
    """Return the index of the last sample of each constant-scan plateau.

    Parameters
    ----------
    scan_values : numpy.ndarray
        The scan (DAC) values.

    Returns
    -------
    numpy.ndarray
        Indices of the last sample of each plateau.
    """
    if len(scan_values) == 0:
        return np.array([], dtype=int)
    change = np.empty(len(scan_values), dtype=bool)
    change[0] = True
    change[1:] = scan_values[1:] != scan_values[:-1]
    last_idx = np.where(change)[0] - 1
    last_idx = np.append(last_idx, len(scan_values) - 1)
    return last_idx[last_idx >= 0]


# ----------------------------------------------------------------------
def read_scan_data(path_data: str, files: List[str],
                   verbose: bool) -> Tuple[np.ndarray, np.ndarray, str]:
    """Read and aggregate the scan data of the given files.

    For each scan step the last frame is kept, and the error is averaged over
    the pixels of the frame.

    Parameters
    ----------
    path_data : str
        Path to the scan data directory.
    files : list of str
        Scan file names to read.
    verbose : bool
        If True, print the name of each read file.

    Returns
    -------
    scan : numpy.ndarray
        Aggregated scan (DAC) values.
    error : numpy.ndarray
        Aggregated error values.
    scan_type : str
        The scan type ('Feedback' or 'Offset').
    """
    all_scan = []
    all_error = []
    scan_type = None

    for file in files:
        if verbose:
            print("Reading data from file ", file)
        x_name, _, file_scan, file_error = rddt.read_scan(
            os.path.join(path_data, file))
        scan_type = x_name
        # average the error over the pixels of each frame
        file_error = np.mean(file_error[:, :], axis=0)

        # keep the last frame of each scan step
        last_idx = last_of_each_plateau(file_scan)
        all_scan.append(file_scan[last_idx])
        all_error.append(file_error[last_idx])

    scan = np.concatenate(all_scan)
    error = np.concatenate(all_error)

    # remove duplicate scan values, keeping the first occurrence
    unique_values, unique_i = np.unique(scan, return_index=True)
    scan = scan[unique_i]
    error = error[unique_i]

    return scan, error, scan_type


# ----------------------------------------------------------------------
def get_plot_config(scan_type: str, col: int, bxl: int, dmx_model: str,
                    board_id: int, fw_version: int,
                    session_name: str) -> Dict[str, object]:
    """Build the plotting parameters for the given scan type and column.

    Parameters
    ----------
    scan_type : str
        The scan type ('Feedback' or 'Offset').
    col : int
        The DEMUX column number.
    bxl : int
        The boxcar length.
    dmx_model : str
        The DEMUX model.
    board_id : int
        The board identifier.
    fw_version : int
        The firmware version.
    session_name : str
        The session name.

    Returns
    -------
    dict
        A dictionary containing the x label, the x/y limits and the plot
        titles and file names.
    """
    config = {"xtit": None, "xlim": None, "ylim": None,
              "data_plot_filename": None, "dnl_plot_filename": None,
              "inl_plot_filename": None, "fig_suptitle": None,
              "fig_suptitle_dnl": None, "fig_suptitle_inl": None}

    base_suptitle = ('(' + dmx_model + " {0:}".format(board_id)
                     + ", Firmware version: {0:}, ".format(fw_version)
                     + session_name + ')')

    if scan_type == FEEDBACK:
        config["xtit"] = "Feedback values (ADU)"
        config["xlim"] = [-cst.fsrDACFdbkADU / 2, cst.fsrDACFdbkADU / 2]
        config["ylim"] = [-cst.fsrADCErrorADU / 2 - YMARGIN,
                          cst.fsrADCErrorADU / 2 + YMARGIN]
        type_text1 = "fdbkAndError"
        type_text2 = "Feedback + Error"
    elif scan_type == OFFSET:
        config["xtit"] = "Offset values (ADU)"
        config["xlim"] = [0, cst.fsrDACOfcoCoarseADU]
        config["ylim"] = [0 - YMARGIN, cst.fsrADCErrorADU / 2 + YMARGIN]
        type_text1 = "ofcoAndError"
        type_text2 = "Offset + Error"
    else:
        raise ValueError("Unknown scan type: {0!r}".format(scan_type))

    config["data_plot_filename"] = (
        type_text1+'_col{0:}_bxl{1:}.png'.format(col, bxl))
    config["dnl_plot_filename"] = (
        type_text1+'DNL_col{0:}_bxl{1:}.png'.format(col, bxl))
    config["inl_plot_filename"] = (
        type_text1+'INL_col{0:}_bxl{1:}.png'.format(col, bxl))
    config["fig_suptitle"] = (
        type_text2 + ' non linearity for column {0:}\n'.format(col)
        + base_suptitle)
    config["fig_suptitle_dnl"] = (
        type_text2 + ' DNL for column {0:}\n'.format(col)
        + base_suptitle)
    config["fig_suptitle_inl"] = (
        type_text2 + ' INL for column {0:} (Scan - linear fit)\n'
        .format(col) + base_suptitle)

    return config


# ----------------------------------------------------------------------
def plot_scan_data(scan: np.ndarray, error: np.ndarray, fit: np.ndarray,
                   i_ok: np.ndarray, coeffs: np.ndarray,
                   config: Dict[str, object], path_plot: str) -> None:
    """Plot the output versus input non-linearity data.

    Parameters
    ----------
    scan : numpy.ndarray
        Sorted scan (DAC) values.
    error : numpy.ndarray
        Sorted error values.
    fit : numpy.ndarray
        Linear fit of the error values.
    i_ok : numpy.ndarray
        Indices of the non-saturated data.
    coeffs : numpy.ndarray
        Coefficients of the linear fit.
    config : dict
        Plot configuration (see get_plot_config).
    path_plot : str
        Path to the plot directory.

    Returns
    -------
    None
    """
    fig0 = plt.figure(figsize=(9, 6))
    plot_full_filename = os.path.join(path_plot,
                                      config["data_plot_filename"])
    fig0.suptitle(config["fig_suptitle"], fontsize=14)
    ax0 = fig0.add_subplot(1, 1, 1)  # output vs input

    ax0.scatter(scan, error, s=DOT_SIZE, c='k', label='Scan')
    if coeffs[1] > 0:
        sign_str = ' + '
    else:
        sign_str = ' - '
    lbl = ('Linear fit (Y = {0:.4} X'.format(coeffs[0]) + sign_str
           + '{0:.4})'.format(abs(coeffs[1])))
    ax0.plot(scan[i_ok], fit, ':', color=COLORS[0], linewidth=1, label=lbl)
    ax0.set_xlim(config["xlim"])
    ax0.set_ylim(config["ylim"])
    ax0.set_xlabel(config["xtit"])
    ax0.set_ylabel("Error values (ADU)")
    ax0.legend(loc='upper left')
    set_grid(ax0, True, True, XMAJ, XMIN, XMAJ, XMIN)

    fig0.tight_layout()
    plt.savefig(plot_full_filename, dpi=300, bbox_inches='tight')
    print("Linearity data plotted in file "
          + config["data_plot_filename"])


# ----------------------------------------------------------------------
def plot_inl(scan: np.ndarray, deviation_lsb: np.ndarray, i_ok: np.ndarray,
             config: Dict[str, object], path_plot: str) -> None:
    """Plot the Integral Non-Linearity (INL) results.

    Parameters
    ----------
    scan : numpy.ndarray
        Sorted scan (DAC) values.
    deviation_lsb : numpy.ndarray
        Deviation from the linear fit in LSB.
    i_ok : numpy.ndarray
        Indices of the non-saturated data.
    config : dict
        Plot configuration (see get_plot_config).
    path_plot : str
        Path to the plot directory.

    Returns
    -------
    None
    """
    fig1 = plt.figure(figsize=(9, 6))
    plot_full_filename = os.path.join(path_plot,
                                      config["inl_plot_filename"])
    fig1.suptitle(config["fig_suptitle_inl"], fontsize=14)
    ax1 = fig1.add_subplot(1, 1, 1)  # non linearity

    val1 = max(np.abs(deviation_lsb))
    val2 = val1 * 100 / cst.fsrADCErrorADU
    lbl = 'The INL is {0:2.1f} LSB or {1:.2} %)'.format(val1, val2)
    ax1.scatter(scan[i_ok], deviation_lsb, s=DOT_SIZE, color=COLORS[0],
                label=lbl)
    ax1.set_xlim(config["xlim"])
    ax1.set_xlabel(config["xtit"])
    ax1.set_ylabel("Non linearity (ADC LSB)")

    ax1.set_ylim([-1 * YLIM_INL, YLIM_INL])
    ax1.legend(loc='upper left')

    # second y axis for percentage units
    ax11 = ax1.twinx()
    ylims = ax1.get_ylim()
    ylims11 = [ylims[0] * 100 / cst.fsrADCErrorADU,
               ylims[1] * 100 / cst.fsrADCErrorADU]
    ax11.set_ylim(ylims11)
    ax11.set_ylabel("Non Linearity (% of ADC FSR)")
    ax11.plot(config["xlim"], [REQ_LEVEL, REQ_LEVEL], color='r',
              linewidth=2)
    ax11.plot(config["xlim"], [-REQ_LEVEL, -REQ_LEVEL], color='r',
              linewidth=2)

    ymaj_inl_adu = YMAJ_INL_PC * cst.fsrADCErrorADU / 100
    ymin_inl_adu = YMIN_INL_PC * cst.fsrADCErrorADU / 100
    set_grid(ax1, True, True, XMAJ, XMIN, ymaj_inl_adu, ymin_inl_adu)

    fig1.tight_layout()
    plt.savefig(plot_full_filename, dpi=300, bbox_inches='tight')
    print("INL results plotted in file " + config["inl_plot_filename"])


# ----------------------------------------------------------------------
def plot_dnl(scan: np.ndarray, dnl: np.ndarray, i_ok: np.ndarray,
             config: Dict[str, object], path_plot: str) -> None:
    """Plot the Differential Non-Linearity (DNL) results.

    Parameters
    ----------
    scan : numpy.ndarray
        Sorted scan (DAC) values.
    dnl : numpy.ndarray
        Differential non-linearity values.
    i_ok : numpy.ndarray
        Indices of the non-saturated data.
    config : dict
        Plot configuration (see get_plot_config).
    path_plot : str
        Path to the plot directory.

    Returns
    -------
    None
    """
    fig2 = plt.figure(figsize=(9, 6))
    plot_full_filename = os.path.join(path_plot,
                                      config["dnl_plot_filename"])
    fig2.suptitle(config["fig_suptitle_dnl"], fontsize=14)
    ax2 = fig2.add_subplot(1, 1, 1)

    val1 = max(np.abs(dnl))
    val2 = val1 * 100 / cst.fsrADCErrorADU
    lbl = 'The DNL is {0:2.1f} LSB or {1:.2} %)'.format(val1, val2)
    ax2.scatter(scan[i_ok][:-1], dnl, s=DOT_SIZE, color=COLORS[0], label=lbl)
    ax2.plot(config["xlim"], [0, 0], '-k', linewidth=0.5)
    ax2.set_xlim(config["xlim"])

    ax2.set_ylim([-1 * YLIM_DNL, YLIM_DNL])
    ax2.set_xlabel(config["xtit"])
    ax2.set_ylabel("Non linearity (ADC LSB)")
    ax2.legend(loc='upper left')

    ax22 = ax2.twinx()
    ylims = ax2.get_ylim()
    ylims22 = [ylims[0] * 100 / cst.fsrADCErrorADU,
               ylims[1] * 100 / cst.fsrADCErrorADU]
    ax22.set_ylim(ylims22)
    ax22.set_ylabel("Non Linearity (% of ADC FSR)")
    ax22.plot(config["xlim"], [REQ_LEVEL, REQ_LEVEL], color='r',
              linewidth=2)
    ax22.plot(config["xlim"], [-REQ_LEVEL, -REQ_LEVEL], color='r',
              linewidth=2)

    ymaj_dnl_adu = YMAJ_DNL_PC * cst.fsrADCErrorADU / 100
    ymin_dnl_adu = YMIN_DNL_PC * cst.fsrADCErrorADU / 100
    set_grid(ax2, True, True, XMAJ, XMIN, ymaj_dnl_adu, ymin_dnl_adu)

    fig2.tight_layout()
    plt.savefig(plot_full_filename, dpi=300, bbox_inches='tight')
    print("DNL results plotted in file " + config["dnl_plot_filename"])


# ----------------------------------------------------------------------
def nonlinearity(verbose: bool = False) -> None:
    """Characterise the non-linearity of DMX signals from scan data.

    This is applicable to feedback + error signals or offset + error signals.
    The type of scan is detected from the name of the scan files. The result
    is a set of plots saved in the /PLOTS directory of the session. The
    function processes the scans of the different columns.

    Parameters
    ----------
    verbose : bool, optional
        If True, print the name of each read file. Default is False.

    Returns
    -------
    None
    """
    # session configuration
    session_name, bxl, path_data, path_hk, path_plot = get_session_info()

    # looking for DEMUX identifiers (board, model, firmware)
    dmx_model, board_id, fw_version = rddt.read_fwVersion_dmxModel(path_hk)

    # checking the number of scan files
    files = get_scan_files(path_data)
    if len(files) == 0:
        raise ValueError("Error, no scan files found in this session !")

    # checking the scan type (feedback or offset)
    scan_type_ini = rddt.read_scan_type(os.path.join(path_data, files[0]))
    if scan_type_ini == FEEDBACK:
        scan_type_text = "FEEDBACK"
    elif scan_type_ini == OFFSET:
        scan_type_text = "OFFSET"
    else:
        raise ValueError("Unknown scan type: {0!r}".format(scan_type_ini))

    print("/----------------------------------------------------------")
    print("/ Non linearity test: ERROR + " + scan_type_text)
    print("/ Test session name: " + session_name)
    print("/----------------------------------------------------------")
    print("/ DEMUX model:       " + dmx_model + " {0:}".format(board_id))
    print("/ Firmware version:  {0:}".format(fw_version))
    print("/ Box car length:    {0:} samples".format(bxl))
    print("/----------------------------------------------------------\n")

    for col in range(cst.nColPerDemux):
        # searching the scan files of this column
        files = get_scan_files(path_data, col)
        if len(files) == 0:
            raise ValueError(
                "Error, no scan files found for column {0:}!".format(col))

        print("Found {0:} HDF5 file(s) for column {1:}"
              .format(len(files), col))

        # reading and aggregating the scan data
        scan, error, scan_type = read_scan_data(path_data, files, verbose)

        # sorting the data with respect to the DAC values
        sorted_i = np.argsort(scan)
        scan = scan[sorted_i]
        error = error[sorted_i]

        # checking if all DAC values are used
        expected_array = np.arange(scan.min(), scan.max() + 1)
        if not (expected_array == scan).all():
            print("   Error, values are missing in the scan!")

        # ignoring saturations
        i_ok = np.where(
            np.abs(error) < 2 ** (cst.dmxNbBitsADCError - 1) - 1)[0]

        # linear fit of the data (full range)
        coeffs = np.polyfit(scan[i_ok], error[i_ok], 1)
        fit = coeffs[1] + coeffs[0] * scan[i_ok]
        deviation_lsb = error[i_ok] - fit

        # computing DNL and INL
        lsb_ideal = ((error[i_ok].max() - error[i_ok].min())
                     / len(error[i_ok]))
        dnl = ((error[i_ok][1:] - error[i_ok][:-1]) / lsb_ideal - 1)

        # building the plot configuration
        config = get_plot_config(scan_type, col, bxl, dmx_model, board_id,
                                 fw_version, session_name)

        # doing the plots
        plot_scan_data(scan, error, fit, i_ok, coeffs, config, path_plot)
        plot_inl(scan, deviation_lsb, i_ok, config, path_plot)
        plot_dnl(scan, dnl, i_ok, config, path_plot)

        print("/---------------")


# -------------------------------------------------------------------------------------
