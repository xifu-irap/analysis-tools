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
#  cutoffFrequency.py
#
# ---------------------------------------------------------------------------------

# imports
import os
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

import constants as cst
import general_tools as gt
import readData as rddt


def get_config() -> Dict[str, str]:
    """Extract a configuration dictionary from the current session.

    The session name is derived from the directory path and is used to
    extract configuration details such as the signal name.

    Returns
    -------
    dict
        A dictionary containing the configuration:
        - session_name (str): the session name.
        - dir_path (str): the base directory path.
        - signal_name (str): a substring representing the signal.
    """
    # data directory
    dir_path = os.path.join("..", "..")

    # looking for the session name
    session_name = os.path.basename(os.path.realpath(dir_path))

    # looking for the signal name
    signal_name = session_name[:4]

    config = {"session_name": session_name,
              "dir_path": dir_path,
              "signal_name": signal_name,
              }

    return config


def riseDetect(signal: np.ndarray, threshold_pc: float = 5) -> int:
    """Detect the first rising edge of a two-level noisy signal.

    Parameters
    ----------
    signal : array-like
        Two-level signal (with noise).
    threshold_pc : float, optional
        Percentage of the signal amplitude used as the detection
        threshold. Default is 5.

    Returns
    -------
    int
        Index of the first detected rising edge.
    """
    signal = np.array(signal)

    # estimation of the two levels
    threshold = (signal.min()
                 + threshold_pc * (signal.max() - signal.min()) / 100)

    # first search for a phase where the signal is at the low state
    bas = np.where(signal < threshold)[0][0] + 5  # +5 is a margin
    # search for the following phase where the signal is at the high state
    haut = np.where(signal[bas:] > threshold)[0][0] + bas

    return haut


def exponential_response(t: np.ndarray, U0: float,
                         tau: float) -> np.ndarray:
    """Compute a first-order exponential response.

    Parameters
    ----------
    t : numpy.ndarray
        Time values.
    U0 : float
        Asymptotic amplitude.
    tau : float
        Time constant.

    Returns
    -------
    numpy.ndarray
        The exponential response.
    """
    return U0 * (1 - np.exp(-t / tau))


def double_exponential_response(t: np.ndarray, U0: float, tau1: float,
                                tau2: float) -> np.ndarray:
    """Compute a double exponential response.

    Parameters
    ----------
    t : numpy.ndarray
        Time values.
    U0 : float
        Asymptotic amplitude.
    tau1 : float
        First time constant.
    tau2 : float
        Second time constant.

    Returns
    -------
    numpy.ndarray
        The double exponential response.
    """
    return (U0 * (1 - (tau1 * np.exp(-t / tau1)
                       - tau2 * np.exp(-t / tau2)) / (tau1 - tau2)))


def fit_time_constant(t_data: np.ndarray, v_data: np.ndarray,
                      fc_feedback: float) -> Tuple[float, float, float]:
    """Fit the data to find the time constants.

    Parameters
    ----------
    t_data : numpy.ndarray
        Time values of the data to fit.
    v_data : numpy.ndarray
        Voltage values of the data to fit.
    fc_feedback : float
        Cutoff frequency of the feedback output.

    Returns
    -------
    tuple of float
        The fitted parameters U0, tau1 and tau2.
    """
    fc1_bound1 = fc_feedback - 40e6
    fc1_bound2 = fc_feedback + 60e6
    tau1_bound1 = 1 / (2 * np.pi * fc1_bound1)
    tau1_bound2 = 1 / (2 * np.pi * fc1_bound2)
    fc2_bound1 = 5e6
    fc2_bound2 = 30e6
    tau2_bound1 = 1 / (2 * np.pi * fc2_bound1)
    tau2_bound2 = 1 / (2 * np.pi * fc2_bound2)
    # fit the data to find the parameters U0 and tau
    popt, pcov = curve_fit(double_exponential_response, t_data, v_data,
                           bounds=([0, tau1_bound2, tau2_bound2],
                                   [3, tau1_bound1, tau2_bound1]))
    U0, tau1, tau2 = popt
    print("  >>> ", popt)
    return U0, tau1, tau2


def measure_cutoffFreq(tconf: Dict[str, str], fc_feedback: float) -> None:
    """Measure and plot the cutoff frequency from dump files.

    Parameters
    ----------
    tconf : dict
        Test configuration dictionary (see get_config).
    fc_feedback : float
        Cutoff frequency of the feedback output.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If no dump files are found.
    """
    # test configuration data
    dir_path = tconf["dir_path"]

    # data directory
    path_data = os.path.join(dir_path, cst.dataDirName)

    print("Processing dump files from directory " + path_data)

    # creation of a directory for the plot files
    path_plot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(path_plot)

    # plot labels
    xlabel = r'Time (ns)'
    ylabel = r'ADC data (V)'

    # getting data
    files = [f for f in os.listdir(path_data)
             if os.path.isfile(os.path.join(path_data, f))
             and f[-3:] == ".h5" and f[:5] == "dump_"]

    if len(files) == 0:
        raise ValueError('Error, no dump files found')

    lenDumps = 2 * cst.nPixPerCol * cst.nSamplesPerRow

    accuDumps = np.zeros((cst.nColPerDemux, lenDumps))
    for index, file in enumerate(np.sort(files)):
        print(file)
        colDumps, errors = rddt.read_dump_from_hdf5(
            os.path.join(path_data, file))
        accuDumps += colDumps
    accuDumps = (accuDumps / len(files)
                 * (cst.fsrADCErrorV / cst.fsrADCErrorADU))

    t = np.arange(lenDumps) / cst.fSamp

    # doing plots
    for col in range(cst.nColPerDemux):
        plotFileName = os.path.join(
            path_plot, 'cutoffFreq_col{0:}.png'.format(col))

        v = accuDumps[col, :]

        tsize = 20
        it1 = riseDetect(v, 5) - 1
        it2 = it1 + tsize

        # computation of the fit
        vfit = v[it1:it2] - v[it1]
        tfit = t[it1:it2] - t[it1]
        U0, tau1, tau2 = fit_time_constant(tfit, vfit, fc_feedback)
        fc1 = 1 / (2 * np.pi * tau1 * 1e6)
        fc2 = 1 / (2 * np.pi * tau2 * 1e6)
        print("Time constant: {0:6.3f} ns".format(tau2 * 1e9))
        print("Cutoff frequency: {0:6.1f} MHz".format(fc2))

        fig = plt.figure(figsize=(8, 6))
        fc1_limit = 120  # in MHz
        if fc1 > fc1_limit:
            suptitle = ("FDBK (fc > {0:6.1f} MHz) + ERROR step response "
                        "(column {1:})".format(fc1_limit, col))
        else:
            suptitle = ("FDBK (fc={0:6.1f} MHz) + ERROR step response "
                        "(column {1:})".format(fc1, col))
        title = tconf["session_name"]
        fig.suptitle(suptitle, fontsize=12)

        ax1 = fig.add_subplot(2, 1, 1)  # global plot
        ax1.set_title(title, fontsize=10)
        ax2 = fig.add_subplot(2, 1, 2)  # zoomed plot

        ax1.plot(t * 1e9, v)
        ax1.plot(t[it1: it2] * 1e9, v[it1: it2], color='r')
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)
        ax1.grid()

        lbl1 = 'ADC Data'
        ax2.plot(t[it1 - 10:it2] * 1e9, v[it1 - 10:it2])
        ax2.plot(t[it1:it2] * 1e9, v[it1:it2], color='r', label=lbl1)
        lbl2 = ('FDBK + ERROR fit --> ERROR cutoff frequency is '
                '{0:6.1f} MHz)'.format(fc2))
        # building a higher resolution time array to plot the fit
        hr_ratio = 10
        thr = np.arange(hr_ratio * tsize) / (cst.fSamp * hr_ratio)
        ax2.plot((thr + t[it1]) * 1e9,
                 double_exponential_response(thr, U0, tau1, tau2) + v[it1],
                 '--', color='k', label=lbl2)
        ax2.set_xlabel(xlabel)
        ax2.set_ylabel(ylabel)
        ax2.legend(loc='best')
        ax2.grid()

        fig.tight_layout()

        plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
        print("results plotted in file " + plotFileName)


# -------------------------------------------------------------------------------------

def cutoffFreq(verbose: bool = True) -> None:
    """Measure the cutoff frequency of dump data.

    Parameters
    ----------
    verbose : bool, optional
        If True, print detailed information during execution.
        Default is True.

    Returns
    -------
    None
    """
    config = get_config()

    fc_feedback = 80e6  # cutoff frequency of feedback output

    # data directory
    hk_path = os.path.join(config["dir_path"], cst.hkDirName)

    # looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(hk_path)

    if verbose:
        print("/----------------------------------------------------------")
        print("/ ERROR band shape test:     ")
        print("/ Test session name:   " + config["session_name"])
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:    {0:}".format(fwVersion))
        print("/ Signal:              " + config["signal_name"])
        print("/----------------------------------------------------------\n")

    measure_cutoffFreq(config, fc_feedback)

# -------------------------------------------------------------------------------------
