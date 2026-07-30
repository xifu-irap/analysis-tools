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
#  thermal_analysis.py
#
# ---------------------------------------------------------------------------------

import os

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gt
import readData as rddt


def get_config():
    """
    Extracts and returns a configuration dictionary from the current session name and directory structure.

    This function analyzes the session name derived from the directory path to extract specific configuration details,
    such as boxcar length, feedback settings, and coarse-grained parameter.

    Returns:
        dict: A dictionary containing the configuration extracted from the session name:
            - session_name (str): The name of the session extracted from the directory path.
            - dir_path (str): The base directory path.
            - setup (str): A substring representing the session setup.
            - rate (str): Placeholder for rate information (currently empty).
            - bxl (int): The boxcar length extracted from the session name.
            - fdbk (str): The feedback setting identifier extracted from the session name.
            - ofco (str): The coarse-grained parameter extracted from the session name.
    """

    # Data directory
    dirpath = os.path.join("..", "..")
    pathHk = os.path.join(dirpath, cst.hkDirName)

    # Looking for the session name and test configuration : "ERRO_ONLY" or "ERRO_FDBK" or "ERRO_OFCO"
    session_name = os.path.basename(os.path.realpath(dirpath))

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(pathHk)

    config = {"session_name": session_name,
              "dir_path": dirpath,
              "dmxModel": dmxModel,
              "boardId": boardId,
              "fwVersion": fwVersion
              }

    return config


def calculate_allan_variance(data, fs, num_blocks=5):
    """
    Calcule la variance d'Allan et les barres d'erreur basées sur des sous-ensembles de données.

    :param data: Tableau 1D contenant les données temporelles.
    :param fs: Fréquence d'échantillonnage (en Hz).
    :param num_blocks: Nombre de sous-ensembles pour calculer les barres d'erreur.
    :return: tau (intervalles de temps), variance d'Allan et erreur estimée.
    """
    N = len(data)
    max_m = N // 2
    tau_values = []
    allan_vars = []
    allan_errors = []

    for m in range(2, max_m):
        tau = m / fs
        tau_values.append(tau)

        # Découper les données en blocs (num_blocks)
        block_size = N // num_blocks
        block_avgs = []
        for block in range(num_blocks):
            subset_data = data[block * block_size:(block + 1) * block_size]
            avg_data = np.array([np.mean(subset_data[i:i + m]) for i in range(0, len(subset_data) - m, m)])
            block_avgs.append(0.5 * np.mean(np.diff(avg_data) ** 2))

        # Moyenne des variances sur les sous-ensembles
        allan_var = np.mean(block_avgs)
        allan_vars.append(allan_var)

        # Calcul de l'erreur standard (erreur-type)
        allan_error = np.std(block_avgs) / np.sqrt(num_blocks)
        allan_errors.append(allan_error)

    return tau_values, allan_vars, allan_errors


def analyseThermalAccuracy(config, hk_name, date_name='Date(EGSE)', start=0, verbose=True):
    """
    Analyse la precision des mesures thermiques par calcul de la variance d'Allan.
    Les résultats sont tracés dans une figure.

    :param path: chemin vers les fichiers de données.
    :param dmxModel: string, indicates the model of the DEMUX module (dm_dmx0, dm_dmx2)
    :start: index pour démarrer la prise en compte des données (par défaut 0)
    """
    import os
    import pandas as pd

    # Data directory
    dirpath = config["dir_path"]
    session_name = config["session_name"]

    pathHk = os.path.join(dirpath, cst.hkDirName)
    pathPlot = os.path.join(dirpath, cst.plotDirName)
    plotFileName = 'thermal_stability_' + config["module"] + '_' + hk_name[:-4] + '.png'
    plotFullFileName = os.path.join(pathPlot, plotFileName)

    title = 'Thermistor stability (' + config["module"] + ', ' + hk_name[:-4] + ')\n' + session_name

    gt.createdir(pathPlot)

    files = [f for f in os.listdir(pathHk) \
             if os.path.isfile(os.path.join(pathHk, f)) \
             and f[:8] == "Hks_" + config["module"] \
             and f[-4:] == ".csv"]

    # Checking number of files
    if len(files) == 0:
        raise ValueError("Error, no HK files found in this session for module ", config["module"])

    fig = plt.figure(figsize=(7, 9))

    # Lecture du fichier CSV
    df = pd.read_csv(os.path.join(pathHk, files[0]), sep=';', encoding='cp1252')

    hk = df[hk_name]
    time_csv = df[date_name]

    time = gt.convert_time(time_csv)
    # average time step
    timeStep = (time - np.roll(time,1))[1:].mean()
    # sampling frequency
    fs = 1/timeStep
    print("Sampling frequency is {0:4.2f} Hz".format(fs))

    time -= time[0]

    # Calcul de la variance d'Allan
    tau_values, allan_vars, allan_errors = calculate_allan_variance(hk[start:], fs, 100)

    # Tracé des résultats

    xreq = 1e3
    yreq = 20e-3
    ymin = 1e-3
    tmin = 0.8
    tmax = 2000
    tmaxLineThick = 2

    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(time, hk, 'k')
    ax1.plot(time[start:], hk[start:], 'r')
    ax1.set_title(title)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Temperature (°C)")
    ax1.grid(True, which="both", ls="--")

    yLabel = 'Allan standard deviation (°C)'
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.errorbar(tau_values, np.sqrt(allan_vars), yerr=allan_errors, fmt='o', capsize=5)
    ax2.loglog(tau_values, np.sqrt(allan_vars))
    ax2.loglog([tmin, tmax], [yreq, yreq], '--', color='red')
    ax2.loglog([xreq, xreq], [1e-9, 1e3], '--', color='k', linewidth = tmaxLineThick)
    ax2.set_xlabel('Duration of averaging (s)')
    ax2.set_ylabel(yLabel)
    ax2.set_xlim([tmin, tmax])
    stdMax = 2*np.sqrt(allan_vars[0])
    ax2.set_ylim([ymin, min(10e-1, stdMax)])
    ax2.grid(True, which="both", ls="--")

    fig.tight_layout()
    plt.savefig(plotFullFileName, dpi=300, bbox_inches='tight')

    print("plot saved in file ", plotFullFileName)

# ---------------------------------------------------------------------------------

def thermal_analysis(verbose=True):
    """
    Analyzes the thermal behavior of the DEMUX.

    Args:
        verbose (bool): Determines if detailed output should be printed during
            execution.
    """

    config = get_config()
    config["module"] = "DMXA"

    if verbose:
        print("/----------------------------------------------------------")
        print("/ Thermal stability test: ")
        print("/ Test session name:      " + config["session_name"])
        print("/----------------------------------------------------------")
        print("/ DEMUX model:            " + config["dmxModel"] + " {0:}".format(config["boardId"]))
        print("/ Firmware version:       {0:}".format(config["fwVersion"]))
        print("/----------------------------------------------------------\n")

    analyseThermalAccuracy(config, "TEMP_MAX(°C)", start=0, verbose=verbose)
    analyseThermalAccuracy(config, "TEMP_AVE(°C)", start=0, verbose=verbose)

# -------------------------------------------------------------------------------------

# TODO: Add the model / module name on the plots


