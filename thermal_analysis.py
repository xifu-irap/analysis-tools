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

import numpy as np
import matplotlib.pyplot as plt
import glob

import general_tools as gentools

path = "/Users/laurent/Data/20.001/test_hk_temp_2/"
#path = "/Users/laurent/Data/20.001/test_hk_temp_3/"


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

    for m in range(1, max_m):
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


def analyseThermalAccuracy(path, start=0):
    """
    Analyse la precision des mesures thermiques par calcul de la variance d'Allan.
    Les résultats sont tracés dans une figure.

    :param path: chemin vers les fichier de données.
    :start: index pour démarrer la prise en compte des données (par défaut 0)
    """

    files = glob.glob(path + "*.csv")
    plotfilename = path+'tempSensorAccuracy.png'

    nFiles = len(files)
    fig = plt.figure(figsize=(6*nFiles, 10))

    for iFile in range(nFiles):

        print("Reading temperature data from file " + files[iFile])
        time, hk_adu = gentools.readHkFromCsv(files[iFile])
        # average time step
        timeStep = (time - np.roll(time,1))[1:].mean()
        # sampling frequency
        fs = 1/timeStep
        print("Sampling frequency is {0:4.2f} Hz".format(fs))

        hk_adu = hk_adu[start:]
        time = time[start:]
        if (files[iFile][-16:] == '-HK_TEMP_MAX.csv'):
            a = 0.02159
            b = -33.8663
            title1 = 'Hk MAX temperature'
            yLabel = 'K rms'
        elif (files[iFile][-16:] == '-HK_TEMP_AVE.csv'):
            a = 0.03150
            b = -54.59583
            title1 = 'Hk AVE temperature'
            yLabel = 'K rms'
        else:
            a = 1
            b = 0
            title1 = 'Unknown HK'
            yLabel = 'ADU rms'

        hk_temp = hk_adu * a + b

        # Calcul de la variance d'Allan
        tau_values, allan_vars, allan_errors = calculate_allan_variance(hk_temp, fs, 100)

        # Tracé des résultats

        xreq = 1e3
        yreq = 20e-3
        ymin = 1e-3
        tmin = 0.8
        tmax = 2000
        tmaxLineThick = 2

        ax1 = fig.add_subplot(2, nFiles, iFile+1)
        ax1.plot(time, hk_temp)
        ax1.set_title(title1)
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Temperature (°C)")

        title2 = 'Allan stdev'
        ax2 = fig.add_subplot(2, nFiles, nFiles+iFile+1)
        ax2.errorbar(tau_values, np.sqrt(allan_vars), yerr=allan_errors, fmt='o', capsize=5)
        ax2.loglog(tau_values, np.sqrt(allan_vars))
        ax2.loglog([tmin, tmax], [yreq, yreq], '--', color='red')
        ax2.loglog([xreq, xreq], [1e-9, 1e3], '--', color='k', linewidth = tmaxLineThick)
        ax2.set_xlabel('Average duration (s)')
        ax2.set_ylabel(yLabel)
        ax2.set_title(title2)
        ax2.set_xlim([tmin, tmax])
        stdMax = 2*np.sqrt(allan_vars[0])
        ax2.set_ylim([ymin, min(10e-1, stdMax)])
        ax2.grid(True, which="both", ls="--")

    fig.tight_layout()
    plt.savefig(plotfilename, dpi=300, bbox_inches='tight')

# ---------------------------------------------------------------------------------

analyseThermalAccuracy(path, 600)
