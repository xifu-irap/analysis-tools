# imports
import os
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

import constants as cst
import general_tools as gt
import readData as rddt


# Détection d'un front de montée
def riseDetect(data, thresholdPercent=5):
    threshold = (data.max() - data.min()) * thresholdPercent / 100 + data.min()

    decal = 30
    i_sup = np.where(data > threshold)[0][0]
    if i_sup == 0:  # Shiftting a bit if the rising edgae was before the beginning
        i_sup = np.where(data[decal:] > threshold)[0]

    i = i_sup - 1

    return (i)


# Définition de la fonction exponentielle
def exponential_response(t, U0, tau):
    return U0 * (1 - np.exp(-t / tau))


# Fonction de fit pour trouver la constante de temps
def fit_time_constant(t_data, v_data):
    # Ajustement des données pour trouver les paramètres U0 et tau
    popt, pcov = curve_fit(exponential_response, t_data, v_data, p0=[max(v_data), t_data[np.argmax(v_data)]])
    U0, tau = popt
    return U0, tau


def cutoffFreq(tconf):
    # Test configuration data
    dir_path = tconf.file_path

    # Data directory
    path_data = os.path.join(dir_path, cst.dataDirName)

    print("Processing dump files from directory " + path_data)

    # Creation of a directory for the plot files
    path_plot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(path_plot)

    # Plot labels
    xlabel = r'Time (ns)'
    ylabel = r'ADC data (V)'

    # Getting data
    files = [f for f in os.listdir(path_data) \
             if os.path.isfile(os.path.join(path_data, f)) \
             and f[-5:] == ".fits" and f[:5] == "dump_"]

    if len(files) == 0:
        raise ValueError('Error, no dump files found')

    lenDumps = 2 * cst.nPixPerCol * cst.nSamplesPerRow

    accuDumps = np.zeros((cst.nColPerDemux, lenDumps))
    for index, file in enumerate(np.sort(files)):
        print(file)
        colDumps, errors = rddt.readDumpFile(os.path.join(path_data, file))
        accuDumps += colDumps
    accuDumps = (accuDumps / len(files)) * (cst.fsrADCErrorV / cst.fsrADCErrorADU)

    t = np.arange(lenDumps) / cst.fSamp

    # Doing plots
    for col in range(cst.nColPerDemux):
        plotFileName = os.path.join(path_plot, 'cutoffFreq_col{0:}.png'.format(col))

        v = accuDumps[col, :]

        tsize = 20
        it1 = riseDetect(v)
        it2 = it1 + tsize

        fig = plt.figure(figsize=(8, 6))
        suptitle = "Error cutoff frequency (column {0:})".format(col)
        title = tconf.testPlanPath
        fig.suptitle(suptitle, fontsize=12)

        ax1 = fig.add_subplot(2, 1, 1)  # global plot
        ax1.set_title(title, fontsize=10)
        ax2 = fig.add_subplot(2, 1, 2)  # global plot

        ax1.plot(t * 1e9, v)
        ax1.plot(t[it1: it2] * 1e9, v[it1: it2], color='r')
        ax1.set_xlabel(xlabel)
        ax1.set_ylabel(ylabel)
        ax1.grid()

        # Calcul du fit
        v = v[it1:it2] - v[it1]
        tfit = t[it1:it2] - t[it1]
        U0, tau = fit_time_constant(tfit, v)
        fc = 1 / (2 * np.pi * tau * 1e6)
        print("Time constant: {0:6.3f} ns".format(tau * 1e9))
        print("Cutoff frequency: {0:6.3f} MHz".format(fc))

        lbl1 = 'ADC Data'
        ax2.plot(tfit * 1e9, v, color='r', label=lbl1)
        lbl2 = "First order fit (fc = {0:6.2f} MHz)".format(fc)
        # Building a higher resolution time array to plot the fit
        hr_ratio = 10
        thr = np.arange(hr_ratio * tsize) / (cst.fSamp * hr_ratio)
        ax2.plot(thr * 1e9, exponential_response(thr, U0, tau), '--', color='k', label=lbl2)
        ax2.set_xlabel(xlabel)
        ax2.set_ylabel(ylabel)
        ax2.legend(loc='best')
        ax2.grid()

        fig.tight_layout()

        plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
        print("results plotted in file " + plotFileName)


# -------------------------------------------------------------------------------------

@dataclass
class TestConfig:
    testPlanPath: str
    session_name: str

    @property
    def file_path(self) -> str:
        return f"{cst.BASE_DATA_PATH}/{self.testPlanPath}/{self.session_name}"


TP27_TURBO45_PATH = "TestPlan27_DM-DMX2_Func_and_Perfs/FW-turbo-45"

list_of_configs = [
    TestConfig(TP27_TURBO45_PATH, "20250618_182150_caracErrorBandShape")
]

# -------------------------------------------------------------------------------------
for conf in list_of_configs[-1:]:
    cutoffFreq(conf)

# -------------------------------------------------------------------------------------
