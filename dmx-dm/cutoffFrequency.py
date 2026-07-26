# imports
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

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
            - signal_name (str): A substring representing the signal.
    """

    # Data directory
    dir_path = os.path.join("..", "..")

    # Looking for the session name
    session_name = os.path.basename(os.path.realpath(dir_path))

    # Looking for the signal name
    signal_name = session_name[:4]

    config = {"session_name": session_name,
              "dir_path": dir_path,
              "signal_name": signal_name,
              }

    return config


# Détection d'un front de montée
def riseDetect_old(data):
    threshold = (data.max() + data.min()) / 2

    i_sup = np.where(data > threshold)[0][0]
    i_inf = np.where(data < threshold)[0][0]

    # Shifting the data if the dump is at the high value at the beginning
    if i_sup == 0:
        data = np.roll(data, -1 * (i_inf + 1))

        i_sup = np.where(data[decal:] > threshold)[0]

    i = i_sup - 1

    return (i)

    # Détection d'un front de montée


def riseDetect(signal, threshold_pc=5):
    """
    Détecte le premier front montant d'un signal bruité à deux niveaux.

    Paramètres
    ----------
    signal : array-like
        Signal à deux niveaux (avec bruit).
    hysteresis : float
        Fraction d’hystérésis pour éviter les faux fronts.
        Par exemple 0.1 = 10% entre le niveau bas et le niveau haut.

    Retour
    ------
    int ou None
        Indice du premier front montant détecté.
    """
    signal = np.array(signal)

    # Estimation des deux niveaux
    threshold = signal.min() + threshold_pc * (signal.max() - signal.min()) / 100

    # On commence par rechercher une phase où le signal est à l'état bas
    bas = np.where(signal < threshold)[0][0] + 5  # +5 is a margin
    # On recherche la phase qui suit où le signal est à l'état haut
    haut = np.where(signal[bas:] > threshold)[0][0] + bas

    return haut

# Définition de la fonction exponentielle
def exponential_response(t, U0, tau):
    return U0 * (1 - np.exp(-t / tau))


# Fonction de fit pour trouver la constante de temps
def fit_time_constant(t_data, v_data):
    # Ajustement des données pour trouver les paramètres U0 et tau
    popt, pcov = curve_fit(exponential_response, t_data, v_data, p0=[max(v_data), t_data[np.argmax(v_data)]])
    U0, tau = popt
    return U0, tau


def measure_cutoffFreq(tconf):
    # Test configuration data
    dir_path = tconf["dir_path"]

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
             and f[-3:] == ".h5" and f[:5] == "dump_"]

    if len(files) == 0:
        raise ValueError('Error, no dump files found')

    lenDumps = 2 * cst.nPixPerCol * cst.nSamplesPerRow

    accuDumps = np.zeros((cst.nColPerDemux, lenDumps))
    for index, file in enumerate(np.sort(files)):
        print(file)
        colDumps, errors = rddt.read_dump_from_hdf5(os.path.join(path_data, file))
        accuDumps += colDumps
    accuDumps = (accuDumps / len(files)) * (cst.fsrADCErrorV / cst.fsrADCErrorADU)

    t = np.arange(lenDumps) / cst.fSamp

    # Doing plots
    for col in range(cst.nColPerDemux):
        plotFileName = os.path.join(path_plot, 'cutoffFreq_col{0:}.png'.format(col))

        v = accuDumps[col, :]

        tsize = 20
        it1 = riseDetect(v, 3)
        it2 = it1 + tsize

        fig = plt.figure(figsize=(8, 6))
        suptitle = "Error cutoff frequency (column {0:})".format(col)
        title = tconf.testPlanPath + '        ' + os.path.basename(dir_path)
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
        vfit = v[it1:it2] - v[it1]
        tfit = t[it1:it2] - t[it1]
        U0, tau = fit_time_constant(tfit, vfit)
        fc = 1 / (2 * np.pi * tau * 1e6)
        print("Time constant: {0:6.3f} ns".format(tau * 1e9))
        print("Cutoff frequency: {0:6.3f} MHz".format(fc))

        lbl1 = 'ADC Data'
        ax2.plot(t[it1 - 10:it2] * 1e9, v[it1 - 10:it2])
        ax2.plot(t[it1:it2] * 1e9, v[it1:it2], color='r', label=lbl1)
        lbl2 = "First order fit (fc = {0:6.2f} MHz)".format(fc)
        # Building a higher resolution time array to plot the fit
        hr_ratio = 10
        thr = np.arange(hr_ratio * tsize) / (cst.fSamp * hr_ratio)
        ax2.plot((thr + t[it1]) * 1e9, exponential_response(thr, U0, tau) + v[it1], '--', color='k', label=lbl2)
        ax2.set_xlabel(xlabel)
        ax2.set_ylabel(ylabel)
        ax2.legend(loc='best')
        ax2.grid()

        fig.tight_layout()

        plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
        print("results plotted in file " + plotFileName)


# -------------------------------------------------------------------------------------

def cutoffFreq(verbose=True):
    """
    Measure the cuttOff frequency of Dump data

    Args:
        verbose (bool): Determines if detailed output should be printed during
            execution. If True, provides additional information about the noise
            testing process and file processing steps. Defaults to True.
    """

    config = get_config()

    # Data directory
    hk_path = os.path.join(config["dir_path"], cst.hkDirName)

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(hk_path)

    if verbose:
        print("/----------------------------------------------------------")
        print("/ Band shape test:     " + config["setup"])
        print("/ Test session name:   " + config["session_name"])
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:    {0:}".format(fwVersion))
        print("/ Signal:              " + config["signal_name"])
        print("/----------------------------------------------------------\n")

    measure_cutoffFreq(config)

# -------------------------------------------------------------------------------------
