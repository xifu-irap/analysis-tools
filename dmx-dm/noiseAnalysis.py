# imports
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import zipfile
import constants as cst
import general_tools as gt
import noiseAnalysisTools as nat
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter


# Plotting noise spectral density for one column
def plot_col_spectrum(dir_path, col_id, win_name, acq_mode, enob=11.5):
    """
    Plots a column spectrum based on the provided input parameters.

    Processes scientific data files from a specified directory path corresponding to
    a specific column, applies windowing functions, and computes power spectrum
    data. Converts the spectrum to an appropriate scale, computes equivalent
    noise bandwidth, signal-to-noise ratio (SNR), and the noise floor. Optionally
    fits 1/f noise if the acquisition mode is 'error'. Generates and saves plots
    in the specified directory.

    Args:
        dir_path (str): The directory path containing data files for processing.
        col_id (int): The column identifier for the science data to process.
        win_name (str): The window function to apply (e.g., "blackman").
        acq_mode (str): The acquisition mode of the data, either 'dump' or 'error'.
        enob (float, optional): The Effective Number Of Bits (ENOB) for the
            computation of SNR and noise floor. Defaults to 11.5.
    """

    signal = dir_path[-9:]
    xlabel_lin = r'Frequencies (MHz)'
    xlabel_log = r'Frequencies (Hz)'
    ylabel = r'Error signal (nV / $\sqrt{Hz}$)'
    xlims_erro = [1, cst.fRow/2]
    xlims_dump = [1e5, cst.fSamp/2]
    ylims_erro = [10, 1e4]
    ylims_dump = [10, 100]

    # Selection of a noise model
    path_models = 'noise_models'
    model = True # the initial hypothesis is "a model exists"
    if signal == 'erro-only' or signal == 'erro-100o':
        model_filename = os.path.join(path_models, "erro-only.txt")
    elif signal == 'erro-fdbk':
        if col_id == 0 or col_id == 3:
            model_filename = os.path.join(path_models, "erro-fdbk-awaxe.txt")
        else:
            model_filename = os.path.join(path_models, "erro-fdbk-rhf200.txt")
    elif signal == 'erro-ofco':
        model_filename = os.path.join(path_models, "erro-ofco.txt")
    else:
        model = False # file type is unknown, no model exists

    # Data directory
    path_data = os.path.join(dir_path, cst.dataDirName)

    # Session name
    session_name = os.path.basename(dir_path)

    # Creation of a directory for the plot files
    path_plot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(path_plot)

    # Processing science files
    if acq_mode == 'dump':
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum = nat.power_spectrum_from_dumps(path_data, col_id, npts, win_name)
        plot_file_name = 'noise_'+signal+'_dumps_c{0:}'.format(col_id)
        fs = cst.fSamp
        xlims = xlims_dump
        ylims = ylims_dump

    elif acq_mode == 'error':
        npts = 2**23
        xf, power_spectrum = nat.power_spectrum_from_error(path_data, col_id, npts, win_name)
        plot_file_name = 'noise_'+signal+'_error_c{0:}'.format(col_id)
        fs = cst.fRow
        xlims = xlims_erro
        ylims = ylims_erro

    xlims1 = [0, fs / 2 / 1e6]
    plot_full_file_name = os.path.join(path_plot, plot_file_name)

    # converting the spectrum to V/sqrt(Hz)
    rbw = xf[1]
    spectrum = np.sqrt(power_spectrum/rbw)

    # Computation of the SNR equivalent to the requested ENOB
    snr_db = 6.02*enob + 1.76
    snr = 10**(snr_db/20)

    # Computation of the noise floor per sqrt(Hz) corresponding to the ENOB
    noise_floor = cst.fsrADCErrorV / (snr * np.sqrt(fs/2))

    # Measuring the noise floor at high frequencies
    noise_floor_hf = spectrum[-100:-2].mean()

    if model:
        # Getting theoretical noise data from 0 to 125 MHz
        f_theo, noise_theo = gt.read_two_vectors_from_file(model_filename)
        # Aliasing the noise at 62.5 MHz
        f_theo, noise_theo = nat.undersampling_with_aliasing(f_theo, noise_theo, cst.fSamp, 2*cst.fSamp)

    # Fit of 1/f behavior
    if acq_mode == 'error':
        xf_start = 3 # to avoid DC perturbations
        f_stop = 1e3 # max frequency of 1/f area
        xf_stop = np.where(xf < f_stop)[0][-1]
        a = nat.fit_one_over_f(xf[xf_start:xf_stop], spectrum[xf_start:xf_stop])

        # Computing theoretical noise data after aliasing at fRow
        if model:
            f_theo, noise_theo = nat.undersampling_with_aliasing(f_theo, noise_theo, cst.fRow, cst.fSamp)

    # Doing plot
    fig = plt.figure(figsize=(8, 10))
    title = session_name + "  col {0:}".format(col_id)
    ax1 = fig.add_subplot(2, 1, 1)

    lbl1 = "Spectrum"
    ax1.semilogy(xf/1e6, spectrum*1e9, label=lbl1)
    lbl2 = "Expected noise floor for ENOB={0:}\nand bandwidth = {1:} MHz".format(enob, fs/2/1e6)
    ax1.semilogy(xlims1, [noise_floor*1e9, noise_floor*1e9], ':', color='purple', label=lbl2)
    #lbl3 = r'{0:3.1f}'.format(noise_floor_hf*1e9)+r' nV/$\sqrt{Hz}$'
    #ax1.semilogy([xf[0]/1e6, xf[-1]/1e6], [noise_floor_hf, noise_floor_hf], '--', color='red', label=lbl3)
    if model:
        lbl4 = 'Model (from datasheets)'
        ax1.semilogy(f_theo/1e6, noise_theo*1e9, '--', color='r', label=lbl4)

    ax1.set_xlim([xlims[0]/1e6, xlims[1]/1e6])
    ax1.set_ylim(ylims)
    ax1.set_title(title)
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel_lin)

    if acq_mode == 'dump':
        # Désactiver la notation scientifique pour l'axe Y
        ax1.yaxis.set_major_formatter(ScalarFormatter())  # Utiliser ScalarFormatter
        ax1.yaxis.set_minor_formatter(ScalarFormatter())  # Idem pour les ticks mineurs
        ax1.ticklabel_format(style='plain', axis='y')  # Forcer le style "plain" (non scientifique)

        # Configurer les ticks pour afficher uniquement des entiers
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax1.grid(True, which='both', linestyle='--')
    ax1.legend(loc='best', framealpha=1)


    ax2 = fig.add_subplot(2, 1, 2)

    ax2.loglog(xf[1:], spectrum[1:]*1e9, label=lbl1)
    ax2.loglog(xlims, [noise_floor*1e9, noise_floor*1e9], ':', color='purple', label=lbl2)
    #ax2.loglog([xf[0], xf[-1]], [noise_floor_hf*1e9, noise_floor_hf*1e9], '--', color='red', linewidth=0.8, label=lbl3)
    if model:
        ax2.loglog(f_theo, noise_theo*1e9, '--', color='r', label=lbl4)

    if acq_mode == 'error':
        x = np.arange(1e4)
        lbl6 = r'1/f noise ({0:3.1f}'.format(nat.one_over_f(1, a) * 1e6) + r' µV/$\sqrt{Hz}$ at 1 Hz)'
        ax2.loglog(x[1:], nat.one_over_f(x[1:], a)*1e9, '-.', color='k', label=lbl6)

    ax2.set_xlim(xlims)
    ax2.set_ylim(ylims)
    ax2.set_title(title)
    ax2.set_ylabel(ylabel)
    ax2.set_xlabel(xlabel_log)

    if acq_mode == 'dump':
        # Désactiver la notation scientifique pour l'axe Y
        ax2.yaxis.set_major_formatter(ScalarFormatter())  # Utiliser ScalarFormatter
        ax2.yaxis.set_minor_formatter(ScalarFormatter())  # Idem pour les ticks mineurs
        ax2.ticklabel_format(style='plain', axis='y')  # Forcer le style "plain" (non scientifique)

        # Configurer les ticks pour afficher uniquement des entiers
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax2.grid(True, which='both', linestyle='--')
    ax2.legend(loc='best', framealpha=1)

    fig.tight_layout()

    plt.savefig(plot_full_file_name, dpi=300, bbox_inches='tight')
    print("Results plotted in file ", plot_full_file_name)


def process_list_of_dir(dir_list, win_name, process_dump, process_error):
    """
    Processes a list of directories, extracts data from zip files if required, applies specific
    data processing functions, and optionally cleans up extracted files. The function handles
    two main tasks: unzipping and plotting datasets for the specified directory paths, depending
    on the supplied parameters.

    Args:
        dir_list (list[str]): List of directory paths to process.
        win_name (str): Window name or identifier used in the data processing.
        process_dump (bool): Flag indicating whether to process dump data.
        process_error (bool): Flag indicating whether to process error data.
    """

    for p in dir_list:
        print("Processing data from ", p)

        data_from_zip_file = False
        if process_error:
            # Looking for zip files if any
            dataPath = os.path.join(p, cst.dataDirName)
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

        for col in range(cst.nColPerDemux):
            if process_dump:
                plot_col_spectrum(p, col, win_name, "dump")
            if process_error:
                file_name = os.path.join(dataPath, "error_noise_C{0:}.fits".format(col))
                if os.path.isfile(file_name):
                    plot_col_spectrum(p, col, win_name, "error")
                    # Removing fits files if extracted from a zip file
                    if data_from_zip_file:
                        os.remove(file_name)

#-------------------------------------------------------------------------------------

list_of_dir = [
    "/Users/laurent/Data/TestPlan21-perfo/20250113_170000_noise_erro-only",
    "/Users/laurent/Data/TestPlan21-perfo/20250121_110000_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250124_170723_noise_erro-only",
    "/Users/laurent/Data/TestPlan21-perfo/20250124_171221_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_153842_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_175444_noise_erro-only",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_175802_noise_erro-ofco",
    "/Users/laurent/Data/TestPlan21-perfo/20250130_110005_noise_conf-FPAs",
    "/Users/laurent/Data/TestPlan21-perfo/20250207_151952_tst_squid-sim",
    "/Users/laurent/Data/TestPlan21-perfo/20250207_153737_tst_squid-sim",
    "/Users/laurent/Data/TestPlan21-perfo/20250207_160403_noise_erro-100o",
    "/Users/laurent/Data/TestPlan21-perfo/20250207_160741_noise_erro-100o",
    "/Users/laurent/Data/TestPlan21-perfo/20250210_173047_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250210_173229_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250210_175013_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250210_175158_noise_erro-fdbk"
]

win = "none"
#win = "blackman"
dump = True
error = True

#-------------------------------------------------------------------------------------
process_list_of_dir(list_of_dir[-1:], win, dump, error)

#-------------------------------------------------------------------------------------
