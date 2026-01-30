# imports
import os
import zipfile

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

import constants as cst
import general_tools as gt
import noiseAnalysisTools as nat

# NSD definitions
## Wide band noise spectral density in V/sqrt(Hz) with BW = 62.5 MHz
nsd_dict = {"ERRO_ONLY": 25e-9, "ERRO_FDBK": 20e-9, "ERRO_OFCO": 15e-9}


# Plotting noise spectral density for one column
def plot_col_spectrum(dir_path, acq_mode, config, verbose=False, lpf=False):

    """
    Plots a column spectrum based on the provided input parameters.

    Processes scientific data files from a specified directory path corresponding to
    a specific column, applies windowing functions, and computes power spectrum
    data. Converts the spectrum to an appropriate scale, computes equivalent
    noise bandwidth, signal-to-noise ratio (SNR), and the noise floor. Optionally
    fits 1/f noise if the acquisition mode is 'error'. Generates and saves plots
    in the specified directory.

    Args:
        tconf : test configurations.
        acq_mode (str): The acquisition mode of the data, either 'dump' or 'error'.
        verbose (boolean): If True some text is displayed
    """

    # Test configuration data
    nsd = nsd_dict[config]

    # session name
    session_name = os.path.basename(os.path.realpath(dir_path))

    # Data directory
    data_path = os.path.join(dir_path, cst.dataDirName)

    # Creation of a directory for the plot files
    plot_path = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(plot_path)

    # Processing science files
    data_exists = [True, True, True, True]
    plot_file_name = 'noise_' + config + '_' + acq_mode
    if acq_mode == 'DUMP':
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum = nat.power_spectrum_from_dumps(data_path, npts)

    elif acq_mode == 'ERRO':
        npts = 2 ** 23
        xf, power_spectrum, data_exists = nat.power_spectrum_from_error(data_path, npts)


    for col_id in range(cst.nColPerDemux):

        if data_exists[col_id]:

            rbw = xf[1]

            spectrum = np.sqrt(power_spectrum[col_id, :] / rbw)

            col_plot_file_name = plot_file_name + '_c{0}'.format(col_id)
            plot_full_file_name = os.path.join(plot_path, col_plot_file_name)

            # Selection of a noise model
            path_models = 'noise_models'
            if config == 'ERRO_ONLY' and acq_mode == 'DUMP':
                model_filename = os.path.join(path_models, "erro-only_125.txt")
            elif config == 'ERRO_ONLY' and acq_mode == 'ERRO':
                model_filename = os.path.join(path_models, "erro-only_6p25.txt")
            # elif config == 'ERRO-FDBK':
            #    if col_id == 0 or col_id == 3:
            #        model_filename = os.path.join(path_models, "erro-fdbk-awaxe.txt")
            #    else:
            #        model_filename = os.path.join(path_models, "erro-fdbk-rhf200.txt")
            # elif signal == 'ERRO-OFCO':
            #    model_filename = os.path.join(path_models, "erro-ofco.txt")
            else:
                model_filename = ''  # file type is unknown, no model exists

            # Doing the plot
            fig, ax = plt.subplots(1, 1, figsize=(8, 6))
            suptitle = config + ' (' + acq_mode + ' acquisition mode in column {0:})'.format(col_id)

            fig.suptitle(suptitle, fontsize=12)
            ax.set_title(session_name, fontsize=10)

            nat.plot_spectrum(ax, xf, spectrum, acq_mode, model_filename, nsd)

            if lpf != 0:  # comparison with a low pass filter model
                a_dc_estimated = spectrum[4:50].mean()
                popt, pcov = curve_fit(nat.low_pass_filter_1, xf, spectrum, np.array([a_dc_estimated, lpf]))
                a_dc = popt[0]
                fc = popt[1]
                if verbose:
                    print("Estimated cutoff frequency: {0:4.0f} MHz".format(fc / 1e6))
                    print("Estimated DC level: {0:4.1f} nV/rtHz".format(a_dc * 1e9))

                # Plotting the 1st order LPF fit
                lbl = '1st order LPF fit (fc = {0:4.0f} MHz)'.format(fc / 1e6)
                ax.loglog(xf, nat.low_pass_filter_1(xf, a_dc, fc) * 1e9, '--', color="orange", label=lbl)

            ax.legend(loc='lower left')
            ax.grid(True)
            fig.tight_layout()
            plt.savefig(plot_full_file_name, dpi=300, bbox_inches='tight')
            plt.close()
            if verbose:
                print("Results plotted in file ", plot_full_file_name)


def noiseAnalysis(verbose):
    """
    Processes a list of directories, extracts data from zip files if required, applies specific
    data processing functions, and optionally cleans up extracted files. The function handles
    two main tasks: unzipping and plotting datasets for the specified directory paths, depending
    on the supplied parameters.

    Args:

    """

    #    import time
    #    start = time.time()

    # Data directory
    dir_path = os.path.join("..", "..")

    # Data directory
    data_path = os.path.join(dir_path, cst.dataDirName)

    # Looking for the session name and test configuration : "ERRO_ONLY" or "ERRO_FDBK" or "ERRO_OFCO"
    session_name = os.path.basename(os.path.realpath(dir_path))
    tst_config = session_name[22:31]

    if verbose:
        print("/----------------------------------------------------------")
        print("/ Processing noise from session:")
        print("/ ", session_name)
        print("/----------------------------------------------------------\n")

    if verbose:
        print("  Processing DUMP files")
    plot_col_spectrum(dir_path, "DUMP", tst_config, verbose=verbose)

    if verbose:
        print("  Processing ERROR files")

    # Looking for zip files if any
    zipfiles = [f for f in os.listdir(data_path) \
                if os.path.isfile(os.path.join(data_path, f)) \
                and f[-4:] == '.zip']

    # Unzipping zip files if any
    data_from_zip_file = False
    if len(zipfiles) > 0:
        data_from_zip_file = True
        for z in zipfiles:
            print('    Unzipping ERROR files from ', os.path.join(data_path, z))
            with zipfile.ZipFile(os.path.join(data_path, z), 'r') as zip_ref:
                error_file_names = zip_ref.namelist()  # liste des noms dans l’archive
                zip_ref.extractall(data_path)

    plot_col_spectrum(dir_path, "ERRO", tst_config, verbose=verbose)

    # Removing fits files if extracted from a zip file
    if (data_from_zip_file):
        print('    Removing ERROR files... ')
        for file in error_file_names:
            os.remove(os.path.join(data_path, file))


#    duration = time.time() - start
#    print(f"Duration: {duration:.1f} s")


#-------------------------------------------------------------------------------------
