# imports
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

import constants as cst
import general_tools as gt
import noiseAnalysisTools as nat
import readData as rddt

# NSD definitions
## Wide band noise spectral density in V/sqrt(Hz) with BW = 62.5 MHz
nsd_dict = {"ERRO-ONLY": 25e-9, "FDBK-ERRO": 20e-9, "OFCO-ERRO": 15e-9}


# TODO
## Change this function to retrieve the information from a configuration file (xml format)

def get_config(session_name):
    """
    Extracts configuration details from a given session name string.

    This function parses the input session name to extract specific values associated with
    boxcar length (`BXL`), feedback (`FDBK`), and coarse offset calibration (`OFCO`).
    The extracted data is stored in a dictionary along with a setup identifier.

    Args:
        session_name (str): The session name string containing encoded configuration details.

    Returns:
        dict: A dictionary containing the parsed configuration with the following keys:
            - "setup" (str): A substring representing the setup identifier.
            - "bxl" (int): The boxcar length value extracted from the session name.
            - "fdbk" (str): The feedback value extracted from the session name.
            - "ofco" (str): The coarse offset calibration value extracted from the session name.
    """
    # Looking for boxcar length value from session name
    index_bxl = session_name.find("BXL") + len("BXL")
    bxl = int(session_name[index_bxl])

    # Looking for FB0 value
    txt = "_FDBK"
    index_fdbk = session_name.find(txt) + len(txt)
    fdbk = session_name[index_fdbk:].split('_')[0]

    # Looking for OFCO COARSE value
    txt = "_OFCO"
    index_ofco = session_name.find(txt) + len(txt)
    ofco = session_name[index_ofco:]

    config = {"setup": session_name[22:31], "bxl": bxl, "fdbk": fdbk, "ofco": ofco}

    return config


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
    nsd = nsd_dict[config["setup"]]

    # session name
    session_name = os.path.basename(os.path.realpath(dir_path))

    # Data directory
    data_path = os.path.join(dir_path, cst.dataDirName)

    # Creation of a directory for the plot files
    plot_path = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(plot_path)

    # Processing science files
    data_exists = [True, True, True, True]

    plot_file_name = 'noise_' + config["setup"] + '_' + acq_mode + '_BXL' + \
                     "{0:}".format(config["bxl"]) + '_FDBK' + config["fdbk"] + '_OFCO' + config["ofco"]

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
            if config["setup"] == 'ERRO-ONLY' and acq_mode == 'DUMP':
                model_filename = os.path.join(path_models, "erro-only_125.txt")
            elif config["setup"] == 'ERRO-ONLY' and acq_mode == 'ERRO':
                model_filename = os.path.join(path_models, "erro-only_6p25.txt")
            # elif config["setup"] == 'FDBK-ERRO':
            #    if col_id == 0 or col_id == 3:
            #        model_filename = os.path.join(path_models, "erro-fdbk-awaxe.txt")
            #    else:
            #        model_filename = os.path.join(path_models, "erro-fdbk-rhf200.txt")
            # elif signal == 'OFCO-ERRO':
            #    model_filename = os.path.join(path_models, "erro-ofco.txt")
            else:
                model_filename = ''  # file type is unknown, no model exists

            # Doing the plot
            fig, ax = plt.subplots(1, 1, figsize=(8, 6))
            suptitle = config["setup"] + ' (' + acq_mode + ' acquisition mode in column {0:})\n'.format(col_id) \
                       + ' BXL = {0:}'.format(config["bxl"]) + ', FDBK = ' + config["fdbk"] + ', OFCO = ' + config[
                           "ofco"]

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


def noiseAnalysis(verbose=True):
    """
    Processes a list of directories, extracts data from zip files if required, applies specific
    data processing functions, and optionally cleans up extracted files. The function handles
    two main tasks: unzipping and plotting datasets for the specified directory paths, depending
    on the supplied parameters.

    Args:

    """

    # Data directory
    dir_path = os.path.join("..", "..")

    # Data directory
    hk_path = os.path.join(dir_path, cst.hkDirName)

    # Looking for the session name and test configuration : "ERRO_ONLY" or "ERRO_FDBK" or "ERRO_OFCO"
    session_name = os.path.basename(os.path.realpath(dir_path))

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(hk_path)

    config = get_config(session_name)

    if verbose:
        print("/----------------------------------------------------------")
        print("/ Noise test:          " + config["setup"])
        print("/ Test session name:   " + session_name)
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:    {0:}".format(fwVersion))
        print("/ Box car length:      {0:} samples".format(config["bxl"] + 1))
        print("/ Feedback:            " + config["fdbk"])
        print("/ Offset compensation: " + config["ofco"])
        print("/----------------------------------------------------------\n")

    if verbose:
        print("  Processing DUMP files")
    plot_col_spectrum(dir_path, "DUMP", config, verbose=verbose)

    if verbose:
        print("  Processing ERROR files")
    plot_col_spectrum(dir_path, "ERRO", config, verbose=verbose)


#-------------------------------------------------------------------------------------
