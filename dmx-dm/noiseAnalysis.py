# imports
import os

import numpy as np

import constants as cst
import general_tools as gt
import noiseAnalysisTools as nat
import readData as rddt


# TODO
## Change this function to retrieve the information from a configuration file (xml format)

def get_config():
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

    # Data directory
    dir_path = os.path.join("..", "..")

    # Looking for the session name and test configuration : "ERRO_ONLY" or "ERRO_FDBK" or "ERRO_OFCO"
    session_name = os.path.basename(os.path.realpath(dir_path))

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

    config = {"session_name": session_name,
              "dir_path": dir_path,
              "setup": session_name[22:31],
              "rate": '',
              "bxl": bxl,
              "fdbk": fdbk,
              "ofco": ofco
              }

    return config


# Plotting noise spectral density for one column
def plot_col_spectrum(acq_mode, config, verbose=False, lpf=False):
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

    # Data directory
    data_path = os.path.join(config["dir_path"], cst.dataDirName)

    # Processing science files

    if acq_mode == 'DUMP':
        config["rate"] = "FREF"
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum = nat.power_spectrum_from_dumps(data_path, npts, config["rate"])
        plot_acqmode_col_spectrum(xf, power_spectrum, acq_mode, config, [True, True, True, True], lpf, verbose)

    elif acq_mode == 'ERRO':
        config["rate"] = "FROW"
        npts = 2 ** 23
        xf, power_spectrum, data_exists = nat.power_spectrum_from_error(data_path, npts, config["rate"])
        plot_acqmode_col_spectrum(xf, power_spectrum, acq_mode, config, data_exists, lpf, verbose)

        config["rate"] = "FFRAME"
        npts = 2 ** 18
        xf, power_spectrum, data_exists = nat.power_spectrum_from_error(data_path, npts, config["rate"])
        plot_acqmode_col_spectrum(xf, power_spectrum, acq_mode, config, data_exists, lpf, verbose)


def plot_acqmode_col_spectrum(xf, power_spectrum, acq_mode, config, data_exists, lpf=False, verbose=False):
    rbw = xf[1]
    xf = xf[1:]  # Removing f=0 to make log scale plots easier to manage
    spectrum = np.sqrt(power_spectrum[:, 1:] / rbw)

    spectra_path = os.path.join(config["dir_path"], cst.spectraDirname)
    gt.createdir(spectra_path)

    # Creation of a directory for the plot files
    plot_path = os.path.join(config["dir_path"], cst.plotDirName)
    gt.createdir(plot_path)

    if config["setup"] != "ERRO-ONLY":
        substract_error = True
    else:
        substract_error = False

    # Iterates columns; plots spectrum with noise model
    for col_id in range(cst.nColPerDemux):

        if data_exists[col_id]:

            nat.plot_spectrum(xf, spectrum[col_id, :], acq_mode, col_id, config, lpf, verbose)

            if substract_error and acq_mode == "ERRO":
                # Soustraction du spectre du signal d'erreur si disponible
                error_spectra_path = os.path.join(".", cst.errorSpectraDirname)
                error_spectrum_fit_file_name = os.path.join(error_spectra_path,
                                                            "ERRO-ONLY_" + config[
                                                                "rate"] + "_col{0:}_spectrum_fit.txt".format(col_id))

                if os.path.isfile(error_spectrum_fit_file_name):
                    print(">>>", error_spectrum_fit_file_name)
                    f_error_fit, noise_error_fit = gt.read_two_vectors_from_file(error_spectrum_fit_file_name)
                    spectrum_without_error = np.sqrt(spectrum[col_id, :] ** 2 - noise_error_fit ** 2)
                    config["setup"] = config["setup"][:-4] + "ONLY"
                    nat.plot_spectrum(xf, spectrum_without_error, acq_mode, col_id, config, lpf, verbose)


def noiseAnalysis(verbose=True):
    """
    Processes a list of directories, extracts data from zip files if required, applies specific
    data processing functions, and optionally cleans up extracted files. The function handles
    two main tasks: unzipping and plotting datasets for the specified directory paths, depending
    on the supplied parameters.

    Args:

    """

    config = get_config()

    # Data directory
    hk_path = os.path.join(config["dir_path"], cst.hkDirName)

    # Looking for DEMUX identifiers (board, model, firmware)
    dmxModel, boardId, fwVersion = rddt.read_fwVersion_dmxModel(hk_path)


    if verbose:
        print("/----------------------------------------------------------")
        print("/ Noise test:          " + config["setup"])
        print("/ Test session name:   " + config["session_name"])
        print("/----------------------------------------------------------")
        print("/ DEMUX model:         " + dmxModel + " {0:}".format(boardId))
        print("/ Firmware version:    {0:}".format(fwVersion))
        print("/ Box car length:      {0:} samples".format(config["bxl"] + 1))
        print("/ Feedback:            " + config["fdbk"])
        print("/ Offset compensation: " + config["ofco"])
        print("/----------------------------------------------------------\n")

    if verbose:
        print("  Processing DUMP files")
    plot_col_spectrum("DUMP", config, verbose=verbose)

    if verbose:
        print("  Processing ERROR files")
    plot_col_spectrum("ERRO", config, verbose=verbose)


#-------------------------------------------------------------------------------------
