# imports
import os

import numpy as np

import constants as cst
import general_tools as gt
import noiseAnalysisTools as nat
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
def plot_col_spectrum(acq_mode, config, process_frow=False, peak_detect=True, plot_model=False, lpf=0, verbose=False):
    """
    Plots the column spectrum based on the acquisition mode provided. This function processes science
    files and generates power spectra for either 'DUMP' or 'ERRO' acquisition modes. Visualization
    parameters and configurations are set from the provided arguments.

    Args:
        acq_mode (str): Acquisition mode specifying the data processing approach. Supports 'DUMP'
            or 'ERRO' modes.
        config (dict): Configuration dictionary containing all required parameters, including the
            directory path and rate settings.
        process_frow (bool, optional): If True, the data are processed at the sample rate FROW.
            Defaults to False.
        verbose (bool, optional): If True, enables verbose output during processing and plotting.
            Defaults to False.
        lpf (bool, optional): Low-pass filter toggle for modifying the spectrum plot. Defaults to False.

    """
    # Data directory
    data_path = os.path.join(config["dir_path"], cst.dataDirName)

    # Processing science files

    if acq_mode == 'DUMP':
        config["rate"] = "FREF"
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum = nat.power_spectrum_from_dumps(data_path, npts, config["rate"])
        plot_acq_mode_col_spectrum(xf, power_spectrum, acq_mode, config, [True, True, True, True],
                                   plot_model=plot_model, lpf=lpf,
                                   peak_detect=peak_detect, verbose=verbose)

    elif acq_mode == 'ERRO':
        if process_frow:
            config["rate"] = "FROW"
            npts = 2 ** 23
            xf, power_spectrum, data_exists = nat.power_spectrum_from_error(data_path, npts, config["rate"])
            plot_acq_mode_col_spectrum(xf, power_spectrum, acq_mode, config, data_exists, plot_model=plot_model,
                                       lpf=lpf, peak_detect=peak_detect,
                                       verbose=verbose)

        config["rate"] = "FFRAME"
        npts = 2 ** 18
        xf, power_spectrum, data_exists = nat.power_spectrum_from_error(data_path, npts, config["rate"])
        plot_acq_mode_col_spectrum(xf, power_spectrum, acq_mode, config, data_exists, plot_model=plot_model, lpf=lpf,
                                   peak_detect=peak_detect,
                                   verbose=verbose)


def plot_acq_mode_col_spectrum(xf, power_spectrum, acq_mode, config,
                               data_exists, plot_model=False, lpf=0,
                               peak_detect=True, verbose=False):
    """
    Plots the acquisition mode column spectrum based on the given parameters.

    This function processes the acquisition mode spectrum data and visualizes it. Parameters control
    the configuration of the plot as well as the application of data processing steps such as adding
    a low-pass filter or increasing verbosity for debugging purposes.

    Args:
        xf: Frequency data used in plotting.
        power_spectrum: Power spectral density values corresponding to the frequency data.
        acq_mode: Acquisition mode information used for labeling or processing the data.
        config: Configuration settings dictating how the spectrum should be plotted or processed.
        data_exists: Flag indicating whether the necessary data is available for plotting.
        lpf: Optional; Boolean flag to specify if a low-pass filter should be applied to the data.
        verbose: Optional; Boolean flag to enable or disable detailed output for debugging purposes.
    """
    rbw = xf[1]
    xf = xf[1:]  # Removing f=0 to make log scale plots easier to manage
    spectrum = np.sqrt(power_spectrum[:, 1:] / rbw)

    spectra_path = os.path.join(config["dir_path"], cst.spectraDirname)
    gt.createdir(spectra_path)

    # Creation of a directory for the plot files
    plot_path = os.path.join(config["dir_path"], cst.plotDirName)
    gt.createdir(plot_path)

    if config["setup"] == "ERRO-ONLY":
        substract_error = False
    else:
        substract_error = True

    # Iterates columns; plots spectrum with noise model
    for col_id in range(cst.nColPerDemux):

        if data_exists[col_id]:

            nat.plot_spectrum(xf, spectrum[col_id, :], acq_mode, col_id, config, plot_model=plot_model, lpf=lpf,
                              peak_detect=peak_detect, verbose=verbose)

            if substract_error and acq_mode == "ERRO":

                # Soustraction du spectre du signal d'erreur si disponible
                error_spectra_path = os.path.join(".", cst.errorSpectraDirname)
                error_param_fit_file_name = os.path.join(error_spectra_path, "ERRO-ONLY_" + config[
                                                             "rate"] + "_col{0:}_param_fit.txt".format(col_id))
                if os.path.isfile(error_param_fit_file_name):
                    nat.plot_fit_spectrum(xf, col_id, config, verbose)


def noiseAnalysis(process_frow=False, plot_model=False, lpf=0, verbose=True):
    """
    Analyzes noise characteristics based on DEMUX identifiers and configurations,
    processes specific file types, and optionally provides verbose output for
    detailed insights about the analysis steps, file types, and configurations.

    Args:
        process_frow (bool, optional): If True, the ERROR data are processed at the sample rates FFRAME and FROW.
            If False the ERROR data are processed at FFRAME only. Defaults to False.
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
    plot_col_spectrum("DUMP", config, process_frow=False, peak_detect=False, plot_model=plot_model, lpf=lpf,
                      verbose=verbose)

    if verbose:
        print("  Processing ERROR files")
    plot_col_spectrum("ERRO", config, process_frow=process_frow, peak_detect=False, plot_model=plot_model, lpf=lpf,
                      verbose=verbose)


#-------------------------------------------------------------------------------------
