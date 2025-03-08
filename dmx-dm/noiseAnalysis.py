# imports
import numpy as np
import matplotlib.pyplot as plt
import os
import zipfile
import constants as cst
import general_tools as gt
import noiseAnalysisTools as nat


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

    # Data directory
    path_data = os.path.join(dir_path, cst.dataDirName)

    signal = dir_path[-9:]

    # Creation of a directory for the plot files
    path_plot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(path_plot)

    # Processing science files
    if acq_mode == 'dump':
        plot_file_name = 'noise_'+signal+'_dumps_c{0:}'.format(col_id)
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum = nat.power_spectrum_from_dumps(path_data, col_id, npts, win_name)

    elif acq_mode == 'error':
        plot_file_name = 'noise_'+signal+'_error_c{0:}'.format(col_id)
        npts = 2**23
        xf, power_spectrum = nat.power_spectrum_from_error(path_data, col_id, npts, win_name)

    plot_full_file_name = os.path.join(path_plot, plot_file_name)

    # converting the spectrum to V/sqrt(Hz)
    rbw = xf[1]
    spectrum = np.sqrt(power_spectrum/rbw)

    # Selection of a noise model
    signal = dir_path[-9:]
    path_models = 'noise_models'
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
        model_filename = '' # file type is unknown, no model exists

    # Doing the plot
    fig, ax = plt.subplots(2, 1, figsize=(8, 10))
    suptitle = signal + ' (' + acq_mode + ' acquisition mode in column {0:})'.format(col_id)
    title = os.path.basename(dir_path)

    fig.suptitle(suptitle, fontsize=12)
    ax[0].set_title(title, fontsize=10)

    nat.plot_spectrum(ax, xf, spectrum, acq_mode, model_filename, enob)

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
process_list_of_dir(list_of_dir[3:4], win, dump, error)

#-------------------------------------------------------------------------------------
