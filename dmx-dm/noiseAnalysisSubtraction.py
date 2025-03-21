# imports
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import zipfile
import constants as cst
import general_tools as gt
import noiseAnalysisTools as nat

plt.rcParams['agg.path.chunksize']=200

def subtract_noises(dir_path1, col_id1, dir_path2, col_id2, win_name, acq_mode, enob=11.5):
    """
    Subtract noise spectra from two data directories and plot the results.

    This function processes noise data from two specified directories, either in 'dump'
    or 'error' acquisition mode. It handles the extraction of data from zip files if
    necessary, computes the power spectrum for both data sets, and generates a plot of
    the resulting noise spectrum after subtraction.

    Parameters:
    ----------
    dir_path1 : str
        The file path to the first data directory containing the first signal data.

    col_id1 : int
        The column index of the first signal data to be processed.

    dir_path2 : str
        The file path to the second data directory containing the second signal data.

    col_id2 : int
        The column index of the second signal data to be processed.

    win_name : str
        The name of the window function to be applied during the power spectrum calculation.

    acq_mode : str
        The acquisition mode, which can be either 'dump' or 'error'. This determines how
        the power spectrum is computed.

    enob : float, optional
        The Effective Number of Bits for the measurement. Default is 11.5.

    Returns:
    -------
    None
        This function does not return any value. It generates and saves a plot of the
        noise spectrum after subtraction to the specified directory.

    Raises:
    ------
    FileNotFoundError
        If the specified directories do not exist or contain no zip files when in 'error' mode.

    Notes:
    -----
    - The function assumes that constants and helper functions such as `cst`, `nat`,
      and `gt` are defined elsewhere in the code.
    - The resulting plot is saved in the specified plot directory with a filename
      that reflects the processed signals and acquisition mode.
    """

    # unzipping error data
    data_from_zip_file = [False, False]
    if acq_mode == 'error':
        for index, dir_path in enumerate([dir_path1, dir_path2]):
            # Looking for zip files if any
            dataPath = os.path.join(dir_path, cst.dataDirName)
            zipfiles = [f for f in os.listdir(dataPath) \
                        if os.path.isfile(os.path.join(dataPath, f)) \
                        and f[-4:] == '.zip']

            # Unzipping zip files
            if len(zipfiles) == 0:
                print('There is no zip file...')
            else:
                data_from_zip_file[index] = True
                for z in zipfiles:
                    print('Unzipping Error data from ', z)
                    with zipfile.ZipFile(os.path.join(dataPath, z), 'r') as zip_ref:
                        zip_ref.extractall(dataPath)

    signal1 = dir_path1[-9:]
    signal2 = dir_path2[-9:]

    # Data directories
    path_data1 = os.path.join(dir_path1, cst.dataDirName)
    path_data2 = os.path.join(dir_path2, cst.dataDirName)
    print(path_data1, path_data2)

    # Creation of a directory for the plot files
    path_plot = os.path.join(dir_path1, cst.plotDirName)
    gt.createdir(path_plot)

    # Processing science files
    if acq_mode == 'dump':
        plot_file_name = 'noise_' + signal1 + 'C{0}_minus_'.format(col_id1) + signal2 + 'C{0:}_dumps'.format(col_id2)
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum1 = nat.power_spectrum_from_dumps(path_data1, col_id1, npts, win_name)
        xf, power_spectrum2 = nat.power_spectrum_from_dumps(path_data2, col_id2, npts, win_name)

    elif acq_mode == 'error':
        plot_file_name = 'noise_' + signal1 + 'C{0}_minus_'.format(col_id1) + signal2 + 'C{0:}_error'.format(col_id2)
        npts = 2**23
        xf, power_spectrum1 = nat.power_spectrum_from_error(path_data1, col_id1, npts, win_name)
        xf, power_spectrum2 = nat.power_spectrum_from_error(path_data2, col_id2, npts, win_name)
        # Removing data files if extracted from a zip
        prefix = 'error_noise_C'
        suffix = '.fits'
        if data_from_zip_file[0]:
            path1 = os.path.join(path_data1, cst.dataDirName)
            prefix1 = os.path.join(path1, prefix)
            file_names1 = glob.glob(prefix1 + "*" + suffix)
            for file in file_names1:
                os.remove(file)
        if data_from_zip_file[1]:
            path2 = os.path.join(path_data2, cst.dataDirName)
            prefix2 = os.path.join(path2, prefix)
            file_names2 = glob.glob(prefix2 + "*" + suffix)
            for file in file_names2:
                os.remove(file)

    plot_full_file_name = os.path.join(path_plot, plot_file_name)

    # converting the spectrum to V/sqrt(Hz)
    power_spectrum = power_spectrum1 - power_spectrum2
    rbw = xf[1]
    spectrum = np.sqrt(power_spectrum/rbw)

    # Selection of a noise model
    model_filename = ''

    # Doing the plot
    fig, ax = plt.subplots(2, 1, figsize=(8, 10))
    suptitle = signal1 + ' minus ' + signal2 + ' (' + acq_mode + ' acquisition mode in column {0:})'.format(col_id2)
    title = os.path.basename(dir_path1) + '  ' + os.path.basename(dir_path2)

    fig.suptitle(suptitle, fontsize=12)
    ax[0].set_title(title, fontsize=10)

    nat.plot_spectrum(ax, xf, spectrum, acq_mode, model_filename, enob)

    fig.tight_layout()

    plt.savefig(plot_full_file_name, dpi=300, bbox_inches='tight')
    print("Results plotted in file ", plot_full_file_name)


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
    "/Users/laurent/Data/TestPlan21-perfo/20250210_175158_noise_erro-fdbk"
]

win = "none"
#win = "blackman"

#-------------------------------------------------------------------------------------
for col in range(cst.nColPerDemux):
    subtract_noises(list_of_dir[12], col, list_of_dir[2], col, win, 'dump', enob=12.9)
    subtract_noises(list_of_dir[12], col, list_of_dir[2], col, win, 'error', enob=12.9)
    #subtract_noises(list_of_dir[3], col, list_of_dir[2], col, win, 'dump', enob=12.9)
    #subtract_noises(list_of_dir[3], col, list_of_dir[2], col, win, 'error', enob=12.9)
    #subtract_noises(list_of_dir[6], col, list_of_dir[2], col, win, 'dump', enob=11.5)

#-------------------------------------------------------------------------------------
