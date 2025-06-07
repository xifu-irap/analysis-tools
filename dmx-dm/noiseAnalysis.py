# imports
import os
import zipfile
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

import constants as cst
import general_tools as gt
import noiseAnalysisTools as nat


# Plotting noise spectral density for one column
# def plot_col_spectrum(dir_path, col_id, win_name, acq_mode, enob=11.5, lpf=0):
def plot_col_spectrum(tconf, col_id, acq_mode):

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
        col_id (int): The column identifier for the science data to process.
        acq_mode (str): The acquisition mode of the data, either 'dump' or 'error'.
    """

    # Test configuration data
    dir_path = tconf.file_path
    win_name = tconf.win
    enob = tconf.enob
    lpf = tconf.lpf

    # Data directory
    path_data = os.path.join(dir_path, cst.dataDirName)

    signal = dir_path[-9:]

    # Creation of a directory for the plot files
    path_plot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(path_plot)

    # Processing science files
    if acq_mode == 'dump':
        plot_file_name = 'noise_' + signal + '_dumps_c{0:}'.format(col_id)
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum = nat.power_spectrum_from_dumps(path_data, col_id, npts, win_name)

    elif acq_mode == 'error':
        plot_file_name = 'noise_' + signal + '_error_c{0:}'.format(col_id)
        npts = 2 ** 23
        xf, power_spectrum = nat.power_spectrum_from_error(path_data, col_id, npts, win_name)

    plot_full_file_name = os.path.join(path_plot, plot_file_name)

    # converting the spectrum to V/sqrt(Hz)
    rbw = xf[1]
    spectrum = np.sqrt(power_spectrum / rbw)

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
        model_filename = ''  # file type is unknown, no model exists

    # Doing the plot
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    suptitle = signal + ' (' + acq_mode + ' acquisition mode in column {0:})'.format(col_id)
    title = tconf.testPlanPath + '        ' + os.path.basename(dir_path)

    fig.suptitle(suptitle, fontsize=12)
    ax.set_title(title, fontsize=10)

    nat.plot_spectrum(ax, xf, spectrum, acq_mode, model_filename, enob)

    if lpf != 0:  # comparison with a low pass filter model
        a_dc_estimated = spectrum[2:10].mean()
        popt, pcov = curve_fit(nat.low_pass_filter_1, xf, spectrum, np.array([a_dc_estimated, lpf]))
        a_dc = popt[0]
        fc = popt[1]
        print("Estimated cutoff frequency: {0:4.0f} MHz".format(fc / 1e6))
        print("Estimated DC level: {0:4.1f} nV/rtHz".format(a_dc * 1e9))

        # Plotting the 1st order LPF fit
        lbl = '1rst order LPF fit (fc = {0:4.0f} MHz)'.format(fc / 1e6)
        ax.loglog(xf, nat.low_pass_filter_1(xf, a_dc, fc) * 1e9, '--', color="orange", label=lbl)

    ax.legend(loc='upper right')
    ax.grid(True)
    fig.tight_layout()
    plt.savefig(plot_full_file_name, dpi=300, bbox_inches='tight')
    print("Results plotted in file ", plot_full_file_name)


def process_list_of_dir(TestConfigList):
    """
    Processes a list of directories, extracts data from zip files if required, applies specific
    data processing functions, and optionally cleans up extracted files. The function handles
    two main tasks: unzipping and plotting datasets for the specified directory paths, depending
    on the supplied parameters.

    Args:
        TestConfigList : List of test configurations to process.

    """
    for tc in TestConfigList:
        print("Processing data from ", tc.file_path)

        data_from_zip_file = False
        if tc.error:
            # Looking for zip files if any
            dataPath = os.path.join(tc.file_path, cst.dataDirName)
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
            if tc.dump:
                plot_col_spectrum(tc, col, "dump")
            if tc.error:
                file_name = os.path.join(dataPath, "error_noise_C{0:}.fits".format(col))
                if os.path.isfile(file_name):
                    plot_col_spectrum(tc, col, "error")
                    # Removing fits files if extracted from a zip file
                    if data_from_zip_file:
                        os.remove(file_name)


#-------------------------------------------------------------------------------------

BASE_DATA_PATH = "/Users/laurent/Data"
DEFAULT_ENOB = 11.5
DEFAULT_LPF = 0.0
TP21_PATH = "TestPlan21-perfo"
TP27_PATH = "TestPlan27_DM-DMX2_Func_and_Perfs"


# TODO: Move TestConfig class in a common file of the project (general_tools.py)
@dataclass
class TestConfig:
    testPlanPath: str
    session_name: str
    enob: float = DEFAULT_ENOB
    lpf: float = DEFAULT_LPF
    dump: bool = True
    error: bool = True
    win: str = "none"

    @property
    def file_path(self) -> str:
        return f"{BASE_DATA_PATH}/{self.testPlanPath}/{self.session_name}"

list_of_dir = [
    TestConfig(TP21_PATH, "20250113_170000_noise_erro-only"),
    TestConfig(TP21_PATH, "20250121_110000_noise_erro-fdbk"),
    TestConfig(TP21_PATH, "20250124_170723_noise_erro-only"),
    TestConfig(TP21_PATH, "20250124_171221_noise_erro-fdbk"),
    TestConfig(TP21_PATH, "20250127_153842_noise_erro-fdbk"),
    TestConfig(TP21_PATH, "20250127_175444_noise_erro-only"),
    TestConfig(TP21_PATH, "20250127_175802_noise_erro-ofco"),
    TestConfig(TP21_PATH, "20250130_110005_noise_conf-FPAs"),
    TestConfig(TP21_PATH, "20250207_151952_tst_squid-sim"),
    TestConfig(TP21_PATH, "20250207_153737_tst_squid-sim"),
    TestConfig(TP21_PATH, "20250207_160403_noise_erro-100o"),
    TestConfig(TP21_PATH, "20250207_160741_noise_erro-100o"),
    TestConfig(TP21_PATH, "20250210_173047_noise_erro-fdbk"),
    TestConfig(TP21_PATH, "20250210_173229_noise_erro-fdbk"),
    TestConfig(TP21_PATH, "20250210_175013_noise_erro-fdbk"),
    TestConfig(TP21_PATH, "20250210_175158_noise_erro-fdbk"),
    TestConfig(TP27_PATH, "20250514_173402_noise_erro-only"),
    TestConfig(TP27_PATH, "20250605_113009_noise_erro-only", DEFAULT_ENOB, 20e6, True, False, "none")
]

#-------------------------------------------------------------------------------------
process_list_of_dir(list_of_dir[-1:])

#-------------------------------------------------------------------------------------
