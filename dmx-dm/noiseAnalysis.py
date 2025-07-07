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
def plot_col_spectrum(tconf, acq_mode, verbose=True):

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
    dir_path = tconf.file_path
    win_name = tconf.win
    enob = tconf.enob
    lpf = tconf.lpf

    # Data directory
    path_data = os.path.join(dir_path, cst.dataDirName)

    signal = tconf.signal

    # Creation of a directory for the plot files
    path_plot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(path_plot)

    # Processing science files
    if acq_mode == 'dump':
        plot_file_name = 'noise_' + signal + '_dumps'
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum = nat.power_spectrum_from_dumps(path_data, npts, win_name)

    elif acq_mode == 'error':
        plot_file_name = 'noise_' + signal + '_error'
        npts = 2 ** 23
        xf, power_spectrum = nat.power_spectrum_from_error(path_data, npts, win_name)

    # converting the spectrum to V/sqrt(Hz)
    rbw = xf[1]

    for col_id in range(cst.nColPerDemux):

        spectrum = np.sqrt(power_spectrum[col_id, :] / rbw)

        col_plot_file_name = plot_file_name + '_c{0}'.format(col_id)
        plot_full_file_name = os.path.join(path_plot, col_plot_file_name)

        # Selection of a noise model
        signal = tconf.signal
        path_models = 'noise_models'
        if signal == 'erro-only' or signal == 'erro-100o':
            model_filename = os.path.join(path_models, "erro-only.txt")
        elif signal == 'erro-fdbk':
            if col_id == 0 or col_id == 3:
                model_filename = os.path.join(path_models, "erro-fdbk-awaxe.txt")
            else:
                model_filename = os.path.join(path_models, "erro-fdbk-rhf200.txt")
        elif signal == 'erro-ofco':
            model_filename = os.path.join(path_models, "erro-fco.txt")
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
            if verbose:
                print("Estimated cutoff frequency: {0:4.0f} MHz".format(fc / 1e6))
                print("Estimated DC level: {0:4.1f} nV/rtHz".format(a_dc * 1e9))

            # Plotting the 1st order LPF fit
            lbl = '1rst order LPF fit (fc = {0:4.0f} MHz)'.format(fc / 1e6)
            ax.loglog(xf, nat.low_pass_filter_1(xf, a_dc, fc) * 1e9, '--', color="orange", label=lbl)

        ax.legend(loc='upper right')
        ax.grid(True)
        fig.tight_layout()
        plt.savefig(plot_full_file_name, dpi=300, bbox_inches='tight')
        if verbose:
            print("Results plotted in file ", plot_full_file_name)


def process_list_of_configs(TestConfigList):
    """
    Processes a list of directories, extracts data from zip files if required, applies specific
    data processing functions, and optionally cleans up extracted files. The function handles
    two main tasks: unzipping and plotting datasets for the specified directory paths, depending
    on the supplied parameters.

    Args:
        TestConfigList : List of test configurations to process.

    """

    #    import time
    #    start = time.time()

    for tc in TestConfigList:
        print("Processing data from ", tc.file_path)

        # Processing data from dump files
        if tc.dump:
            plot_col_spectrum(tc, "dump")

        # Processing data from Error files
        if tc.error:

            # Looking for zip files if any
            data_from_zip_file = False
            dataPath = os.path.join(tc.file_path, cst.dataDirName)
            zipfiles = [f for f in os.listdir(dataPath) \
                        if os.path.isfile(os.path.join(dataPath, f)) \
                        and f[-4:] == '.zip']

            # Unzipping zip files if any
            if len(zipfiles) > 0:
                data_from_zip_file = True
                for z in zipfiles:
                    print('Unzipping Error data from ', z)
                    with zipfile.ZipFile(os.path.join(dataPath, z), 'r') as zip_ref:
                        error_file_names = zip_ref.namelist()  # liste des noms dans l’archive
                        zip_ref.extractall(dataPath)

            # Processing the 4 error files
            plot_col_spectrum(tc, "error")

            # Removing fits files if extracted from a zip file
            if data_from_zip_file:
                print("Removing error files")
                for file in error_file_names:
                    os.remove(os.path.join(dataPath, file))


#    duration = time.time() - start
#    print(f"Duration: {duration:.1f} s")


#-------------------------------------------------------------------------------------

DEFAULT_ENOB = 11.5
DEFAULT_LPF = 0.0


@dataclass
class TestConfig:
    testPlanPath: str
    session_name: str
    signal: str = ""
    enob: float = DEFAULT_ENOB
    lpf: float = DEFAULT_LPF
    dump: bool = True
    error: bool = True
    win: str = "none"

    @property
    def file_path(self) -> str:
        return f"{cst.BASE_DATA_PATH}/{self.testPlanPath}/{self.session_name}"


TP27_TURBO45_PATH = "TestPlan27_DM-DMX2_Func_and_Perfs/FW-turbo-45"
TP27_TURBO38_PATH = "TestPlan27_DM-DMX2_Func_and_Perfs/FW-turbo-38"
TP27_TURBO39_PATH = "TestPlan27_DM-DMX2_Func_and_Perfs/FW-turbo-39"
TP27_STR36_PATH = "TestPlan27_DM-DMX2_Func_and_Perfs/FW-str-36"
TP27_STRFAST36_PATH = "TestPlan27_DM-DMX2_Func_and_Perfs/FW-str-fast-36"
TP27_FAST36_PATH = "TestPlan27_DM-DMX2_Func_and_Perfs/FW-fast-36"
TP27_TURBO36_PATH = "TestPlan27_DM-DMX2_Func_and_Perfs/FW-turbo-36"
TP27_STRFAST30_PATH = "TestPlan27_DM-DMX2_Func_and_Perfs/FW-str-fast-30"

list_of_configs = [
    TestConfig(cst.TP21_PATH, "20250113_170000_noise_erro-only", "erro-only"),
    TestConfig(cst.TP21_PATH, "20250121_110000_noise_erro-fdbk", "erro-fdbk"),
    TestConfig(cst.TP21_PATH, "20250124_170723_noise_erro-only", "erro-only"),
    TestConfig(cst.TP21_PATH, "20250124_171221_noise_erro-fdbk", "erro-fdbk"),
    TestConfig(cst.TP21_PATH, "20250127_153842_noise_erro-fdbk", "erro-fdbk"),
    TestConfig(cst.TP21_PATH, "20250127_175444_noise_erro-only", "erro-only"),
    TestConfig(cst.TP21_PATH, "20250127_175802_noise_erro-ofco", "erro-ofco"),
    TestConfig(cst.TP21_PATH, "20250130_110005_noise_conf-FPAs"),
    TestConfig(cst.TP21_PATH, "20250207_151952_tst_squid-sim"),
    TestConfig(cst.TP21_PATH, "20250207_153737_tst_squid-sim"),
    TestConfig(cst.TP21_PATH, "20250207_160403_noise_erro-100o", "erro-only"),
    TestConfig(cst.TP21_PATH, "20250207_160741_noise_erro-100o", "erro-only"),
    TestConfig(cst.TP21_PATH, "20250210_173047_noise_erro-fdbk", "erro-fdbk"),
    TestConfig(cst.TP21_PATH, "20250210_173229_noise_erro-fdbk", "erro-fdbk"),
    TestConfig(cst.TP21_PATH, "20250210_175013_noise_erro-fdbk", "erro-fdbk"),
    TestConfig(cst.TP21_PATH, "20250210_175158_noise_erro-fdbk", "erro-fdbk"),
    TestConfig(cst.TP27_PATH, "20250514_173402_noise_erro-only", "erro-only"),
    TestConfig(cst.TP27_PATH, "20250605_113009_noise_erro-only", "erro-only", DEFAULT_ENOB, 20e6, True, True, "none"),
    TestConfig(cst.TP27_PATH, "20250610_150520_noise_erro-only", "erro-only"),
    TestConfig(TP27_TURBO45_PATH, "20250618_172726_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_TURBO45_PATH, "20250618_171649_noise_erro-fdbk_Fb-0_5V", "erro-fdbk"),
    TestConfig(TP27_TURBO45_PATH, "20250618_171815_noise_erro-fdbk_Fb0V", "erro-fdbk"),
    TestConfig(TP27_TURBO45_PATH, "20250618_171941_noise_erro-fdbk_Fb0_5V", "erro-fdbk"),
    TestConfig(TP27_TURBO45_PATH, "20250618_175720_noise_erro-ofco_Fb0V", "erro-ofco"),
    TestConfig(TP27_STR36_PATH, "20250630_161836_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_STRFAST36_PATH, "20250630_165513_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_TURBO36_PATH, "20250630_172635_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_TURBO36_PATH, "20250701_111622_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_TURBO36_PATH, "20250701_112036_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_TURBO36_PATH, "20250701_112437_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_STR36_PATH, "20250701_113636_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_STR36_PATH, "20250701_114030_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_STRFAST36_PATH, "20250701_114557_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_STRFAST36_PATH, "20250701_115007_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_FAST36_PATH, "20250701_115527_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_FAST36_PATH, "20250701_115858_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_TURBO38_PATH, "20250701_141934_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_TURBO39_PATH, "20250701_145529_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_TURBO36_PATH, "20250701_151718_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_STRFAST30_PATH, "20250701_163206_noise_erro-fdbk_Fb0V", "erro-only"),
    TestConfig(TP27_STRFAST30_PATH, "20250701_163611_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_STRFAST36_PATH, "20250701_170300_noise_erro-only_Fb0V", "erro-only"),
    TestConfig(TP27_FAST36_PATH, "20250701_171143_noise_erro-only_Fb0V", "erro-only")
]

#-------------------------------------------------------------------------------------
process_list_of_configs(list_of_configs[-1:])

#-------------------------------------------------------------------------------------
