# imports
import os
import zipfile

import matplotlib.pyplot as plt
import numpy as np

import constants as cst
import general_tools as gt
import noiseAnalysisTools as nat


# Plotting cross-noise spectral density for one column
def plot_col_cross_spectrum(dir_path, col_ref, win_name):
    """
    Plots the cross-spectrum between two columns based on the provided input parameters.

    Processes scientific data files from a specified directory path corresponding to
    two specific columns, applies windowing functions,

    Args:
        dir_path(str): The directory path containing data files for processing.
        col_ref(int): The column identifier reference column.
        win_name(str): The window function to apply (e.g., "blackman").
    """

    # Data directory
    path_data = os.path.join(dir_path, cst.dataDirName)

    signal = dir_path[-9:]

    # Creation of a directory for the plot files
    path_plot = os.path.join(dir_path, cst.plotDirName)
    gt.createdir(path_plot)

    plot_file_name = 'Xnoise_'+signal+'_error_c{0:}.png'.format(col_ref)
    npts = 2**23
    xf, cross_power_spectrum = nat.cross_power_spectrum_from_error(path_data, col_ref, npts, win_name)

    plot_full_file_name = os.path.join(path_plot, plot_file_name)

    # converting the spectrum to V/sqrt(Hz)
    rbw = xf[1]
    cross_power_spectrum = np.sqrt(cross_power_spectrum/rbw)

    # Doing the plot
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)

    sup_title = 'Spectrum of cross correlation with col {0:} ('.format(col_ref) + signal + ')'
    fig.suptitle(sup_title, fontsize=12)

    title = os.path.basename(dir_path)
    ax.set_title(title, fontsize=10)

    nat.plot_cross_spectrum(ax, xf, cross_power_spectrum, col_ref)

    fig.tight_layout()
    plt.savefig(plot_full_file_name, dpi=300, bbox_inches='tight')
    print("    -> Results plotted in file ", plot_full_file_name)


def process_list_of_dir(dir_list, col_ref, win_name):
    """
    Traite une liste de répertoires et leurs données associées.
    
    Args:
        dir_list (list[str]): Liste des chemins de répertoires à traiter
        col_ref (int): Colonne de référence pour le traitement
        win_name (str): Nom de la fenêtre pour le traitement
    """
    ZIP_EXTENSION = '.zip'
    ERROR_FILE_PATTERN = "error_noise_C{0:}.fits"

    def extract_zip_files(data_path):
        """Extrait tous les fichiers zip dans le répertoire donné."""
        zip_files = [f for f in os.listdir(data_path) 
                    if os.path.isfile(os.path.join(data_path, f)) 
                    and f.endswith(ZIP_EXTENSION)]
        
        if not zip_files:
            return False
            
        for zip_file in zip_files:
            print(f'    Extraction des données depuis {zip_file}...', end=' ')
            with zipfile.ZipFile(os.path.join(data_path, zip_file), 'r') as zip_ref:
                zip_ref.extractall(data_path)
            print("Terminé")
        return True

    for directory in dir_list:
        print(f"Traitement des données depuis {directory}")
        data_path = os.path.join(directory, cst.dataDirName)
        
        has_extracted_files = extract_zip_files(data_path)
        
        plot_col_cross_spectrum(directory, col_ref, win_name)

        if has_extracted_files:
            for col in range(cst.nColPerDemux):
                error_file = os.path.join(data_path, ERROR_FILE_PATTERN.format(col))
                os.remove(error_file)

#-------------------------------------------------------------------------------------

list_of_dir = [
    "/Users/laurent/Data/TestPlan21-perfo/20250113_170000_noise_erro-only",
    "/Users/laurent/Data/TestPlan21-perfo/20250121_110000_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250124_170723_noise_erro-only",
    "/Users/laurent/Data/TestPlan21-perfo/20250124_171221_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_153842_noise_erro-fdbk",
    "/Users/laurent/Data/TestPlan21-perfo/20250127_175444_noise_erro-only_autreC3",
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
col_ref = 0
process_list_of_dir(list_of_dir[5:6], col_ref, win)

#-------------------------------------------------------------------------------------