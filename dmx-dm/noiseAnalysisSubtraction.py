# imports
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import zipfile
import constants as cst
import general_tools as gt
import noiseAnalysisTools as nat
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter

plt.rcParams['agg.path.chunksize']=200

# Plotting noise spectral density for one column
def subtract_noises(dir_path1, col_id1, dir_path2, col_id2, win_name, acq_mode, enob=11.5):

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
    xlabel_lin = r'Frequencies (MHz)'
    xlabel_log = r'Frequencies (Hz)'
    ylabel = r'Signal (nV / $\sqrt{Hz}$)'
    xlims_erro = [1, cst.fRow/2]
    xlims_dump = [1e5, cst.fSamp/2]
    ylims_erro = [1, 1e4]
    ylims_dump = [1, 100]

    # Data directories
    path_data1 = os.path.join(dir_path1, cst.dataDirName)
    path_data2 = os.path.join(dir_path2, cst.dataDirName)
    print(path_data1, path_data2)

    # Creation of a directory for the plot files
    path_plot = os.path.join(dir_path1, cst.plotDirName)
    gt.createdir(path_plot)

    # Processing science files
    if acq_mode == 'dump':
        npts = 2 * cst.nSamplesPerRow * cst.muxFactor
        xf, power_spectrum1 = nat.power_spectrum_from_dumps(path_data1, col_id1, npts, win_name)
        xf, power_spectrum2 = nat.power_spectrum_from_dumps(path_data2, col_id2, npts, win_name)
        power_spectrum = power_spectrum1 - power_spectrum2
        plot_file_name = 'Dump_' + signal1 + '_c{0:}-'.format(col_id1) + signal2 + '_c{0:}'.format(col_id2)
        fs = cst.fSamp
        xlims = xlims_dump
        ylims = ylims_dump

    elif acq_mode == 'error':
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

            power_spectrum = power_spectrum1 - power_spectrum2
        plot_file_name = 'Error_' + signal1 + '_c{0:}-'.format(col_id1) + signal2 + '_c{0:}'.format(col_id2)
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

    # Fit of 1/f behavior
    if acq_mode == 'error':
        xf_start = 3 # to avoid DC perturbations
        f_stop = 1e3 # max frequency of 1/f area
        xf_stop = np.where(xf < f_stop)[0][-1]
        a = nat.fit_one_over_f(xf[xf_start:xf_stop], spectrum[xf_start:xf_stop])

    # Doing plot
    fig, ax = plt.subplots(2, 1, figsize=(8, 10))

    title = signal1 + ' minus ' + signal2 + ' (' + acq_mode + ' acquisition mode in column {0:})'.format(col_id2)
    fig.suptitle(title, fontsize=12)

    lbl1 = "Spectrum"
    ax[0].semilogy(xf/1e6, spectrum*1e9, label=lbl1)
    lbl2 = "Expected noise floor for ENOB={0:}\nand bandwidth = {1:} MHz".format(enob, fs/2/1e6)
    ax[0].semilogy(xlims1, [noise_floor*1e9, noise_floor*1e9], ':', color='purple', label=lbl2)
    #lbl3 = r'{0:3.1f}'.format(noise_floor_hf*1e9)+r' nV/$\sqrt{Hz}$'
    #ax[0].semilogy([xf[0]/1e6, xf[-1]/1e6], [noise_floor_hf, noise_floor_hf], '--', color='red', label=lbl3)


    ax[0].set_xlim([xlims[0]/1e6, xlims[1]/1e6])
    ax[0].set_ylim(ylims)
    test_name1 = path_data1.split('/')[-2]
    test_name2 = path_data2.split('/')[-2]
    ax[0].set_title(test_name1 + '  ' + test_name2, fontsize=10)
    ax[0].set_ylabel(ylabel)
    ax[0].set_xlabel(xlabel_lin)

    if acq_mode == 'dump':
        # Désactiver la notation scientifique pour l'axe Y
        ax[0].yaxis.set_major_formatter(ScalarFormatter())  # Utiliser ScalarFormatter
        ax[0].yaxis.set_minor_formatter(ScalarFormatter())  # Idem pour les ticks mineurs
        ax[0].ticklabel_format(style='plain', axis='y')  # Forcer le style "plain" (non scientifique)

        # Configurer les ticks pour afficher uniquement des entiers
        ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax[0].grid(True, which='both', linestyle='--')
    ax[0].legend(loc='best', framealpha=1)


    ax[1].loglog(xf[1:], spectrum[1:]*1e9, label=lbl1)
    ax[1].loglog(xlims, [noise_floor*1e9, noise_floor*1e9], ':', color='purple', label=lbl2)
    #ax2.loglog([xf[0], xf[-1]], [noise_floor_hf*1e9, noise_floor_hf*1e9], '--', color='red', linewidth=0.8, label=lbl3)

    if acq_mode == 'error':
        x = np.arange(1e4)
        lbl6 = r'1/f noise ({0:3.1f}'.format(nat.one_over_f(1, a) * 1e6) + r' µV/$\sqrt{Hz}$ at 1 Hz)'
        ax[1].loglog(x[1:], nat.one_over_f(x[1:], a) * 1e9, '-.', color='k', label=lbl6)

    ax[1].set_xlim(xlims)
    ax[1].set_ylim(ylims)
    ax[1].set_ylabel(ylabel)
    ax[1].set_xlabel(xlabel_log)

    if acq_mode == 'dump':
        # Désactiver la notation scientifique pour l'axe Y
        ax[1].yaxis.set_major_formatter(ScalarFormatter())  # Utiliser ScalarFormatter
        ax[1].yaxis.set_minor_formatter(ScalarFormatter())  # Idem pour les ticks mineurs
        ax[1].ticklabel_format(style='plain', axis='y')  # Forcer le style "plain" (non scientifique)

        # Configurer les ticks pour afficher uniquement des entiers
        ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    ax[1].grid(True, which='both', linestyle='--')
    ax[1].legend(loc='best', framealpha=1)

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
    #subtract_noises(list_of_dir[12], col, list_of_dir[2], col, win, 'dump', enob=12.9)
    #subtract_noises(list_of_dir[12], col, list_of_dir[2], col, win, 'error', enob=12.9)
    subtract_noises(list_of_dir[3], col, list_of_dir[2], col, win, 'dump', enob=12.9)
    subtract_noises(list_of_dir[3], col, list_of_dir[2], col, win, 'error', enob=12.9)
    subtract_noises(list_of_dir[6], col, list_of_dir[2], col, win, 'dump', enob=11.5)

#-------------------------------------------------------------------------------------
