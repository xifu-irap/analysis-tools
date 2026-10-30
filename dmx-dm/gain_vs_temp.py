# ---------------------------------------------------------------------------------
# !/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright (C) 2021-2030 Laurent Ravera, IRAP Toulouse.
#  This file is part of the ATHENA X-IFU DRE data analysis tools software.
#
#  analysis-tools is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# ---------------------------------------------------------------------------------
#
#  laurent.ravera@cnrs.fr
#  gain_vs_temp.py
#
# ---------------------------------------------------------------------------------

# imports
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator, ScalarFormatter

import constants as cst
import readData as rddt
import general_tools as gt

colors = ['#FF0000', '#0000FF', "#29E329", '#FFA500', '#800080', '#008B8B', '#8B0000',
 '#696969', '#A52A2A', '#000000', '#4682B4', '#556B2F', '#4B0082', '#B8860B',
 '#C0C0C0', '#CD5C5C', '#8B008B', '#FFD700', '#00008B']


def _darken(color, factor=0.35):
    rgb = np.asarray(mcolors.to_rgb(color))
    return tuple(rgb * (1.0 - factor))


def _lighten(color, factor=0.5):
    rgb = np.asarray(mcolors.to_rgb(color))
    return tuple(rgb + (1.0 - rgb) * factor)

def get_temperatures(session, verbose=True):

    hk_path = os.path.join(session, cst.hkDirName)

    header1 = 'TEMP_MAX(°C)'
    header2 = 'PT104_channel2Value()'
    header3 = 'TEMP_AVE(°C)'
    header4 = 'PT104_channel1Value()'

    filename_dmx_list = [f for f in os.listdir(hk_path)
                        if os.path.isfile(os.path.join(hk_path, f))
                        and f.startswith('Hks_DMXA')]
    if len(filename_dmx_list) != 1:
        print("ERROR: wrong number of DMX hk files: {0:}".format(
            len(filename_dmx_list)))
        return
    filename_dmx = os.path.join(hk_path, filename_dmx_list[0])

    if verbose:
        print("   Getting DEMUX temperatures from file ", filename_dmx)
    hk1 = np.array(rddt.read_hk_name_from_csv(filename_dmx, header1)).astype(float)
    t_max = hk1.mean()
    hk3 = np.array(rddt.read_hk_name_from_csv(filename_dmx, header3)).astype(float)
    t_avg = hk3.mean()


    filename_pt104_list = [f for f in os.listdir(hk_path)
                          if os.path.isfile(os.path.join(hk_path, f))
                          and f.startswith('Hks_Pt104_')]
    if len(filename_pt104_list) != 1:
        print("ERROR: wrong number of Pt100 hk files: {0:}".format(
            len(filename_pt104_list)))
        return
    filename_pt104 = os.path.join(hk_path, filename_pt104_list[0])
    
    if verbose:
        print("   Getting PT100 temperatures from file ", filename_pt104)
    hk2 = np.array(rddt.read_hk_name_from_csv(filename_pt104, header2))
    t_pt100_2 = hk2.mean()
    hk4 = np.array(rddt.read_hk_name_from_csv(filename_pt104, header4))
    t_pt100_1 = hk4.mean()
    
    return t_avg, t_max, t_pt100_1, t_pt100_2
    


def get_gain_adu(session, signal_adu, verbose=True):

    data_path = os.path.join(session, cst.dataDirName)
    
    if verbose:
        print("   Getting DEMUX data from session ", session)       

    erro_adu = np.zeros((cst.nColPerDemux))
    for col in range (cst.nColPerDemux):
        erro_adu[col] = rddt.read_col_science_from_dir(
            data_path = data_path,
            col_id = col,
            flatten = True,
            remove_dc = False,
            verbose = False)[0].mean()
        
    gain_adu = erro_adu / signal_adu
    
    return gain_adu


def get_data(signal_adu, force_read_raw_data=False, verbose=False):
    
    """
    """

    path = "../.."

    if verbose:
        print("/----------------------------------------------------------")
        print("/ Gain vs temperature characterisation")
        print("/----------------------------------------------------------\n")

    if verbose:
        print("Searching sessions")
    
    savefile = 'data.npz'
    if os.path.isfile(os.path.join(path, savefile)):
        print("Retrieving data from npz file ", savefile)
        data = np.load(os.path.join(path, savefile))
        t_dmx_avg = data['t_dmx_avg']
        t_dmx_max = data['t_dmx_max']
        t_pt100_foot = data['t_pt100_foot']
        t_pt100_top = data['t_pt100_top']
        gain_adu = data['gain_adu']

    else:
        sessions = [
            d
            for d in os.listdir(path)
            if os.path.isdir(os.path.join(path, d))
            and d[16:28]=="GAIN_VS_TEMP"
        ]
        
        nb_sessions = len(sessions)
        if verbose:
            print("Found {0:} sessions".format(nb_sessions))
        gain_adu = np.zeros((cst.nColPerDemux, nb_sessions))
        t_dmx_avg = np.zeros((nb_sessions))
        t_dmx_max = np.zeros((nb_sessions))
        t_pt100_foot = np.zeros((nb_sessions))
        t_pt100_top = np.zeros((nb_sessions))
    
        for i, session in enumerate(np.sort(sessions)):
            if verbose:
                print("Session {0:}/{1:} --->  ".format(1+i, nb_sessions), session)
            t_dmx_avg[i], t_dmx_max[i], t_pt100_foot[i], t_pt100_top[i] = get_temperatures(
                os.path.join(path, session))
            gain_adu[:, i] = get_gain_adu(os.path.join(path, session), signal_adu)
        
        np.savez_compressed(os.path.join(path, savefile), t_dmx_avg=t_dmx_avg, t_dmx_max=t_dmx_max,
                            t_pt100_foot=t_pt100_foot, t_pt100_top=t_pt100_top, gain_adu=gain_adu)

    return(t_dmx_avg, t_dmx_max, t_pt100_foot, t_pt100_top, gain_adu)


def plot_gain_vs_temp(temperature, gain, signal, deg_poly_fit = 2, up_down=True, yresiduals_autoscale=True, verbose=True):

    path = "../.."
    plot_path = os.path.join(path, cst.plotDirName)
    gt.createdir(plot_path)

    plot_filename = 'gain_vs_temp.png'
    full_plot_filename = os.path.join(plot_path, plot_filename)

    dt = 1
    vec_temp = np.arange(temperature.min()-dt, temperature.max()+dt, 0.1)  # Temperature array for the fits

    nb_sessions = gain.shape[1]
    half = nb_sessions // 2

    fig1 = plt.figure(figsize=(8, 11), constrained_layout=True)
    gs = fig1.add_gridspec(cst.nColPerDemux + 1, 1,
                           height_ratios=[2.5] + [1.1] * cst.nColPerDemux)
    ax1 = fig1.add_subplot(gs[0])
    suptitle = "Gain versus temperature"
    fig1.suptitle(suptitle, fontsize=12)

    ax_res = [fig1.add_subplot(gs[i + 1], sharex=ax1) for i in range(cst.nColPerDemux)]

    ymax_res = np.zeros(cst.nColPerDemux)
    for col in range(cst.nColPerDemux):
        if up_down:
            ax1.scatter(temperature[:half], gain[col, :half], marker='.',
                        color=_darken(colors[col]))
            ax1.scatter(temperature[half:], gain[col, half:], marker='.',
                        color=_lighten(colors[col]))
        else:
            ax1.scatter(temperature, gain[col, :], marker='.', color=colors[col])
        p = np.poly1d(np.polyfit(temperature, gain[col, :], deg_poly_fit))

        # Building the fit equation for the legend, skipping terms with tiny coefficients
        terms = []
        for deg in range(deg_poly_fit, -1, -1):
            coeff = p.coeffs[deg_poly_fit - deg]
            if abs(coeff) < 1e-8:
                continue
            if deg == 0:
                term = '{0:.3g}'.format(abs(coeff))
            elif deg == 1:
                term = '{0:.3g} $\\times 10^{{-6}}$ T'.format(abs(coeff) * 1e6)
            else:
                term = '{0:.3g} T^{1}'.format(abs(coeff), deg)
            if coeff < 0:
                terms.append('- ' + term)
            else:
                terms.append('+ ' + term)
        if terms:
            terms[0] = terms[0].lstrip('+').strip()
        fit_txt = 'Gain col({0:}) = '.format(col) + ' '.join(terms)
        ax1.plot(vec_temp, p(vec_temp), ':', c=colors[col], label=fit_txt)

        # Residuals: measurement minus fit, scaled by 1e6
        residual = (gain[col, :] - p(temperature)) * 1e6
        if up_down:
            ax_res[col].plot(temperature[:half], residual[:half], marker='.',
                             linestyle='none', c=_darken(colors[col]))
            ax_res[col].plot(temperature[half:], residual[half:], marker='.',
                             linestyle='none', c=_lighten(colors[col]))
        else:
            ax_res[col].plot(temperature, residual, marker='.', linestyle='none',
                             c=colors[col])
        ax_res[col].axhline(0, color='k', linewidth=0.5)
        ax_res[col].set_title(r'Col {0:} gain - Fit ($\times 10^{{6}}$)'.format(col))
        ax_res[col].grid(True, linestyle='--')
        ax_res[col].yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        ax_res[col].ticklabel_format(style='plain', axis='y')

        # Keeping the largest residual magnitude for this column
        ymax_res[col] = np.max(np.abs(residual))

    if yresiduals_autoscale:
        # Autoscale each residual plot individually, centered on 0
        for col in range(cst.nColPerDemux):
            ax_res[col].set_ylim(-ymax_res[col], ymax_res[col])
    else:
        # Same Y range on all residual plots, centered on 0
        ymax_global = np.max(ymax_res)
        for col in range(cst.nColPerDemux):
            ax_res[col].set_ylim(-ymax_global, ymax_global)

    ax1.set_ylabel("DEMUX " + signal + " + Error gain (ADU/ADU)")
    ax1.legend(loc='best')
    ax1.xaxis.set_major_locator(MultipleLocator(2))
    ax_res[-1].set_xlabel("Temperature measured on DEMUX fixation (°C)")

    plt.savefig(full_plot_filename, dpi=300, bbox_inches='tight')
    if verbose:
        print("Plot saved in file ", full_plot_filename)


def plot_temp_vs_reftemp(t_dmx_avg, t_dmx_max, t_pt100_foot, t_pt100_top, up_down=True, verbose=False):

    path = "../.."
    plot_path = os.path.join(path, cst.plotDirName)
    gt.createdir(plot_path)

    nb_sessions = t_dmx_avg.shape[0]
    half = nb_sessions // 2

    # Temperature ratios with respect to the reference temperature pt104_foot
    ratios = [
        ('pt104 top / pt104 foot', t_pt100_top / t_pt100_foot),
        ('dmx avg / pt104 foot', t_dmx_avg / t_pt100_foot),
        ('dmx max / pt104 foot', t_dmx_max / t_pt100_foot),
    ]
    ref_temp = t_pt100_foot

    fig2, ax_ratio = plt.subplots(len(ratios), 1, figsize=(8, 9),
                                  constrained_layout=True, sharex=True)
    for i, (label, ratio) in enumerate(ratios):
        if up_down:
            ax_ratio[i].scatter(ref_temp[:half], ratio[:half], marker='.',
                                color=_lighten(colors[0]))
            ax_ratio[i].scatter(ref_temp[half:], ratio[half:], marker='.',
                                color=_darken(colors[0]))
        else:
            ax_ratio[i].scatter(ref_temp, ratio, marker='.', color=colors[0])

        # Quadratic fit of the ratio versus reference temperature
        p = np.poly1d(np.polyfit(ref_temp, ratio, 2))
        vec_fit = np.linspace(ref_temp.min(), ref_temp.max(), 100)

        # Building the fit equation for the legend: coefficients in engineering
        # notation (power of 10 is a multiple of 3), powers as LaTeX exponents,
        # no explicit T^1
        terms = []
        for deg in range(2, -1, -1):
            coeff = p.coeffs[2 - deg]
            if abs(coeff) < 1e-8:
                continue

            # Engineering notation: exponent multiple of 3
            expo = 3 * int(np.floor(np.log10(abs(coeff))) // 3)
            mantissa = abs(coeff) / 10 ** expo
            mantissa = float('{0:.3g}'.format(mantissa))
            if expo == 0:
                pow_txt = '{0:g}'.format(mantissa)
            else:
                pow_txt = '{0:g} $\\times 10^{{{1}}}$'.format(mantissa, expo)

            if deg == 0:
                term = pow_txt
            elif deg == 1:
                term = pow_txt + ' $T$'
            else:
                term = pow_txt + ' $T^{{{0}}}$'.format(deg)
            if coeff < 0:
                terms.append('- ' + term)
            else:
                terms.append('+ ' + term)
        if terms:
            terms[0] = terms[0].lstrip('+').strip()
        fit_txt = '{0:} = '.format(label) + ' '.join(terms)
        ax_ratio[i].plot(vec_fit, p(vec_fit), ':', c=colors[1], label=fit_txt)

        ax_ratio[i].set_ylabel(label)
        ax_ratio[i].legend(loc='best')
        ax_ratio[i].grid(True, linestyle='--')
    ax_ratio[-1].set_xlabel("pt104 foot temperature (°C)")
    ax_ratio[-1].xaxis.set_major_locator(MultipleLocator(2))

    ratio_filename = 'temp_ratios_vs_temp.png'
    full_ratio_filename = os.path.join(plot_path, ratio_filename)
    plt.savefig(full_ratio_filename, dpi=300, bbox_inches='tight')
    if verbose:
        print("Plot saved in file ", full_ratio_filename)

#-------------------------------------------------------------------------------------

def thermal_behaviour(verbose=False):
    
    path = "../.."
    master_session = os.path.basename(os.path.realpath(path))

    # Looking for signal name and signal value in ADU
    fdbk_name = "FDBK"
    ofco_name = "OFCO"
    adu_ncar = 4
    if master_session.find(fdbk_name) != -1:
        signal = fdbk_name
    else:
        signal = ofco_name
    i_adu = master_session.find(signal) + len(signal)
    signal_adu = int(master_session[i_adu: i_adu + adu_ncar])

    # Reading the data
    t_dmx_avg, t_dmx_max, t_pt100_foot, t_pt100_top, gain_adu = get_data(signal_adu, verbose=verbose)

    # Making plots
    plot_gain_vs_temp(t_dmx_avg, gain_adu, signal, deg_poly_fit = 1, up_down=True, yresiduals_autoscale=True, verbose=verbose)        
    plot_temp_vs_reftemp(t_dmx_avg, t_dmx_max, t_pt100_foot, t_pt100_top, verbose=verbose)
    