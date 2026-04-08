#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module name: dre_noise_models.py
Description: This module contains the noise models of the ATHENA XIFU-DRE.
Author: Yann PAROT / IRAP
Date: 09/03/2026

Copyright (C) 2021-2030 Yann PAROT, IRAP Toulouse.
This file is part of the ATHENA X-IFU DRE data analysis tools software.

analysis-tools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

---------------------------------------------------------------------------------

yann.parot@cnrs.fr

---------------------------------------------------------------------------------
"""

# --------------standard import----------------#

# --------------third party import----------------#

import types

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import LogLocator
from scipy import constants


# --------------local import----------------#

# --------------Constants----------------#

# --------------Functions----------------#

def celsius2kelvin(tc):
    """
    Converts a degrees Celsius to degrees Kelvin

    Parameters
    ----------
    tc : float
        Temperature in degrees Celsius.

    Returns
    -------
    tk : float
        Temperature in degrees Kelvin.

    """
    tk = tc + constants.zero_Celsius
    return tk


def thermal_noise(r, tk):
    """
    Returns the thermal voltage noise for a resistor R at TK degrees Kelvin

    Parameters
    ----------
    r : float
        Resistance (Ohms).
    tk : float
        Degrees Kelvin.

    Returns
    -------
    en : float
        Voltage noise (V/√Hz).

    """
    en = np.sqrt(4 * constants.Boltzmann * r * tk)
    return en


def shot_noise(i):
    """
    Returns the shot noise for a current source.

    Parameters
    ----------
    i : float
        Current (A).

    Returns
    -------
    inn : float
        shot noise (A/√Hz).

    """
    inn = np.sqrt(2 * constants.elementary_charge * i)
    return inn


def quantization_noise(vpp, n, fs):
    """
    Returns the quantization noise for a DAC/ADC of Vpp dynamic range, N bits and a sampling frequency Fs

    Parameters
    ----------
    vpp : float
        Input/Output analog dynamic (V).
    n : int
        Number of bits.
    fs : float
        Sampling frequency (Hz).

    Returns
    -------
    e_q : float
        Quantization noise (V/√Hz).
    rms_q : float
        RMS value of the quantization noise (Vrms).
    lsb : float
        Size of the Least Significant Bite in V.

    """
    lsb = vpp / 2 ** n
    rms_q = lsb / np.sqrt(12)
    e_q = rms_q / np.sqrt(fs / 2)
    return e_q, rms_q, lsb


def aliasing(noise_func, analog_bw, fs, n):
    """
    Returns the n-point aliased waveform for an analog noise, taken into account on a window [0:analog_bw] sampled at fs

    Parameters
    ----------
    noise_func : function of 1 parameter
        Noise function en(f) describing the noise with f the frequency in Hz.
    analog_bw : float
        Highest frequency to take into account for the analog signal.
    fs : float
        Sampling frequency in Hz.
    n : int
        Number of points in the aliased vector [0 : Fs/2].

    Returns
    -------
    f_num : array[float]
        Aliased frequency vector.
    e_fold : float
        Aliased noise.

    """

    if not isinstance(noise_func, types.FunctionType):
        print("noise_func n'est pas une fonction")
        return

    f_nyq = fs / 2

    if analog_bw < f_nyq:
        print("La bande analogique est inférieure à la fréquence de Nyquist")
        return

    # freqs numériques
    f_num = np.unique(np.round(np.logspace(0, np.ceil(np.log10(f_nyq)), n)))  # arrondi au Hz et unique
    f_num = f_num[f_num < f_nyq]
    f_num = np.append(f_num, f_nyq)  # on est sûr de pas dépasser la fréquence de nyquist et d'avoir un point dessus
    # on ajoute les repliement "impaire" Fs-f
    f_num = np.unique(np.sort(np.append(f_num, fs - f_num)))

    # Fréquence analogique (Hz)
    # ce sont les fréquences numériques périodisée sur toute la bande
    k_alias_max = np.ceil(analog_bw / fs)

    # Construction des alias
    k_alias = np.arange(0, k_alias_max + 1)[:, None]  # vecteur colonne
    f_analogique = f_num + k_alias * fs  # ajoute tous les multiples de Fs
    f_analogique = f_analogique[(f_analogique <= analog_bw)]  # ne garde que les fréquences inférieure

    # Tri final
    f_analogique = np.unique(
        np.sort(f_analogique))  # on s'assure de l'unicité (cas où on a une freq à 0 et Fs par exemple)

    # ==== CALCUL ====

    # repliement fréquentiel
    f_modulo = np.mod(f_analogique, fs)
    f_repliee = np.where(f_modulo > f_nyq, fs - f_modulo, f_modulo)

    # indices des fréquences identiques
    indices = np.searchsorted(f_num, f_repliee)

    # Somme des contributions
    e_fold_squared = np.zeros_like(f_analogique)
    e_vect = noise_func(2 * np.pi * f_analogique)
    np.add.at(e_fold_squared, indices, e_vect ** 2)

    # On prend la NSD mais que jusqu'à Fnyq
    ind_nyq = np.where(f_analogique == f_nyq)[0][0]  ## existe vu la définition de f_analogique
    f_num = f_analogique[:ind_nyq]
    e_fold = np.sqrt(e_fold_squared[:ind_nyq])

    return f_num, e_fold


def export_2vectors_to_file(x, y, filename, verbose=False):
    if len(x) != len(y):
        raise ValueError("The vectors have different sizes")
    with open(filename, 'w') as file:
        for i in range(len(x)):
            file.write('{0:}  {1:}\n'.format(x[i], y[i]))
    if verbose:
        print("Data saved in ", filename)


def adc_chain_noise_analog(tc, zin, ein, fs):
    """
    Returns the noise of the DRE ADC chain considering an output impedance of the previous stage equal to zin.
    The noise is "analog": no aliasing is performed
    quantization noise is also returned

    Parameters
    ----------
    tc : float
        temperature in degrees Celsius.
    zin : float
        Previous stage output impedance (Ohms).
        real number (resistor)
        Noiseless
    ein : function
        Previous stage output noise (ADC chain input noise source)
        Is function of the angular frequency w (rad/s)
    fs : float
        ADC sampling frequency (Hz)

    Returns
    -------
    en : function of w (angular frequency)
        analog noise voltage.
    eq : float
        quantization noise voltage (constant defined on [0 : fs/2].

    """
    # -----------Constantes------------#

    ##Composants passifs##

    # Résistances
    R2 = 100
    R3 = 10e3
    R4 = 10e3
    R5 = 15
    R6 = 15

    # Capacités
    C1 = 15e-12
    C2 = 270e-12

    ##Composants actifs##

    # RHF200
    G_RHF200 = 1  # RHF200 is set for a gain of 1. The RHF200 BW is way larger than many other BW. It is not taken into account
    E_RHF200 = 10e-9  # Gain = 1.

    # ADC
    VPP = 2
    N = 14
    ADC_THERMAL_NOISE_RMS = 1.25  # LSB rms
    ADC_THERMAL_NOISE_BW = 150e6 / 2
    ADC_INPUT_BW = 650e6

    # -----------Eléments de bruits calculés------------#

    # Bruits thermiques des résistances
    er2 = thermal_noise(R2, celsius2kelvin(tc))
    er3 = thermal_noise(R3, celsius2kelvin(tc))
    er4 = thermal_noise(R4, celsius2kelvin(tc))
    er5 = thermal_noise(R5, celsius2kelvin(tc))
    er6 = thermal_noise(R6, celsius2kelvin(tc))

    # Bruit de quantification et LSB de l'ADC
    eq, _, lsb = quantization_noise(VPP, N, fs)

    # Bruit thermique
    e_adc_thermal_nyquist = ADC_THERMAL_NOISE_RMS * lsb / np.sqrt(
        ADC_THERMAL_NOISE_BW)  # etalage sur Nyquist vis à vis de la mesure datasheet
    thermal_equi_bw = np.pi / 2 * ADC_INPUT_BW
    e_adc_th = e_adc_thermal_nyquist * np.sqrt(
        ADC_THERMAL_NOISE_BW / thermal_equi_bw)  # étalage sur toute la bande vis à vis repliement

    # -----------Fonctions de transfert------------#

    # Input ADC
    tau_c_adc = 1 / (2 * np.pi * ADC_INPUT_BW)
    h_adc = lambda w: 1 / (1 + 1j * tau_c_adc * w)

    # Gain pour R5 ou R6
    h_r5_r6 = lambda w: 1 / (1 + 1j * (R5 + R6) * C2 * w) * h_adc(w)

    # Gain entrée RHF200
    h_rhf_200 = lambda w: G_RHF200 * h_r5_r6(w)

    # Gain appliqué à l'entrée (zin)
    zequi_1 = lambda w: (R3 + R4) * (1 + 1j * R2 * C1 * w) / (
            1 + 1j * C1 * (R2 + R3 + R4) * w)  # impédance équivalente vue après zin
    h_zin = lambda w: zequi_1(w) / (zin + zequi_1(w)) * h_rhf_200(w)

    # Gain de R2
    zequi_2 = 1 / (1 / zin + 1 / (R3 + R4))
    h_r2 = lambda w: 1j * zequi_2 * C1 * w / (1 + 1j * (zequi_2 + R2) * C1 * w) * h_rhf_200(w)

    # Gain de R3 ou R4
    zequi_3 = lambda w: zin * (1 + 1j * R2 * C1 * w) / (1 + 1j * C1 * (R2 + zin) * w)
    h_r3_r4 = lambda w: 1 / (1 + (R3 + R4) / zequi_3(w)) * h_rhf_200(w)

    # -----------Bruits ramenés en entrée ADC------------#
    ezin_out = lambda w: np.abs(h_zin(w)) * ein(w)
    er2_out = lambda w: np.abs(h_r2(w)) * er2
    er3_out = lambda w: np.abs(h_r3_r4(w)) * er3
    er4_out = lambda w: np.abs(h_r3_r4(w)) * er4
    er5_out = lambda w: np.abs(h_r5_r6(w)) * er5
    er6_out = lambda w: np.abs(h_r5_r6(w)) * er6
    erhf_out = lambda w: np.abs(h_rhf_200(w)) * E_RHF200
    e_adc_th_out = lambda w: np.abs(h_adc(w)) * e_adc_th

    # -----------Bruit total------------#
    en = lambda w: np.sqrt(
        ezin_out(w) ** 2 + er2_out(w) ** 2 + er3_out(w) ** 2 + er4_out(w) ** 2 + er5_out(w) ** 2 + er6_out(
            w) ** 2 + erhf_out(w) ** 2 + e_adc_th_out(w) ** 2)

    return en, eq


def awaxe_fb_chain_noise(tc):
    """
    Returns the noise of the AWAXE version of the FB chain as a model composed of a noise source en and the output equivalent impedance zout (real & noiseless)

    Parameters
    ----------
    tc : float
        Temperature in degrees Celsius.

    Returns
    -------
    en : function of w (angular frequency)
            analog noise voltage.
    zout : float
        Output impedance in Ohms (real).

    """
    ##Composants passifs##

    # Résistances
    R1a = 50
    R1b = 50

    ##Composants actifs##

    # DAC
    IDAC = 20e-3

    # -----------Eléments de bruits calculés------------#

    # Bruits thermiques des résistances
    er1a = thermal_noise(R1a, celsius2kelvin(tc))
    er1b = thermal_noise(R1a, celsius2kelvin(tc))

    # Bruit de shot du DAC
    i_shot_dac = shot_noise(IDAC)

    # Bruit du shot du DAC combiné au miroir de courant AWAXE
    i_shot_tot = np.sqrt(2) * i_shot_dac

    # -----------Fonctions de transfert------------#
    # Miroir de courant
    h_current_mirror = R1a  # Vrai seulement si R1a = R1b / Bande passante non connue

    # -----------Bruit total------------#
    en = lambda w: np.sqrt(er1a ** 2 + er1b ** 2 + (h_current_mirror * i_shot_tot) ** 2)

    # -----------Impédanced de sortie------------#
    zout = R1a + R1b

    return en, zout


def rhf200_fb_chain_noise(tc):
    """
    Returns the noise of the RHF200 version of the FB chain as a model composed of a noise source en and the output equivalent impedance zout (real & noiseless)

    Parameters
    ----------
    tc : float
        Temperature in degrees Celsius.

    Returns
    -------
    en : function of w (angular frequency)
            analog noise voltage.
    zout : float
        Output impedance in Ohms (real).

    """
    ##Composants passifs##

    # Résistances
    R0a = 50
    R0b = 50
    R1a = 50
    R1b = 50

    ##Composants actifs##

    # DAC
    IDAC = 20e-3

    # RHF200
    G_RHF200 = 1  # RHF200 is set for a gain of 1. The RHF200 BW is way larger than many other BW. It is not taken into account
    E_RHF200 = 10e-9  # Gain = 1.

    # -----------Eléments de bruits calculés------------#

    # Bruits thermiques des résistances
    er0a = thermal_noise(R0a, celsius2kelvin(tc))
    er0b = thermal_noise(R0b, celsius2kelvin(tc))
    er1a = thermal_noise(R1a, celsius2kelvin(tc))
    er1b = thermal_noise(R1b, celsius2kelvin(tc))

    # Bruit de shot du DAC
    i_shot_dac = shot_noise(IDAC)

    # -----------Fonctions de transfert------------#

    # R0a / R0b
    h_r0 = G_RHF200

    # DAC
    h_dac = R0a * G_RHF200

    # -----------Bruit total------------#
    en = lambda w: np.sqrt(
        (er0a * h_r0) ** 2 + (er0b * h_r0) ** 2 + (i_shot_dac * h_dac) ** 2 + E_RHF200 ** 2 + er1a ** 2 + er1b ** 2)

    # -----------Impédance de sortie------------#
    zout = R1a + R1b

    return en, zout


def ofco_chain_noise(tc, dac_coarse, dac_fine, dac_fast):
    """
    Returns the noise of the OFCO chain as a model composed of a noise source en and the output equivalent impedance zout (real & noiseless)

    Parameters
    ----------
    tc : float
        Temperature in degrees Celsius.
    dac_coarse : int 0 to 4095
        DAC coarse setting.
    dac_fine : int 0 to 4095
        DAC fine setting.
    dac_fast : int 0 to 7
        DAC fast setting (R2R).

    Returns
    -------
    en : TYPE
        DESCRIPTION.
    zout : TYPE
        DESCRIPTION.

    """
    # Vérification des valeurs d'entrées
    if dac_coarse > 4095:
        dac_coarse = 4095

    if dac_fine > 4095:
        dac_fine = 4095

    if dac_fast > 4095:
        dac_fast = 4095

    ##Composants passifs##

    # Résistances
    R1 = 2e3
    R2 = 2e3
    R3 = 1e3
    R4 = 2e3
    R5 = 2e3
    R6 = 50
    R7 = 1e3
    R8 = 1e3
    R9 = 1e3
    R10 = 499
    R11 = 499
    R12 = 499
    R13 = 1e3
    R14 = 249
    R15 = 169
    R16a = 49.9
    R16b = 49.9

    # Capacités
    C1 = 3.3e-9
    C2 = 68e-12

    ##Composants actifs##

    # ADA4084
    # --bruit
    E_ADA4084_1KHZ = 3.9e-9  # V/rtHz1kHz
    E_ADA4084_1HZ = 9.1e-9
    I_ADA4084 = 0.55e-12  # pA/rtHz
    # E_BF_ADA4084 = 0.14e-6  #Vpp sur 0.1 à 10Hz utilisation courbe datasheet à la place: E_ADA4084_1HZ
    # --fonction de transfert
    # GBP = 15.9e6 #Produit gain bande en Hz
    AVO_DB = 104  # dB
    UGC = 9.9e6  # Unity Gain Crossover en Hz
    FC3 = 2e6
    FC2 = 300e3
    AMP_FC2_DB = 32  # dB

    # ISL21010CHF333
    # --bruit
    # Bandes
    BAND1_LOW, BAND1_HIGH = 0.1, 10
    BAND2_LOW, BAND2_HIGH = 10, 1000
    # Mesures
    NOISE_BAND1_PP = 95e-6
    NOISE_BAND2_RMS = 40e-6
    # --Impédance de sortie

    # RHF200
    E_RHF200 = 15e-9  # for gain 2
    G_RHF200 = 2  # RHF200 is set for a gain of 1. The RHF200 BW is way larger than many other BW. It is not taken into account

    # -----------Eléments de bruits calculés------------#

    # Bruits thermiques des résistances
    er1 = thermal_noise(R1, celsius2kelvin(tc))
    er2 = thermal_noise(R2, celsius2kelvin(tc))
    er3 = thermal_noise(R3, celsius2kelvin(tc))
    er4 = thermal_noise(R4, celsius2kelvin(tc))
    er5 = thermal_noise(R5, celsius2kelvin(tc))
    er6 = thermal_noise(R6, celsius2kelvin(tc))
    # er7 = thermal_noise(R7,celsius2kelvin(tc)) equivalent noise of R2R ladder use instead
    # er8 = thermal_noise(R8,celsius2kelvin(tc))
    # er9 = thermal_noise(R9,celsius2kelvin(tc))
    # er10 = thermal_noise(R10,celsius2kelvin(tc))
    # er11 = thermal_noise(R11,celsius2kelvin(tc))
    # er12 = thermal_noise(R12,celsius2kelvin(tc))
    er13 = thermal_noise(R13, celsius2kelvin(tc))
    er14 = thermal_noise(R14, celsius2kelvin(tc))
    er15 = thermal_noise(R15, celsius2kelvin(tc))
    er16a = thermal_noise(R16a, celsius2kelvin(tc))
    er16b = thermal_noise(R16b, celsius2kelvin(tc))

    # Modèle de bruit en tension ADA4084
    k_ada4084 = (E_ADA4084_1HZ) ** 2 - E_ADA4084_1KHZ ** 2  # Coefficient du bruit 1/f
    # f_corner_ada4084 = k_ada4084 / (E_ADA4084_1KHZ**2)        # Fréquence de coin (corner frequency)
    e_ada4084 = lambda w: np.sqrt(E_ADA4084_1KHZ ** 2 + k_ada4084 / (w / (2 * np.pi)))  # Bruit total (1/f + blanc)

    # Bruit en sortie du suiveur en fonction des résistances feedback, entrées et capa
    e_ada4084_coarse = lambda w: np.sqrt(
        e_ada4084(w) ** 2 + (I_ADA4084 ** 2) * (R2 ** 2 + R1 ** 2 / (1 + 1j * R1 * C1 * w) ** 2))
    e_ada4084_fine = lambda w: np.sqrt(e_ada4084(w) ** 2 + (I_ADA4084 ** 2) * (R5 ** 2 + R4 ** 2))

    # Modèle de bruit ISL21010CHF333
    # Conversion pp → RMS (approximation) En vrai cela dépend du nombre de points: Vpp = 2*sigma * sqrt(2*ln(N)) et non Vpp = 6 * sigma
    noise_band1_rms = NOISE_BAND1_PP / 6.6  # trouvé dans la littérature comme approche prudente pour ce cas.
    p1 = noise_band1_rms ** 2
    p2 = NOISE_BAND2_RMS ** 2
    a1 = np.log(BAND1_HIGH / BAND1_LOW)
    # a2 = np.log(BAND2_HIGH / BAND2_LOW)
    b1 = BAND1_HIGH - BAND1_LOW
    b2 = BAND2_HIGH - BAND2_LOW
    # Résolution analytique
    e_white_sq_isl = (p2 - p1) / (b2 - b1)
    k_isl = (p1 - b1 * e_white_sq_isl) / a1
    e_white_isl = np.sqrt(e_white_sq_isl)
    # f_corner_isl = k_isl / (e_white_isl**2)          # Fréquence de coin (corner frequency)
    e_isl = lambda w: np.sqrt(e_white_isl ** 2 + k_isl / (w / (2 * np.pi)))  # Bruit total (1/f + blanc)

    # -----------Modélisation impédances---------#
    # Impédance de sortie ISL21010CHF333
    R1_MODEL_ISL = 0.5
    C1_MODEL_ISL = 800e-6
    L1_MODEL_ISL = 30e-6
    L2_MODEL_ISL = 4e-6
    R2_MODEL_ISL = 1.5
    C2_MODEL_ISL = 1.9e-6
    L3_MODEL_ISL = 500e-9
    R3_MODEL_ISL = 0.1
    z1_model_isl = lambda w: R1_MODEL_ISL + 1 / (1j * C1_MODEL_ISL * w)
    z2_model_isl = lambda w: z1_model_isl(w) * 1j * L1_MODEL_ISL * w / (
            z1_model_isl(w) + 1j * L1_MODEL_ISL * w) + 1j * L2_MODEL_ISL * w
    z3_model_isl = lambda w: R2_MODEL_ISL + 1 / (1j * C2_MODEL_ISL * w)
    z_out_isl = lambda w: z2_model_isl(w) * z3_model_isl(w) / (
            z2_model_isl(w) + z3_model_isl(w)) + 1j * L3_MODEL_ISL * w + R3_MODEL_ISL

    # Capacité de découplage 300904106475JC
    R0_CD1 = 1.905e-3
    R1_CD1 = 3.597e-3
    R2_CD1 = 71.83e-3
    R3_CD1 = 72.15e-3
    R4_CD1 = 1e9
    L1_CD1 = 1.685e-9
    L2_CD1 = 1.246e-9
    L3_CD1 = 11.99e-9
    C1_CD1 = 3.807e-6
    C2_CD1 = 799e-9
    C3_CD1 = 94e-9
    z1_cd1 = lambda w: 1j * L1_CD1 * w + R1_CD1 + 1 / (1j * C1_CD1 * w)
    z2_cd1 = lambda w: 1j * L2_CD1 * w + R2_CD1 + 1 / (1j * C2_CD1 * w)
    z3_cd1 = lambda w: 1j * L3_CD1 * w + R3_CD1 + 1 / (1j * C3_CD1 * w)
    z_cd1 = lambda w: R0_CD1 + 1 / (1 / z1_cd1(w) + 1 / z2_cd1(w) + 1 / z3_cd1(w) + 1 / R4_CD1)

    # Capacité de découplage 300904102104JC
    R0_CD2 = 100e-6
    R1_CD2 = 17.54e-3
    R2_CD2 = 2.449
    R3_CD2 = 1.385
    R4_CD2 = 1e9
    L1_CD2 = 1.041e-9
    L2_CD2 = 5.32e-9
    L3_CD2 = 21.96e-9
    C1_CD2 = 83e-9
    C2_CD2 = 15e-9
    C3_CD2 = 1.9e-9
    z1_cd2 = lambda w: 1j * L1_CD2 * w + R1_CD2 + 1 / (1j * C1_CD2 * w)
    z2_cd2 = lambda w: 1j * L2_CD2 * w + R2_CD2 + 1 / (1j * C2_CD2 * w)
    z3_cd2 = lambda w: 1j * L3_CD2 * w + R3_CD2 + 1 / (1j * C3_CD2 * w)

    z_cd2 = lambda w: R0_CD2 + 1 / (1 / z1_cd2(w) + 1 / z2_cd2(w) + 1 / z3_cd2(w) + 1 / R4_CD2)

    # Capacité de découplage équivalente à 2 CD1 parallèles, parallèles à 2 Cd2  parallèles
    z_cd = lambda w: 1 / (2 / (z_cd1(w)) + 2 / (z_cd2(w)))

    # Résistance équivalente R2R (On considère R6 négligeable)
    r2r_eq_b0 = R13 * R9 / (R13 + R9) + R12
    r2r_eq_b1 = r2r_eq_b0 * R8 / (r2r_eq_b0 + R8) + R11
    r2r_eq_b2 = r2r_eq_b1 * R7 / (r2r_eq_b1 + R7) + R10

    # -----------Fonctions de transfert------------#

    # ADA4084
    avo = 10 ** (AVO_DB / 20)
    amp_fc2 = 10 ** (AMP_FC2_DB / 20)
    # calcul des fréquences des poles et zeros
    # dernier zero
    tau4 = 1 / (2 * np.pi * UGC)
    # pole
    tau3 = 1 / (2 * np.pi * FC3)
    # zero
    tau2 = 1 / (2 * np.pi * FC2)
    # zero
    fc1 = (amp_fc2 * FC2) / avo
    tau1 = 1 / (2 * np.pi * fc1)
    # fonction de transfert en BO
    h_ada4084_bo = lambda w: avo * (1 + 1j * tau3 * w) ** 2 / (
            (1 + 1j * tau1 * w) * (1 + 1j * tau2 * w) * (1 + 1j * tau4 * w))
    # Fonction de transfert en suiveur (gain unitaire)
    h_ada4084_bf = lambda w: h_ada4084_bo(w) / (1 + h_ada4084_bo(w))

    # ISL21010CHF333 vers DAC121S101
    h_isl_dac = lambda w: z_cd(w) / (z_cd(w) + z_out_isl(w))

    # fonction de transfert des dac fine et coarse
    h_dac_slow = lambda code: code / 4096

    # suiveur coarse
    z_r2r = lambda w: r2r_eq_b2 * R14 / (r2r_eq_b2 + R14 + 1j * r2r_eq_b2 * R14 * C2 * w)
    h_out_suiveur_coarse = lambda w: z_r2r(w) / (z_r2r(w) + R3)

    # DAC coarse
    h_dac_coarse_out = lambda w: 1 / (1 + 1j * R1 * C1 * w) * h_ada4084_bf(w) * h_out_suiveur_coarse(w)

    # R14
    r_vue_r14 = r2r_eq_b2 * R3 / (r2r_eq_b2 + R3)
    h_r14 = lambda w: 1 / (1 + R14 / r_vue_r14 + 1j * R14 * C2 * w)

    # R2R
    # Pour la fonction de transfert du R2R, on considère que R7,R8,R9,R13 = 2R et R10,R11,R12 = R pour alléger le calcul.
    h_r2r = lambda code: code / 8
    h_r2r_out = lambda w: h_out_suiveur_coarse(w)  # approximation: on présente 1k comme le cas du suiveur coarse

    # Le bruit du réseau R2R peut se simplifier aussi avec les égalités sur les résistances précédentes. Il faut lui appliquer H_R2R_out
    e_r2r = er13

    # Le MUX est considéré sans bruit.

    # DAC fine
    h_dac_fine_out = lambda w, code: h_ada4084_bf(w) * h_r2r(code) * h_r2r_out(w)

    # -----------Bruits ramenés en sortie------------#
    # ramenés en entrée RHF200
    e_dac_fine_out = lambda w: np.abs(e_isl(w) * h_isl_dac(w) * h_dac_slow(dac_fine) * h_dac_fine_out(w, dac_fast))
    e_r4_out = lambda w: er4 * np.abs(h_dac_fine_out(w, dac_fast))
    e_r5_out = lambda w: er5 * h_r2r(dac_fast) * np.abs(h_r2r_out(w))
    e_r6_out = lambda w: er6 * h_r2r(dac_fast) * np.abs(h_r2r_out(w))
    e_r2r_out = lambda w: e_r2r * np.abs(h_r2r_out(w))
    e_dac_coarse_out = lambda w: np.abs(e_isl(w) * h_isl_dac(w) * h_dac_slow(dac_coarse) * h_dac_coarse_out(w))
    e_r1_out = lambda w: er1 * np.abs(h_dac_coarse_out(w))
    e_r2_out = lambda w: er2 * np.abs(h_out_suiveur_coarse(w))
    e_r3_out = lambda w: er3 * np.abs(h_out_suiveur_coarse(w))
    e_r14_out = lambda w: er14 * np.abs(h_r14(w))
    e_ada4084_coarse_out = lambda w: np.abs(e_ada4084_coarse(w) * h_out_suiveur_coarse(w))
    e_ada4084_fine_out = lambda w: np.abs(e_ada4084_fine(w) * h_r2r(dac_fast) * h_r2r_out(w))

    # total en entrée RHF200
    e_in_rhf200 = lambda w: np.sqrt(
        e_dac_fine_out(w) ** 2 + e_r4_out(w) ** 2 + e_r5_out(w) ** 2 + e_r6_out(w) ** 2 + e_r2r_out(
            w) ** 2 + e_dac_coarse_out(w) ** 2 + e_r1_out(w) ** 2
        + e_r2_out(w) ** 2 + e_r3_out(w) ** 2 + e_r14_out(w) ** 2 + e_ada4084_coarse_out(w) ** 2 + e_ada4084_fine_out(
            w) ** 2)

    # -----------Bruit total------------#
    en = lambda w: np.sqrt((G_RHF200 * (e_in_rhf200(w) + er15)) ** 2 + er16a ** 2 + er16b ** 2 + E_RHF200 ** 2)

    # -----------Impédance de sortie------------#
    zout = R16a + R16b

    return en, zout


def plot_loglog_with_ticks(x, y, xfigsize=10, yfigsize=6):
    """
    Parameters
    ----------
    x : see matplotlib loglog
        x data.
    y : see matplotlib loglog
        y data.
    xfigsize : float
        see matplotlib.pyplot.figure.(size in inches)
    yfigsize : float
        see matplotlib.pyplot.figure.(size in inches)

    Returns
    -------
    fig : TYPE
        DESCRIPTION.
    ax : TYPE
        DESCRIPTION.

    """
    fig, ax = plt.subplots(figsize=(xfigsize, yfigsize))
    ax.loglog(x, y)
    # grille principale
    ax.grid(which='major', linestyle='-', linewidth=0.6, alpha=0.8)
    # grille secondaire (optionnelle, plus discrète)
    ax.grid(which='minor', linestyle=':', linewidth=0.4, alpha=0.6)
    # --- FORCER les graduations secondaires (minor ticks) ---
    # Major ticks : puissances de 10 (10^n)
    ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
    # Minor ticks : entre 2,3,...,9 * 10^n
    # Utiliser subs = np.arange(2,10) (valeurs 2..9). Certains anciens matplotlib
    # acceptent subs as list/tuple; ici on passe un ndarray.
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    # Ajuster l'apparence des ticks
    ax.tick_params(which='major', length=8, width=1)
    ax.tick_params(which='minor', length=4, width=0.8, color='gray')

    return fig, ax


# --------------Classes----------------#


# --------------Exemples----------------#
if __name__ == "__main__":

    do_plots = False

    F_ANALOG = 100e9
    F_REF = 125e6
    TMUX = 34

    # Température de calcul
    TC = 20  # °C

    # Fréquence d'échantillonage et nom de fichier de sauvegarde
    # N_NUM is the number of points after aliasing
    FS, N_NUM, data_filename_ext = F_REF, 400, '_Fref.txt'
    # FS, N_NUM, data_filename_ext = F_REF / (20), 200, '_Frow.txt'
    # FS, N_NUM, data_filename_ext = F_REF / (20 * TMUX), 50, '_Fframe.txt'

    # Fréquences de calcul
    f = np.logspace(0, 10, 1000)  # 1 Hz à 10 GHz
    w = 2 * np.pi * f

    # Cas ADC avec 100 Ohms en entrée ###################################################################################
    data_filename = 'mod-erro-only' + data_filename_ext
    x_lbl = "Fréquence (Hz)"
    y_lbl = "Densité spectrale de bruit (V/√Hz)"

    zin = 100
    ein = lambda w: thermal_noise(zin, celsius2kelvin(TC))  # Le bruit en entrée doit être une fonction de w
    en, eq = adc_chain_noise_analog(TC, zin, ein, FS)
    tit = "Bruit analogique en entrée ADC"
    tit_ext = " après échantillonnage"

    if do_plots:
        _, ax = plot_loglog_with_ticks(f, en(w))
        ax.set_title(tit)
        ax.set_xlabel(x_lbl)
        ax.set_ylabel(y_lbl)
        plt.tight_layout()
        plt.show()

    # Repliement et bruit de quantification
    f_num, e_fold = aliasing(en, F_ANALOG, FS, N_NUM)
    e_aliased_total = np.sqrt(np.ones_like(f_num) * eq ** 2 + e_fold ** 2)

    if do_plots:
        _, ax = plot_loglog_with_ticks(f_num, e_aliased_total)
        ax.set_title(tit + tit_ext)
        ax.set_xlabel(x_lbl)
        ax.set_ylabel(y_lbl)
        plt.tight_layout()
        plt.show()

    export_2vectors_to_file(f_num, e_aliased_total, data_filename)

    # Cas de la chaine FB awaxe rebouclée sur ADC #######################################################################
    data_filename = 'mod-erro-fdbk-awaxe' + data_filename_ext
    ein, zin = awaxe_fb_chain_noise(TC)
    en, eq = adc_chain_noise_analog(TC, zin, ein, FS)
    tit = "Bruit analogique en entrée ADC - FB awaxe rebouclé"

    if do_plots:
        _, ax = plot_loglog_with_ticks(f, en(w))
        ax.set_title(tit)
        ax.set_xlabel(x_lbl)
        ax.set_ylabel(y_lbl)
        plt.tight_layout()
        plt.show()

    # Repliement et bruit de quantification
    f_num, e_fold = aliasing(en, F_ANALOG, FS, N_NUM)
    e_aliased_total = np.sqrt(np.ones_like(f_num) * eq ** 2 + e_fold ** 2)

    if do_plots:
        _, ax = plot_loglog_with_ticks(f_num, e_aliased_total)
        ax.set_title(tit + tit_ext)
        ax.set_xlabel(x_lbl)
        ax.set_ylabel(y_lbl)
        plt.tight_layout()
        plt.show()

    export_2vectors_to_file(f_num, e_aliased_total, data_filename, verbose=True)

    # Cas de la chaine FB RHF200 rebouclée sur ADC ######################################################################
    data_filename = 'mod-erro-fdbk-rhf200' + data_filename_ext
    ein, zin = rhf200_fb_chain_noise(TC)
    en, eq = adc_chain_noise_analog(TC, zin, ein, FS)
    tit = "Bruit analogique en entrée ADC - FB RHF200 rebouclé"

    if do_plots:
        _, ax = plot_loglog_with_ticks(f, en(w))
        ax.set_title(tit)
        ax.set_xlabel(x_lbl)
        ax.set_ylabel(y_lbl)
        plt.tight_layout()
        plt.show()

    # Repliement et bruit de quantification
    f_num, e_fold = aliasing(en, F_ANALOG, FS, N_NUM)
    e_aliased_total = np.sqrt(np.ones_like(f_num) * eq ** 2 + e_fold ** 2)

    if do_plots:
        _, ax = plot_loglog_with_ticks(f_num, e_aliased_total)
        ax.set_title(tit + tit_ext)
        ax.set_xlabel(x_lbl)
        ax.set_ylabel(y_lbl)
        plt.tight_layout()
        plt.show()

    export_2vectors_to_file(f_num, e_aliased_total, data_filename, verbose=True)

    # Cas de l'OFCO avec DAC coarse réglé à 0 ########################################################################
    ofco_coarse = 0
    data_filename = 'mod-erro-ofco_{0:}'.format(ofco_coarse) + data_filename_ext
    ein, zin = ofco_chain_noise(TC, ofco_coarse, 0, 0)
    en, eq = adc_chain_noise_analog(TC, zin, ein, FS)
    tit = "Bruit analogique en entrée ADC - OFCO rebouclé"

    if do_plots:
        _, ax = plot_loglog_with_ticks(f, en(w))
        ax.set_title(tit)
        ax.set_xlabel(x_lbl)
        ax.set_ylabel(y_lbl)
        plt.tight_layout()
        plt.show()

    # Repliement et bruit de quantification
    f_num, e_fold = aliasing(en, F_ANALOG, FS, N_NUM)
    e_aliased_total = np.sqrt(np.ones_like(f_num) * eq ** 2 + e_fold ** 2)

    if do_plots:
        _, ax = plot_loglog_with_ticks(f_num, e_aliased_total)
        ax.set_title(tit + tit_ext)
        ax.set_xlabel(x_lbl)
        ax.set_ylabel(y_lbl)
        plt.tight_layout()
        plt.show()

    export_2vectors_to_file(f_num, e_aliased_total, data_filename, verbose=True)

    # Cas de l'OFCO avec DAC coarse réglé à 1023 ########################################################################
    ofco_coarse = 1023
    data_filename = 'mod-erro-ofco_{0:}'.format(ofco_coarse) + data_filename_ext
    ein, zin = ofco_chain_noise(TC, ofco_coarse, 0, 0)
    en, eq = adc_chain_noise_analog(TC, zin, ein, FS)
    tit = "Bruit analogique en entrée ADC - OFCO rebouclé"

    if do_plots:
        _, ax = plot_loglog_with_ticks(f, en(w))
        ax.set_title(tit)
        ax.set_xlabel(x_lbl)
        ax.set_ylabel(y_lbl)
        plt.tight_layout()
        plt.show()

    # Repliement et bruit de quantification
    f_num, e_fold = aliasing(en, F_ANALOG, FS, N_NUM)
    e_aliased_total = np.sqrt(np.ones_like(f_num) * eq ** 2 + e_fold ** 2)

    if do_plots:
        _, ax = plot_loglog_with_ticks(f_num, e_aliased_total)
        ax.set_title(tit + tit_ext)
        ax.set_xlabel(x_lbl)
        ax.set_ylabel(y_lbl)
        plt.tight_layout()
        plt.show()

    export_2vectors_to_file(f_num, e_aliased_total, data_filename, verbose=True)
