#!/usr/bin/env python
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
#  laurent.ravera@irap.omp.eu
#  converter.py
#

###############################################################################
def switch_bin2hexa_digit(digit):
    r"""
    This function translates a binary string to an hexadecimal string
    """
    return(hex(int(digit,2))[2:].upper())

###############################################################################
def twos_comp(val, bits):
    r"""
    This function computes the 2's complement of an int value
    """
    if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
        val = val - (1 << bits)        # compute negative value
    return val                         # return positive value as is

###############################################################################
def twos_comp_bin2int(binary_string):
    r"""
    This function converts a two's complement binary value (string) to an integer
    """
    base=2
    return(twos_comp(int(binary_string,base), len(binary_string)))

###############################################################################
def twos_comp_hex2int(hex_string):
    r"""
    This function converts a two's complement hex value (string) to an integer
    """
    base=16
    return(twos_comp(int(hex_string,base), 4*len(hex_string)))

###############################################################################
def switch_bin2hexa(bin):
    r"""
    This function makes the conversion from binary ascii to hexadecimal ascci.
    """
    nbits=len(bin)
    offset=0
    hexa=''
    while offset<nbits:
        hexa=hexa+switch_bin2hexa_digit(bin[offset:offset+4])
        offset+=4
    return(hexa)

###############################################################################
def switch_natbin2dec(bin):
    r"""
    This function makes the conversion from natural binary ascii to decimal value.
    """
    return(int(bin,2))

###############################################################################
def dec2cad(dec, n):
    r"""
    This function computes the 2s complement binary value of an integer.
    
    Parameters
    ----------
    dec : number
        The decimal input integer to be converted
    n : number
        The number of bits of the output string

    Output
    ------
    bin_str : string
        The binary value in a string format

    Raises
    ------
    ValueError
        When n is too small to do the conversion

    See Also
    --------
    num : string to number conversion

    Examples
    --------
    >>> dec2cad(7, 8)
    '00000111'
    >>> dec2cad(-7, 8)
    '11111001'
    >>> dec2cad(7, 2)
    ValueError: Requested size is too small for this value

    """
    bin_str = ''
    if dec >= 0:
        while n > 0:
            bin_str = str(dec % 2) + bin_str
            dec >>= 1
            n += -1
    else:
        dec = abs(dec + 1)
        while n > 0:
            rev_str = ('1', '0')
            bin_str = rev_str[dec % 2] + bin_str
            dec >>= 1
            n += -1
    if dec > 0:
        raise ValueError('Requested size is too small for this value')
    return bin_str

###############################################################################
def dec2natbin(dec, n):
    r"""
    This function computes the natural binary value of a positive integer.
    
    Parameters
    ----------
    dec : number
        The decimal input integer to be converted
    n : number
        The number of bits of the output string

    Output
    ------
    bin_str : string
        The binary value in a string format

    Raises
    ------
    ValueError
        When n is too small to do the conversion

    See Also
    --------
    num : string to number conversion

    """
    bin_str = bin(dec)[2:]
    l=len(bin_str)
    if l<=n:
        bin_str = '0'*(n-l) + bin_str
    else:
        raise ValueError('Requested size is too small for this value')
    return bin_str

###############################################################################
def signed_hexa_2_dec(hexstr,nbits):
    r"""
    This function computes the decimal value of a 2's complement signed hexa string.
    
    Parameters
    ----------
    hexstr : string
        The hexadecimal value to be converted
    nbits : number
        The number of bits of the data

    Output
    ------
    value : integer
        The converted value

    Examples
    --------
    >>> signed_hexa_2_dec('FFFE',16)
    -2
    >>> signed_hexa_2_dec('7FFF',16)
    32767
    >>> signed_hexa_2_dec('7F',8)
    127
    >>> signed_hexa_2_dec('FF',8)
    -1
    """
    value = int(hexstr,16)
    if value & (1 << (nbits-1)):
        value -= 1 << nbits
    return value

###############################################################################
def dec_to_signed_hexa(signeddec,nbits):
    r"""
    This function computes the hexadecimal string of a signed decimal value.
    
    Parameters
    ----------
    signeddec : decimal
        The decimal value to be converted
    nbits : number
        The number of bits of the data

    Output
    ------
    value : string
        The converted value

    Examples
    --------
    >>> dec_to_signed_hexa(2,8)                                                           
    '0x2'
    >>> dec_to_signed_hexa(2,16)                                                          
    '0x2'
    >>> dec_to_signed_hexa(-2,16)                                                         
    '0xfffe'
    """

    return(hex(signeddec & (2**nbits-1))[2:])

###############################################################################
