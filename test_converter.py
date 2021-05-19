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

import converter

###############################################################################
def test_switch_bin2hexa_digit():
    assert converter.switch_bin2hexa_digit("0000") == "0"
    assert converter.switch_bin2hexa_digit("0001") == "1"
    assert converter.switch_bin2hexa_digit("0010") == "2"
    assert converter.switch_bin2hexa_digit("0011") == "3"
    assert converter.switch_bin2hexa_digit("0100") == "4"
    assert converter.switch_bin2hexa_digit("0101") == "5"
    assert converter.switch_bin2hexa_digit("0110") == "6"
    assert converter.switch_bin2hexa_digit("0111") == "7"
    assert converter.switch_bin2hexa_digit("1000") == "8"
    assert converter.switch_bin2hexa_digit("1001") == "9"
    assert converter.switch_bin2hexa_digit("1010") == "A"
    assert converter.switch_bin2hexa_digit("1011") == "B"
    assert converter.switch_bin2hexa_digit("1100") == "C"
    assert converter.switch_bin2hexa_digit("1101") == "D"
    assert converter.switch_bin2hexa_digit("1110") == "E"
    assert converter.switch_bin2hexa_digit("1111") == "F"
    assert converter.switch_bin2hexa_digit("10000") == "10"

###############################################################################
def test_twos_comp():
    assert converter.twos_comp(0,4) == 0
    assert converter.twos_comp(4,4) == 4
    assert converter.twos_comp(4,3) == -4
    assert converter.twos_comp(15,4) == -1
    assert converter.twos_comp(7,4) == 7
    
###############################################################################
def test_twos_comp_bin2int():
    assert converter.twos_comp_bin2int("000000") == 0
    assert converter.twos_comp_bin2int("111111") == -1
    assert converter.twos_comp_bin2int("011111") == 31
    assert converter.twos_comp_bin2int("100000") == -32
    assert converter.twos_comp_bin2int("1") == -1

###############################################################################
def test_twos_comp_hex2int():
    assert converter.twos_comp_hex2int("0") == 0
    assert converter.twos_comp_hex2int("FF") == -1
    assert converter.twos_comp_hex2int("1F") == 31
    assert converter.twos_comp_hex2int("A0") == -96
    assert converter.twos_comp_hex2int("1") == 1

###############################################################################
def test_switch_natbin2dec():
    assert converter.switch_natbin2dec("0000") == 0
    assert converter.switch_natbin2dec("0001") == 1
    assert converter.switch_natbin2dec("0010") == 2
    assert converter.switch_natbin2dec("0011") == 3
    assert converter.switch_natbin2dec("0100") == 4
    assert converter.switch_natbin2dec("0101") == 5
    assert converter.switch_natbin2dec("0110") == 6
    assert converter.switch_natbin2dec("0111") == 7
    assert converter.switch_natbin2dec("1000") == 8
    assert converter.switch_natbin2dec("1001") == 9
    assert converter.switch_natbin2dec("1010") == 10
    assert converter.switch_natbin2dec("1011") == 11
    assert converter.switch_natbin2dec("1100") == 12
    assert converter.switch_natbin2dec("1101") == 13
    assert converter.switch_natbin2dec("1110") == 14
    assert converter.switch_natbin2dec("1111") == 15

###############################################################################
def test_dec2cad():
    assert converter.dec2cad(7, 8) == '00000111'
    assert converter.dec2cad(-7, 8) == '11111001'

###############################################################################
def test_dec2natbin():
    assert converter.dec2natbin(7, 8) == '00000111'
    assert converter.dec2natbin(0, 8) == '00000000'
    assert converter.dec2natbin(3, 2) == '11'
    
###############################################################################
def test_signed_hexa_2_dec():
    assert converter.signed_hexa_2_dec('FFFE',16) == -2
    assert converter.signed_hexa_2_dec('7FFF',16) == 32767
    assert converter.signed_hexa_2_dec('7F',8) == 127
    assert converter.signed_hexa_2_dec('FF',8) == -1

###############################################################################
def test_dec_to_signed_hexa():
    assert converter.dec_to_signed_hexa(2,8) == '2'
    assert converter.dec_to_signed_hexa(2,16) == '2'
    assert converter.dec_to_signed_hexa(-2,16) == 'fffe'
    assert converter.dec_to_signed_hexa(-32768,16) == '8000'

###############################################################################
