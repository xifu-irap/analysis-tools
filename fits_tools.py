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
#  fits_tools.py
#

#---------------------------------------------------------------------------------    
from astropy.io import fits


#---------------------------------------------------------------------------------    
VERSION = "FAKE-1"
ADC_DUMP_EXTENSION = "ADC_DUMP"

#---------------------------------------------------------------------------------    
def write_adc_dump_fits_file(file_name, dump):
    """This function saves an ADC dump into a fits file. 

    Args:
        file_name (string): Name of the fits file.
        dump (np_array): dump dataset
    """
    # Empty primary HDU with header only    
    hdr = fits.Header()
    hdr['BITSMP'] = (16, 'Number of bits per sample')
    empty_primary = fits.PrimaryHDU(header=hdr)
        
    # Secondary HDU (contains the fake dump)
    col = fits.Column(name='AdcValue', format='I', array=dump)
    cols = fits.ColDefs([col])
    hdu = fits.BinTableHDU.from_columns(cols)
    
    hdu.header['EXTNAME']=ADC_DUMP_EXTENSION
    hdu.header['VERSION']=VERSION
    hdu.header['INSTRUME']='XIFU'
    hdu.header['ORIGIN']='IRAP'
    hdu.header['FILETYPE']='AdcDump'
    
    hdul = fits.HDUList([empty_primary, hdu])
    hdul.writeto(file_name, overwrite=True)


#---------------------------------------------------------------------------------    
def read_adc_dump_fits_file(file_name, verbose=False):
    """This function reads the content of an ADC dump fits file.
    
    Args:
        file_name (string): Name of the fits file.
        verbose (boolean): verbosity level. Defaults to False.    
    """
    hdul = fits.open(file_name)
    hdul.info()
    header = hdul['ADC_DUMP'].header
    if header['EXTNAME'] != ADC_DUMP_EXTENSION:
        raise NameError("Error ! This is not an ADC dump fits file. hdu extension is: ", header['EXTNAME'])
    if verbose:
        print("-------------------------------------------------")
        print(" Reading ADC dump from file :", file_name)
        print(" File format version:", header['VERSION'])
        print("-------------------------------------------------")
    
    data = hdul['ADC_DUMP'].data
    
    return(data['AdcValue'])
    

#---------------------------------------------------------------------------------    
