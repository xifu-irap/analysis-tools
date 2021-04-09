#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, csv
import numpy as np

# -----------------------------------------------------------------------
def checkdir(dirname):
    r"""
        This function checks if a directory exists. If not it creates it.

        Parameters:
        -----------
        dirname: String
        Name of the directory to be verified / created.

        Returns
        -------
        Nothing

        """
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    return()

# -----------------------------------------------------------------------
def get_csv(filename):
    r"""
        This function reads a dictionnary from a csv file.

        Parameters:
        -----------
        filename: string
        The name of the csv file

        Returns
        -------
        dictionnary: dictionnary

        """

    dictionnary={}

    if not os.path.exists(os.path.join(filename)):
        print("File "+filename+" not found.")
    else:
        with open(filename, newline='') as csvfile:
            dict_reader = csv.reader(csvfile, delimiter='=', quotechar='|')
            for row in dict_reader:
                try:    # for numbers
                    dictionnary[row[0]]=float(row[1].replace(',','.'))
                except: # for strings
                    dictionnary[row[0]]=row[1]
    return(dictionnary)

# -----------------------------------------------------------------------
def do_spectrum(x, npts, window="none"):
    r"""
        This function computes the spectrum of the input vector.
        If the input vector is long enough, several computations are averaged.
        A Blackman window is applied before the rfft.

        Parameters:
        -----------
        x: numpy array
        input vector

        npts: number
        Number of values to be used in the rfft.

        window: string
        indicates if a window shall be applied ("blackman" or "none") 
        (default is "none")

        Returns
        -------
        spectrum: numpy array
        computed spectrum.

        """
    from numpy.fft import rfft
    from scipy import signal

    if window=="blakman":
        w=signal.blackman(npts)
    else:
        w=np.ones(npts)

    if len(x)<npts:
        raise ValueError("Not enough values in input vector to compute spectra.")
    else:
        if is_even(npts):
            spectrum=np.zeros(int(npts/2)+1)
        else:
            spectrum=np.zeros(int((npts+1)/2))
        nslices=int(len(x)/npts)
        for slice in range(nslices):
            spectrum+=abs(rfft(x[slice*npts:(slice+1)*npts]*w))
    return(spectrum/nslices)

# -----------------------------------------------------------------------
def is_even(x):
    r"""
        This function checks if a number is even.

        Parameters:
        -----------
        x: number
        The value to be verified

        Returns
        -------
        result: boolean
        True if the input is even, False if the input is odd.

        """
    return(x%2==0)

# -----------------------------------------------------------------------
class configuration:
    def __init__(self, csv_file):
        self.config = get_csv(csv_file+".csv")

# -----------------------------------------------------------------------
