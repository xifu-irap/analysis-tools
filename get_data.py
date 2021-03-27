import numpy as np

# Some general text definitions
plotting_message="Plotting dump data..."
time_label="Time (ms)"
frequency_label="Frequency (Hz)"
spectral_power_density_label=r"Power spectral density (dB/$\sqrt{Hz}$)"

# -----------------------------------------------------------------------------
def dumptype_from_header(header_dumptype_field):
    r"""
        This function convertsbthe dumptype id into a human readable string

        Parameters
        ----------
        header_dumptype_field : unit16
        The dumptype field from the file header

        Returns
        -------
        dumptype_int : number
        Dumptype id.

        dumptype_str : string
        Description of the dumptype.
        """

    header24=int(header_dumptype_field/2**12)
    header23=int((header_dumptype_field-header24*2**12)/2**8)
    dumptype_int=int((header_dumptype_field-header24*2**12-header23*2**8)/2**4)
 
    if dumptype_int==0:
        dumptype_str = "DACFB1, SCIENCE"
    elif dumptype_int==1:
        dumptype_str = "ERROR5MHZ, SCIENCE"
    elif dumptype_int==2:
        dumptype_str = "DACFB2, SCIENCE"
    elif dumptype_int==4:
        dumptype_str = "DACFB1, DACFB2"        
    elif dumptype_int==5:
        dumptype_str = "ERROR40MHZ"        
    elif dumptype_int==6:
        dumptype_str = "DACFB1"        
    elif dumptype_int==8:
        dumptype_str = "SCIENCE PACKETIZED"        
    elif dumptype_int==9:
        dumptype_str = "ERROR5MHz, DACFB1"                
    elif dumptype_int==15 or dumptype_int==10:
        dumptype_str = "COUNTER32"                
    else:
        raise ValueError('Wrong dump type')

    return(dumptype_int, dumptype_str)

# -----------------------------------------------------------------------------
def read_dumptype(dumpfilename):
    r"""
        This function reads the dumptype of a dumpfile

        Parameters
        ----------
        dumpfilename : string
        The name of the dump file (with the path and the extension)

        Returns
        -------
        dumptype_int : number
        the dumptype id.
        
        """
    fdat=open(dumpfilename, 'rb')
    data=np.fromfile(fdat, dtype='<h')
    fdat.close()
        
    DADA=-9510      # 0xDADA interpreted as int16
    if data[0] != DADA:
        raise ValueError('Problem with file format!')
    dumptype_int, dumptype_str=dumptype_from_header(data[1].astype('uint16'))
    print("# Dumptype is:", dumptype_str)
    return(dumptype_int)

# -----------------------------------------------------------------------------
def read_dump(dumpfilename, config):
    r"""
        This function reads data from a DRE adc dump file, and returns 1 array
        (format <h).
        The file header (first 32-bit) word is removed from the data

        Parameters
        ----------
        dumpfilename : string
        The name of the dump file (with the path and the extension)

        config : dictionnary
        contains different informations such as path and directory names

        Returns
        -------
        data : array type
        contains the values of the file (format int16).

        dumptype_int : number
        Dumptype id.

        dumptype_str : string
        Description of the dumptype.

        """
    
    fdat=open(dumpfilename, 'rb')
    data=np.fromfile(fdat, dtype='<h')
    fdat.close()
        
    DADA=-9510      # 0xDADA interpreted as int16
    if data[0] != DADA:
        raise ValueError('Problem with file format!')

    dumptype_int, dumptype_str=dumptype_from_header(data[1].astype('uint16'))
    print("# Dumpfile name: ", dumpfilename)
    print("# Dumptype     : ", dumptype_str)

    # removing header
    data=data[2:]

    # decommutation of data according to dumptype
    ## 5MHz dump data
    if dumptype_int == 0 or dumptype_int == 1 or dumptype_int == 2 or dumptype_int == 4 or dumptype_int == 9:
        data = np.resize(data, (len(data)//2, 2))

    ## ADC dump data
    if dumptype_int == 5:
        # input values are in the wrong order
        data2 = np.resize(data, (len(data)//2, 2))
        data[0::2]=data2[:,1]
        data[1::2]=data2[:,0]

    ## science data
    if dumptype_int == 8:
        # merging 16-bit words in 32-bit words
        data = np.resize(data, (len(data)//2, 2)).astype('uint16')
        data = data[:, 1]*2.**16 + data[:, 0]

        # the first header has been removed by read_dump function. I need it here. I add a fake one.
        data = np.append(np.zeros(1), data)

        # reordering data
        npix = int(config["npix"])
        data = np.transpose(np.resize(data, (len(data)//(npix+2), npix+2)))

        # removing periodic headers
        data = data[1:,:]

    ## 32-bit counter data
    if dumptype_int == 15:
        data = np.resize(data, (len(data)//2, 2)).astype('uint16')
        data = data[:, 1]*2.**16 + data[:, 0]

    return(data, dumptype_int, dumptype_str)

# -----------------------------------------------------------------------------
