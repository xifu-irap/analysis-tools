from astropy.io import fits
import numpy as np
import os


def readEventRecords(fits_file):
    """Reads DEMUX science data (error or science) from a fits file.

    Parameters:
        fits_file (string): Name of the fits file (includes the path).

    Returns:
        firstDf (int): number of the Data Frame corresponding to the beginning of the fits file.

        data (array): event records. Each element corresponds to a record, it contains:
            - the data frame number of the beginning of the record
            - the data
    """
    with fits.open(fits_file) as hdul:
        # Lecture des données et de l'en-têtes
        firstDf = hdul[1].header['FIRST_DF']
        data = hdul[1].data
        return firstDf, data

# Chemin vers le fichier FITS
#file_path = '~laurent/Data/TestPlan21-perfo/20241121_eventRecords/pulses_20241121-164729.fits'
#firstDf, data = readEventRecords(file_path)
#print("First Data frame: ", firstDf)
#print("Data: ")
#print(data)

def readScienceFile(fits_file):
    """Reads DEMUX science data (error or science) from a fits file.

    Parameters:
        fits_file (string): Name of the fits file (includes the path).

    Returns:
        error (array): The error values (one value per pixel and per step)

        ctrl (array): The control words (one value per step)
    """
    with fits.open(fits_file) as hdul:

        # Checking the extension name 'pixels_data'
        if 'pixels_data' in hdul:
            data_ext = hdul['pixels_data']
        else:
            raise ValueError("'pixels_data' extension not found in the FITS file.")

        # Extraire les données en tant que tableau NumPy
        data = data_ext.data
        dataArray = np.array(data.tolist()).T


        # -1 parce que les dernières données sont à 0 (this is corrected in XIFU STUDIO)
        #return dataArray[1:,:-1], dataArray[0,:-1]
        return dataArray[1:,:], dataArray[0,:]


def readScienceDataOneCol(data_path, col_id, remove_dc=True):
    """Reads DEMUX science data for one column from a fits file.

    Parameters:
        data_path (string): Path to the data files.
        col_id (int): Column ID (0 to 1359)
        remove_dc (boolean): If True the dc of the data is removed (default is True)

    Returns:
        col_data (numpy array): The data
    """
    files = [f for f in os.listdir(data_path) \
             if os.path.isfile(os.path.join(data_path, f)) \
             and f[:6] == 'error_' and f[-6:] == '{0:}.fits'.format(col_id)]

    if len(files) == 0:
        raise ValueError('No fits file for column {0:}'.format((col_id)))

    file_name = files[0]
    file_name_with_path = os.path.join(data_path, file_name)

    print("    Reading ERROR data from ", file_name_with_path, ".... ", end = '')

    col_data, _ = readScienceFile(file_name_with_path)
    # Flattening the array to have data at Frow
    # Dividing by 4 because Error data are in S(16,2) format
    col_data = col_data.flatten('F') / 4

    # removing DC
    if remove_dc:
        col_data -= col_data.mean()
    print("Done")

    return col_data


def readDumpFile(fits_file):
    """Reads DEMUX dump data from a fits file.

    Parameters
        fits_file (sring): Name of the fits file (includes the path).

    Returns:
        dump (array): the dump data (4 x 1360 values)

        adc_error (array): the ADC error data (1360 bytes)
    """
    with fits.open(fits_file) as hdul:

        # Checking the extension name 'Dumps'
        if 'Dumps' in hdul:
            data_ext = hdul['Dumps']
        else:
            raise ValueError("'Dumps' extension not found in the FITS file.")

        # Extraire les données en tant que tableau NumPy
        data = data_ext.data
        dump = np.array(data.tolist())[0]

        # Return columns and errors
        return dump[0:-1,:], dump[-1,:]


def readScan(fits_file):
    """Reads DEMUX scan data from a fits file.

    Parameters
        fits_file (string): Name of the fits file (includes the path).

    Returns:
        xName (string): the name of the signal on the X axis

        ctrl (array): array with the control words

        xValues (array): array with the x values

        error (array): array with the error values (one value per pixel and per step)
    """
    with fits.open(fits_file) as hdul:

        # Checking the extension name 'pixels_data'
        if 'pixels_data' in hdul:
            data_ext = hdul['pixels_data']
        else:
            raise ValueError("'pixels_data' extension not found in the FITS file.")

        xName = data_ext.header['TTYPE2']

        # Extraire les données en tant que tableau NumPy
        data = data_ext.data
        scan = np.array(data.tolist())
        errors = scan[:,2:].T

        # Return xName, CTRL, xValues and error per pixels
        return xName, scan[:,0], scan[:,1], errors

def save_dump_txt(fits_file, col):
    """Saves the data of a column from a dump to a text file.

    Parameters
        fits_file (string): the name of the dump fits_file

        col (integer): the col id
    """
    data, _ = readDumpFile(fits_file)
    np.savetxt(fits_file[0:-5] + '_col_{0:}.txt'.format(col), data[col,:])

file = "/users/laurent/Data/TestPlan26-DMX2/debug_com_ADC/dump_20250523-153003.fits"
save_dump_txt(file, 1)
