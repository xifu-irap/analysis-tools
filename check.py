
# -----------------------------------------------------------------------
# Imports
import os
import numpy as np
import matplotlib
matplotlib.use('agg') # to avoid errors when used from Questa
import matplotlib.pyplot as plt
import converter

# -----------------------------------------------------------------------
def check_mk_dir(dirname):
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
    dir_already_exists=True
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
        dir_already_exists=False
    return(dir_already_exists)


# -----------------------------------------------------------------------
def get_data_from_ascii_file(filename):
    r"""
        This function gets testbench data from a ascii file

        Parameters
        ----------
        filename : string
        The name of the data file (this includes the path)

        Returns
        -------
        data : array
        The data in an array of integers.
        
        """
    if os.path.isfile(filename):
        data=[]
        f=open(filename, 'r')
        for line in f:
            if line[0] != '#' and line[0] != ' ':
                data.append(converter.twos_comp_bin2int(line.strip()))
        f.close()
    else:
        raise ValueError("Error, file ", filename, "does not exists.")
    return(np.array(data))

# -----------------------------------------------------------------------
def demux(column_data, mux_factor=34, offset=0):
    r"""
        This function demultiplexes data. i.e. it extracts 1 over mux_factor values.

        Parameters
        ----------
        column_data : array
        array of column values

        mux_factor : int
        multiplexing factor (default=34)

        offset : int
        offset that can be applied to correct for processing time (default=0)

        Returns
        -------
        pix_data : array
        array containing original data re-aranged in mux_factor lines 
        and n_frames columns
        
        """
    nframes=int((len(column_data)-offset)/mux_factor)
    column_data=column_data[offset:nframes*mux_factor+offset]
    return(column_data.reshape((mux_factor,nframes),order='F'))

# -----------------------------------------------------------------------
def mosaic(dirname, science_shift, overwrite=False, mux_factor=34):
    r"""
        This function plots the data for a specific pixel.

        Parameters
        ----------
        dirname : string
        path to the data

        science_shift : int
        shift to be applied to align the science data with the TES data
        (compensates the processing time)

        overwrite : boolean
        if True the previous plots are overwritten (default=False)

        mux_factor : int
        multiplexing factor (default=34)

        Returns
        -------
        Nothing
        
        """

    # getting data from files
    data_tes=demux(get_data_from_ascii_file(os.path.join(dirname, 'vtes.log')))
    data_error=demux(get_data_from_ascii_file(os.path.join(dirname, 'error.log')))
    data_feedback=demux(get_data_from_ascii_file(os.path.join(dirname, 'feedback.log')))
    data_science=demux(get_data_from_ascii_file(os.path.join(dirname, 'science.log')), offset=science_shift)

    xmin=0
    xmax=len(data_tes[0])
    vect_x=np.arange(xmax)
    extension='_pulse' # type of plots

    plotdirname=os.path.join(os.path.join(dirname,'..'),'PLOTS')
    plotdirname_already_exists=check_mk_dir(plotdirname)

    x_zoom_max=20
    major_ticks = np.arange(0, xmax, 200)
    minor_ticks = np.arange(0, xmax, 50)
 
    if not plotdirname_already_exists or overwrite:

        for pixel in range(mux_factor):

            if data_tes[pixel].max()>0 and abs(data_error[pixel]).max()>0 \
                and data_feedback[pixel].max()>0 and data_science[pixel].max()>0:

                for plot_index in range(2):

                    plotfilename=os.path.join(plotdirname, "pix-{0:1d}".format(pixel).zfill(2)+extension+".png")

                    fig = plt.figure(figsize=(12, 12))
                    xtitle = "Samples"

                    # Plotting TES data
                    ax1 = fig.add_subplot(3, 3, 1)
                    ax1.step(vect_x, data_tes[pixel], where='mid')
                    ax1.plot(vect_x, data_tes[pixel], '.k')
                    ax1.set_title('TES data (pix {0:2d})'.format(pixel))
                    ax1.set_xlabel(xtitle)
                    ax1.set_xticks(major_ticks)
                    ax1.set_xticks(minor_ticks, minor=True)
                    ax1.grid(which='minor', alpha=0.2)
                    ax1.grid(which='major', alpha=0.5)
                    #ax1.grid(color='k', linestyle=':', linewidth=0.5)
                    ax1.set_xlim(xmin,xmax)

                    # Plotting error data
                    ax2 = fig.add_subplot(3, 3, 2)
                    ax2.step(vect_x, data_error[pixel], where='mid')
                    ax2.plot(vect_x, data_error[pixel], '.k')
                    ax2.set_title('Error data (pix {0:2d})'.format(pixel))
                    ax2.set_xlabel(xtitle)
                    ax2.set_xticks(major_ticks)
                    ax2.set_xticks(minor_ticks, minor=True)
                    ax2.grid(which='minor', alpha=0.2)
                    ax2.grid(which='major', alpha=0.5)
                    ax2.set_xlim(xmin,xmax)

                    # Plotting feedback data
                    ax3 = fig.add_subplot(3, 3, 4)
                    ax3.step(vect_x, data_feedback[pixel], where='mid')
                    ax3.plot(vect_x, data_feedback[pixel], '.k')
                    ax3.set_title('Feedback data (pix {0:2d})'.format(pixel))
                    ax3.set_xlabel(xtitle)
                    ax3.set_xticks(major_ticks)
                    ax3.set_xticks(minor_ticks, minor=True)
                    ax3.grid(which='minor', alpha=0.2)
                    ax3.grid(which='major', alpha=0.5)
                    ax3.set_xlim(xmin,xmax)

                    # Plotting science data
                    ax4 = fig.add_subplot(3, 3, 5)
                    ax4.step(vect_x, data_science[pixel], where='mid')
                    ax4.plot(vect_x, data_science[pixel], '.k')
                    ax4.set_title('Science data (pix {0:2d})'.format(pixel))
                    ax4.set_xlabel(xtitle)
                    ax4.set_xticks(major_ticks)
                    ax4.set_xticks(minor_ticks, minor=True)
                    ax4.grid(which='minor', alpha=0.2)
                    ax4.grid(which='major', alpha=0.5)
                    ax4.set_xlim(xmin,xmax)

                    # Plotting FB + error data
                    ax5 = fig.add_subplot(3, 3, 6)
                    ax5.step(vect_x, data_feedback[pixel]+data_error[pixel], where='mid')
                    ax5.plot(vect_x, data_feedback[pixel]+data_error[pixel], '.k')
                    ax5.set_title('Feedback+error (pix {0:2d})'.format(pixel))
                    ax5.set_xlabel(xtitle)
                    ax5.set_xticks(major_ticks)
                    ax5.set_xticks(minor_ticks, minor=True)
                    ax5.grid(which='minor', alpha=0.2)
                    ax5.grid(which='major', alpha=0.5)
                    ax5.set_xlim(xmin,xmax)

                    # Plotting the comparison science vs tes
                    ax6 = fig.add_subplot(3, 3, 7)
                    ax6.step(vect_x, data_tes[pixel]/data_tes[pixel][-1], where='mid')
                    ax6.plot(vect_x, data_tes[pixel]/data_tes[pixel][-1], '.k')
                    ax6.step(vect_x, data_science[pixel]/data_science[pixel][-1], where='mid')
                    ax6.plot(vect_x, data_science[pixel]/data_science[pixel][-1], '.r')
                    ax6.set_title('norm tes & norm science (pix {0:2d})'.format(pixel))
                    ax6.set_xlabel(xtitle)
                    ax6.set_xticks(major_ticks)
                    ax6.set_xticks(minor_ticks, minor=True)
                    ax6.grid(which='minor', alpha=0.2)
                    ax6.grid(which='major', alpha=0.5)
                    ax6.set_xlim(xmin,xmax)

                    # Plotting the comparison science vs tes
                    ax7 = fig.add_subplot(3, 3, 8)
                    ax7.step(vect_x[1:], data_tes[pixel][:-1]/data_tes[pixel][-1]-data_science[pixel][1:]/data_science[pixel][-1], where='mid')
                    ax7.plot(vect_x[1:], data_tes[pixel][:-1]/data_tes[pixel][-1]-data_science[pixel][1:]/data_science[pixel][-1], '.k')
                    ax7.set_title('norm tes - norm shifted science (pix {0:2d})'.format(pixel))
                    ax7.set_xlabel(xtitle)
                    ax7.set_xticks(major_ticks)
                    ax7.set_xticks(minor_ticks, minor=True)
                    ax7.grid(which='minor', alpha=0.2)
                    ax7.grid(which='major', alpha=0.5)
                    ax7.set_xlim(xmin,xmax)

                    fig.tight_layout()
                    #plt.show()
                    plt.savefig(plotfilename, bbox_inches='tight')
                    plt.close()

                    xmin=0
                    xmax=x_zoom_max
                    major_ticks = np.arange(0, x_zoom_max, 4)
                    minor_ticks = np.arange(0, x_zoom_max, 1)
                    extension='_start' # type of plots (zoom or not)



# -----------------------------------------------------------------------

science_shift=1

name = '../DATA'
#name='/Users/laurent/MyCore_xifu/20_DRE/70_Architecture/00_Firmware/TDM_simulations/2020-10-06_12-28-13/DATA'
mosaic(name, science_shift, overwrite=True, mux_factor=34)
