
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

    xmin=20
    xmax=len(data_tes[0])-1
    x_zoom_min=0
    x_zoom_max=20
    extension='_pulse' # type of plots

    plotdirname=os.path.join(os.path.join(dirname,'..'),'PLOTS')
    plotdirname_already_exists=check_mk_dir(plotdirname)

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
                    vect_x=np.arange(xmax-xmin)

                    # Plotting TES data
                    ax1 = fig.add_subplot(3, 3, 1)
                    ax1.step(vect_x, data_tes[pixel][xmin:xmax], where='mid')
                    ax1.plot(vect_x, data_tes[pixel][xmin:xmax], '.k')
                    ax1.set_title('TES input (pix {0:2d})'.format(pixel))
                    ax1.set_xlabel(xtitle)
                    ax1.set_xticks(major_ticks)
                    ax1.set_xticks(minor_ticks, minor=True)
                    ax1.grid(which='minor', alpha=0.2)
                    ax1.grid(which='major', alpha=0.5)
                    #ax1.grid(color='k', linestyle=':', linewidth=0.5)
                    ax1.set_xlim(xmin,xmax)

                    # Plotting error data
                    ax2 = fig.add_subplot(3, 3, 2)
                    ax2.step(vect_x, data_error[pixel][xmin:xmax], where='mid')
                    ax2.plot(vect_x, data_error[pixel][xmin:xmax], '.k')
                    ax2.set_title('Error (pix {0:2d})'.format(pixel))
                    ax2.set_xlabel(xtitle)
                    ax2.set_xticks(major_ticks)
                    ax2.set_xticks(minor_ticks, minor=True)
                    ax2.grid(which='minor', alpha=0.2)
                    ax2.grid(which='major', alpha=0.5)
                    ax2.set_xlim(xmin,xmax)

                    # Plotting feedback data
                    ax3 = fig.add_subplot(3, 3, 3)
                    ax3.step(vect_x, data_feedback[pixel][xmin:xmax], where='mid')
                    ax3.plot(vect_x, data_feedback[pixel][xmin:xmax], '.k')
                    ax3.set_title('Feedback (pix {0:2d})'.format(pixel))
                    ax3.set_xlabel(xtitle)
                    ax3.set_xticks(major_ticks)
                    ax3.set_xticks(minor_ticks, minor=True)
                    ax3.grid(which='minor', alpha=0.2)
                    ax3.grid(which='major', alpha=0.5)
                    ax3.set_xlim(xmin,xmax)

                    # Plotting science data
                    ax4 = fig.add_subplot(3, 3, 4)
                    ax4.step(vect_x, data_science[pixel][xmin:xmax], where='mid')
                    ax4.plot(vect_x, data_science[pixel][xmin:xmax], '.k')
                    ax4.set_title('Science output (pix {0:2d})'.format(pixel))
                    ax4.set_xlabel(xtitle)
                    ax4.set_xticks(major_ticks)
                    ax4.set_xticks(minor_ticks, minor=True)
                    ax4.grid(which='minor', alpha=0.2)
                    ax4.grid(which='major', alpha=0.5)
                    ax4.set_xlim(xmin,xmax)

                    # Plotting FB + error data
                    ax5 = fig.add_subplot(3, 3, 5)
                    ax5.step(vect_x, data_feedback[pixel][xmin:xmax]+data_error[pixel][xmin:xmax], where='mid')
                    ax5.plot(vect_x, data_feedback[pixel][xmin:xmax]+data_error[pixel][xmin:xmax], '.k')
                    ax5.set_title('Feedback+error (pix {0:2d})'.format(pixel))
                    ax5.set_xlabel(xtitle)
                    ax5.set_xticks(major_ticks)
                    ax5.set_xticks(minor_ticks, minor=True)
                    ax5.grid(which='minor', alpha=0.2)
                    ax5.grid(which='major', alpha=0.5)
                    ax5.set_xlim(xmin,xmax)

                    # Checking science: Plotting Science - (FB + error data) 
                    ax6 = fig.add_subplot(3, 3, 6)
                    ax6.step(vect_x, 100*(data_science[pixel][xmin+1:xmax+1]-\
                        (data_feedback[pixel][xmin:xmax]+data_error[pixel][xmin:xmax]))/ \
                        data_science[pixel][xmin:xmax].max(), where='mid')
                    ax6.plot(vect_x, 100*(data_science[pixel][xmin+1:xmax+1]-\
                        (data_feedback[pixel][xmin:xmax]+data_error[pixel][xmin:xmax]))/ \
                        data_science[pixel][xmin:xmax].max(), '.k')
                    ax6.set_title('Science-(Feedback+error) (%, pix {0:2d})'.format(pixel))
                    ax6.set_xlabel(xtitle)
                    ax6.set_xticks(major_ticks)
                    ax6.set_xticks(minor_ticks, minor=True)
                    ax6.grid(which='minor', alpha=0.2)
                    ax6.grid(which='major', alpha=0.5)
                    ax6.set_xlim(xmin,xmax)

                    # Plotting the comparison science vs tes
                    ax7 = fig.add_subplot(3, 3, 7)
                    ax7.step(vect_x, data_tes[pixel][xmin:xmax]/data_tes[pixel][-1], where='mid')
                    ax7.plot(vect_x, data_tes[pixel][xmin:xmax]/data_tes[pixel][-1], '.k', label='TES input (norm.)')
                    ax7.step(vect_x, data_science[pixel][xmin:xmax]/data_science[pixel][-1], where='mid')
                    ax7.plot(vect_x, data_science[pixel][xmin:xmax]/data_science[pixel][-1], '.r', label='Science output (norm.)')
                    ax7.set_title('TES input & science output (pix {0:2d})'.format(pixel))
                    ax7.legend(loc='best')
                    ax7.set_xlabel(xtitle)
                    ax7.set_xticks(major_ticks)
                    ax7.set_xticks(minor_ticks, minor=True)
                    ax7.grid(which='minor', alpha=0.2)
                    ax7.grid(which='major', alpha=0.5)
                    ax7.set_xlim(xmin,xmax)

                    # Plotting the comparison science vs tes
                    ax7 = fig.add_subplot(3, 3, 8)
                    diff_pc=100.*(data_tes[pixel][xmin:xmax-1]/data_tes[pixel][-1]-data_science[pixel][xmin+1:xmax]/data_science[pixel][-1])
                    ax7.step(vect_x[:-1], diff_pc, where='mid')
                    ax7.plot(vect_x[:-1], diff_pc, '.k')
                    ax7.set_title('norm tes - norm shifted science (%, pix {0:2d})'.format(pixel))
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

                    xmin=x_zoom_min
                    xmax=x_zoom_max
                    major_ticks = np.arange(0, x_zoom_max, 4)
                    minor_ticks = np.arange(0, x_zoom_max, 1)
                    extension='_start' # type of plots (zoom or not)

# -----------------------------------------------------------------------

science_shift=1

name = '../DATA'
mosaic(name, science_shift, overwrite=True, mux_factor=34)
