
# imports
import numpy as np
import matplotlib.pyplot as plt

# custom imports
import constants as cst
import general_tools as gtl

#----------------------------------------------------------------
def get_tm_mode(value):
    """
    This function returns the name of the tm mode corresponding to a 
    given 8-bit control word

    Args:
        value (string): hexadecimal 8-bit control word expressed as a string

    Returns:
        string: The name of the tm mode
    """    
    return [key for key, val in cst.tm_dict.items() if val == value][0]


#----------------------------------------------------------------
class packet:
    """
    Class corresponding to a science TM packet
    A packet is defined with:
        A filename (includes the path)
        A time in ns
        An array witth the tm_mode
        A number of frame
        A multiplexing factor mux
        An array with the data containing one sub-array per column:
            the colomun 0 of size ((nb of frames, nb of pixels))
            the colomun 1 of size ((nb of frames, nb of pixels))
            the colomun 2 of size ((nb of frames, nb of pixels))
            the colomun 3 of size ((nb of frames, nb of pixels))        
    """
    def __init__(self, fullfilename, time, tm_mode, nframes, mux, col0, col1, col2, col3):
        self.fullfilename=fullfilename
        self.time=time
        self.tm_mode=[tm_mode]
        self.nframes=nframes
        self.mux=mux
        self.data=[col0, col1, col2, col3]

    #------------------------------------------------------------
    def print_pixel(self, pixel):
        """
        Method to print the packet content for a specific pixel

        Args:
            pixel (integer): pixel index whose data will be printed
        """
        if self.tm_mode != cst.sci_name:
            print("Impossible to do this print in this tm mode")            
        else:
            print("--------------")
            print("Time: ", self.time)
            for frame in range(self.nframes):
                print("Telemetry mode: ", self.tm_mode[frame])
            print("Nframes: ", self.nframes)
            print("Multiplexing factor: ", self.mux)
            for column in range(cst.nb_col):
                print("\nb_cols {0:d}, pixel {1:d}:".format(column, pixel))
                for frame in range(self.nframes):
                    print('{0:X} '.format(self.data[column][frame, pixel]), end='')
            print("\n")

    #------------------------------------------------------------
    def print_frame(self, frame):
        """
        Method to print the packet content for a specific frame

        Args:
            frame (integer): frame index whose data will be printed
        """
        print("--------------")
        print("Time: ", self.time)
        print("Telemetry mode: ", self.tm_mode[frame])
        print("Nframes: ", self.nframes)
        print("Multiplexing factor: ", self.mux)
        for column in range(cst.nb_col):
            print("\nb_cols {0:d}, frame {1:d}:".format(column, frame))
            for pix in range(self.mux):
                print('{0:X} '.format(self.data[column][frame, pix]), end='')
        print("\n")


    #------------------------------------------------------------
    def packet_append(self, lines_chunk):
        """
        Method to append the packet content with a scd file slice

        Args:
            lines_chunk (array of strings): string lines corresponding 
        to a frame
        """
        p = import_packet(lines_chunk)
        self.tm_mode=np.append(self.tm_mode, p.tm_mode)
        self.nframes=self.nframes+p.nframes
        self.data[0]=np.append(self.data[0], p.data[0], axis=0)
        self.data[1]=np.append(self.data[1], p.data[1], axis=0)
        self.data[2]=np.append(self.data[2], p.data[2], axis=0)
        self.data[3]=np.append(self.data[3], p.data[3], axis=0)
        

    #------------------------------------------------------------
    def plot(self, column=4, pixel=34, zoom=0, xy_limits=[0, 0, 0, 0]):
        """
        Method to plot the content of an ADC dump file

        Args:
            column (integer): column index whose data will be ploted
            if 4 a multiplot will be made with all the columns
            (default = 4)
            pixel (integer): pixel index whose data will be plotted
            if 34 a multiplot will be made with all the pixels
            (default = 34)
            zoom (integer): number of first samples to be plotted.
            If 0, no zoom is applied (defaut = 0)
        """

        filename = self.fullfilename.split('/')[-1]

        if (self.tm_mode[0] == cst.dmp_name):
            plot_attributes = {
                "ticks" : 'o-',
                "xlabel" : "Samples",
                "title" : filename + ", dump file"
            }
            plot_tm_adc_dmp(self.data, self.fullfilename, plot_attributes, column, zoom)
        elif (self.tm_mode[0] == cst.adc_name):
            plot_attributes = {
                "ticks" : '.',
                "xlabel" : "Samples",
                "title" : filename + ", ADC file"
            }
            plot_tm_adc_dmp(self.data, self.fullfilename, plot_attributes, column, zoom)
        elif (self.tm_mode[0] == cst.sci_name):
            plot_attributes = {
                "ticks" : ' ',
                "xlabel1" : "Time (µs)",
                "xlabel2" : "Frames",
                "xy_limits" : xy_limits,
                "title" : filename + ", science file"
            }
            plot_tm_sci(self.data, self.fullfilename, plot_attributes, column, pixel)
        else:            
            raise ValueError("Wrong file header!")

    #------------------------------------------------------------
    def min(self):
        """This function returns the minimum value of each column
        """        
        return(np.array([min(self.data[0][0][:]), min(self.data[1][0][:]), min(self.data[2][0][:]), min(self.data[3][0][:])]))             

    #------------------------------------------------------------
    def max(self):
        """This function returns the maximum value of each column
        """        
        return(np.array([max(self.data[0][0][:]), max(self.data[1][0][:]), max(self.data[2][0][:]), max(self.data[3][0][:])]))             

    #------------------------------------------------------------
    def mean(self):
        """This function returns the average value of each column
        """        
        return(np.array([np.mean(self.data[0][0][:]), np.mean(self.data[1][0][:]), np.mean(self.data[2][0][:]), np.mean(self.data[3][0][:])]))             

    #------------------------------------------------------------
    def edge_detect_old(self, col=0):
        """This function detects the edges of the signal and returns the sample number

        Args:
            col (int, optional): index of the column to consider. Defaults to 0.
        """        
        trigger = self.mean()
        x_higher = np.nonzero(self.data[col][0] > trigger[col])[0]
        if x_higher[0] == 0:    # case of a rising edge
            x_edge = np.nonzero(self.data[col][0] < trigger[col])[0][0]
        else:                   # case of a falling edge
            x_edge = x_higher[0]
        return(x_edge)
    
    #------------------------------------------------------------
    def edge_detect_0(self, col=0):
        """This function detects a pulse upward or downward in a signal and 
           it measures the position of the rising or falling edge.
           
        Args:
            col (int, optional): index of the column to consider. Defaults to 0.
        """        
        trigger = self.mean()[col]
        ratio_high = 100 * (trigger - self.min()[col]) / (self.max()[col] - self.min()[col])
        pulse_high = ratio_high < 50 # Signal is low with a short pulse upward
        
        x_higher = np.nonzero(self.data[col][0] > trigger)[0]

        if pulse_high:  # case of a pulse upward
            print("  The signal is low with a pulse upward", end="")
            x_edge = x_higher[0]
        else:           # case of a pulse downward
            print("  The signal is high with a pulse downward", end="")
            x_edge = np.nonzero(self.data[col][0] < trigger)[0][0]
        print("(the signal is ~{0:2.1f}% low and ~{1:2.1f}% high)".format(100-ratio_high, ratio_high))

        return(x_edge)
    
    #------------------------------------------------------------
    def edge_detect(self, col=0):
        """This function detects a pulse upward or downward in a signal and 
           it measures the position of the rising or falling edge.
           
        Args:
            col (int, optional): index of the column to consider. Defaults to 0.
        """        
        trigger = self.mean()[col]
        data=self.data[col][0]

        ratio_high = 100 * (trigger - self.min()[col]) / (self.max()[col] - self.min()[col])
        pulse_high = ratio_high < 50 # Signal is low with a short pulse upward
        
        if pulse_high:  # case of a pulse upward
            print("  The signal is low with a pulse upward", end="")
            x_edges = np.flatnonzero((data[:-1] < trigger) & (data[1:] > trigger))+1
        else:           # case of a pulse downward
            print("  The signal is high with a pulse downward", end="")
            x_edges = np.flatnonzero((data[:-1] > trigger) & (data[1:] < trigger))+1
        print("(the signal is ~{0:2.1f}% low and ~{1:2.1f}% high)".format(100-ratio_high, ratio_high))

        if len(x_edges) > 2 :
            print(x_edges)
            raise ValueError('  --> The number of Edges in the signal is incorrect (>2)!')
        
        frame_length = cst.nb_pix * cst.nb_samples_per_col
        x_edge =  x_edges[0] - frame_length        

        return(x_edge)

#------------------------------------------------------------
def dac_vs_adc_delay_plot(p1, p2, p3, dac_name, zoom=0):
    """
    Method to plot the dump files to characterize the DAC versus ADC delay

    Args:
        zoom (integer): number of first samples to be plotted.
        If 0, no zoom is applied (defaut = 0)
    """

    if (p1.tm_mode[0] != cst.dmp_name) \
        or (p2.tm_mode[0] != cst.dmp_name) \
        or (p3.tm_mode[0] != cst.dmp_name) :
        raise ValueError("Wrong file header!")
    else:
        start_test_ref=8    # Beginning of the test reference in the file name (i.e. UT_5018)
        filename_err = p1.fullfilename.split('/')[-1][start_test_ref:]
        filename_dac = p2.fullfilename.split('/')[-1][start_test_ref:]
        filename_dac_corr = p3.fullfilename.split('/')[-1][start_test_ref:]
        xlabel = "Samples"
        ticks1 = '.-'
        ticks2 = 'o-'
        ticks3 = 'x-'
        col1 = 'k'
        col2 = 'red'
        col3 = 'green'
        title = "Expertise " + dac_name + " delay (Error dump: " + filename_err + \
            ", " + dac_name + " dump: " + filename_dac + ", " + dac_name + " delayed dump: " + filename_dac_corr + ")"
        if zoom == 0:
            last_sample = len(p1.data[0][0,:])
        else:
            last_sample = zoom

        fig = plt.figure(figsize=(13, 13))
        for col in range(cst.nb_col):
            ax1 = fig.add_subplot(4, 1, 1+col)
            ax1.plot(p1.data[col][0,:last_sample], ticks1, color=col1)
            if col == 0:
                ax1.set_title(title)
            elif col == cst.nb_col-1 :
                ax1.set_xlabel(xlabel)
            ax1.set_ylabel("FAS -> ADC, col {0:d} (ADU)".format(col))
            
            ax2 = ax1.twinx()
            mksz = 3
            ax2.plot(p2.data[col][0,:last_sample], ticks2, color=col2, markersize=mksz, label="FB, no delay correction")
            ax2.plot(p3.data[col][0,:last_sample], ticks3, color=col3, markersize=mksz, label="FB, with delay correction")
            ax2.set_ylabel("DAC -> ADC, col {0:d} (ADU)".format(col))
            ax2.tick_params(axis ='y', labelcolor = col2) 

            ax2.legend(loc="best")

        fig.tight_layout()
        plot_file_name = p2.fullfilename.split('.')[0]+ ' ' + dac_name + '.png'
        plt.savefig(plot_file_name, dpi=300, bbox_inches='tight')
        plt.close()

#----------------------------------------------------------------
def import_packet(lines_chunk, verbose=False):
    """
    This function translates a slice / frame axpressed as string lines
    into a packet object

    Args:
        lines_chunk (array of strings): string lines corresponding 
            to a frame
        verbose (boolean): verbosity level. If true, some informations
            are provided (default = False)
    """
    line=lines_chunk[1]
    time_str=line.split(':')[1]
    time=int(time_str.split(' ')[1])

    line=lines_chunk[2]
    header=line.split(':')[1][1:-1]
    tm_mode = get_tm_mode(header)

    nframes=1

    if verbose:
        print("Importing packet ...")
        print("   packet type: ", tm_mode)
        
    col0,mux=read_column(lines_chunk[4])
    col1,_=read_column(lines_chunk[5])
    col2,_=read_column(lines_chunk[6])
    col3,_=read_column(lines_chunk[7])
    
    if tm_mode == cst.dmp_name:
        col0[0,:] = gtl.to_signed((col0[0,:]),14)
        col1[0,:] = gtl.to_signed((col1[0,:]),14)
        col2[0,:] = gtl.to_signed((col2[0,:]),14)
        col3[0,:] = gtl.to_signed((col3[0,:]),14)
    elif tm_mode == cst.adc_name or tm_mode == cst.dtv_name:
        col0[0,:] = gtl.to_signed((col0[0,:]),16)/4
        col1[0,:] = gtl.to_signed((col1[0,:]),16)/4
        col2[0,:] = gtl.to_signed((col2[0,:]),16)/4
        col3[0,:] = gtl.to_signed((col3[0,:]),16)/4
    elif tm_mode == cst.sci_name:
        col0[0,:] = gtl.to_signed((col0[0,:]),16)
        col1[0,:] = gtl.to_signed((col1[0,:]),16)
        col2[0,:] = gtl.to_signed((col2[0,:]),16)
        col3[0,:] = gtl.to_signed((col3[0,:]),16)
           
    p=packet('', time, tm_mode, nframes, mux, col0, col1, col2, col3)

    return(p)


#----------------------------------------------------------------
def read_column(line):
    """
    This function translates a line of the tm file into an array of
    individual values (string hexa)

    Args:
        line (string): one line of the tm data file
    """
    col_str=line.split(':')[1]
    col_array_str=np.array(col_str.split(',')[:-1])
    mux=len(col_array_str)
    col=np.zeros((1,mux)).astype(int)
    for pix in range (mux):
        col[0,pix]=int(col_array_str[pix], 16)
    return(col,mux)


#----------------------------------------------------------------
def read_tm(fullfilename, verbose=False):
    """
    This function imports a data packet from a science data file
    into a packet object

    Args:
        fullfilename (string): name (with path) of the tm data file
        verbose (boolean): specify the verbosity level. default to 0
    """
    f = open(fullfilename,"r")
    filename = fullfilename.split('/')[-1]

    lines = f.readlines()
    
    nblines = len(lines)
    nbpackets = int((nblines - cst.nb_lines_end_tmfile) / cst.nb_lines_packet)
    
    if verbose:
        print("The file ", filename, " contains {0:3d} packets".format(nbpackets))

    p=import_packet(lines[0:cst.nb_lines_packet]) # reading data of first packet
    p.fullfilename = fullfilename

    for ipacket in range(nbpackets-1): # appending the data of other packets
        p.packet_append(lines[(ipacket+1)*cst.nb_lines_packet:(ipacket+2)*cst.nb_lines_packet])
    
    return(p)


#----------------------------------------------------------------

def split_scd(fullfilename, verbose=False):
    """
    This function splits a data packet according to the tm type
    into a packet object

    Args:
        fullfilename (string): name (with path) of the tm data file
        verbose (boolean): specify the verbosity level. default to 0        
    """

    filename = fullfilename.split('/')[-1]

    f = open(fullfilename, "r")
    lines = f.readlines()
    
    nb_lines = len(lines)
    nb_packets = int((nb_lines - cst.nb_lines_end_tmfile) / cst.nb_lines_packet)
    if verbose:
        print("Splitting file ", filename, "...")
        print("   This file contains {0:3d} packets".format(nb_packets))

    packet_index = 0
    split_index = 0
    line_index = 0
    files_list=[]
    while packet_index < nb_packets:
        current_nb_packet = 1
        split_filename = fullfilename+"{0:d}".format(split_index)

        # Checking tm mode of the next packet
        current_tm_mode = get_tm_mode(lines[2 + packet_index * cst.nb_lines_packet].split(':')[1][1:-1])
        
        # Comparing the tm mode of the packets to the tm mode of the first packet
        while packet_index < nb_packets-1 \
            and get_tm_mode(lines[2 + (packet_index+1) * cst.nb_lines_packet].split(':')[1][1:-1]) == current_tm_mode \
            and current_tm_mode != cst.dmp_name:
            packet_index += 1
            current_nb_packet += 1
        if verbose:
            print("     Found {0:3d} packet(s) with the tm mode ".format(current_nb_packet), current_tm_mode, " => ", split_filename)
        
        # Saving the packets with the same tm mode in a file
        # if all the packets in the initial file have the same tm mode this step is skipped
        all_the_packets_of_the_file_have_same_tm_mode = (current_nb_packet == nb_packets and split_index == 0)
        if not all_the_packets_of_the_file_have_same_tm_mode: 
            files_list.append(split_filename)
            f_split = open(split_filename, "w")
            for l in range (current_nb_packet*cst.nb_lines_packet):
                f_split.write(lines[line_index+l])
            # printing end of file
            for l in range (cst.nb_lines_end_tmfile):
                f_split.write(lines[-cst.nb_lines_end_tmfile + l])
            f_split.close()

        packet_index += 1
        split_index += 1
        line_index += current_nb_packet * cst.nb_lines_packet

    return(files_list)

#----------------------------------------------------------------

def plot_tm_adc_dmp(data, fullfilename, plt_at, column=cst.nb_col, zoom=0):
    
    length = len(data[0][0,:]) * len(data[0][:,0])
    if zoom == 0:
        last_sample = length    # plot range
        extension = ""          # filename extension
        mksz = 1                # markersie on the plot
    else:
        last_sample = min(zoom, length)
        extension="_zoom"
        mksz = 4

    if column == 4:
        fig = plt.figure(figsize=(12, 12))
        for col in range(cst.nb_col):
            ax = fig.add_subplot(4, 1, 1+col)
            
            d = np.reshape(data[col],(len(data[col][:,0])*len(data[col][0,:])))
            ax.plot(d[:last_sample], plt_at["ticks"], markersize=mksz, label='Column {0:d}, '.format(col))

            ax.grid()
            ax.set_ylabel("ADU")
            ax.legend(loc="upper right")

            if col == 0:
                ax.set_title(plt_at["title"])
            elif col == cst.nb_col-1 :
                ax.set_xlabel(plt_at["xlabel"])

        fig.tight_layout()
        plot_file_name = fullfilename.split('.')[0]+extension+'_4cols.png'
            
    else:
        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_subplot(1, 1, 1)
        
        d = np.reshape(data[column],(len(data[column][:,0])*len(data[column][0,:])))
        ax.plot(d[:last_sample],plt_at["ticks"], markersize=mksz)
        
        ax.grid()
        ax.set_ylabel("ADU")
        ax.set_title('Column {0:d}, '.format(column) + plt_at["title"])
        ax.set_xlabel(plt_at["xlabel"])

        fig.tight_layout()
        plot_file_name = fullfilename.split('.')[0]+extension+'_col{0:d}.png'.format(column)
    
    plt.savefig(plot_file_name, dpi=300, bbox_inches='tight')
    plt.close()
    print("  Plot saved in ", plot_file_name)

#----------------------------------------------------------------

def plot_tm_sci(data, fullfilename, plt_at, column=cst.nb_col, pixel=cst.nb_pix):
    
    # Here we do one plot per column 
    if column == cst.nb_col:
        fig = plt.figure(figsize=(12, 12))
        plot_file_name = fullfilename.split('.')[0]
        for col in range(column):
            ax = fig.add_subplot(4, 1, 1+col)
            plot_file_name = column_plot_tm_sci(ax, data, plot_file_name, plt_at, col, pixel)                                            
    # Here we do the plot for a single column 
    else:
        fig = plt.figure(figsize=(8, 5))
        plot_file_name = fullfilename.split('.')[0] + '_col{0:d}'.format(column)
        ax = fig.add_subplot(1, 1, 1)
        plot_file_name = column_plot_tm_sci(ax, data, plot_file_name, plt_at, column, pixel)                 

    fig.tight_layout()    
    plt.savefig(plot_file_name, dpi=300, bbox_inches='tight')
    plt.close()

    print("  Plot saved in ", plot_file_name)

#----------------------------------------------------------------

def column_plot_tm_sci(ax, data, plot_file_name, plt_at, col, pixel):
    char_size = 8
    
    l = len(data[0][:,0])
    time = np.arange(l)*1e6/cst.f_frame
    # Here we do the plot for all the pixels 
    if pixel == cst.nb_pix:
        cmap = plt.get_cmap('hsv')
        colors = cmap(np.linspace(0, 1.0, cst.nb_pix))
        pixels = np.arange(cst.nb_pix)
        for pix, color in zip(pixels, colors):
            ax.plot(time, data[col][:,pix], marker=plt_at["ticks"], linewidth=0.8, color=color)
        ax.set_title(plt_at["title"] + ', Column {0:d}'.format(col))
        ax.set_xlabel(plt_at["xlabel1"])
        # Doing x zoom if requested
        if plt_at["xy_limits"][0] != 0 or plt_at["xy_limits"][1] != 0:
            ax.set_xlim(plt_at["xy_limits"][0], plt_at["xy_limits"][1])
            plot_file_name = plot_file_name+'_xzoom'
        # Doing y zoom if requested
        if plt_at["xy_limits"][2] != 0 or plt_at["xy_limits"][3] != 0:
            ax.set_ylim(plt_at["xy_limits"][2], plt_at["xy_limits"][3])
            plot_file_name = plot_file_name+'_yzoom'
        ax2 = ax.twiny()
        ax2.plot(np.arange(l), np.zeros(l), ' ') # fake plot
        ax2.set_xlim(0, l)
        ax2.set_xlabel(plt_at["xlabel2"])
        plot_file_name = plot_file_name+'.png'
    # Here we do the plot for a single pixel 
    else:
        ax.plot(time, data[col][:,pixel], marker=plt_at["ticks"])
        ax.set_title(plt_at["title"] + ', Column {0:d}, pixel {1:d}'.format(col, pixel))
        plot_file_name = plot_file_name+'_pix{0:d}.png'.format(pixel)
        ax.set_xlabel(plt_at["xlabel1"])
        ax2 = ax.twiny()
        ax2.plot(np.arange(l), np.zeros(l), ' ') # fake plot
        ax2.set_xlim(0, l)
        ax2.set_xlabel(plt_at["xlabel2"])

        for item in ([ax.xaxis.label, ax.yaxis.label, ax2.xaxis.label]):
            item.set_weight('bold')
            item.set_fontsize(char_size)
        for item in (ax.get_xticklabels() + ax.get_yticklabels() + ax2.get_xticklabels()):
            item.set_fontsize(char_size)
   
    ax.set_xlim(time[0], time[-1])
    ax.grid()
    ax.set_ylabel("ADU")
    return(plot_file_name)

#----------------------------------------------------------------
