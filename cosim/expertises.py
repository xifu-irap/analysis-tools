
#imports
import numpy as np
import scd_tools, general_tools, constants
import matplotlib.pyplot as plt


#----------------------------------------------------------------
def sampling_delay(path, filename_err, col=0):
    """This function characterize the delay of a DAC signal
    with respect to the address line signal

    Args:
        path (string): path to the data files
        filename_err (string): name of the error dump tm file
        col (int, optional): index of the column to be characterized. Defaults to 0.
    """

    print(constants.draw_line)
    print("Characterizing the sampling delay")
    print("    Pixel 0 has a specific TES setting point.")
    print(constants.measure_file_text, filename_err)
    fullfilename_err = path+filename_err
    files_err = scd_tools.split_scd(fullfilename_err)
    p_err = scd_tools.read_scd(files_err[1])
    p_err.plot(column=col)
   
    x_edge_err = p_err.edge_detect(col=col)
    print("\n The pixel 0 starts at sample # {0:d} after the synchronisation".format(x_edge_err))

    sample_id = 16
    delay = sample_id + x_edge_err
    if delay > constants.sampling_max_delay:
        raise ValueError(" >> The sampling is too high, it cannot be corrected")
    else:
        print(" Considering we want the sample # {0:d}/20 of each pixel, the delay shall be set to {1:d}".format(sample_id, delay))
        print("\n CY_SAMPLING_DELAY = " + hex(delay) + "\n")    
    
#----------------------------------------------------------------
def dac_delay(path, filename_err, filename_dac, dac_delay_dict, col=0):
    """This function characterize the delay of a DAC signal
    with respect to the address line signal

    Args:
        path (string): path to the data files
        filename_err (string): name of the error dump tm file
        filename_dac (string): name of the feedback dump tm file
        dac_delay_dict (dictionnary): with signal_name (string) and parameter_name (string)
        col (int, optional): index of the column to be characterized. Defaults to 0.
    """

    print(constants.draw_line)
    print("Characterizing the delay of the", dac_delay_dict["signal_name"], "signal with respect to the error signal")

    # Processing delay reference data
    print("  Measurement 1:")
    print("    Pixel 0 has a specific TES setting point")
    print(constants.measure_file_text, filename_err)
    fullfilename_err = path+filename_err
    files_err = scd_tools.split_scd(fullfilename_err)
    p_err = scd_tools.read_scd(files_err[1])
    x_edge_err = p_err.edge_detect(col=col)
    print("    The pixel 0 starts at sample # {0:d} after the synchronisation".format(x_edge_err))

    # Processing data generated from a DRE DAC (no delay correction)
    print("\n  Measurement 2:")
    print("    In open loop, pixel 0 has a specific FB0")
    print(constants.measure_file_text, filename_dac)
    fullfilename_dac = path+filename_dac
    files_dac = scd_tools.split_scd(fullfilename_dac)
    p_dac = scd_tools.read_scd(files_dac[0])
    x_edge_dac = p_dac.edge_detect(col=col)
    print("    The pixel 0 starts at sample # {0:d} after the synchronisation".format(x_edge_dac))

    # Computing delay correction
    delay = x_edge_err - x_edge_dac
    nbits_format = 10
    print("\n  -> The " + dac_delay_dict["signal_name"] + " signal shall be delayed by {0:d} samples".format(delay))
    print("\n " + dac_delay_dict["parameter_name"] + " = " + general_tools.dec_to_signed_hexa(delay, nbits_format) + "\n")    

    # Processing data generated from a DRE DAC (after delay correction)
    p_dac_corrected = scd_tools.read_scd(files_dac[1])
    x_edge_dac_corrected = p_dac_corrected.edge_detect(col=col)
    print("    The edge of the", dac_delay_dict["signal_name"], "signal after delay correction is on sample {0:d}\n".format(x_edge_dac_corrected))
    
    scd_tools.dac_vs_adc_delay_plot(p_err, p_dac, p_dac_corrected, dac_delay_dict["signal_name"], zoom=200)
        
#----------------------------------------------------------------
def amp_sq_vphi(path, filename, setup, col=0):
    """This function plots the V/Phi curve from an AMP SQUID expertise measurement

    Args:
        path (string): path to the data files
        filename_err (string): name of the error tm file
        setup (dictionnary): contains the setting of the expertise measurement
        col (int, optional): index of the column to be characterized. Defaults to 0.
    """

    # Defining plot text
    title = "AMP SQUID $V(\Phi)$ curve"
    xlabel = "Expected DRE AMP SQUID offset voltage (V)"
    ylabel = "DRE ADC value (ADU)"
    label1 = "Data from pattern"
    label2 = "Data from shifted pattern"
    
    print(constants.draw_line)
    print("Characterizing the V/Phi curve of the AMP SQUID on column {0:d}".format(col))
    print(" TM error file:     ", filename)
    fullfilename = path+filename

    # Getting data from file
    p = scd_tools.read_scd(fullfilename, verbose=True)

    #   - In a frame pixels 0-16 use an offset of 0, pixels 17-33 use an offset of FSR/2
    #   - The few first data needs to be ignored (settling time of the slow DAC)
    
    # Counting the number of steps (data valids)
    modes, counts = np.unique(p.tm_mode, return_counts=True)
    sorted_modes = dict(zip(modes, counts))
    n_steps = sorted_modes['datavld']
    print('The V/Phi curve contains {0:d} steps'.format(n_steps))
    if n_steps != setup["n_steps"]:
        raise ValueError("The number of steps is incorrect!")
    
    # Making X array
    nbits_amp_sq = 12
    nbits_pattern = 14
    correction_factor = nbits_amp_sq - nbits_pattern
    v_offset_low = 1/6 * (setup["start_value"] + np.arange(n_steps)*setup["step_size"]) * constants.amp_sq_quantum * 2**(correction_factor)
    v_offset_high = v_offset_low + 1/6 * ( setup["shift_OFFSET"])
    
    # Making Y array
    pixel = 14
    delta = int(p.mux / 2)
    
    counter = 0
    while p.tm_mode[counter] != 'datavld':
        counter+=1
    
    scan_low = p.data[col][counter:counter+setup["frames_per_steps"]*n_steps, pixel]
    data_low = scan_low[0::setup["frames_per_steps"]]
    scan_high = p.data[col][counter:counter+setup["frames_per_steps"]*n_steps, pixel+delta]
    data_high = scan_high[0::setup["frames_per_steps"]]

    # Doing the plot
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(1, 1, 1)
    
    mksz = 4
    ax.plot(v_offset_low, data_low, '3r', markersize=mksz, label=label1)
    ax.plot(v_offset_high, data_high, '4b', markersize=mksz, label=label2)

    ax.set_xlim(0,1)

    ax.grid()
    ax.set_title(title, color='blue')
    ax.title.set_weight('bold')
    ax.title.set_fontsize(16)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc='best')

    for item in ([ax.xaxis.label, ax.yaxis.label]):
        item.set_weight('bold')
        item.set_fontsize(12)
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(12)
    
    fig.tight_layout()
    plot_file_name = p.fullfilename.split('.')[0]+'.png'
    plt.savefig(plot_file_name, dpi=300, bbox_inches='tight')
   
#----------------------------------------------------------------
