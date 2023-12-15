
# TDM constants
nb_col = 4
nb_pix = 34
f_master = 125e6
nb_samples_per_col = 20
f_row = f_master / nb_samples_per_col
f_frame = f_row / nb_pix

# Hardware constants
nb_bits_adc = 14
amp_sq_quantum = 3/0xFFF # 3V

# Delays
sampling_max_delay = 19

# DEMUX science data modes
dmp_name = "   dump"
sci_name = "science"
dtv_name = "datavld"
adc_name = "    ADC"
tm_dict = {dmp_name: "C2",
           sci_name: "C8",
           dtv_name: "CA",
           adc_name: "E2"
}

# packets constants
nb_lines_packet = 8
nb_lines_end_tmfile = 3

# strings
draw_line = "\n------------------------------------------------"   
measure_file_text = "  We measure the following dump file:     "                