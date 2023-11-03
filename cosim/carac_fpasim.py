
#imports
import scd_tools as stl

path='/Users/laurent/Data/TestPlan_10_cosim/'


#----------------------------------------------------------------
# Processing selection  |
#------------------------

do_carac_pulse_shape = True
do_carac_sqa_tf = True
do_carac_sqm_tf_1phi = True
do_carac_sqm_tf_8phi = True

def plot_tm_curve(fullfilename):
    print("  Processing file...")
    p = stl.read_tm(fullfilename, verbose=True)
    mini, maxi = p.data[0].min(), p.data[0].max()
    print("  The characteristic ranges from {0:d} to {1:d}".format(mini, maxi))
    print("  The full scale range for this signal is -8192 / 8191 (14 bits signed).")
    p.plot(column=0)

#----------------------------------------------------------------

# pulse shape data file
filename_pulse='UT_5025/DRE_DMX_UT_5025_scd'

# AMP SQUID TF data file
filename_sqa='UT_5026/DRE_DMX_UT_5026_scd'

# MUX SQUID TF data file
filename_sqm1='UT_5027/DRE_DMX_UT_5027_scd'

# AMP SQUID TF data file
filename_sqm8='UT_5028/DRE_DMX_UT_5028_scd'

if do_carac_pulse_shape:

    print("------------------------------------------------------------------------")
    print("Characterisation of the FPAsim pulse shape")
    fullfilename = path+filename_pulse
    plot_tm_curve(fullfilename)

if do_carac_sqa_tf:

    print("------------------------------------------------------------------------")
    print("Characterisation of the FPAsim AMP SQUID transfer function")
    fullfilename = path+filename_sqa
    plot_tm_curve(fullfilename)

if do_carac_sqm_tf_1phi:

    print("------------------------------------------------------------------------")
    print("Characterisation of the FPAsim MUX SQUID transfer function (1 Phi0)")
    fullfilename = path+filename_sqm1
    plot_tm_curve(fullfilename)

if do_carac_sqm_tf_8phi:

    print("------------------------------------------------------------------------")
    print("Characterisation of the FPAsim MUX SQUID transfer function (8 Phi0)")
    fullfilename = path+filename_sqm8
    plot_tm_curve(fullfilename)

#----------------------------------------------------------------
