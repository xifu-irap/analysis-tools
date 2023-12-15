
#imports
import expertises as exp
import constants as cst

path='/Users/laurent/Data/TestPlan_10_cosim/'


#----------------------------------------------------------------
# Processing selection  |
#------------------------

print('--------')
do_exp_syn_delay = False
do_exp_seq_delay = True
do_exp_mux_delay = False
do_exp_amp_delay = False
do_exp_amp_chara = False
do_exp_mux_chara = False        

#----------------------------------------------------------------

# ADC dump tm file
filename_ER='UT_5017/DRE_DMX_UT_5017_scd' 

# SYNC delay dump tm file
filename_SYNC_ref='UT_5031/DRE_DMX_UT_5031a_scd' # No sync delay correction
filename_SYNC_corb='UT_5031/DRE_DMX_UT_5031b_scd' # With sync delay correction = 1
filename_SYNC_corc='UT_5031/DRE_DMX_UT_5031c_scd' # With sync delay correction = 2
filename_SYNC_cord='UT_5031/DRE_DMX_UT_5031d_scd' # With sync delay correction = 3
filename_SYNC_core='UT_5031/DRE_DMX_UT_5031e_scd' # With sync delay correction = 4
filename_SYNC_corf='UT_5031/DRE_DMX_UT_5031f_scd' # With sync delay correction = 5
filename_SYNC_corg='UT_5031/DRE_DMX_UT_5031g_scd' # With sync delay correction = 6

# ERROR delay dump tm file
filename_ERROR_ref='UT_5032/DRE_DMX_UT_5032a_scd' # No sequence delay correction
error_delay_ref_dict={"config": "fpasim", "g_RA_DELAY": 6, "g_ERROR_DELAY": 0}
filename_ERROR_cor='UT_5032/DRE_DMX_UT_5032b_scd' # With sequence delay correction = 2
error_delay_cor_dict={"config": "fpasim", "g_RA_DELAY": 6, "g_ERROR_DELAY": 2}

# Feedback dump tm file
filename_FB='UT_5018/DRE_DMX_UT_5018_scd'
mux_delay_dict={"signal_name": "MUX SQUID feedback", "parameter_name": "CY_MUX_SQ_FB_DELAY"}

# Offset dump tm file
filename_OF='UT_5019/DRE_DMX_UT_5019_scd'
amp_delay_dict={"signal_name": "AMP SQUID offset", "parameter_name": "CY_AMP_SQ_OFFSET_MUX_DELAY"}

# amp squid characteristic adc tm file
filename_amp_charac='UT_5021_20230929/DRE_DMX_UT_5021_scd'
# test setup definition
OFFSET_LSB = 0xFFF * cst.amp_sq_quantum
OFFSET_FINE = [1, 1, 1]
shift = OFFSET_LSB * (0.5 * OFFSET_FINE[2] + 0.25 * OFFSET_FINE[1] + 0.125 * OFFSET_FINE[0])
amp_sq_VPhi_setup={"start_value": 0/4, "frames_per_steps": 2, "n_steps": 150, "step_size": 0x195/4, "shift_OFFSET": shift}

if do_exp_syn_delay:
    exp.sync_delay(path, filename_SYNC_ref, reference=True)
    exp.sync_delay(path, filename_SYNC_corb)
    exp.sync_delay(path, filename_SYNC_corc)
    exp.sync_delay(path, filename_SYNC_cord)
    exp.sync_delay(path, filename_SYNC_core)
    exp.sync_delay(path, filename_SYNC_corf)
    exp.sync_delay(path, filename_SYNC_corg)

if do_exp_seq_delay:
    exp.sequence_delay(path, filename_ERROR_ref, error_delay_ref_dict)
    exp.sequence_delay(path, filename_ERROR_cor, error_delay_cor_dict)

if do_exp_mux_delay:
    exp.dac_delay(path, filename_ER, filename_FB, mux_delay_dict)

if do_exp_amp_delay:
    exp.dac_delay(path, filename_ER, filename_OF, amp_delay_dict)

if do_exp_amp_chara:
    exp.amp_sq_vphi(path, filename_amp_charac, amp_sq_VPhi_setup, col=0)

#----------------------------------------------------------------
