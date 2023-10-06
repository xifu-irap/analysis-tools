
#imports
import expertises, constants

path='/Users/laurent/Data/TestPlan_10_cosim/'


#----------------------------------------------------------------
# Processing selection  |
#------------------------

print('--------')
do_exp_smp_delay = True
do_exp_mux_delay = True
do_exp_amp_delay = True
do_exp_amp_chara = True
do_exp_mux_chara = True        

#----------------------------------------------------------------

# ADC dump tm file
filename_ER='UT_5017/DRE_DMX_UT_5017_scd'

# Feedback dump tm file
filename_FB='UT_5018/DRE_DMX_UT_5018_scd'
mux_delay_dict={"signal_name": "MUX SQUID feedback", "parameter_name": "CY_MUX_SQ_FB_DELAY"}

# Offset dump tm file
filename_OF='UT_5019/DRE_DMX_UT_5019_scd'
amp_delay_dict={"signal_name": "AMP SQUID offset", "parameter_name": "CY_AMP_SQ_OFFSET_MUX_DELAY"}

# amp squid characteristic adc tm file
filename_amp_charac='UT_5021_20230929/DRE_DMX_UT_5021_scd'
# test setup definition
OFFSET_LSB = 0xFFF * constants.amp_sq_quantum
OFFSET_FINE = [1, 1, 1]
shift = OFFSET_LSB * (0.5 * OFFSET_FINE[2] + 0.25 * OFFSET_FINE[1] + 0.125 * OFFSET_FINE[0])
amp_sq_VPhi_setup={"start_value": 0/4, "frames_per_steps": 2, "n_steps": 150, "step_size": 0x195/4, "shift_OFFSET": shift}

if do_exp_smp_delay:
    expertises.sampling_delay(path, filename_ER)

if do_exp_mux_delay:
    expertises.dac_delay(path, filename_ER, filename_FB, mux_delay_dict)

if do_exp_amp_delay:
    expertises.dac_delay(path, filename_ER, filename_OF, amp_delay_dict)

if do_exp_amp_chara:
    expertises.amp_sq_vphi(path, filename_amp_charac, amp_sq_VPhi_setup, col=0)

#----------------------------------------------------------------
