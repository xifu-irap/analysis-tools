import os
import get_data
import general_tools

config = general_tools.get_csv('demux_tools_cfg.csv')

datadirname = os.path.join(os.path.normcase(config['dir_data']))

dumpfilenames = [f for f in os.listdir(datadirname) \
                if os.path.isfile(os.path.join(datadirname, f)) \
                and f[-4:]==".dat"]

for file in dumpfilenames:
    print('\n#---------------------')
    get_data.check_dump(config, file, max_duration=0.000001, spectral=True)

