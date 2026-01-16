# imports
import os

import fdbkDelayAnalysis
import samplingDelayAnalysis


def do_processing(verbose=False):
    # Data directory
    dir_path = os.path.join("..", "..")
    full_session_name = os.path.realpath(dir_path)
    session_name = os.path.basename(full_session_name)
    test_type = session_name[16: 25]
    print(test_type)

    if test_type == "fdbkDelay":
        fdbkDelayAnalysis.fdbkDelayAnalysis(verbose);
    if test_type == "sampDelay":
        samplingDelayAnalysis.samplingDelayAnalysis();


do_processing()
