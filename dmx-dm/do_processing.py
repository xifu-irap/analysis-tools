# imports
import os

import fdbkDelayAnalysis
import nonLinearity
import samplingDelayAnalysis


def do_processing(verbose=False):
    """
    This function launch the data processing of the DEMUX performance tests

    Args:
        verbose: (boolean) if True some informations are displayed during the processing
                default is False

    Returns:

    """
    # Data directory
    dir_path = os.path.join("..", "..")
    full_session_name = os.path.realpath(dir_path)
    session_name = os.path.basename(full_session_name)
    test_type = session_name[16: 30]
    print(test_type)

    match test_type:
        case "FDBK_DELAY____":
            fdbkDelayAnalysis.fdbkDelayAnalysis(verbose)
        case "SAMP_DELAY____":
            samplingDelayAnalysis.samplingDelayAnalysis()
        case "FDBK_ERROR_LIN" | "OFCO_ERROR_LIN":
            nonLinearity.nonlinearity()


do_processing()
