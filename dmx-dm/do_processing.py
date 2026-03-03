# imports
import os

import delayAnalysis
import noiseAnalysis
import nonLinearity
import samplingDelayAnalysis


# import XTalkAnalysis


def do_processing(verbose=True):
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
    test_type = session_name[16: 31]

    match test_type:
        case "FDBK-DELAY-----" | "OFCOMUX-DELAY--" | "OFCODAC-DELAY--":
            delayAnalysis.delayAnalysis(test_type, verbose)
        case "SAMP-DELAY-----":
            samplingDelayAnalysis.samplingDelayAnalysis(verbose)
        case "FDBK-ERROR-LIN-" | "OFCO-ERROR-LIN-":
            nonLinearity.nonlinearity()
        case "NOISE-ERRO-ONLY" | "NOISE-FDBK-ERRO" | "NOISE-OFCO-ERRO":
            noiseAnalysis.noiseAnalysis(verbose)


#        case "XTALK-PERP-FDBK" | "XTALK-PERP-OFCO":
#            XTalkAnalysis.xtalkAnalysis(verbose)


do_processing()
