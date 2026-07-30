# ---------------------------------------------------------------------------------
# !/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright (C) 2021-2030 Laurent Ravera, IRAP Toulouse.
#  This file is part of the ATHENA X-IFU DRE data analysis tools software.
#
#  analysis-tools is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# ---------------------------------------------------------------------------------
#
#  laurent.ravera@irap.omp.eu
#  do_analysis.py
#
# ---------------------------------------------------------------------------------

# imports
import os

import XTalkAnalysis
import cutoffFrequency
import delayAnalysis
import delayRangeAnalysis
import noiseAnalysis
import nonLinearity
import ofcoCoarseAnalysis
import ofcoFineAnalysis
import samplingDelayAnalysis
import thermal_analysis


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
        case "FDBK-DELAY-----":
            delayRangeAnalysis.delayRangeAnalysis(test_type, verbose)
            delayAnalysis.delayAnalysis(test_type, 10, verbose)
        case "OFCOMUX-DELAY--":
            delayRangeAnalysis.delayRangeAnalysis(test_type, verbose)
            delayAnalysis.delayAnalysis(test_type, 10, verbose)
        case "SAMP-DELAY-----":
            samplingDelayAnalysis.samplingDelayAnalysis(verbose)
        case "FDBK-ERROR-LIN-" | "OFCO-ERROR-LIN-":
            nonLinearity.nonlinearity()
        case "NOISE-ERRO-ONLY" | "NOISE-FDBK-ERRO" | "NOISE-OFCO-ERRO":
            process_frow = False
            plot_model = False
            noiseAnalysis.noiseAnalysis(process_frow=process_frow, plot_model=plot_model, lpf=0, verbose=verbose)
        case "XTALK-PERP-FDBK" | "XTALK-PERP-OFCO":
            XTalkAnalysis.xtalkAnalysis(verbose)
        case "ERRO_BANDSHAPE-":
            cutoffFrequency.cutoffFreq(verbose)
        case "OFCO-FINE-HRESO" | "OFCO-FINE-HRANG":
            ofcoFineAnalysis.ofcoFineResoAndRange(test_type, verbose)
        case "OFCO-COAR-TPT--":
            ofcoCoarseAnalysis.ofcoCoarseReso(verbose)
        case "THERMAL_CHARAC-":
            thermal_analysis.thermal_analysis(verbose)


do_processing()
