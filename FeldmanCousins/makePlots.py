import os
import logging, sys
import copy
import ROOT
import PlotUtils
import numpy as np
np.set_printoptions(precision=4)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)
from scipy import optimize, integrate

import argparse
ccnueroot = os.environ.get('CCNUEROOT')

import math
import psutil
#import multiprocessing
import time
#import threading
#nthreads = 4
from array import array

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from Tools.Histogram import *
from Tools.PlotTools import *

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
#MNVPLOTTER.chi2_use_overflow_err = True
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.04
#MNVPLOTTER.extra_top_margin = -.035# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46 
MNVPLOTTER.mc_bkgd_line_color = 46
MNVPLOTTER.legend_n_columns = 3
#MNVPLOTTER.mc_line_width = 0
MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  
MNVPLOTTER.axis_maximum = 500 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .03
MNVPLOTTER.legend_offset_x           = .15

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

def CompareResults(histogram):
    dataHist = histogram.GetDataHistogram()
    mcHist = histogram.GetMCHistogram()

    data = np.array(dataHist)[1:-1]
    mc = np.array(mcHist)[1:-1]

    universes = histogram.GetFluxUniverses()
    invCov = histogram.GetInverseCovarianceMatrix(sansFlux=True) 

    A = histogram.GetAMatrix()
    V = invCov    
    C = data - mc
    I = np.identity(len(universes))

    L = 2 * A @ V @ C
    Q = A @ V @ A.T + I
    solution = np.linalg.inv(Q) @ L/2

    np.savetxt("ryan_quadraticTerm.csv",Q,delimiter=',')
    np.savetxt("ryan_linearTerm.csv",L,delimiter=',')
    np.savetxt("ryan_invQuadraticTerm.csv",np.linalg.inv(Q),delimiter=',')
    np.savetxt("ryan_invCVcovariance.csv",invCov,delimiter=',')
    np.savetxt("ryan_Cvector.csv",C,delimiter=',')
    np.savetxt("ryan_CVcovariance.csv",histogram.GetCovarianceMatrix(sansFlux=True),delimiter=',')
    np.savetxt("ryan_CVcovariance_NoFlux.csv",histogram.GetCovarianceMatrix(sansFlux=True),delimiter=',')
    np.savetxt("ryan_Fluxcovariance.csv",histogram.GetCovarianceMatrix(sansFlux=True) - histogram.GetCovarianceMatrix(sansFlux=True),delimiter=',')
    np.savetxt("ryan_Amatrix.csv",A,delimiter=',')

if __name__ == "__main__":
    BLUEARC = "/exp/minerva/data/users/{}/surfaces".format(os.environ["USER"])
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--grid",
                        action="store_true",
                        default = False,
                        help = "Run macro on grid, Input/Output path must be updated to avoid direct access to BlueArc"
    )
    parser.add_argument("-o", "--output",
                        dest = "output_dir",
                        help="Use alternate location for output file.",
                        default=BLUEARC
    )
    parser.add_argument("-m", "--delta_m",
                        dest = "delta_m",
                        help="Delta m^2 value to probe.",
                        type=float,
                        default=0
    )
    parser.add_argument("-Ue4", "--U_e4",
                        dest = "U_e4",
                        help="U_e4 parameter to probe.",
                        type=float,
                        default=0
    )
    parser.add_argument("-Umu4", "--U_mu4",
                        dest = "U_mu4",
                        help="U_mu4 parameter to probe.",
                        type=float,
                        default=0
    )
    parser.add_argument("-U", "--U_tau4",
                        dest = "U_tau4",
                        help="U_tau4 parameter to probe.",
                        type=float,
                        default=0
    )
    parser.add_argument("--pseudodata",
                        dest = "pseudodata",
                        default = False,
                        action="store_true",
    )

    args = parser.parse_args()
    outdir_surface = args.output_dir
    delta_m = args.delta_m
    runongrid = args.grid
    U_e4 = args.U_e4
    U_mu4 = args.U_mu4
    U_tau4 = args.U_tau4

    filename = ""
    filename = "NuE_stitched_hists.root"
    
    file_path = "{}/FeldmanCousins/{}".format(ccnueroot,filename)

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)
    sample_histogram.SetPlottingStyle()

    n4 = {
            'm':7,
            'ue4':.14,
            'umu4':.003,
            'utau4':0.5
        }

    OscillateHistogram(sample_histogram, n4['m'], n4['ue4'], n4['umu4'], n4['utau4'])

    invCov=sample_histogram.GetInverseCovarianceMatrix(sansFlux=True)
    nullSolution,nullPen = FluxSolution(sample_histogram,invCov=invCov)

    sample_histogram.PlotSamples(fluxSolution=nullSolution,plotName="AllSamples")
    PlotOscillationEffects(sample_histogram,n4,"Neutrino4",plotSamples=True)
    PlotOscillationRatios(sample_histogram,n4,"Neutrino4")
    PlotFluxMarginalizationEffects(sample_histogram,n4,"Neutrino4")
    PlotSampleMarginalizationEffects(sample_histogram)
