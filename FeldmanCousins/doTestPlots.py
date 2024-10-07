import os
import logging, sys
import ROOT
import PlotUtils
import numpy as np
np.random.seed(0)
np.set_printoptions(precision=1)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)
from scipy import optimize, integrate
import argparse
ccnueroot = os.environ.get('CCNUEROOT')

import math
import psutil
import multiprocessing
import threading
nthreads = 4
from array import array

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from fit_tools.FitTools import *

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
MNVPLOTTER.mc_line_width = 0
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

rowstodo = 100
hist_dict = {"fhc mu sel":  6, "fhc e sel":  31, "rhc nueel":  57, "rhc mu sel":  63, "rhc e sel":  89}

if __name__ == "__main__":
    BLUEARC = "/exp/exp/minerva/data/users/{}/surfaces".format(os.environ["USER"])
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
                        default=0
    )
    parser.add_argument("-U", "--U_tau4",
                        dest = "U_tau4",
                        help="U_tau4 parameter to probe.",
                        default=0
    )
    args = parser.parse_args()
    outdir_surface = args.output_dir
    delta_m = int(args.delta_m)
    runongrid = args.grid
    U_tau4 = float(args.U_tau4)

    stitched_data = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('data_stitched')
    stitched_mc = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched')
    
    stitched_nueTemp = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('LE_template_nue')
    stitched_numuTemp = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('LE_template_numu')

    stitched_nue_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_nue')
    stitched_numu_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_numu')
    stitched_nutau_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_nutau')
    stitched_nueselection_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_nueselection')

    templates = {
            "nue":stitched_nueTemp,
            "numu":stitched_numuTemp,
            "nue_energy":stitched_nue_energy,
            "numu_energy":stitched_numu_energy,
            "nutau_energy":stitched_nutau_energy,
            "nueselection_energy":stitched_nueselection_energy,
            }

    #doFit(stitched_mc,templates,stitched_data,makePlot=True)
    #exit()

    for m in [10,50]:
        for U_e4 in [0.15]:
            for U_mu4 in [.183,.244,.306,.367,.428,.489]:
                for U_tau4 in [0]:
                    universe = GetOscillatedHistogram(stitched_data, templates, m, U_e4, U_mu4, U_tau4, doSys=True)
                    doFit(stitched_mc,templates,universe, makePlot = True, plotArgs = [m,U_e4,U_mu4,U_tau4])
    exit()

    dodelta = False
    if not runongrid: # surface plot
        for m in range(1,100,3):
            logging.info("initializing multiprocess threads")
            t1 = multiprocessing.Process(target=makeChi2Surface, args=(dodelta,m), name='t1')
            t2 = multiprocessing.Process(target=makeChi2Surface, args=(dodelta,m+1), name='t2')
            t3 = multiprocessing.Process(target=makeChi2Surface, args=(dodelta,m+2), name='t3')

            logging.info("starting multiprocess threads")
            t1.start()
            t2.start()
            t3.start()

            logging.info("joining multiprocess threads")
            t1.join()
            t2.join()
            t3.join()
        exit()
    else:
        m = delta_m / 2
        makeChi2Surface(dodelta,m,U_tau4)
