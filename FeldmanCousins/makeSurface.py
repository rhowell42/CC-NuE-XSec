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
from Tools.FitTools import *

#logging.basicConfig(stream=sys.stdout, level=logging.INFO)

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
    parser.add_argument("-U", "--U_tau4",
                        dest = "U_tau4",
                        help="U_tau4 parameter to probe.",
                        type=float,
                        default=0
    )
    parser.add_argument("-Ue4", "--U_e4",
                        dest = "U_e4",
                        help="U_e4 parameter to probe.",
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
    U_tau4 = args.U_tau4
    U_e4 = args.U_e4
    pseudodata = args.pseudodata

    filename = ""
    if pseudodata:
        filename = "NuE_stitched_hists_pseudo.root"
    else:
        filename = "NuE_stitched_hists.root"

    stitched_data = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('data_stitched')
    stitched_mc = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched')
    
    stitched_nueTemp = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('LE_template_nue')
    stitched_numuTemp = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('LE_template_numu')
    stitched_swapTemp = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('LE_template_swap')

    stitched_nue_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_nue')
    stitched_numu_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_numu')
    stitched_nutau_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_nutau')
    stitched_nueselection_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_nueselection')
    stitched_ratio_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_ratio')
    stitched_fhc_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_fhc')
    stitched_swap_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_swap')

    templates = {
            "nue":stitched_nueTemp,
            "numu":stitched_numuTemp,
            "swap":stitched_swapTemp,
            "nue_energy":stitched_nue_energy,
            "numu_energy":stitched_numu_energy,
            "nutau_energy":stitched_nutau_energy,
            "fhc_energy":stitched_nutau_energy,
            "swap_energy":stitched_swap_energy,
            "ratio_energy":stitched_ratio_energy,
            "nueselection_energy":stitched_nueselection_energy,
    }

    dodelta = False
    if not runongrid: # surface plot
        m_toloop = np.logspace(0,2,60)
        for m in m_toloop:
            makeChi2Surface(templates,stitched_data,stitched_mc,outdir_surface,dodelta,m,U_e4,U_tau4)
            #logging.info("completed delta m2 = {:.2f} run".format(m))
            #t1 = multiprocessing.Process(target=makeChi2Surface, args=(templates,stitched_data,stitched_mc,outdir_surface,dodelta,dm,U_tau4), name='t1')
            #t2 = multiprocessing.Process(target=makeChi2Surface, args=(templates,stitched_data,stitched_mc,outdir_surface,dodelta,dm+0.5,U_tau4), name='t2')
            #t3 = multiprocessing.Process(target=makeChi2Surface, args=(templates,stitched_data,stitched_mc,outdir_surface,dodelta,dm+1,U_tau4), name='t3')
            #t4 = multiprocessing.Process(target=makeChi2Surface, args=(templates,stitched_data,stitched_mc,outdir_surface,dodelta,dm+1.5,U_tau4), name='t4')

            #logging.info("starting multiprocess threads")
            #t1.start()
            #t2.start()
            #t3.start()
            #t4.start()

            #logging.info("joining multiprocess threads")
            #t1.join()
            #t2.join()
            #t3.join()
            #t4.join()
        exit()
    else:
        print("******* {} *******".format(U_tau4))
        makeChi2Surface(templates,stitched_data,stitched_mc,outdir_surface,dodelta,delta_m,U_e4,U_tau4)
