import os
import logging, sys
import ROOT
import PlotUtils

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
import numpy as np
np.set_printoptions(precision=5)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)
from scipy import optimize, integrate
import argparse
import subprocess
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
from Tools.PlotTools import *

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

def MakeSurface(histogram,outdir,dodeltachi2=False,deltam=1,U_e4=0,U_tau4=0,makePlot=False):
    U_mu4s = 0.41*np.logspace(-5,0,100)
    U_mu4s[0] = 0
    asimov_surface    = np.zeros(np.shape(U_mu4s)[0],dtype='f')
    data_surface      = np.zeros(np.shape(U_mu4s)[0],dtype='f')
    data_penalties    = np.zeros(np.shape(U_mu4s)[0],dtype='f')
    asimov_penalties  = np.zeros(np.shape(U_mu4s)[0],dtype='f')
    fits              = np.zeros(np.shape(U_mu4s)[0],dtype='f')
    count = 0

    for i in range(U_mu4s.shape[0]):
        count+=1
        U_mu4 = U_mu4s[i]

        OscillateHistogram(histogram, deltam, U_e4, U_mu4, U_tau4)

        chi2_data,data_penalty = Chi2DataMC(histogram,invCov=histogram.GetInverseCovarianceMatrix(),marginalize=True,useOsc=True)
        chi2_asimov,asimov_penalty = Chi2DataMC(histogram,invCov=histogram.GetInverseCovarianceMatrix(),marginalize=True,useOsc=True,usePseudo=True)

        if makePlot:
            res={"m":deltam,"ue4":U_e4,"umu4":U_mu4,"utau4":0}
            PlotOscillationEffects(histogram,res,"asimov",False,useMarg=False,usePseudo=True)
            PlotOscillationEffects(histogram,res,"asimov_marg",True,useMarg=True,usePseudo=True)

        data_surface[i] = chi2_data
        data_penalties[i] = data_penalty
        asimov_surface[i] = chi2_asimov
        asimov_penalties[i] = asimov_penalty
        logging.info("{:.2f}% done with chi2s. Current data, asimov chi2s = {:.4f}, {:.4f}".format(100*count/(U_mu4s.shape[0]),data_surface[i],asimov_surface[i]))

    np.save('{}/chi2_surface_data_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),data_surface)
    np.save('{}/chi2_surface_pseudodata_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),asimov_surface)
    np.save('{}/chi2_penalty_data_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),data_penalties)
    np.save('{}/chi2_penalty_pseudodata_m_{}_Ue4_{}.dat'.format(outdir,deltam,U_e4),asimov_penalties)

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
                        help="index of U_e4 parameter to probe.",
                        type=int,
                        default=0
    )
    parser.add_argument("--pseudodata",
                        dest = "pseudodata",
                        default = False,
                        action="store_true",
    )
    parser.add_argument("--dodeltachi2",
                        dest = "dodelta",
                        default = False,
                        action="store_true",
    )

    args = parser.parse_args()
    outdir = args.output_dir
    delta_m = args.delta_m
    runongrid = args.grid
    U_tau4 = args.U_tau4
    i_U_e4 = args.U_e4
    dodelta = args.dodelta

    filename = "NuE_stitched_hists.root"
    file_path = "{}/FeldmanCousins/{}".format(ccnueroot,filename)

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    U_e4s = 0.15*np.logspace(-4,0,100)
    U_e4s[0] = 0
    U_e4 = U_e4s[i_U_e4]

    if not runongrid: # surface plot
        m_toloop = np.logspace(-1,2,100)
        for m in m_toloop:
            MakeSurface(sample_histogram,outdir,dodelta,m,U_e4,U_tau4)
    else:
        MakeSurface(sample_histogram,outdir,dodelta,delta_m,U_e4,U_tau4,makePlot=False)
