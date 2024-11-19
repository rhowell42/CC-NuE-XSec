import os
import logging, sys
import ROOT
import PlotUtils
import numpy as np
from root_numpy import matrix
np.set_printoptions(precision=1)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)
from scipy import optimize, integrate
from scipy.stats import multivariate_normal
import argparse
ccnueroot = os.environ.get('CCNUEROOT')
process = os.environ.get('PROCESS')

import math
from array import array

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from Tools.FitTools import *
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

def ThrowSystematics(stitched_mc,stitched_data,n_samples=50):
    includeStatError = False
    errorAsFraction  = True
    useOnlyShapeErrors = False

    covMatrixTmp  = stitched_mc.GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)
    covMatrixTmp += stitched_data.GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)
    for i in range(stitched_mc.GetNbinsX()+1):
        for j in range(i,stitched_mc.GetNbinsX()+1):
            covMatrixTmp[i][j] = covMatrixTmp[i][j] * stitched_mc.GetBinContent(i) * stitched_mc.GetBinContent(j)
            covMatrixTmp[j][i] = covMatrixTmp[i][j]

    covMatrix = np.asarray(matrix(covMatrixTmp))[1:-1,1:-1] # convert to numpy array, exclude over/underflow bins
    pred_vals = np.array(stitched_mc)[1:-1] # store MC bin contents excluding over/underflow bins

    sys_throws = multivariate_normal.rvs(mean=pred_vals,cov=covMatrix,size=n_samples)
    return(sys_throws)

def ThrowPoissons(lambdas):
    throws = []
    for lam in lambdas:
        while lam[lam<0].any():
            logging.error("Negative value in sys throws, rerunning this throw...")
            lam = ThrowSystematics(stitched_mc,stitched_data,1) 
        throw = np.random.poisson(lam)
        throws.append(throw)
    throws = np.array(throws)
    return(throws)

def FitToyExperiments(stitched_mc,stitched_data,experiments,templates):
    results = []
    for toy in experiments:
        weights = stitched_data.Clone().GetCVHistoWithStatError()
        for i in range(1,weights.GetNbinsX()+1):
            weight = stitched_data.GetBinContent(i) / toy[i-1] if toy[i-1] != 0 else stitched_data.GetBinContent(i)
            weights.SetBinContent(i,weight)
            weights.SetBinError(i,0)

        histogram = stitched_data.Clone()
        histogram.DivideSingle(histogram,weights)

        chi2_fit, res = doFit(histogram, templates, stitched_mc)
        chi2_mod = Chi2DataMC(histogram,stitched_mc)
        if chi2_mod == -1:
            Nbins = stitched_mc.GetNbinsX()
            #get the covariance matrix
            covMatrix = np.zeros(shape=[Nbins,Nbins],dtype='f')
            useOnlyShapeErrors = False
            includeStatError   = True
            errorAsFraction    = False
            covMatrixTmp  =   stitched_mc.GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)
            covMatrixTmp += histogram.GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)

            for i in range(0,Nbins):
                for j in range(0,Nbins):
                    covMatrix[i][j] = covMatrixTmp[i+1][j+1]

            errorMatrix = np.linalg.inv(covMatrix)
            np.savetxt("{}/err_inverse_{}.csv".format(outdir_surface,len(results)),errorMatrix,delimiter=',')
            np.savetxt("{}/err_{}.csv".format(outdir_surface,len(results)),covMatrix,delimiter=',')
            for name in stitched_mc.GetVertErrorBandNames():
                print(name)
                for n in range(stitched_mc.GetVertErrorBand(name).GetNHists()):
                    mc_vals = []
                    as_vals = []
                    for m in range(1,stitched_mc.GetNbinsX()+1):
                        mc_vals.append(stitched_mc.GetVertErrorBand(name).GetHist(n).GetBinContent(m))
                        as_vals.append(histogram.GetVertErrorBand(name).GetHist(n).GetBinContent(m))
                    print(mc_vals)
                    print(as_vals)
            DataMCCVPlot(histogram,stitched_mc,"{}/pseudo_experiment_{}.png".format(outdir_surface,len(results)))
            exit()

        elif chi2_fit > chi2_mod:
            logging.error("Negative chi2")
            chi2_fit = chi2_mod

        if chi2_fit == -1:
            logging.error("Fit returned -1 Chi2")
        elif chi2_fit < 0:
            logging.error("Fit returned negative Chi2")

        results.append(chi2_mod-chi2_fit)

    results = np.array(results)
    return(results)

def SplitList(alist,wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] for i in range(wanted_parts) ]

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
    stitched_swap_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_swap')

    templates = {
            "nue":stitched_nueTemp,
            "numu":stitched_numuTemp,
            "swap":stitched_swapTemp,
            "nue_energy":stitched_nue_energy,
            "numu_energy":stitched_numu_energy,
            "nutau_energy":stitched_nutau_energy,
            "swap_energy":stitched_swap_energy,
            "ratio_energy":stitched_ratio_energy,
            "nueselection_energy":stitched_nueselection_energy,
    }

    sys_throws  = ThrowSystematics(stitched_mc,stitched_data,50)
    experiments = ThrowPoissons(sys_throws)

    dchi2s = FitToyExperiments(stitched_mc,stitched_data,experiments,templates)
    np.savetxt('{}/sample_dchi2s_{}.csv'.format(outdir_surface,process),dchi2s,delimiter=',')
