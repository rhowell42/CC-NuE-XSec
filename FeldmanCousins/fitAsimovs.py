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

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

def ThrowSystematics(histogram,n_samples=50):
    includeStatError = False
    errorAsFraction  = True
    useOnlyShapeErrors = False

    covMatrixTmp  = histogram.GetMCHistogram().GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)
    covMatrixTmp += histogram.GetPseudoHistogram().GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)

    mc_hist = histogram.GetMCHistogram()
    for i in range(mc_hist.GetNbinsX()+1):
        for j in range(i,mc_hist.GetNbinsX()+1):
            covMatrixTmp[i][j] = covMatrixTmp[i][j] * mc_hist.GetBinContent(i) * mc_hist.GetBinContent(j)
            covMatrixTmp[j][i] = covMatrixTmp[i][j]

    covMatrix = np.asarray(matrix(covMatrixTmp))[1:-1,1:-1] # convert to numpy array, exclude over/underflow bins
    pred_vals = np.array(mc_hist)[1:-1] # store MC bin contents excluding over/underflow bins

    sys_throws = multivariate_normal.rvs(mean=pred_vals,cov=covMatrix,size=n_samples)
    return(sys_throws)

def ThrowPoissons(lambdas,histogram):
    throws = []
    for lam in lambdas:
        while lam[lam<0].any():
            logging.error("Negative value in sys throws, rerunning this throw...")
            lam = ThrowSystematics(histogram,1) 
        throw = np.random.poisson(lam)
        throws.append(throw)
    throws = np.array(throws)
    return(throws)

def FitToyExperiments(histogram,experiments):
    stitched_data = histogram.GetDataHistogram()
    results = []
    for toy in experiments:
        weights = stitched_data.Clone().GetCVHistoWithStatError()
        for i in range(1,weights.GetNbinsX()+1):
            weight = stitched_data.GetBinContent(i) / toy[i-1] if toy[i-1] != 0 else stitched_data.GetBinContent(i)
            weights.SetBinContent(i,weight)
            weights.SetBinError(i,0)

        data_histogram = stitched_data.Clone()
        data_histogram.DivideSingle(data_histogram,weights)
        histogram.SetDataHistogram(data_histogram)

        chi2_fit, res = DoFit(histogram)
        chi2_mod = Chi2DataMC(histogram.GetDataHistogram(),histogram.GetMCHistogram(),histogram.GetDataCov(),histogram.GetMCCov())

        results.append(chi2_mod-chi2_fit)
        histogram.SetDataHistogram(stitched_data)

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

    filename = ""
    filename = "NuE_stitched_hists.root"

    file_path = "{}/FeldmanCousins/{}".format(ccnueroot,filename)

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    sys_throws  = ThrowSystematics(sample_histogram,50)
    experiments = ThrowPoissons(sys_throws,sample_histogram)

    dchi2s = FitToyExperiments(sample_histogram,experiments)
    np.savetxt('{}/sample_dchi2s_{}.csv'.format(outdir_surface,process),dchi2s,delimiter=',')
