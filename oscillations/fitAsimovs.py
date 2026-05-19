import os
import logging, sys
import ROOT
import PlotUtils
import numpy as np
np.set_printoptions(precision=1)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)
from scipy import optimize, integrate, linalg
from scipy.stats import multivariate_normal
import argparse
ccnueroot = os.environ.get('CCNUEROOT')
process = os.environ.get('PROCESS')

import math
from array import array

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from tools.Fitters import *
from tools.StitchedHistogram import *

from config.AnalysisConfig import AnalysisConfig

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

# Kevin test
# throw first from the flux covariance using cholesky decomp
# throw rest of the systematics around that result
# throw poisson from that result
# do rest as normal

def ThrowSystematics(histogram,throwFlux=False,useDataSubMCCov=True,n_samples=50):
    includeStatError = False
    errorAsFraction  = True
    useOnlyShapeErrors = False

    pred_vals = np.array(histogram.GetMCHistogram())[1:-1] # store MC bin contents excluding over/underflow bins

    if useDataSubMCCov:
        covMatrix = histogram.GetCovarianceMatrix()
    else:
        covMatrixTmp  = histogram.GetMCHistogram().GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)
        covMatrixTmp += histogram.GetPseudoHistogram().GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)

        mc_hist = histogram.GetMCHistogram()
        for i in range(mc_hist.GetNbinsX()+1):
            for j in range(i,mc_hist.GetNbinsX()+1):
                covMatrixTmp[i][j] = covMatrixTmp[i][j] * mc_hist.GetBinContent(i) * mc_hist.GetBinContent(j)
                covMatrixTmp[j][i] = covMatrixTmp[i][j]

        covMatrix = np.asarray(matrix(covMatrixTmp))[1:-1,1:-1] # convert to numpy array, exclude over/underflow bins

    if throwFlux:
        # Isolate the flux covariance matrix
        flux_cov_matrix = histogram.GetCovarianceMatrix() - histogram.GetCovarianceMatrix(sansFlux=True)

        # Get flux Cholesky factor (Lower triangular matrix 'L')
        L_flux = linalg.cholesky(flux_cov_matrix, lower=True)

        # Get the non-flux systematic covariance and its Cholesky factor
        sans_flux_matrix = histogram.GetCovarianceMatrix(sansFlux=True)
        L_sans = linalg.cholesky(sans_flux_matrix, lower=True)

        num_bins = len(pred_vals)

        # Generate standard normal random numbers for both steps
        z_flux = np.random.normal(0, 1, size=(num_bins, n_samples))
        z_sans = np.random.normal(0, 1, size=(num_bins, n_samples))

        # Correlated flux throws around the predicted mean
        flux_throws = pred_vals[:, np.newaxis] + (L_flux @ z_flux)

        # Correlated systematic throws around the flux throws
        sys_throws = flux_throws + (L_sans @ z_sans)

        # Ensure no expected bin counts are negative before Poisson sampling
        sys_throws = np.clip(sys_throws, a_min=0, a_max=None)

        # Transpose to get the standard shape: (n_samples, num_bins)
        sys_throws = sys_throws.T
    else:
        sys_throws = multivariate_normal.rvs(mean=pred_vals,cov=covMatrix,size=n_samples) # randomly sample covariance matrix around null hypothesis

    return(sys_throws)

def ThrowPoissons(lambdas,histogram):
    throws = []
    for lam in lambdas:
        while lam[lam<0].any():
            logging.error("Negative value in sys throws, rerunning this throw...")
            lam = ThrowSystematics(histogram, throwFlux = True, n_samples=1) 
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
        
        stat = Statistics(sample_histogram,lam=AnalysisConfig.lambdaValue,exclude=AnalysisConfig.exclude)

        chi2_null,penalty = stat.Chi2DataMC(marginalize=True)

        fitter = OscillationFitter(sample_histogram,lam=AnalysisConfig.lambdaValue,exclude=AnalysisConfig.exclude)
        chi2_fit,res = fitter.DoFit()

        dchi2 = chi2_null - chi2_fit
        results.append(dchi2)
        histogram.SetDataHistogram(stitched_data)

    results = np.array(results)
    return(results)

def SplitList(alist,wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] for i in range(wanted_parts) ]

if __name__ == "__main__":
    filename = "rootfiles/NuE_stitched_hists.root"

    file_path = "{}/oscillations/{}".format(ccnueroot,filename)

    sample_histogram = StitchedHistogram("sample")
    sample_histogram.Load(file_path)

    sys_throws  = ThrowSystematics(sample_histogram,throwFlux=True,n_samples=50)
    experiments = ThrowPoissons(sys_throws,sample_histogram)

    dchi2s = FitToyExperiments(sample_histogram,experiments)
    np.savetxt('{}/sample_dchi2s_{}.csv'.format(AnalysisConfig.output_dir,process),dchi2s,delimiter=',')
