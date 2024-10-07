import os
import sys
import ROOT
import PlotUtils
import math
import copy
from array import array
from collections import OrderedDict

from tools.PlotLibrary import HistHolder
from config.AnalysisConfig import AnalysisConfig
from config import BackgroundFitConfig
from tools import Utilities,PlotTools

ROOT.TH1.AddDirectory(False)

def ExtractTHXDFromMnvHXD(mnvHXD,errorband,i):
    pass

def PushErrorBandToMnvHXD(mnvHXD,hist_list,name):
    pass

def HistogramLinearlizer(reco_data,migration,mig_reco=None,mig_truth=None):
    #sanity check
    try:
        oneD = ValidateInput1D(reco_data,migration,mig_reco,mig_truth)
    except ValueError as e:
        oneD = False
    twoD = not oneD and ValidateInput2D(reco_data,migration,mig_reco,mig_truth)
    
    migration_tmp = EnlargeHistogram2D(iput) if oneD else migration
    mig_reco_tmp = LinearlizeHistogram(mig_reco) if mig_reco is not None else None
    mig_truth_tmp = LinearlizeHistogram(mig_truth) if mig_truth is not None else None
    reco_data_tmp = LinearlizeHistogram(reco_data)
    response = RooUnfold.RooUnfoldResponse(mig_reco,mig_truth,migration_tmp)
    return response,reco_data_tmp


def LinearlizeHistogram(iput):
    n= iput.GetSize()
    new_hist = ROOT.TH1D("{}Clone",format(iput.GetName()),iput.GetName(),n,0,n)
    for i in range(n):
        new_hist.SetBinContent(i+1,iput.GetBinContent(i))
        new_hist.SetBinError(i+1,iput.GetBinError(i))
    return new_hist

def EnlargeHistogram2D(iput):
    nx = iput.GetNbinX()+2
    ny = iput.GetNbinY()+2
    new_hist = ROOT.TH2D("{}Clone",format(iput.GetName()),iput.GetName(),nx,0,nx,ny,0,ny)
    for i in range(nx):
        for j in range(ny):
            new_hist.SetBinContent(i+1,j+1,iput.GetBinContent(i,j))
            new_hist.SetBinError(i+1,j+1,iput.GetBinError(i,j))
    return new_hist

def ValidateInput1D(reco_data,migration,mig_reco,mig_truth):
    if not reco_data:
        raise ValueError("reco data empty")
    if not migration:
        raise ValueError("migration empty")
    if not migration.GetNbinsX() == reco_data.GetNbinsX():
        raise ValueError("number of reco bins not match: migration and reco data")
    if mig_reco and migration.GetNbinsX() != mig_reco.GetNbinsX():
        raise ValueError("number of reco bins not match: migration and mig reco")
    if mig_truth and migration.GetNbinsY() != mig_truth.GetNbinsX():
        raise ValueError("number of reco bins not match: migration and mig truth")
    return True


def ValidateInput2D(reco_data,migration,mig_reco,mig_truth):
    if not reco_data:
        raise ValueError("reco data empty")
    if not migration:
        raise ValueError("migration empty")
    if not migration.GetNbinsX() == reco_data.GetSize():
        raise ValueError("number of reco bins not match: migration and reco data")
    if mig_reco and migration.GetNbinsX() != mig_reco.GetSize():
        raise ValueError("number of reco bins not match: migration and mig reco")
    if mig_truth and migration.GetNbinsY() != mig_truth.GetSize():
        raise ValueError("number of reco bins not match: migration and mig truth")
    return True

def Approach1(mig,mc_reco,mc_truth,fakedata_reco,fakedata_bkg,fakedata_truth):
    tmp =fakedata_reco.Clone()
    tmp.Add(fakedata_bkg,-1)
    options = DefaultOptions.copy()
    chi2s = []
    for i in range(NStat):
        toUnfold = StatShiftedReco(tmp)
        unfolded,cov = MyUnfoldWrapper(mig,mc_reco,mc_truth,toUnfold,options)
        DrawUnfoldingResult(unfolded,fakedata_truth,cov,"approach1{}".format(i))
        chi2 = CalcChi2(unfolded,fakedata_truth,cov)
        chi2s.append(chi2)
    return chi2s

def MyDelinerizer(hist_sample,hist_input):
    return_hist=hist_sample.Clone()
    nd = hist_sample.GetSize()
    for i in range(nd):
        return_hist.SetBinContent(i,hist_input.GetBinContent(i+1))
        return_hist.SetBinError(i,hist_input.GetBinError(i+1))
    return return_hist

def CalcChi2(h1,h2,cov):
    l1 = h1.GetSize()
    l2 = h2.GetSize()

    v1 = ROOT.TMatrixD(1,l1)
    v2 = ROOT.TMatrixD(l2,1)

    if not (l1 == l2 == cov.GetNcols() == cov.GetNrows()):
        print((l1,l2,cov.GetNcols(),cov.GetNrows()))
        raise ValueError("cov matrix/histograms dont have same dimention to calculate correlation.")

    for i in range(l1):
        v1[0][i] = h1.GetBinContent(i) - h2.GetBinContent(i)
        v2[i][0] = h1.GetBinContent(i) - h2.GetBinContent(i)

    SVD =ROOT.TDecompSVD(cov)
    tmp = SVD.Invert()

    tmp = tmp*v2
    tmp = v1*tmp
    return tmp[0][0]


def SubtractPoissonDist(h1,h2):
    for i in range(h1.GetSize()):
        mu = h1.GetBinContent(i)-h2.GetBinContent(i)
        sig = h1.GetBinContent(i)+h2.GetBinContent(i)
        h1.SetBinContent(i,mu)
        h1.SetBinError(i,math.sqrt(sig))

def StatShiftedReco(hist):
    hist2 = hist.Clone()
    # for i in range(hist.GetSize()):
    #     mu = RandomGenerator.Poisson(hist.GetBinContent(i))
    #     #print mu,hist.GetBinContent(i)
    #     hist2.SetBinContent(i,mu)
    #     hist2.SetBinError(i,math.sqrt(mu))
    return hist2


if __name__ == '__main__':
    f = ROOT.TFile.Open('/exp/minerva/data/users/hsu/nu_e/kin_dist_mcme1A-D_nx_col4_Sarah.root')
    # mig = ROOT.TH2D(f.Get("Eavail_migration"))
    # mc_reco = ROOT.TH1D(f.Get("Eavail"))
    # mc_truth = ROOT.TH1D(f.Get("Eavail_migration").ProjectionY())
    # fakedata_reco = ROOT.TH1D(f.Get("Eavail"))
    # fakedata_truth = ROOT.TH1D(f.Get("Eavail_migration").ProjectionY())
    # fakedata_signal = ROOT.TH1D(f.Get("Eavail_migration").ProjectionX())

    mig = ROOT.TH2D(f.Get("Eavail_q3_migration"))
    mc_reco = ROOT.TH2D(f.Get("Eavail_q3"))
    mc_truth = ROOT.TH2D(f.Get("Eavail_q3_migration_truth"))
    fakedata_reco = ROOT.TH2D(f.Get("Eavail_q3"))
    fakedata_truth = ROOT.TH2D(f.Get("Eavail_q3_migration_truth"))
    fakedata_signal = ROOT.TH2D(f.Get("Eavail_q3_migration_reco"))
    fakedata_bkg = fakedata_reco.Clone()
    fakedata_bkg.Add(fakedata_signal,-1)

    methods = [Approach1]
    for i in methods:
        chi2s = i(mig,mc_reco,mc_truth,fakedata_reco,fakedata_bkg,fakedata_truth)
        print(chi2s)
        print(sum(chi2s)/NStat)
