import os
import sys
import ROOT
from ROOT import *
import PlotUtils
from functools import partial
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities
from tools.PlotLibrary import HistHolder
from config.PlotConfig import ELECTRON_ENERGY_BINNING, UPSTREAM_INLINE_ENERGY_BINS
from PlotUtils.HistWrapper import HistWrapper
import math
import array

channel = "Electron Energy"
REGION = "Signal"
REGIONSig1 = "CCQElike"
REGIONSig2 = "notCCQElike"

fitPath = '/exp/minerva/data/users/rhowell/antinu_e/'

def CalcExcess(mc_hist,totSig_hist,data_hist,channel):
    data_hist = data_hist[channel]
    mc_hist = mc_hist[channel]
      
    excess_hist = data_hist.Clone()
    mctot_hist = mc_hist.Clone()      
   
    mc_hist.Add(totSig_hist,-1.0)     
    excess_hist.Add(mc_hist,-1.0)

    return excess_hist

def GetMCHist(channel,data_file,mc_file):
    hist = HistHolder(channel,mc_file,REGION,True)
    mc_hist = hist.GetHist().Clone()
    return mc_hist

def GetSignal1Hist(channel,data_file,mc_file):
    hist = HistHolder(channel,mc_file,REGIONSig1,True)
    mc_hist = hist.GetHist().Clone()
    return mc_hist

def GetSignal2Hist(channel,data_file,mc_file):
    hist = HistHolder(channel,mc_file,REGIONSig2,True)
    mc_hist = hist.GetHist().Clone()
    return mc_hist

def AddSignal(Sig1_hist,Sig2_hist):
    sig1_hist = Sig1_hist[channel]
    sig2_hist = Sig2_hist[channel]
    totSig = sig1_hist.Clone()
    totSig.Add(sig2_hist)  
    return totSig

def GetDataHist(channel,data_file,mc_file):
    hist = HistHolder(channel,data_file,REGION,True)
    data_hist = hist.GetHist().Clone()
    return data_hist

def GetShapeWeight(num,den):
    local_num = num.Clone()
    local_den = den.Clone()
    local_num.Divide(local_num,local_den)
    return local_num

def GetAvgWeight(weights):
    bins = weights.GetNbinsX()
    j=0
    for i in range(1,bins+1):
      j = j + weights.GetBinContent(i) 
    return j/bins

if __name__ == "__main__":

    data_file = ROOT.TFile.Open('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_dataNewSidebandPi0_nx_new_fspline.root')
    mc_file = ROOT.TFile.Open('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcNewSidebandPi0_nx_new_fspline.root')

    data_pot = Utilities.getPOTFromFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_dataNewSidebandPi0_nx_new_fspline.root')
    mc_pot = Utilities.getPOTFromFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcNewSidebandPi0_nx_new_fspline.root')

    mc_hists = dict(zip([channel],map(partial(GetMCHist,data_file=data_file, mc_file=mc_file), [channel])))
    Sig1_hist = dict(zip([channel],map(partial(GetSignal1Hist,data_file=data_file, mc_file=mc_file), [channel])))
    Sig2_hist = dict(zip([channel],map(partial(GetSignal2Hist,data_file=data_file, mc_file=mc_file), [channel])))
    data_hists = dict(zip([channel],map(partial(GetDataHist,data_file=data_file, mc_file=mc_file), [channel])))

    mc_hists[channel].Scale(data_pot/mc_pot)
    Sig1_hist[channel].Scale(data_pot/mc_pot)
    Sig2_hist[channel].Scale(data_pot/mc_pot)
    
    totSig = AddSignal(Sig1_hist,Sig2_hist)
    excess_hist = CalcExcess(mc_hists,totSig,data_hists,channel) 

    weights = GetShapeWeight(excess_hist,totSig)
    avgWeight = GetAvgWeight(weights) 
    print "This is the average weight for Signal(CCQElike,notCCQElike:)", avgWeight 
     
