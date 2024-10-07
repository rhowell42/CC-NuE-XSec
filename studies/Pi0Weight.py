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

channel = "Psi*Ee vs Ee"
REGION = "High_Psi_and_dEdX"
REGIONPi0 = "NCPi0"


fitPath = '/exp/minerva/data/users/rhowell/antinu_e/'

def CalcExcess(mc_hist,pi0_hist,data_hist,channel):
    data_hist = data_hist[channel]
    mc_hist = mc_hist[channel]
    pi0_hist = pi0_hist[channel]
 
    excess_hist = data_hist.Clone()
    mctot_hist = mc_hist.Clone()      
    
    mc_hist.Add(pi0_hist,-1.0) 
    excess_hist.Add(mc_hist,-1.0)
   
    return excess_hist

def GetMCHist(channel,data_file,mc_file):
    hist = HistHolder(channel,mc_file,REGION,True)
    mc_hist = hist.GetHist().Clone()
    return mc_hist

def GetPi0Hist(channel,data_file,mc_file):
    hist = HistHolder(channel,mc_file,REGIONPi0,True)
    mc_hist = hist.GetHist().Clone()
    return mc_hist

def GetDataHist(channel,data_file,mc_file):
    hist = HistHolder(channel,data_file,REGION,True)
    data_hist = hist.GetHist().Clone()
    return data_hist

def GetShapeWeight(num,den):
    local_num = num.Clone() 
    local_den = den.Clone()
   
    local_num.Divide(local_num,local_den)
  
    return local_num

def GetAvgWeight(excess,pi0):
    weights = []
    ex_tot = []
    pi0_tot = []
    
    binPsiEe = excess.GetNbinsX()
    binEe = excess.GetNbinsY()
    for i in range(binEe+1):
        pi0_count = 0
        excess_count = 0
        for j in range(binPsiEe+1): 
            ex = excess.GetBinContent(j, i)
            pi0b = pi0.GetBinContent(j, i)
            pi0_count = pi0_count + pi0b
            excess_count = excess_count + ex
        if pi0_count != 0:
            weight = excess_count/pi0_count
        else:
            weight = 0
        weights.append(weight)
        ex_tot.append(excess_count)
        pi0_tot.append(pi0_count)
    return weights


if __name__ == "__main__":

    data_file = ROOT.TFile.Open('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_datatest_nx_new_fspline.root')
    mc_file = ROOT.TFile.Open('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mctest_nx_new_fspline.root')

    data_pot = Utilities.getPOTFromFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_datatest_nx_new_fspline.root')
    mc_pot = Utilities.getPOTFromFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mctest_nx_new_fspline.root')

    mc_hists = dict(zip([channel],map(partial(GetMCHist,data_file=data_file, mc_file=mc_file), [channel])))
    pi0_hist = dict(zip([channel],map(partial(GetPi0Hist,data_file=data_file, mc_file=mc_file), [channel]))) 
    data_hists = dict(zip([channel],map(partial(GetDataHist,data_file=data_file, mc_file=mc_file), [channel])))

    mc_hists[channel].Scale(data_pot/mc_pot)
    pi0_hist[channel].Scale(data_pot/mc_pot)
    
    excess_hist = CalcExcess(mc_hists,pi0_hist,data_hists,channel)  
    weight_hist = GetShapeWeight(excess_hist,pi0_hist[channel])
    calcWeights = GetAvgWeight(excess_hist,pi0_hist[channel])
    #avgWeight = GetAvgWeight(weights)
    #print "This is the average weight for Signal(CCQElike,notCCQElike:)", avgWeight 
    print calcWeights     
