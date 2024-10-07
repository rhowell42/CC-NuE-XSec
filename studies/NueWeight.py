#!/usr/bin/env python

import os, ROOT, sys, uuid

from tools import KinematicsCalculator
from tools import TruthTools
import PlotUtils
from PlotUtils import FluxReweighter
import random
from PlotUtils.HistWrapper import HistWrapper
from tools.PlotLibrary import HistHolder


#fitPath  = '{0}/data/Reweight'.format(os.environ['MPARAMFILESROOT'])
fitPath  = '/exp/minerva/data/users/rhowell/antinu_e'

file_FHC  =  'kin_dist_mcFHCEstimator_nx_new_fspline.root'.format(fitPath)
file_RHC = 'kin_dist_mcRHCEstimator_nx_new_fspline.root'.format(fitPath)


def GetFHCHist(): #get FHC estimator (Ee+Eavail)
    fFHC = ROOT.TFile.Open("{}/{}".format(fitPath,file_FHC))   
    histFHCHolder = HistHolder("True Neutrino Energy", fFHC, "Nue", True) 
    histFHC = histFHCHolder.GetHist().Clone() 
    histFHC.SetDirectory(0) 
    return histFHC
 
def GetFHCBarHist(): #get FHC estimator (Ee+Eavail)
    fFHC = ROOT.TFile.Open("{}/{}".format(fitPath,file_FHC))
    histFHCHolder = HistHolder("True Neutrino Energy", fFHC, "NueBar", True)
    histFHC = histFHCHolder.GetHist().Clone()
    histFHC.SetDirectory(0)
    return histFHC

def GetRHCHist(): #get RHC estimator (Ee+Eavail)
    fRHC = ROOT.TFile.Open("{}/{}".format(fitPath,file_RHC))  
    histRHCHolder = HistHolder("True Neutrino Energy", fRHC, "Nue", True)
    histRHC = histRHCHolder.GetHist().Clone()
    histRHC.SetDirectory(0)
    return histRHC

def GetRHCBarHist(): #get RHC estimator (Ee+Eavail)
    fRHC = ROOT.TFile.Open("{}/{}".format(fitPath,file_RHC))
    histRHCHolder = HistHolder("True Neutrino Energy", fRHC, "NueBar", True)
    histRHC = histRHCHolder.GetHist().Clone()
    histRHC.SetDirectory(0)
    return histRHC

def GetNueWeightHist(): #tEnu ratio RHC/FHC nue
    nueFHC = GetFHCHist()
    nueRHC = GetRHCHist()
    local_num = nueRHC.Clone()
    local_den = nueFHC.Clone()
    local_num.Divide(local_num,local_den) 
    return local_num

def GetNueBarWeightHist(): #tEnu ratio FHC/RHC nuebar
    nueFHC = GetFHCBarHist()
    nueRHC = GetRHCBarHist()
    local_den = nueRHC.Clone()
    local_num = nueFHC.Clone()
    local_num.Divide(local_num,local_den)
    return local_num


def GetNueWeight(event): #this should be applied to FHC events
    nueHist = GetNueWeightHist()
    bins = nueHist.GetNbinsX()
    nue = abs(event.mc_incoming) == 12
    Enu = event.mc_incomingE/1e3
    run = event.mc_run
    
    q0 = event.kin_cal.reco_q0
    visE = event.kin_cal.reco_visE
    Ee = event.kin_cal.reco_E_e 

    if run < 120000: #at this point we want to scale EVERYTHING including background 
    #if nue and run < 120000: #if FHC and nue event 
        if Ee is not None: #why am i getting SO many kin_cal None values
            estimator = Ee+visE 
            weight = nueHist.GetBinContent(nueHist.FindBin(estimator))
            #print Ee, estimator, weight, " Ee estimator weight" 
            return weight
        else:
            return 0.0
    else:
        return 0.0

def GetIsNue(event):
    nue = event.mc_incoming == 12
    #print event.mc_incoming
    if nue:
        return 0.0
    else:
        return 1.0

def GetNueBarWeight(event): #this should be applied to RHC events
    nuebarHist = GetNueBarWeightHist()
    bins = nuebarHist.GetNbinsX()
    nue = abs(event.mc_incoming) == 12
    Enu = event.mc_incomingE/1e3
    run = event.mc_run

    q0 = event.kin_cal.reco_q0
    visE = event.kin_cal.reco_visE
    Ee = event.kin_cal.reco_E_e

    if nue and run > 120000: #if RHC and nue event 
        if Ee is not None: #why am i getting SO many kin_cal None values
            estimator = Ee+visE
            weight = nuebarHist.GetBinContent(nuebarHist.FindBin(estimator)) 
            return weight
        else:
            return 0.0
    else:
        return 0.0

