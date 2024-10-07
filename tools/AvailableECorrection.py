#!/usr/bin/env python

import os, ROOT, sys, uuid

from tools.KinematicsCalculator import KinematicsCalculator
import PlotUtils
from PlotUtils import FluxReweighter
import random
from tools.TruthTools import Rebin
from PlotUtils.HistWrapper import HistWrapper
from tools.PlotLibrary import HistHolder
from config.PlotConfig import LOW_RECOIL_BIN_Q0
from ROOT import *
from studies.drawVisECorr import fit_func

fitPath = "{0}/tools/".format(os.environ["CCNUEROOT"])
file_correction = 'outVisE_AuditSys.txt'.format(fitPath)


def MakeList(): #convert text file to list
    #a_file = open("/exp/minerva/app/users/shenry/cmtuser/Minerva_v22r1p1/Ana/CCNuE/macros/Low_Recoil/scripts/outVisE.txt","r")
    a_file = open("{}/{}".format(fitPath,file_correction))
    
    listtxt = a_file.readlines()
    list_of_corrections = [float(i) for i in listtxt]
    
    a_file.close()
    return list_of_corrections

def GetVisECorrection(event):
    list_of_corrections = MakeList() 
    visE = event.recoile_passive_tracker/1e3+event.recoile_passive_ecal/1e3
    if len(list_of_corrections) != (len(LOW_RECOIL_BIN_Q0)-1):
        print("WARNING: # of bins and length of list do not match")

    for i in range(0,len(list_of_corrections)):
        if LOW_RECOIL_BIN_Q0[i] < visE <= LOW_RECOIL_BIN_Q0[i+1]: 
            return list_of_corrections[i] #return proper correction
        
    return 0.0 #return no correction if outside of our scope

def GetSmoothVisECorrection(event):
    list_of_corrections = MakeList()
    xbins = LOW_RECOIL_BIN_Q0
    avg = []
    datlist = []
    for i in range(len(xbins)-1):
        v1 = xbins[i]
        v2 = xbins[i+1]
        v3 = (v1+v2)/2
        avg.append(v3)
    #for value in list_of_corrections:
    #    value = value.rstrip('\n')
    #    datlist.append(double(value))

    h2 = ROOT.TProfile('h2','',len(avg),0.0,1.2)
    #fit_func = TF1("f1", "pol3",len(avg))
    
    for i in range(0,len(avg)):
        h2.Fill(avg[i],list_of_corrections[i])


    visE = event.recoile_passive_tracker/1e3+event.recoile_passive_ecal/1e3
    if len(list_of_corrections) != (len(LOW_RECOIL_BIN_Q0)-1):
        print("WARNING: # of bins and length of list do not match")
    if visE <= max(LOW_RECOIL_BIN_Q0): 
        return fit_func.Eval(visE)
    else:
        return 0.0 
