#!/usr/bin/env python
"""
Contains methods for getting weights from a event.
Includes flux, genie, minerva_genie_tunes(non-res-pion, rpa, and 2p2h)
"""
import os, ROOT, sys, uuid

from array import *
import math

import PlotUtils
from PlotUtils import FluxReweighter

fitPath = '{0}/data/Reweight'.format(os.environ['MPARAMFILESROOT'])
weight_cv_2p2h = PlotUtils.weight_2p2h('{0}/fit-mec-2d-noScaleDown-penalty00300-best-fit'.format(fitPath))
weight_nn_2p2h = PlotUtils.weight_2p2h('{0}/fit-mec-2d-nn-only-noScaleDown-penalty00300-best-fit'.format(fitPath))
weight_np_2p2h = PlotUtils.weight_2p2h('{0}/fit-mec-2d-np-only-noScaleDown-penalty02000-best-fit'.format(fitPath))
weight_qe_2p2h = PlotUtils.weight_2p2h('{0}/fit-qe-gaussian-noScaleDown-penalty02000-best-fit'.format(fitPath))
weight_cv_and_var_RPA = PlotUtils.weightRPA('{0}/outNievesRPAratio-nu12C-20GeV-20170202.root'.format(fitPath))
weight_LowQ2Pi = PlotUtils.weightLowQ2Pi()
flux_reweighter = None

# Minerva genie tune mod on Non resonace pion:
CV_nonres=0.43
deuteriumNonResNorm1sig = 0.04

def GetNonResonantPionWeight(event,var = 0):
    if var not in [-1,0,1]:
        return 1.0

    isGenieNonRes1pi = event.truth_genie_wgt_Rvn1pi[2] < 1.0 or event.truth_genie_wgt_Rvp1pi[2] < 1.0
    if isGenieNonRes1pi:
        return CV_nonres + deuteriumNonResNorm1sig * var
    else:
        return 1

def GetRPAWeight(event,var=0):

    ## Variations
    #0 CV value
    #1 hq2pos
    #2 hq2neg
    #3 lq2pos
    #4 lq2neg

    if not event.mc_intType == 1:  return 1 # Only calculate RPA weight if event is CCQE
    if event.mc_targetZ < 6:       return 1 # Intranuclear correlations don't occur in hydrogen

    q0 = event.kin_cal.true_q0
    q3 = event.kin_cal.true_q3

    if var == 0: return weight_cv_and_var_RPA.getWeight(q0,q3)
    if var == 1: return weight_cv_and_var_RPA.getWeightHighQ2(q0,q3,1)
    if var == 2: return weight_cv_and_var_RPA.getWeightHighQ2(q0,q3,-1)
    if var == 3: return weight_cv_and_var_RPA.getWeightLowQ2(q0,q3,1)
    if var == 4: return weight_cv_and_var_RPA.getWeightLowQ2(q0,q3,-1)

    print("I'm inside getRPAWeight, but I shouldn't have gotten this far!")
    return 

def Get2p2hWeight(event,var=0):
    ## Variations
    #0 CV value
    #1 nn+pp pairs only
    #2 np pair only
    #3 qe 1p1h variation
    #4 wgt = 1

    if var == 4: return 1.0

    # Only proceed to calculate 2p2h weight if event is CCQE (mc_intType == 1) or MEC (mc_intType == 8)
    if not (event.mc_intType == 1 or event.mc_intType == 8): return 1.0

    applyOn2p2h = True if (var == 0 or var == 1 or var ==2) else False
    applyOn1p1h = True if (var == 3) else False

    # Target analysis
    target = event.mc_targetNucleon
    isnnorpp = True if (target-2000000200 == 0 or target-2000000200 == 2) else False
    isnp = True if (target-2000000200 == 1) else False

    ## Handle the various permutations in which a weight shouldn't be applied
    # If CCQE and don't apply 1p1h, don't apply weights
    if event.mc_intType == 1 and not applyOn1p1h: return 1.0
    # If MEC and don't apply 2p2h, don't apply weights
    if event.mc_intType == 8 and not applyOn2p2h: return 1.0
    # If MEC and do apply 1p1h, don't apply weights
    if event.mc_intType == 8 and applyOn1p1h: return 1.0
    # Variation 1 is for nn/pp only interactions
    if var == 1 and not isnnorpp: return 1.0
    # Variation 2 is for np only interactions
    if var == 2 and not isnp: return 1.0
    
    q0 = event.kin_cal.true_q0
    q3 = event.kin_cal.true_q3

    #return self.weight_2p2h.getWeight(q0,q3)

    if var == 0: return weight_cv_2p2h.getWeight(q0,q3)
    if var == 1: return weight_nn_2p2h.getWeight(q0,q3)
    if var == 2: return weight_np_2p2h.getWeight(q0,q3)
    if var == 3: return weight_qe_2p2h.getWeight(q0,q3)

    print("I'm inside get2p2hWeight, but I shouldn't have gotten this far!")
    return


def GetLowQ2PiWeight(event,shift=0):
    #shift
    #0:CVweight
    #-1 : weaker weighting universe
    #+1 : stronger weighting universe
    #channel options:
    #MINOS, JOINT, NU1PI, NUNPI, NUPI0, NUBARPI0
    channel = "JOINT"
    is_ccres = event.mc_intType == 2 and event.mc_current == 1 # CC and RES 
    q2 = event.kin_cal.true_q2
    return weight_LowQ2Pi.getWeight(q2, channel, shift) if is_ccres else 1.0

def GetFluxWeight(event, i_universe = None):
    print("here!")
    exit()
    global flux_reweighter
    if flux_reweighter is None:
        flux_reweighter =  FluxReweighter(-12,False,
                                          PlaylistLookup(event.mc_run),
                                          FluxReweighter.gen2thin,
                                          FluxReweighter.g4numiv6)
    if i_universe is None:
        return flux_reweighter.GetFluxCVWeight(event.mc_incomingE/1e3,event.mc_incoming)
    else :
        return flux_reweighter.GetFluxErrorWeight(event.mc_incomingE/1e3, event.mc_incoming,i_universe)

def PlaylistLookup(run):
    myplist = FluxReweighter.minervame1A
    if run>=110000 and run<=110020:
        myplist = FluxReweighter.minervame1A
    elif run>=111000 and run<=111005:
        myplist = FluxReweighter.minervame1B
    elif run>=111030 and run<=111040:
        myplist = FluxReweighter.minervame1C
    elif run>=111100 and run<=111130:
        myplist = FluxReweighter.minervame1D
    elif run>=111325 and run<=111350:
        myplist = FluxReweighter.minervame1E
    elif run>=111490 and run<=111525:
        myplist = FluxReweighter.minervame1F
    elif run>=110150 and run<=110180:
        myplist = FluxReweighter.minervame1G
    elif run>=113000 and run<=113002:
        myplist = FluxReweighter.minervame1L
    elif run>=113020 and run<=113095:
        myplist = FluxReweighter.minervame1M
    elif run>=113270 and run<=113585:
        myplist = FluxReweighter.minervame1N
    elif run>=113375 and run<=113382:
        myplist = FluxReweighter.minervame1O
    elif run>=112000 and run<=112010:
        myplist = FluxReweighter.minervame1P

    return myplist
