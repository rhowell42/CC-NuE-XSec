import os
import time
import logging, sys
import ROOT
import PlotUtils
import numpy as np

minBinCont = 1

errorbandDict = { #keys are the errorbands that need to be renamed, values are what to rename them to
        "GENIE_D2_MaRES":"GENIE_MaRES",
        "GENIE_D2_NormCCRES":"GENIE_NormCCRES",
        "EnergyScale":"ElectronScale",
        "HCALScale":"elE_HCAL",
        "Other_Response":"response_other",
        "Pion_Response":"response_meson",
        "Proton_Response":"response_p",
        #"EM_Response":"response_em",
        "Crosstalk":"response_xtalk",
        "Numu Coh Scale 1":"numucoh1",
        "Numu Coh Scale 2":"numucoh2",
        "Numu Coh Scale 3":"numucoh3",
        "Numu Coh Scale 4":"numucoh4",
        "Numu Coh Scale 5":"numucoh5",
        "Numu Coh Scale 6":"numucoh6",
        "Numu Scale":"Numu CC Scale",
        "RPAHigh":"RPA_HighQ2",
        "RPALow":"RPA_LowQ2",
        "Target_Mass":"Target_Mass_CH",
        }

def SeparateBeamAngle(hlist):
    for h in hlist:
        hists = h.GetVertErrorBand("beam_angle").GetHists()
        h.PopVertErrorBand("beam_angle")
        xhists = hists[:2]
        yhists = hists[2:]
        h.AddVertErrorBand("BeamAngleX",xhists)
        h.AddVertErrorBand("BeamAngleY",yhists)

def LateralToVertical(hlist):
    for h in hlist:
        for name in h.GetLatErrorBandNames():
            hists = h.GetLatErrorBand(name).GetHists()
            h.PopLatErrorBand(name)
            h.AddVertErrorBand(name,hists)


def GetCovarianceMatrix(mnv_mc,mnv_data,ftag):
    NbinsTotal = mnv_mc.GetNbinsX()
    covMatrix = np.zeros(shape=[NbinsTotal,NbinsTotal],dtype='f')
    includeStatError = False
    errorAsFraction  = False
    useOnlyShapeErrors = False

    covMatrixTmp = mnv_mc.GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)
    for i in range(0,NbinsTotal):
        for j in range(0,NbinsTotal):
            covMatrix[i][j] = (covMatrixTmp[i+1][j+1])#*mnv_data.GetBinContent(i+1)*mnv_data.GetBinContent(j+1)
    np.savetxt("mc_covmatrix_"+ftag+".csv",covMatrix,delimiter=",")

    covMatrixTmp = mnv_data.GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)
    for i in range(0,NbinsTotal):
        for j in range(0,NbinsTotal):
            covMatrix[i][j] = (covMatrixTmp[i+1][j+1])#*mnv_data.GetBinContent(i+1)*mnv_data.GetBinContent(j+1)
    np.savetxt("data_covmatrix_"+ftag+".csv",covMatrix,delimiter=",")

    return(covMatrixTmp,covMatrix)

def StitchThis2D(h_new,h_olds,h_tests,sample=""):
    i_new = 0
    for h in h_olds:
        h_old = h_olds[h]
        h_test = h_tests[h]

        for x in range(1, h_test.GetNbinsX()+1):
            if h_test.GetBinContent(x) <= minBinCont:
                continue
            i_new += 1
            colInt = 0
            for c in range(1,h_old.GetNbinsX()+1):
                if isinstance(h_old,ROOT.TH1D):
                    colInt+=h_old.GetBinContent(c)
                else:
                    colInt+=h_old.GetBinContent(c,x)

            for c in range(1,h_old.GetNbinsX()+1):
                if isinstance(h_old,ROOT.TH1D):
                    bin_c = h_old.GetBinContent(c)
                else:
                    bin_c = h_old.GetBinContent(c,x)
                h_new.SetBinContent(i_new, c, bin_c/colInt if colInt > 0 else 0)

def FillErrorBandfromHist2(h_new,h_olds,mchists=None,isData=False,pseudodata=False):
    offset = 1
    for h in h_olds:
        h_old = h_olds[h]
        Nbins = h_old.GetNbinsX()+1 if not mchists else mchists[h].GetNbinsX()+1

        errorband_names_vert = h_old.GetVertErrorBandNames()
        n_univ = 0
        sys_bc = 0.0

        for error_band in errorband_names_vert:
            if error_band == 'Flux':
                n_univ = 100
            else:
                n_univ = h_old.GetVertErrorBand(error_band).GetNHists()
                if error_band == "GENIE_MaZExpCCQE":
                    n_univ = 100
                elif error_band == "GENIE_MaCCQE":
                    n_univ = 2

            if not h_new.HasVertErrorBand( error_band ) and h_old.HasVertErrorBand( error_band ):
                h_new.AddVertErrorBandAndFillWithCV( error_band, n_univ )
                h_new.GetVertErrorBand(error_band).SetUseSpreadError( h_old.GetVertErrorBand(error_band).GetUseSpreadError())

            for universe in range(0, n_univ):
                bin_offset = offset
                for b in range(1, Nbins):
                    bin_c = h_old.GetBinContent(b)
                    mcbin_c = mchists[h].GetBinContent(b)
                    sys_bc = h_old.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)
                    ratio = sys_bc/bin_c if bin_c != 0 else 0
                    if pseudodata:
                        sys_bc = mcbin_c*ratio

                    if mcbin_c <= minBinCont:
                        bin_offset += -1
                        continue
                    bin_new = b + bin_offset - 1
                    h_new.GetVertErrorBand(error_band).GetHist(universe).SetBinContent( bin_new, sys_bc )
        for i in range(1,Nbins):
            if mchists[h].GetBinContent(i) <= minBinCont:
                continue
            offset+=1

    for error in h_new.GetVertErrorBandNames():
        for n in range(h_new.GetVertErrorBand(error).GetNHists()):
            for i in range(h_new.GetNbinsX()+1):
                if h_new.GetBinContent(i) == 0:
                    h_new.GetVertErrorBand(error).SetBinContent(i,0)
                    h_new.GetVertErrorBand(error).GetHist(n).SetBinContent(i,0)

def StitchThis(h_new,h_olds,mchists=None,data=False,pseudodata=False):
    i_new = 0
    for h in h_olds:
        h_old = h_olds[h]
        for i in range(1,h_old.GetNbinsX()+1):
            bin_c = h_old.GetBinContent(i)
            if pseudodata:
                bin_c = mchists[h].GetBinContent(i)
            
            if mchists[h].GetBinContent(i) <= minBinCont:
                continue
            i_new += 1

            h_new.SetBinContent(i_new,bin_c)
            if data:
                continue

            # no statistical error on elastic scattering special production
            if 'nueel' in h:
                h_new.SetBinError(i_new,0)

def EmptyHist(h):
    for i in range(0,h.GetNbinsX()+1):
        h.SetBinContent(i,0)
        h.SetBinError(i,0)

def UnityHist(h):
    for i in range(0,h.GetNbinsX()+1):
        h.SetBinContent(i,1)
        h.SetBinError(i,0)

def SyncErrorBandsv2(hists):
    rename_bands = True
    if rename_bands:
        for h in hists:
            for name in h.GetVertErrorBandNames():
                if str(name) in errorbandDict.keys():
                    universes = h.GetVertErrorBand(name).GetHists()
                    useSpread = h.GetVertErrorBand(name).GetUseSpreadError()
                    h.PopVertErrorBand(name)
                    h.AddVertErrorBand(errorbandDict[str(name)],universes)
                    h.GetVertErrorBand(errorbandDict[str(name)]).SetUseSpreadError(useSpread)
                    for i in range(h.GetNbinsX()+1):
                        h.GetVertErrorBand(errorbandDict[str(name)]).SetBinContent(i,h.GetBinContent(i))

    for _h1 in range(len(hists)):
        for _h2 in range(len(hists)):
            h1 = hists[_h1]
            h2 = hists[_h2]
            if h1 == h2:
                continue

            for name in h1.GetVertErrorBandNames():
                if name == "Flux":
                    if h1.GetVertErrorBand(name).GetNHists() > 100:
                        h1_hists = h1.GetVertErrorBand(name).GetHists()
                        h1_hists = [h1_hists[i] for i in range(100)]
                        useSpread = h1.GetVertErrorBand(name).GetUseSpreadError()
                        h1.PopVertErrorBand(name)
                        h1.AddVertErrorBand(name,h1_hists)
                        h1.GetVertErrorBand(name).SetUseSpreadError(useSpread)
                        for i in range(h1.GetNbinsX()+1):
                            h1.GetVertErrorBand(name).SetBinContent(i,h1.GetBinContent(i))
                elif(name not in h2.GetVertErrorBandNames()):
                    n_universes = h1.GetVertErrorBand(name).GetNHists()
                    h2.AddVertErrorBandAndFillWithCV(name, n_universes)
                    
            for name in h2.GetVertErrorBandNames():
                if name == "Flux":
                    if h2.GetVertErrorBand(name).GetNHists() > 100:
                        h2_hists = h2.GetVertErrorBand(name).GetHists()
                        h2_hists = [h2_hists[i] for i in range(100)]
                        useSpread = h2.GetVertErrorBand(name).GetUseSpreadError()
                        h2.PopVertErrorBand(name)
                        h2.AddVertErrorBand(name,h2_hists)
                        h2.GetVertErrorBand(name).SetUseSpreadError(useSpread)
                        for i in range(h2.GetNbinsX()+1):
                            h2.GetVertErrorBand(name).SetBinContent(i,h2.GetBinContent(i))
                elif(name not in h1.GetVertErrorBandNames()):
                    n_universes = h2.GetVertErrorBand(name).GetNHists()
                    h1.AddVertErrorBandAndFillWithCV(name, n_universes)
