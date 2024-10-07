import os
from collections import OrderedDict
import argparse
import logging, sys
import ROOT
import PlotUtils
from fit_tools.FitTools import *
from fit_tools.PlotTools import *
import numpy as np

import math
from array import array

from config.SignalDef import SWAP_SIGNAL_DEFINATION, SIGNAL_DEFINATION
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
from tools import Utilities
from tools.PlotLibrary import HistHolder

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
MNVPLOTTER.draw_normalized_to_bin_width=False
#MNVPLOTTER.extra_top_margin = -.033# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46 
MNVPLOTTER.mc_bkgd_line_color = 46
#MNVPLOTTER.mc_line_width = 0

MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.1
MNVPLOTTER.width_xspace_per_letter = .3
MNVPLOTTER.legend_text_size        = .07
MNVPLOTTER.legend_n_columns        = 1
MNVPLOTTER.axis_draw_grid_x = True
MNVPLOTTER.axis_draw_grid_y = True

MNVPLOTTER.error_summary_group_map.clear();
for k,v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v :
        vec.push_back(vs)
    MNVPLOTTER.error_summary_group_map[k]= vec

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

minBinCont = 1
def DrawStacked(TArray,mc_hist,data_hist,U_mu4,U_e4,m,title):
    c1 = ROOT.TCanvas()
    margin = .12
    bottomFraction = .2
    overall = ROOT.TCanvas("Data/MC")
    top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()

    mc_hist.GetXaxis().SetTitle("Bin number")
    mc_hist.GetYaxis().SetTitle("Entries")
    MNVPLOTTER.ApplyAxisStyle(mc_hist,True,True)

    MNVPLOTTER.DrawDataStackedMC(data_hist,TArray,1,"TR","Data",0,0,3015)

    if title == "FHC Elastic Scattering":
        null_hist = h_fhc_nueel_mc.Clone()
    elif title == "RHC Elastic Scattering":
        null_hist = h_rhc_nueel_mc.Clone()
    elif title == "RHC CC Muon Neutrino":
        null_hist = h_rhc_muselection_mc.Clone()
        data_hist.GetXaxis().SetTitle("E_{#nu} Estimator (GeV)")
    elif title == "RHC CC Electron Neutrino":
        null_hist = h_rhc_selection_mc.Clone()
        data_hist.GetXaxis().SetTitle("E_{#nu} Estimator (GeV)")

    ratio =  data_hist.Clone()
    nullRatio =  null_hist.Clone()
    oscRatio =  mc_hist.Clone()

    null_hist.SetLineColor(ROOT.kBlue)
    null_hist.SetLineWidth(3)
    null_hist.SetMarkerStyle(0)
    null_hist.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
    null_hist.GetYaxis().SetTitle("Nevents")
    null_err = null_hist.GetCVHistoWithError()
    null_err.Draw("E2 same")
    null_line = null_hist.Clone()
    null_line.SetFillColor(0)
    null_line.Draw('hist same')

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    osc = mc_hist.GetCVHistoWithError()
    osc.SetLineColor(ROOT.kRed)
    osc.SetLineWidth(3)
    osc.SetMarkerStyle(0)
    osc.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
    osc.Draw("E2 SAME")

    null = mc_hist.GetCVHistoWithError()
    null.SetLineColor(ROOT.kBlue)
    null.SetLineWidth(3)
    null.SetMarkerStyle(0)
    null.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
    null.Draw("E2 SAME")

    oscRatio.Divide(oscRatio, null_hist)
    nullRatio.Divide(nullRatio,null_hist)
    ratio.Divide(ratio, null_hist)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    nullErrors = null_hist.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
        nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
        nullErrors.SetBinContent(whichBin, 1)

    oscRatio.SetTitle("")
    oscRatio.SetLineColor(ROOT.kRed)
    nullRatio.SetLineColor(ROOT.kBlue)
    oscRatio.SetLineWidth(3)
    nullRatio.SetLineWidth(3)
    oscRatio.SetTitleSize(0)

    #Error envelope for the MC
    nullErrors.SetLineWidth(0)
    nullErrors.SetMarkerStyle(0)
    nullErrors.SetFillColorAlpha(ROOT.kBlue + 1, 0.4)

    nullErrors.GetYaxis().SetTitle("#splitline{Ratio to}{Null Hypothesis}")
    nullErrors.GetYaxis().SetLabelSize(.13)
    nullErrors.GetYaxis().SetTitleSize(0.1)
    nullErrors.GetYaxis().SetTitleOffset(0.6)
    nullErrors.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
    nullErrors.GetXaxis().SetTitleSize(0.16)
    nullErrors.GetXaxis().SetTitleOffset(0.9)
    nullErrors.GetXaxis().SetLabelSize(.15)
    nullErrors.GetXaxis().SetTitle("Bin Number")
    nullErrors.SetMinimum(0.5)
    nullErrors.SetMaximum(1.5)

    
    #Draw the data ratios
    oscRatio.SetMinimum(0.5)
    oscRatio.SetMaximum(1.5)
    nullRatio.SetMinimum(0.5)
    nullRatio.SetMaximum(1.5)
    ratio.SetLineWidth(3)
    ratio.SetMarkerSize(1)
    ratio.SetMarkerStyle(20)
    ratio.SetLineColor(ROOT.kBlack)

    #Draw a flat line at 1 for oscRatio of MC to itself
    straightLine = nullErrors.Clone()
    for whichBin in range(0, straightLine.GetXaxis().GetNbins()+1): 
        straightLine.SetBinError(whichBin, 0)
        oscRatio.SetBinError(whichBin,0)
        if oscRatio.GetBinContent(whichBin) == 0:
            oscRatio.SetBinContent(whichBin,1)
        if straightLine.GetBinContent(whichBin) == 0:
            straightLine.SetBinContent(whichBin,1)

    nullErrors.Draw("E2")
    straightLine.SetLineColor(ROOT.kBlue)
    straightLine.SetLineWidth(3)
    straightLine.SetFillColor(0)
    straightLine.Draw("SAME HIST L")
    ratio.Draw("same")
    oscRatio.SetFillColor(0)
    oscRatio.Draw("hist l same")

    overall.Print("fit_m-{}_ue4-{}_umu4-{}_{}_thesis.png".format(m,U_e4,U_mu4,title))

def MakeOscillationPlot(bestParams,mc_hist,data_hist,title):
    if title == "FHC Elastic Scattering":
        numutemplate = h2_fhc_nueel_template_anumu.Clone()
        numutemplate.Add(h2_fhc_nueel_template_numu.Clone())
        nuetemplate = h2_fhc_nueel_template_anue.Clone()
        nuetemplate.Add(h2_fhc_nueel_template_nue.Clone()) 
        nueHist = h_fhc_nueel_mc.Clone()
        swapHist = h_fhc_nueel_mc.Clone()
        numuHist = h_fhc_nueel_mc.Clone()
        nutauHist = h_fhc_nueel_mc.Clone()
        for i in range(h_fhc_muselection_mc.GetNbinsX()+1):
            nueHist.SetBinContent(i,0)
            swapHist.SetBinContent(i,0)
            numuHist.SetBinContent(i,0)
            nutauHist.SetBinContent(i,0)
    elif title == "RHC Elastic Scattering":
        numutemplate = h2_rhc_nueel_template_anumu.Clone()
        numutemplate.Add(h2_rhc_nueel_template_numu.Clone())
        nuetemplate = h2_rhc_nueel_template_anue.Clone()
        nuetemplate.Add(h2_rhc_nueel_template_nue.Clone()) 
        nueHist = h_fhc_nueel_mc.Clone()
        swapHist = h_fhc_nueel_mc.Clone()
        numuHist = h_fhc_nueel_mc.Clone()
        nutauHist = h_fhc_nueel_mc.Clone()
        for i in range(h_fhc_muselection_mc.GetNbinsX()+1):
            nueHist.SetBinContent(i,0)
            swapHist.SetBinContent(i,0)
            numuHist.SetBinContent(i,0)
            nutauHist.SetBinContent(i,0)
    elif title == "RHC CC Muon Neutrino":
        numutemplate = h2_rhc_musel_template.Clone()
        nuetemplate = h2_rhc_musel_template.Clone()
        nueHist = h_fhc_selection_mc.Clone()
        swapHist = h_fhc_nueel_mc.Clone()
        numuHist = h_fhc_selection_mc.Clone()
        nutauHist = h_fhc_selection_mc.Clone()
        for i in range(h_fhc_muselection_mc.GetNbinsX()+1):
            nueHist.SetBinContent(i,0)
            swapHist.SetBinContent(i,0)
            numuHist.SetBinContent(i,0)
            nutauHist.SetBinContent(i,0)
    elif title == "RHC CC Electron Neutrino":
        nuetemplate = h2_rhc_template.Clone()
        swaptemplate = h2_rhc_swap_template.Clone()
        nueHist = h_fhc_selection_mc.Clone()
        swapHist = h_fhc_nueel_mc.Clone()
        numuHist = h_fhc_selection_mc.Clone()
        nutauHist = h_fhc_selection_mc.Clone()
        for i in range(h_fhc_muselection_mc.GetNbinsX()+1):
            nueHist.SetBinContent(i,0)
            swapHist.SetBinContent(i,0)
            numuHist.SetBinContent(i,0)
            nutauHist.SetBinContent(i,0)

    for i in range(len(bestParams["m"])):
        m = bestParams["m"][i]
        U_e4 = bestParams["ue4"][i]
        U_mu4 = bestParams["umu4"][i]
        U_tau4 = bestParams["utau4"][i]
        if title == "FHC Elastic Scattering":
            for q in range(mc_hist.GetNbinsX()+1):
                nue_sin = sin_average(q,m,nuetemplate)
                numu_sin = sin_average(q,m,numutemplate)

                P_ee = 1 - 4*U_e4*(1-U_e4)*nue_sin
                eeCont = P_ee * (fhcnueelnue.GetBinContent(q)+fhcnueelanumu.GetBinContent(q))
                P_mue = 4*U_e4*U_mu4*numu_sin
                mueCont = P_mue * (fhcnueelanumu.GetBinContent(q)+fhcnueelnumu.GetBinContent(q))
                P_mumu = 1 - 4*U_mu4*(1-U_mu4)*numu_sin
                mumuCont = P_mumu * (fhcnueelnumu.GetBinContent(q)+fhcnueelanumu.GetBinContent(q))
                P_mutau = 4*U_tau4*U_mu4*numu_sin
                mutauCont = P_mutau * (fhcnueelanumu.GetBinContent(q)+fhcnueelnumu.GetBinContent(q))
                P_etau = 4*U_e4*U_tau4*nue_sin
                etauCont = P_etau * (fhcnueelnue.GetBinContent(q)+fhcnueelanue.GetBinContent(q))

                nueHist.SetBinContent(q,eeCont)
                numuHist.SetBinContent(q,mumuCont)
                swapHist.SetBinContent(q,mueCont)
                nutauHist.SetBinContent(q,etauCont+mutauCont)
            mc_hists = [nueHist,numuHist,swapHist,nutauHist]
            titles = ["remaining #nu_{e}","remaining #nu_{#mu}","appeared #nu_{e}","appeared #nu_{#tau}"]
            color = [ROOT.kRed,ROOT.kBlue,ROOT.kTeal,ROOT.kGray]
            TArray = ROOT.TObjArray()
            for i in range(len(mc_hists)):
                if color is not None:
                    mc_hists[i].SetFillColor(color[i])
                if titles is not None:
                    mc_hists[i].SetTitle(titles[i])
                TArray.Add(mc_hists[i])
            mc_hist = nueHist.Clone()
            mc_hist.Add(numuHist)
            mc_hist.Add(nutauHist)
            DrawStacked(TArray,mc_hist,data_hist,U_mu4,U_e4,m,title)
        elif title == "RHC Elastic Scattering":
            for q in range(mc_hist.GetNbinsX()+1):
                nue_sin = sin_average(q,m,nuetemplate)
                numu_sin = sin_average(q,m,numutemplate)

                P_ee = 1 - 4*U_e4*(1-U_e4)*nue_sin
                eeCont = P_ee * (rhcnueelnue.GetBinContent(q)+rhcnueelanue.GetBinContent(q))
                P_mue = 4*U_e4*U_mu4*numu_sin
                mueCont = P_mue * (rhcnueelanumu.GetBinContent(q)+rhcnueelnumu.GetBinContent(q))
                P_mumu = 1 - 4*U_mu4*(1-U_mu4)*numu_sin
                mumuCont = P_mumu * (rhcnueelnumu.GetBinContent(q)+rhcnueelanumu.GetBinContent(q))
                P_mutau = 4*U_tau4*U_mu4*numu_sin
                mutauCont = P_mutau * (rhcnueelanumu.GetBinContent(q)+rhcnueelnumu.GetBinContent(q))
                P_etau = 4*U_e4*U_tau4*nue_sin
                etauCont = P_etau * (rhcnueelnue.GetBinContent(q)+rhcnueelanue.GetBinContent(q))

                nueHist.SetBinContent(q,eeCont)
                numuHist.SetBinContent(q,mumuCont)
                swapHist.SetBinContent(q,mueCont)
                nutauHist.SetBinContent(q,etauCont+mutauCont)
            mc_hists = [nueHist,numuHist,swapHist,nutauHist]
            titles = ["remaining #nu_{e}","remaining #nu_{#mu}","appeared #nu_{e}","appeared #nu_{#tau}"]
            color = [ROOT.kRed,ROOT.kBlue,ROOT.kTeal,ROOT.kGray]
            TArray = ROOT.TObjArray()
            for i in range(len(mc_hists)):
                if color is not None:
                    mc_hists[i].SetFillColor(color[i])
                if titles is not None:
                    mc_hists[i].SetTitle(titles[i])
                TArray.Add(mc_hists[i])
            mc_hist = nueHist.Clone()
            mc_hist.Add(numuHist)
            mc_hist.Add(nutauHist)
            DrawStacked(TArray,mc_hist,data_hist,U_mu4,U_e4,m,title)
        elif title == "RHC CC Muon Neutrino":
            for q in range(mc_hist.GetNbinsX()+1):
                numu_sin = sin_average(q,m,numutemplate)

                P_mumu = 1 - 4*U_mu4*(1-U_mu4)*numu_sin
                mumuCont = P_mumu * (h_rhc_muselection_mc.GetBinContent(q))
                numuHist.SetBinContent(q,mumuCont)

            mc_hists = [numuHist]
            titles = ["remaining #nu_{#mu}"]
            color = [ROOT.kBlue]
            TArray = ROOT.TObjArray()
            for i in range(len(mc_hists)):
                if color is not None:
                    mc_hists[i].SetFillColor(color[i])
                if title is not None:
                    mc_hists[i].SetTitle(titles[i])
                TArray.Add(mc_hists[i])
            mc_hist = numuHist.Clone()
            DrawStacked(TArray,mc_hist,data_hist,U_mu4,U_e4,m,title)
        elif title == "RHC CC Electron Neutrino":
            for q in range(mc_hist.GetNbinsX()+1):
                nue_sin = sin_average(q,m,nuetemplate)
                numu_sin = sin_average(q,m,swaptemplate)

                P_ee = 1 - 4*U_e4*(1-U_e4)*nue_sin
                eeCont = P_ee * (h_rhc_selection_mc.GetBinContent(q))
                P_mue = 4*U_e4*U_mu4*numu_sin
                mueCont = P_mue * (h_rhc_swap_selection.GetBinContent(q))

                nueHist.SetBinContent(q,eeCont)
                numuHist.SetBinContent(q,mueCont)
            mc_hists = [nueHist,numuHist]
            titles = ["remaining #nu_{e}","appeared #nu_{e}"]
            color = [ROOT.kRed,ROOT.kTeal]
            TArray = ROOT.TObjArray()
            for i in range(len(mc_hists)):
                if color is not None:
                    mc_hists[i].SetFillColor(color[i])
                if titles is not None:
                    mc_hists[i].SetTitle(titles[i])
                TArray.Add(mc_hists[i])
            mc_hist = numuHist.Clone()
            mc_hist.Add(nueHist)
            DrawStacked(TArray,mc_hist,data_hist,U_mu4,U_e4,m,title)

def ReSumHists(h_list):
    if len(h_list) == 0:
        return None
    hnew = h_list[0].Clone()
    for i in range(1,len(h_list)):
        hnew.Add(h_list[i])
    return hnew

def SyncErrorBandsv2(hists):
    for h1 in hists:
        for h2 in hists:
            if h1 == h2:
                continue
            for name in h1.GetVertErrorBandNames():
                if name == "Flux":
                    universes = 100
                elif(name not in h2.GetVertErrorBandNames()):
                    universes = h1.GetVertErrorBand(name).GetNHists()
                    h2.AddVertErrorBandAndFillWithCV(name, universes)

            for name in h1.GetLatErrorBandNames():
                if name == "Flux":
                    universes = 100
                elif(name not in h2.GetLatErrorBandNames()):
                    universes = h1.GetLatErrorBand(name).GetNHists()
                    h2.AddLatErrorBandAndFillWithCV(name, universes)

            for name in h2.GetVertErrorBandNames():
                if name == "Flux":
                    universes = 100
                elif(name not in h1.GetVertErrorBandNames()):
                    universes = h2.GetVertErrorBand(name).GetNHists()
                    h1.AddVertErrorBandAndFillWithCV(name, universes)

            for name in h2.GetLatErrorBandNames():
                if name == "Flux":
                    universes = 100
                elif(name not in h1.GetLatErrorBandNames()):
                    universes = h2.GetLatErrorBand(name).GetNHists()
                    h1.AddLatErrorBandAndFillWithCV(name, universes)

def FillErrorBandfromHist(h_new,h_olds,mchists = None,pseudodata = False, isData=False,sample=""):
    offset = 1
    for h in h_olds:
        h_old = h_olds[h]
        Nbins = h_old.GetNbinsX()+1 if not mchists else mchists[h].GetNbinsX()+1

        if sample == 'nue':
            if 'muselection' in h or 'imd' in h:
                for i in range(1,Nbins):
                    if mchists[h].GetBinContent(i) <= minBinCont:
                        continue
                    offset+=1
                continue
            if "nueel" in h and "fhc" in h:
                h_old = fhcnueelnue.Clone()
                h_old.Add(fhcnueelanue)
            elif "nueel" in h and "rhc" in h:
                h_old = rhcnueelnue.Clone()
                h_old.Add(rhcnueelanue)
            elif "ratio" in h and "fhc" in h:
                h_old = h_fhc_selection_mc
            elif "ratio" in h and "rhc" in h:
                h_old = h_rhc_selection_mc

        elif sample == 'numu':
            if '_selection' in h:
                for i in range(1,Nbins):
                    if mchists[h].GetBinContent(i) <= minBinCont:
                        continue
                    offset+=1
                continue
            if "nueel" in h and "fhc" in h:
                h_old = fhcnueelnumu.Clone()
                h_old.Add(fhcnueelanumu)
            elif "nueel" in h and "rhc" in h:
                h_old = rhcnueelnumu.Clone()
                h_old.Add(rhcnueelanumu)
            elif "ratio" in h and "fhc" in h:
                h_old = h_fhc_muselection_mc
            elif "ratio" in h and "rhc" in h:
                h_old = h_rhc_muselection_mc

        elif sample == 'swap':
            if 'fhc_selection' == h or 'fhc_ratio' == h:
                for i in range(1,Nbins):
                    if mchists[h].GetBinContent(i) <= minBinCont:
                        continue
                h_old = h_fhc_swap_selection.Clone()
            elif 'rhc_selection' == h or 'rhc_ratio' == h:
                for i in range(1,Nbins):
                    if mchists[h].GetBinContent(i) <= minBinCont:
                        continue
                h_old = h_rhc_swap_selection.Clone()
            else:
                for i in range(1,Nbins):
                    if mchists[h].GetBinContent(i) <= minBinCont:
                        continue
                    offset+=1
                continue
        else:
            continue

        errorband_names_vert = h_old.GetVertErrorBandNames()
        errorband_names_lat = h_old.GetLatErrorBandNames()
        n_univ = 0
        sys_bc = 0.0

        for error_band in errorband_names_vert:
            #print("looping over vert error ", error_band)
            if error_band == 'Flux' and isData:            #hacked for the flux since fhc files have 100 and rhc have 1000 universes
                n_univ = 100
            elif error_band == 'Flux':
                n_univ = 100
            else:
                n_univ = h_old.GetVertErrorBand(error_band).GetNHists()
                if error_band == "GENIE_MaZExpCCQE":
                    n_univ = 100
                elif error_band == "GENIE_MaCCQE":
                    n_univ = 2

            if error_band == "Muon_Energy_MINERvA" or error_band == "Muon_Energy_MINOS" or error_band == "Muon_Energy_Resolution":
                new_band = "v_{}".format(error_band)
            else:
                new_band = error_band

            if not h_new.HasVertErrorBand( new_band ) and h_old.HasVertErrorBand( error_band ):
                h_new.AddVertErrorBandAndFillWithCV( new_band, n_univ )
                h_new.GetVertErrorBand(new_band).SetUseSpreadError( h_old.GetVertErrorBand(error_band).GetUseSpreadError())

            for universe in range(0, n_univ):
                bin_offset = offset
                for b in range(1, Nbins):
                    bin_c = h_old.GetBinContent(b)
                    mcbin_c = mchists[h].GetBinContent(b)
                    #mcbin_c = h_old.GetBinContent(b)
                    sys_bc = h_old.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)
                    ratio = sys_bc/bin_c if bin_c != 0 else 0
                    if mcbin_c <= minBinCont:
                        bin_offset += -1
                        continue
                    bin_new = b + bin_offset - 1
                    if pseudodata:
                        sys_bc = mcbin_c*ratio
                    h_new.GetVertErrorBand(new_band).GetHist(universe).SetBinContent( bin_new, sys_bc )

        for error_band in errorband_names_lat:
            if error_band == "Muon_Energy_MINERvA" or error_band == "Muon_Energy_MINOS" or error_band == "Muon_Energy_Resolution":
                new_band = "h_{}".format(error_band)
            else:
                new_band = error_band
            #print("looping over lat error ", error_band)
            if not h_new.HasLatErrorBand(new_band) and h_old.HasLatErrorBand(error_band):
                n_univ = h_old.GetLatErrorBand(error_band).GetNHists()
                #print 'old hist has errorband ', error_band, ' with ', n_univ, ' universes'
            else:
                n_univ = h_old.GetLatErrorBand(error_band).GetNHists()

            if not h_new.HasLatErrorBand( new_band ):
                h_new.AddLatErrorBandAndFillWithCV( new_band , n_univ )
                h_new.GetLatErrorBand(new_band).SetUseSpreadError( h_old.GetLatErrorBand(error_band).GetUseSpreadError())
            for universe in range(0, n_univ):
                bin_offset = offset
                for b in range(1, Nbins):  #bin 0 and 1 have no entries. Bin seven is the last bin
                    bin_c = h_old.GetBinContent(b)
                    mcbin_c = mchists[h].GetBinContent(b)
                    #mcbin_c = h_old.GetBinContent(b)
                    sys_bc = h_old.GetLatErrorBand(error_band).GetHist(universe).GetBinContent(b)
                    ratio = sys_bc/bin_c if bin_c !=0 else 0
                    if mcbin_c <= minBinCont:
                        bin_offset += -1
                        continue
                    bin_new = b  + bin_offset - 1
                    if pseudodata:
                        sys_bc = mcbin_c*ratio
                    h_new.GetLatErrorBand(new_band).GetHist(universe).SetBinContent( bin_new, sys_bc )

        for i in range(1,Nbins):
            if mchists[h].GetBinContent(i) <= minBinCont:
                continue
            offset+=1

def FillErrorBandfromHist2(h_new,h_olds,mchists = None, pseudodata = False, isData=False):
    offset = 1
    for h in h_olds:
        h_old = h_olds[h]
        Nbins = h_old.GetNbinsX()+1 if not mchists else mchists[h].GetNbinsX()+1

        errorband_names_vert = h_old.GetVertErrorBandNames()
        errorband_names_lat = h_old.GetLatErrorBandNames()
        n_univ = 0
        sys_bc = 0.0

        for error_band in errorband_names_vert:
            #print("looping over vert error ", error_band)
            if error_band == 'Flux' and isData:# and h_old.GetVertErrorBand(error_band).GetNHists() < 1000:            #hacked for the flux since fhc files have 100 and rhc have 1000 universes
                n_univ = 100
            elif error_band == 'Flux':
                n_univ = 100
            else:
                n_univ = h_old.GetVertErrorBand(error_band).GetNHists()
                if error_band == "GENIE_MaZExpCCQE":
                    n_univ = 100
                elif error_band == "GENIE_MaCCQE":
                    n_univ = 2

            if error_band == "Muon_Energy_MINERvA" or error_band == "Muon_Energy_MINOS" or error_band == "Muon_Energy_Resolution":
                new_band = "v_{}".format(error_band)
            else:
                new_band = error_band

            if not h_new.HasVertErrorBand( new_band ) and h_old.HasVertErrorBand( error_band ):
                h_new.AddVertErrorBandAndFillWithCV( new_band, n_univ )
                h_new.GetVertErrorBand(new_band).SetUseSpreadError( h_old.GetVertErrorBand(error_band).GetUseSpreadError())

            for universe in range(0, n_univ):
                bin_offset = offset
                for b in range(1, Nbins):
                    bin_c = h_old.GetBinContent(b)
                    mcbin_c = mchists[h].GetBinContent(b)
                    sys_bc = h_old.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)
                    ratio = sys_bc/bin_c if bin_c != 0 else 0
                    if mcbin_c <= minBinCont:
                        bin_offset += -1
                        continue
                    bin_new = b + bin_offset - 1
                    if pseudodata:
                        sys_bc = mcbin_c*ratio
                    h_new.GetVertErrorBand(new_band).GetHist(universe).SetBinContent( bin_new, sys_bc )

        for error_band in errorband_names_lat:
            if error_band == "Muon_Energy_MINERvA" or error_band == "Muon_Energy_MINOS" or error_band == "Muon_Energy_Resolution":
                new_band = "h_{}".format(error_band)
            else:
                new_band = error_band
            #print("looping over lat error ", error_band)
            if not h_new.HasLatErrorBand(new_band) and h_old.HasLatErrorBand(error_band):
                n_univ = h_old.GetLatErrorBand(error_band).GetNHists()
                #print 'old hist has errorband ', error_band, ' with ', n_univ, ' universes'
            else:
                continue
            if not h_new.HasLatErrorBand( new_band ):
                h_new.AddLatErrorBandAndFillWithCV( new_band , n_univ )
                h_new.GetLatErrorBand(new_band).SetUseSpreadError( h_old.GetLatErrorBand(error_band).GetUseSpreadError())
            for universe in range(0, n_univ):
                bin_offset = offset
                for b in range(1, Nbins):  #bin 0 and 1 have no entries. Bin seven is the last bin
                    bin_c = h_old.GetBinContent(b)
                    mcbin_c = mchists[h].GetBinContent(b)
                    sys_bc = h_old.GetLatErrorBand(error_band).GetHist(universe).GetBinContent(b)
                    ratio = sys_bc/bin_c if bin_c !=0 else 0
                    if mcbin_c <= minBinCont:
                        bin_offset += -1
                        continue
                    bin_new = b  + bin_offset - 1
                    if pseudodata:
                        sys_bc = mcbin_c*ratio
                    h_new.GetLatErrorBand(new_band).GetHist(universe).SetBinContent( bin_new, sys_bc )
        for i in range(1,Nbins):
            if mchists[h].GetBinContent(i) <= minBinCont:
                continue
            offset+=1

def StitchThis(h_new,h_olds,sample="",mchists=None,pseudodata=False):
    i_new = 0
    for h in h_olds:
        h_old = h_olds[h]
        for i in range(1,h_old.GetNbinsX()+1):
            if pseudodata:
                bin_c = mchists[h].GetBinContent(i)
            else:
                bin_c = h_old.GetBinContent(i)
            
            if mchists[h].GetBinContent(i) <= minBinCont:
                continue
            i_new += 1

            h_new.SetBinContent(i_new,bin_c)
            if sample=="data":
                continue

            # no statistical error on elastic scattering special production
            if 'nueel' in h:
                h_new.SetBinError(i_new,0)

            # handle nue vs numu content of elastic scattering histograms
            if h == "fhc_nueel":
                tot = fhcnueelnumu.GetBinContent(i)+fhcnueelanumu.GetBinContent(i)+fhcnueelnue.GetBinContent(i)+fhcnueelanue.GetBinContent(i)
                fraction = 1
                if sample == "numu":
                    fraction = (fhcnueelanumu.GetBinContent(i)+fhcnueelnumu.GetBinContent(i))/tot
                elif sample == "nue":
                    fraction = (fhcnueelanue.GetBinContent(i)+fhcnueelnue.GetBinContent(i))/tot
                h_new.SetBinContent(i_new,bin_c*fraction)
            elif h == "rhc_nueel":
                tot = rhcnueelnumu.GetBinContent(i)+rhcnueelanumu.GetBinContent(i)+rhcnueelnue.GetBinContent(i)+rhcnueelanue.GetBinContent(i)
                fraction = 1
                if sample == "numu":
                    fraction = (rhcnueelanumu.GetBinContent(i)+rhcnueelnumu.GetBinContent(i))/tot
                elif sample == "nue":
                    fraction = (rhcnueelanue.GetBinContent(i)+rhcnueelnue.GetBinContent(i))/tot
                h_new.SetBinContent(i_new,bin_c*fraction)

            # handle nue vs numu content of selection histograms
            if sample == "nue" and ("_muselection" in h or "imd" in h):
                h_new.SetBinContent(i_new,0)
            elif sample == "nue" and "rhc" in h and ("_selection" in h or "_ratio" in h):
                h_new.SetBinContent(i_new,h_rhc_selection_mc.GetBinContent(i))
            elif sample == "nue" and "fhc" in h and ("_selection" in h or "_ratio" in h):
                h_new.SetBinContent(i_new,h_fhc_selection_mc.GetBinContent(i))

            if sample == "numu" and "_selection" in h:
                h_new.SetBinContent(i_new,0)
            elif sample == "numu" and "rhc" in h and ("_muselection" in h or "_ratio" in h):
                h_new.SetBinContent(i_new,h_rhc_muselection_mc.GetBinContent(i))
            elif sample == "numu" and "fhc" in h and ("_muselection" in h or "_ratio" in h):
                h_new.SetBinContent(i_new,h_fhc_muselection_mc.GetBinContent(i))

            # handle ID histogram for elastic scattering histograms
            if sample == "nutau" and "nueel" in h:
                h_new.SetBinContent(i_new,1)
            elif sample == "nutau":
                h_new.SetBinContent(i_new,0)

            # handle ID histogram from CC nue selection histograms
            if sample == "nueselection" and "_selection" in h:
                h_new.SetBinContent(i_new,0)
            elif sample == "nueselection":
                h_new.SetBinContent(i_new,1)

            if sample == 'ratio' and 'ratio' in h:
                h_new.SetBinContent(i_new,1)
            elif sample == 'ratio':
                h_new.SetBinContent(i_new,0)

            if sample == 'swap' and 'rhc_ratio' in h:
                h_new.SetBinContent(i_new,h_rhc_swap_selection.GetBinContent(i))
            elif sample == 'swap' and 'fhc_ratio' in h:
                h_new.SetBinContent(i_new,h_fhc_swap_selection.GetBinContent(i))
            elif sample == 'swap' and 'fhc_selection' in h:
                h_new.SetBinContent(i_new,h_fhc_swap_selection.GetBinContent(i))
            elif sample == 'swap' and 'rhc_selection' in h:
                h_new.SetBinContent(i_new,h_rhc_swap_selection.GetBinContent(i))
            elif sample == 'swap':
                h_new.SetBinContent(i_new,0)

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

def GetCovarianceMatrix(mnv_mc,mnv_data):
    NbinsTotal = mnv_mc.GetNbinsX()
    covMatrix = np.zeros(shape=[NbinsTotal,NbinsTotal],dtype='f')
    includeStatError = True
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

def MakeThesisPlot(mc_hist,data_hist,title):
    c1 = ROOT.TCanvas()
    margin = .12
    bottomFraction = .2
    overall = ROOT.TCanvas(title)
    top = ROOT.TPad(title+"DATAMC", title+"DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad(title+"Ratio", title+"Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()

    mc_hist.SetTitle(title)
    data_hist.SetTitle(title)
    
    if "Elastic" in title:
        mc_hist.GetXaxis().SetTitle("Electron Energy (GeV)")
    mc_hist.GetYaxis().SetTitle("Entries")
    MNVPLOTTER.ApplyAxisStyle(mc_hist,True,True)

    ratio =  data_hist.Clone()
    MNVPLOTTER.DrawDataMC(data_hist,mc_hist,1,"TR")
    #MNVPLOTTER.DrawDataMC(data_hist,mc_hist,1,"TR")
    mc = mc_hist.GetCVHistoWithError()
    mc.SetLineColor(ROOT.kRed)
    mc.SetLineWidth(3)
    mc.SetMarkerStyle(0)
    mc.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
    mc.Draw("E2 SAME")
    
    ratio.Divide(ratio, mc_hist)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    #Now fill mcRatio with 1 for bin content and fractional error
    mcRatio = mc_hist.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    if "Elastic" in title:
        mcRatio.GetXaxis().SetTitle("Electron Energy (GeV)")
    for whichBin in range(1, mcRatio.GetXaxis().GetNbins()+1): 
        mcRatio.SetBinError(whichBin, max(mcRatio.GetBinContent(whichBin), 1e-9))
        mcRatio.SetBinContent(whichBin, 1)

    ratio.SetTitle("")
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(3)
    ratio.SetTitleSize(0)

    if "Elastic" in title:
        ratio.GetXaxis().SetTitle("Electron Energy (GeV)")
    else:
        ratio.GetXaxis().SetTitle("E_{#nu} Estimator (GeV)")
    ratio.GetYaxis().SetTitle("Data / MC")
    ratio.GetYaxis().SetLabelSize(.15)
    ratio.GetYaxis().SetTitleSize(0.16)
    ratio.GetYaxis().SetTitleOffset(0.4)
    ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.

    ratio.GetXaxis().SetTitleSize(0.16)
    ratio.GetXaxis().SetTitleOffset(0.9)
    ratio.GetXaxis().SetLabelSize(.15)

    ratio.SetMinimum(0.5)
    ratio.SetMaximum(1.5)
    ratio.Draw()

    #Error envelope for the MC
    mcRatio.SetLineColor(ROOT.kRed)
    mcRatio.SetLineWidth(3)
    mcRatio.SetMarkerStyle(0)
    mcRatio.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
    mcRatio.Draw("E2 SAME")

    #Draw a flat line at 1 for ratio of MC to itself
    straightLine = mcRatio.Clone()
    straightLine.SetFillStyle(0)
    straightLine.Draw("HIST L SAME")
    top.cd()

    overall.Print("{}_thesis.png".format(title))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bin_width",
                        dest = "binwidth",
                        default = False,
                        action="store_true",
    )
    parser.add_argument("--thesis_plots",
                        dest = "thesis_plots",
                        default = False,
                        action="store_true",
    )
    parser.add_argument("--pseudodata",
                        dest = "pseudodata",
                        default = False,
                        action="store_true",
    )
    parser.add_argument("--exclude",
                        dest ="exclude",
                        default = "None",
    )
    parser.add_argument("--ratio",
                        dest ="ratio",
                        default = False,
                        action = "store_true",
    )

    args = parser.parse_args()
    binwidthScale = args.binwidth
    pseudodata = args.pseudodata
    exclude = str(args.exclude).lower()
    doratio = args.ratio
    thesis_plots = args.thesis_plots
    if thesis_plots:
        binwidthScale = True
    if doratio:
        ftag = "ratio"
    else:
        ftag = "stitched"


    type_path_map = {'data':'/exp/minerva/data/users/rhowell/nu_e/kin_dist_dataFHC_Selection_100Univ_thesis_MAD.root','mc':'/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Selection_100Univ_thesis_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    fhc_selection_mc = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/bkgfit_FHC_Selection_100Univ_N4_tune_thesis_MAD.root")
    fhc_selection_mcHold = HistHolder("Predicted MC",fhc_selection_mc,"Signal",True,mc_pot,standPOT)
    h_fhc_selection_data = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/bkgfit_FHC_Selection_100Univ_N4_tune_thesis_MAD.root").Get("EN4_data_bkgSubbed")
    fhc_sel = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Selection_100Univ_Genie_thesis_MAD.root")
    fhc_sel_template = HistHolder("Reco Energy vs L/E",fhc_sel,"Signal",True,mc_pot,standPOT)


    type_path_map = {'mc':'/exp/minerva/data/users/rhowell/nu_e_swap/kin_dist_mcFHC_Selection_100Univ_thesis_swap_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    print("fhc_swap: ",mc_pot,data_pot)
    fhc_sel_swap = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e_swap/kin_dist_mcFHC_Selection_100Univ_thesis_swap_MAD.root")
    fhc_swap_sel_template = HistHolder("Reco Energy vs L/E",fhc_sel_swap,"Signal",True,mc_pot,standPOT)
    fhc_swap_selection_mc = HistHolder("Biased Neutrino Energy",fhc_sel_swap,"Signal",True,mc_pot,standPOT)
    
    type_path_map = {'data':'/exp/minerva/data/users/rhowell/antinu_e/kin_dist_dataRHC_Selection_1000Univ_thesis_MAD.root','mc':'/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Selection_1000Univ_thesis_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_1000Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    rhc_selection_mc = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/bkgfit_RHC_Selection_1000Univ_N4_tune_thesis_MAD.root")
    rhc_selection_mcHold = HistHolder("Predicted MC",rhc_selection_mc,"Signal",True,mc_pot,standPOT)
    h_rhc_selection_data = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/bkgfit_RHC_Selection_1000Univ_N4_tune_thesis_MAD.root").Get("EN4_data_bkgSubbed")
    rhc_sel = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Selection_Genie_thesis_MAD.root")
    rhc_sel_template = HistHolder("Reco Energy vs L/E",rhc_sel,"Signal",True,mc_pot,standPOT)

    type_path_map = {'mc':'/exp/minerva/data/users/rhowell/antinu_e_swap/kin_dist_mcRHC_Selection_1000Univ_thesis_swap_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_1000Univ",type_path_map,"MAD",True)
    print("rhc_swap: ",mc_pot,data_pot)
    rhc_sel_swap = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e_swap/kin_dist_mcRHC_Selection_1000Univ_thesis_swap_MAD.root")
    rhc_swap_sel_template = HistHolder("Reco Energy vs L/E",rhc_sel_swap,"Signal",True,mc_pot,standPOT)
    rhc_swap_selection_mc = HistHolder("Biased Neutrino Energy",rhc_sel_swap,"Signal",True,mc_pot,standPOT)

    ### NuMu Selection ###
    cates = ["CCQE","CCDelta","CCDIS","CC2p2h","CCOther","CCWrongSign"]

    type_path_map = {'data':'/exp/minerva/data/users/rhowell/nu_mu/kin_dist_dataFHC_Selection_100Univ_thesis_muon_MAD.root','mc':'/exp/minerva/data/users/rhowell/nu_mu/kin_dist_mcFHC_Selection_100Univ_thesis_muon_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    fhc_muselection_mcHold = HistHolder("Biased Neutrino Energy",mc_file,"Signal",True,mc_pot,standPOT)
    fhc_muselection_dataHold = HistHolder("Biased Neutrino Energy",data_file,"Signal",False,data_pot,standPOT)
    fhc_musel_template = HistHolder("Reco Energy vs L/E",mc_file,"Signal",True,mc_pot,standPOT)

    type_path_map = {'data':'/exp/minerva/data/users/rhowell/antinu_mu/kin_dist_dataRHC_Selection_1000Univ_thesis_muon_MAD.root','mc':'/exp/minerva/data/users/rhowell/antinu_mu/kin_dist_mcRHC_Selection_1000Univ_thesis_muon_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_1000Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    rhc_muselection_mcHold = HistHolder("Biased Neutrino Energy",mc_file,"Signal",True,mc_pot,standPOT)
    rhc_muselection_dataHold = HistHolder("Biased Neutrino Energy",data_file,"Signal",False,data_pot,standPOT)
    rhc_musel_template = HistHolder("Reco Energy vs L/E",mc_file,"Signal",True,mc_pot,standPOT)

    h2_fhc_nueel_template_nue = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_nue')
    h2_fhc_nueel_template_numu = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_numu')
    h2_fhc_nueel_template_anue = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinue')
    h2_fhc_nueel_template_anumu = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinumu')

    h_fhc_imd_template = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_mu/FHC_imd_scattering_template.root').Get('imd_scattering_template')
    h_rhc_imd_template = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_mu/RHC_imd_scattering_template.root').Get('imd_scattering_template')

    h2_rhc_nueel_template_nue = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_nue')
    h2_rhc_nueel_template_numu = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_numu')
    h2_rhc_nueel_template_anue = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinue')
    h2_rhc_nueel_template_anumu = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinumu')

    h_fhc_swap_selection = fhc_swap_selection_mc.hists["Total"].Clone("sigTotal") 
    h_rhc_swap_selection = rhc_swap_selection_mc.hists["Total"].Clone("sigTotal") 

    h_fhc_muselection_mc = fhc_muselection_mcHold.GetHist()
    h_rhc_muselection_mc = rhc_muselection_mcHold.GetHist()
    h_fhc_muselection_data = fhc_muselection_dataHold.GetHist()
    h_rhc_muselection_data = rhc_muselection_dataHold.GetHist()
    h2_fhc_musel_template = fhc_musel_template.GetHist()
    h2_rhc_musel_template = rhc_musel_template.GetHist()

    h2_fhc_template = fhc_sel_template.hists["Total"].Clone("sigTotal")
    h2_rhc_template = rhc_sel_template.hists["Total"].Clone("sigTotal")
    h2_fhc_swap_template = fhc_swap_sel_template.hists["Total"].Clone("sigTotal")
    h2_rhc_swap_template = rhc_swap_sel_template.hists["Total"].Clone("sigTotal")

    fhc_muselection_mcHold.POTScale(binwidthScale)
    rhc_muselection_mcHold.POTScale(binwidthScale)
    fhc_muselection_dataHold.POTScale(binwidthScale)
    rhc_muselection_dataHold.POTScale(binwidthScale)
    fhc_swap_selection_mc.POTScale(binwidthScale)
    rhc_swap_selection_mc.POTScale(binwidthScale)

    h_fhc_swap_selection.Reset()
    h_rhc_swap_selection.Reset()

    h2_fhc_template.Reset()
    h2_rhc_template.Reset()
    h2_fhc_swap_template.Reset()
    h2_rhc_swap_template.Reset()

    h_fhc_muselection_mc.Reset()
    h_rhc_muselection_mc.Reset()
    h2_fhc_musel_template.Reset()
    h2_rhc_musel_template.Reset()

    for group in fhc_sel_template.hists:
        if group == "Total":
            continue
        elif group in SIGNAL_DEFINATION:
            if fhc_sel_template.hists[group]:
                h2_fhc_template.Add(fhc_sel_template.hists[group])

    for group in rhc_sel_template.hists:
        if group == "Total":
            continue
        elif group in SIGNAL_DEFINATION:
            if rhc_sel_template.hists[group]:
                h2_rhc_template.Add(rhc_sel_template.hists[group])

    for group in fhc_swap_sel_template.hists:
        if group == "Total":
            continue
        elif group in SWAP_SIGNAL_DEFINATION:
            if fhc_swap_sel_template.hists[group]:
                h_fhc_swap_selection.Add(fhc_swap_selection_mc.hists[group])
                h2_fhc_swap_template.Add(fhc_swap_sel_template.hists[group])

    for group in rhc_swap_sel_template.hists:
        if group == "Total":
            continue
        elif group in SWAP_SIGNAL_DEFINATION:
            if rhc_swap_sel_template.hists[group]:
                h_rhc_swap_selection.Add(rhc_swap_selection_mc.hists[group])
                h2_rhc_swap_template.Add(rhc_swap_sel_template.hists[group])

    for group in fhc_muselection_mcHold.hists:
        if group == "Total":
            continue
        elif group in cates:
            h_fhc_muselection_mc.Add(fhc_muselection_mcHold.hists[group])
            h2_fhc_musel_template.Add(fhc_musel_template.hists[group])

    for group in fhc_muselection_mcHold.hists:
        if group == "Total":
            continue
        elif group in cates:
            h_rhc_muselection_mc.Add(rhc_muselection_mcHold.hists[group])
            h2_rhc_musel_template.Add(rhc_musel_template.hists[group])

    h_fhc_selection_mc = fhc_selection_mcHold.GetHist()
    h_rhc_selection_mc = rhc_selection_mcHold.GetHist()

    #get FHC electron energy hist
    h_fhc_nueel_mc = ROOT.TFile('/exp/minerva/data/users/rhowell/nueel/DeepikaMc_LauraData_Mnv.root').Get('mnv_eff_cor_mc')
    h_fhc_nueel_data = ROOT.TFile('/exp/minerva/data/users/rhowell/nueel/DeepikaMc_LauraData_Mnv.root').Get('mnv_eff_cor_data')

    #get RHC electron energy hist
    h_rhc_nueel_mc = ROOT.TFile("/exp/minerva/data/users/rhowell/nueel/electronE_bkgsub_FullSetFinalv2_nobinnorm.root").Get('mc_fluxonly')   #.Get('mc_bkgsub_effcor')
    h_rhc_nueel_data =  ROOT.TFile("/exp/minerva/data/users/rhowell/nueel/electronE_bkgsub_FullSetFinalv2_nobinnorm.root").Get('data_bkgsub_effcor')

    #get IMD histograms
    h_fhc_imd_mc = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("FHC_MC")
    h_fhc_imd_data = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("FHC_Data")

    h_rhc_imd_mc = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("RHC_MC")
    h_rhc_imd_data = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("RHC_Data")

    fhcnueelnue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_nue')
    fhcnueelnumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_numu')
    fhcnueelanue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_anue')
    fhcnueelanumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_anumu')
    rhcnueelnue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_nue')
    rhcnueelnumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_numu')
    rhcnueelanue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_anue')
    rhcnueelanumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_anumu')

    #Make all hist have the same errorbands. Filled with CV if they don't have it.
    SyncErrorBandsv2([h_fhc_nueel_mc,h_rhc_nueel_mc,h_fhc_selection_mc,h_rhc_selection_mc,h_fhc_muselection_mc,h_rhc_muselection_mc])
    SyncErrorBandsv2([h_fhc_nueel_data,h_rhc_nueel_data,h_fhc_selection_data,h_rhc_selection_data,h_fhc_muselection_data,h_rhc_muselection_data])
    
    mc_histsToStitch = OrderedDict()
    mc_histsToStitch['fhc_nueel'] = h_fhc_nueel_mc
    mc_histsToStitch['fhc_imd'] = h_fhc_imd_mc
    mc_histsToStitch['fhc_muselection']=h_fhc_muselection_mc
    mc_histsToStitch['fhc_selection']=h_fhc_selection_mc
    mc_histsToStitch['rhc_nueel'] = h_rhc_nueel_mc
    mc_histsToStitch['rhc_imd'] = h_rhc_imd_mc
    mc_histsToStitch['rhc_muselection']=h_rhc_muselection_mc
    mc_histsToStitch['rhc_selection']=h_rhc_selection_mc

    nue_templates = OrderedDict()
    h2_fhc_nueel_template_nue.Add(h2_fhc_nueel_template_anue)
    nue_templates['fhc_nueel'] = h2_fhc_nueel_template_nue.Clone()
    nue_templates['fhc_imd'] = h_fhc_imd_template.Clone() # placeholder
    nue_templates['fhc_imd'].Reset()
    nue_templates['fhc_muselection']=h2_fhc_musel_template.Clone() # placeholder
    nue_templates['fhc_muselection'].Reset()
    nue_templates['fhc_selection']=h2_fhc_template.Clone()
    h2_rhc_nueel_template_anue.Add(h2_rhc_nueel_template_nue)
    nue_templates['rhc_nueel'] = h2_rhc_nueel_template_anue.Clone()
    nue_templates['rhc_imd'] = h_rhc_imd_template.Clone() #  placeholder
    nue_templates['rhc_imd'].Reset()
    nue_templates['rhc_muselection']=h2_rhc_musel_template.Clone()# placeholder
    nue_templates['rhc_muselection'].Reset()
    nue_templates['rhc_selection']=h2_rhc_template.Clone()
    
    numu_templates = OrderedDict()
    h2_fhc_nueel_template_numu.Add(h2_fhc_nueel_template_anumu)
    numu_templates['fhc_nueel'] = h2_fhc_nueel_template_numu.Clone()
    numu_templates['fhc_imd'] = h_fhc_imd_template.Clone() # placeholder
    numu_templates['fhc_muselection']=h2_fhc_musel_template.Clone()
    numu_templates['fhc_selection']=h2_fhc_swap_template.Clone()
    numu_templates['fhc_selection'].Reset()
    h2_rhc_nueel_template_anumu.Add(h2_rhc_nueel_template_numu)
    numu_templates['rhc_nueel'] = h2_rhc_nueel_template_anumu.Clone()
    numu_templates['rhc_imd'] = h_rhc_imd_template.Clone() # placeholder
    numu_templates['rhc_muselection']=h2_rhc_musel_template.Clone()
    numu_templates['rhc_selection']=h2_rhc_swap_template.Clone()
    numu_templates['rhc_selection'].Reset()

    swap_templates = OrderedDict()
    h2_fhc_nueel_template_numu.Add(h2_fhc_nueel_template_anumu)
    swap_templates['fhc_nueel'] = h2_fhc_nueel_template_numu.Clone()
    swap_templates['fhc_imd'] = h_fhc_imd_template.Clone() # placeholder
    swap_templates['fhc_muselection']=h2_fhc_musel_template.Clone()
    swap_templates['fhc_selection']=h2_fhc_swap_template.Clone()
    h2_rhc_nueel_template_anumu.Add(h2_rhc_nueel_template_numu)
    swap_templates['rhc_nueel'] = h2_rhc_nueel_template_anumu.Clone()
    swap_templates['rhc_imd'] = h_rhc_imd_template.Clone() # placeholder
    swap_templates['rhc_muselection']=h2_rhc_musel_template.Clone()
    swap_templates['rhc_selection']=h2_rhc_swap_template.Clone()
    for h in swap_templates.keys():
        if "_selection" not in h:
            swap_templates[h].Reset()

    data_histsToStitch = OrderedDict()
    data_histsToStitch['fhc_nueel'] = h_fhc_nueel_data
    data_histsToStitch['fhc_imd'] = h_fhc_imd_data
    data_histsToStitch['fhc_muselection']=h_fhc_muselection_data
    data_histsToStitch['fhc_selection']=h_fhc_selection_data
    data_histsToStitch['rhc_nueel'] = h_rhc_nueel_data
    data_histsToStitch['rhc_imd'] = h_rhc_imd_data
    data_histsToStitch['rhc_muselection']=h_rhc_muselection_data
    data_histsToStitch['rhc_selection']=h_rhc_selection_data

    if thesis_plots:
        # Create histograms
        h_fhc_selection_mc.Scale(1,"width")
        h_rhc_selection_mc.Scale(1,"width")
        h_fhc_selection_data.Scale(1,"width")
        h_rhc_selection_data.Scale(1,"width")

        h_fhc_nueel_mc.Scale(2.0,"width")
        h_rhc_nueel_mc.Scale(2.0,"width")
        h_fhc_nueel_data.Scale(2.0,"width")
        h_rhc_nueel_data.Scale(2.0,"width")

        fhcnueelnue.Scale(2.0,"width") 
        fhcnueelnumu.Scale(2.0,"width") 
        fhcnueelanue.Scale(2.0,"width") 
        fhcnueelanumu.Scale(2.0,"width") 

        rhcnueelnue.Scale(2.0,"width") 
        rhcnueelnumu.Scale(2.0,"width") 
        rhcnueelanue.Scale(2.0,"width") 
        rhcnueelanumu.Scale(2.0,"width")

        mc_histograms = [h_fhc_nueel_mc]
        mc_histograms.extend([h_rhc_nueel_mc,h_rhc_muselection_mc,h_rhc_selection_mc])
        data_histograms = [h_fhc_nueel_data]
        data_histograms.extend([h_rhc_nueel_data,h_rhc_muselection_data,h_rhc_selection_data])
        titles = ["FHC Elastic Scattering","RHC Elastic Scattering","RHC CC Muon Neutrino","RHC CC Electron Neutrino"]

        for i in range(len(mc_histograms)):
            mc_hist = mc_histograms[i]
            data_hist = data_histograms[i]
            title = titles[i]
            MakeThesisPlot(mc_hist,data_hist,title)

            #bestParams = {"ue4":[0.062,0.079,0.10,0.14,0.14,0.07,0.03,0.03],"umu4":[0,0.001,0.002,0,0.001,0.005,0.031,0.073],"m":[18,17,16,11,11,20,15,10]}
            #bestParams = {"m":[7,7,7],"ue4":[0.12,.12,.12],"umu4":[0,0.005,0.01]}  # RAA
            #bestParams = {"m":[3,2.5,2.5],"ue4":[0.1,.054,0.015],"umu4":[.013,.031,.073],'utau4':[0,0,0]} # best fits
            #bestParams = {"m":[7,4,6,10,80,10,10,10,10],"ue4":[.12,.03,.03,0.03,.03,.02,.1,.03,.03],"umu4":[0,.031,.031,0.031,0.031,.005,.005,.031,.073],'utau4':[0,0,0,0,0,0,0,0.66,0.66]} # m effect
            #bestParams = {"m":[10,10],"ue4":[.03,.03],"umu4":[.031,.073]} # tau effect
            bestParams = {"m":[7],"ue4":[.12],"umu4":[.0],'utau4':[0.0]} # tau effect
            stacked = MakeOscillationPlot(bestParams,mc_hist,data_hist,title)
        exit()
        

    if binwidthScale == True:
        canvas = ROOT.TCanvas("canvas", "Histograms", 1200, 600)
        canvas.Divide(4, 3,0.01,0.00001)
        # Create histograms
        h_fhc_selection_mc.Scale(1,"width")
        h_rhc_selection_mc.Scale(1,"width")
        h_fhc_selection_data.Scale(1,"width")
        h_rhc_selection_data.Scale(1,"width")

        h_fhc_nueel_mc.Scale(2.0,"width")
        h_rhc_nueel_mc.Scale(2.0,"width")
        h_fhc_nueel_data.Scale(2.0,"width")
        h_rhc_nueel_data.Scale(2.0,"width")
        fhchist = h_fhc_selection_mc.Clone()
        rhchist = h_rhc_selection_mc.Clone()
        fhc_CV = h_fhc_selection_mc.GetCVHistoWithStatError()
        rhc_CV = h_rhc_selection_mc.GetCVHistoWithStatError()
        fhchist.DivideSingle(fhchist,fhc_CV)
        rhchist.DivideSingle(rhchist,rhc_CV)

        c1 = ROOT.TCanvas()
        err = "Flux"
        fhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("FHC_"+err+".png")
        rhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("RHC_"+err+".png")

        err = "elE_ECAL"
        fhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("FHC_"+err+".png")
        rhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("RHC_"+err+".png")
        err = "elE_HCAL"
        fhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("FHC_"+err+".png")
        rhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("RHC_"+err+".png")
        err = "eltheta"
        fhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("FHC_"+err+".png")
        rhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("RHC_"+err+".png")
        err = "response_em"
        fhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("FHC_"+err+".png")
        rhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("RHC_"+err+".png")
        err = "Leakage_Uncertainty"
        fhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("FHC_"+err+".png")
        rhchist.GetVertErrorBand(err).DrawAll("hist",True)
        c1.Print("RHC_"+err+".png")
        exit()

        mc_histograms = [h_fhc_nueel_mc,h_fhc_muselection_mc,h_fhc_selection_mc]
        mc_histograms.extend([h_rhc_nueel_mc,h_rhc_muselection_mc,h_rhc_selection_mc])
        data_histograms = [h_fhc_nueel_data,h_fhc_muselection_data,h_fhc_selection_data]
        data_histograms.extend([h_rhc_nueel_data,h_rhc_muselection_data,h_rhc_selection_data])
        for h in range(len(mc_histograms)):
            mc_histograms[h].PopVertErrorBand("LowQ2Pi")
            data_histograms[h].PopVertErrorBand("LowQ2Pi")
            hist = mc_histograms[h]
            hist.GetXaxis().SetTitle("E_{#nu} Estimator (GeV)")
            data_histograms[h].GetXaxis().SetTitle("E_{#nu} Estimator (GeV)")
            if pseudodata:
                for q in range(0,hist.GetNbinsX()+1):
                    data_histograms[h].SetBinContent(q,hist.GetBinContent(q))

        # Plot histograms on the canvas
        h_i = 0
        for i in range(12):
            pad = canvas.cd(i+1)
            if i < 8:# and i != 4 and i != 7:
                pad.SetBottomMargin(0)
            if i > 3:
                pad.SetTopMargin(0)
            if i not in [0,4,8]:
                pad.SetLeftMargin(0.08)

            pad.SetRightMargin(0.001)
            ROOT.TGaxis.SetExponentOffset(0.04,-0.1,"y")
            MNVPLOTTER.legend_text_size=0.07
            if i == 3:
                mcRatio = h_fhc_muselection_mc.Clone()
                mcRatio.Divide(mcRatio,h_fhc_selection_mc)
                dataRatio = h_fhc_muselection_data.Clone()
                dataRatio.Divide(dataRatio,h_fhc_selection_data)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                mc = mcRatio.GetCVHistoWithError()
                mc.SetLineColor(ROOT.kRed)
                mc.SetLineWidth(3)
                mc.SetMarkerStyle(0)
                mc.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
                mc.DrawClone("E2 SAME")
                MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                hist.GetYaxis().SetRangeUser(0,0.35)
            elif i == 7:
                mcRatio = h_rhc_muselection_mc.Clone()
                mcRatio.Divide(mcRatio,h_rhc_selection_mc)
                dataRatio = h_rhc_muselection_data.Clone()
                dataRatio.Divide(dataRatio,h_rhc_selection_data)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                mc = mcRatio.GetCVHistoWithError()
                mc.SetLineColor(ROOT.kRed)
                mc.SetLineWidth(3)
                mc.SetMarkerStyle(0)
                mc.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
                mc.DrawClone("E2 SAME")
                MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                hist.GetYaxis().SetRangeUser(0,0.35)
            elif i == 8:
                mcRatio = h_fhc_nueel_mc.Clone()
                mcRatio.Divide(mcRatio,h_rhc_nueel_mc)
                dataRatio = h_fhc_nueel_data.Clone()
                dataRatio.Divide(dataRatio,h_rhc_nueel_data)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                mc = mcRatio.GetCVHistoWithError()
                mc.SetLineColor(ROOT.kRed)
                mc.SetLineWidth(3)
                mc.SetMarkerStyle(0)
                mc.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
                mc.DrawClone("E2 SAME")
                MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                hist.GetYaxis().SetRangeUser(0,0.35)
            elif i == 9:
                mcRatio = h_fhc_muselection_mc.Clone()
                mcRatio.Divide(mcRatio,h_rhc_muselection_mc)
                dataRatio = h_fhc_muselection_data.Clone()
                dataRatio.Divide(dataRatio,h_rhc_muselection_data)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                hist.GetYaxis().SetRangeUser(0,0.35)
            elif i == 10:
                mcRatio = h_fhc_selection_mc.Clone()
                mcRatio.Divide(mcRatio,h_rhc_selection_mc)
                dataRatio = h_fhc_selection_data.Clone()
                dataRatio.Divide(dataRatio,h_rhc_selection_data)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                hist.GetYaxis().SetRangeUser(0,0.35)
            elif i != 8 and h_i < len(data_histograms):
                MNVPLOTTER.DrawDataMC(data_histograms[h_i],mc_histograms[h_i],1,"TR")
                mc = mc_histograms[h_i].GetCVHistoWithError()
                mc.SetLineColor(ROOT.kRed)
                mc.SetLineWidth(3)
                mc.SetMarkerStyle(0)
                mc.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
                mc.DrawClone("E2 SAME")
                h_i+=1
            elif i == 11:
                fhc_mcRatio = h_fhc_muselection_mc.Clone()
                fhc_mcRatio.Divide(fhc_mcRatio,h_fhc_selection_mc)
                fhc_dataRatio = h_fhc_muselection_data.Clone()
                fhc_dataRatio.Divide(fhc_dataRatio,h_fhc_selection_data)

                rhc_mcRatio = h_rhc_muselection_mc.Clone()
                rhc_mcRatio.Divide(rhc_mcRatio,h_rhc_selection_mc)
                rhc_dataRatio = h_rhc_muselection_data.Clone()
                rhc_dataRatio.Divide(rhc_dataRatio,h_rhc_selection_data)
                mcRatio = fhc_mcRatio.Clone()
                mcRatio.Divide(mcRatio,rhc_mcRatio)
                dataRatio = fhc_dataRatio.Clone()
                dataRatio.Divide(dataRatio,rhc_dataRatio)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                hist.GetYaxis().SetRangeUser(0,0.35)
        # Show the canvas
        canvas.Print("grid_plots_errors.png")

        mc_histograms = [h_fhc_nueel_mc,h_fhc_muselection_mc,h_fhc_selection_mc]
        mc_histograms.extend([h_rhc_nueel_mc,h_rhc_muselection_mc,h_rhc_selection_mc])
        data_histograms = [h_fhc_nueel_data,h_fhc_muselection_data,h_fhc_selection_data]
        data_histograms.extend([h_rhc_nueel_data,h_rhc_muselection_data,h_rhc_selection_data])
        for h in range(len(mc_histograms)):
            mc_histograms[h].PopVertErrorBand("LowQ2Pi")
            data_histograms[h].PopVertErrorBand("LowQ2Pi")
            hist = mc_histograms[h]
            hist.GetXaxis().SetTitle("E_{#nu} Estimator (GeV)")
            data_histograms[h].GetXaxis().SetTitle("E_{#nu} Estimator (GeV)")
            if pseudodata:
                for q in range(0,hist.GetNbinsX()+1):
                    data_histograms[h].SetBinContent(q,hist.GetBinContent(q))

        # Plot histograms on the canvas
        h_i = 0
        for i in range(12):
            pad = canvas.cd(i+1)

            ROOT.TGaxis.SetExponentOffset(0.04,-0.1,"y")
            MNVPLOTTER.legend_text_size=0.07
            if i == 3:
                mcRatio = h_fhc_muselection_mc.Clone()
                mcRatio.Divide(mcRatio,h_fhc_selection_mc)
                dataRatio = h_fhc_muselection_data.Clone()
                dataRatio.Divide(dataRatio,h_fhc_selection_data)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                mc = mcRatio.GetCVHistoWithError()
                mc.SetLineColor(ROOT.kRed)
                mc.SetLineWidth(3)
                mc.SetMarkerStyle(0)
                mc.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
                mc.DrawClone("E2 SAME")
                #MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                #hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                #hist.GetYaxis().SetRangeUser(0,0.35)
            elif i == 7:
                mcRatio = h_rhc_muselection_mc.Clone()
                mcRatio.Divide(mcRatio,h_rhc_selection_mc)
                dataRatio = h_rhc_muselection_data.Clone()
                dataRatio.Divide(dataRatio,h_rhc_selection_data)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                mc = mcRatio.GetCVHistoWithError()
                mc.SetLineColor(ROOT.kRed)
                mc.SetLineWidth(3)
                mc.SetMarkerStyle(0)
                mc.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
                mc.DrawClone("E2 SAME")
                #MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                #hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                #hist.GetYaxis().SetRangeUser(0,0.35)
            elif i == 8:
                mcRatio = h_fhc_nueel_mc.Clone()
                mcRatio.Divide(mcRatio,h_rhc_nueel_mc)
                dataRatio = h_fhc_nueel_data.Clone()
                dataRatio.Divide(dataRatio,h_rhc_nueel_data)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                mc = mcRatio.GetCVHistoWithError()
                mc.SetLineColor(ROOT.kRed)
                mc.SetLineWidth(3)
                mc.SetMarkerStyle(0)
                mc.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
                mc.DrawClone("E2 SAME")
                #MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                #hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                #hist.GetYaxis().SetRangeUser(0,0.35)
            elif i == 9:
                mcRatio = h_fhc_muselection_mc.Clone()
                mcRatio.Divide(mcRatio,h_rhc_muselection_mc)
                dataRatio = h_fhc_muselection_data.Clone()
                dataRatio.Divide(dataRatio,h_rhc_muselection_data)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                #MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                #hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                #hist.GetYaxis().SetRangeUser(0,0.35)
            elif i == 10:
                mcRatio = h_fhc_selection_mc.Clone()
                mcRatio.Divide(mcRatio,h_rhc_selection_mc)
                dataRatio = h_fhc_selection_data.Clone()
                dataRatio.Divide(dataRatio,h_rhc_selection_data)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                #MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                #hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                #hist.GetYaxis().SetRangeUser(0,0.35)
            elif i != 8 and h_i < len(data_histograms):
                MNVPLOTTER.DrawDataMC(data_histograms[h_i],mc_histograms[h_i],1,"TR")
                mc = mc_histograms[h_i].GetCVHistoWithError()
                mc.SetLineColor(ROOT.kRed)
                mc.SetLineWidth(3)
                mc.SetMarkerStyle(0)
                mc.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
                mc.DrawClone("E2 SAME")
                h_i+=1
            elif i == 11:
                fhc_mcRatio = h_fhc_muselection_mc.Clone()
                fhc_mcRatio.Divide(fhc_mcRatio,h_fhc_selection_mc)
                fhc_dataRatio = h_fhc_muselection_data.Clone()
                fhc_dataRatio.Divide(fhc_dataRatio,h_fhc_selection_data)

                rhc_mcRatio = h_rhc_muselection_mc.Clone()
                rhc_mcRatio.Divide(rhc_mcRatio,h_rhc_selection_mc)
                rhc_dataRatio = h_rhc_muselection_data.Clone()
                rhc_dataRatio.Divide(rhc_dataRatio,h_rhc_selection_data)
                mcRatio = fhc_mcRatio.Clone()
                mcRatio.Divide(mcRatio,rhc_mcRatio)
                dataRatio = fhc_dataRatio.Clone()
                dataRatio.Divide(dataRatio,rhc_dataRatio)
                MNVPLOTTER.DrawDataMC(dataRatio,mcRatio,1,"TR")
                #MNVPLOTTER.DrawErrorSummary(mcRatio,"TR",False,True,0)
                #hist = ROOT.gPad.GetListOfPrimitives().FindObject("h_total_err_errSum_5087")
                #hist.GetYaxis().SetRangeUser(0,0.35)
        # Show the canvas
        canvas.Print("grid_plots_hists.png")
        exit()

    
    for h in list(mc_histsToStitch.keys()):
        try:
            if "fhc" in exclude and "selection" in h and "fhc" in h:
                print('deleting '+h)
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
            if "rhc" in exclude and "selection" in h and "rhc" in h:
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
            if "numu" in exclude and "muselection" in h:
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
            if "nue" in exclude and "_selection" in h:
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
            if "elastic" in exclude and "nueel" in h:
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
            if "imd" in exclude and "imd" in h:
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
        except:
            continue

    if doratio:
        for h in list(mc_histsToStitch.keys()):
            if "selection" in h:
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
        if "fhc" not in exclude:
            toClone = h_fhc_muselection_mc.Clone()
            toClone.Divide(toClone,h_fhc_selection_mc)
            mc_histsToStitch['fhc_ratio']=toClone
            toClone = h_fhc_muselection_data.Clone()
            toClone.Divide(toClone,h_fhc_selection_data)
            data_histsToStitch['fhc_ratio']=toClone
            nue_templates['fhc_ratio'] = h2_fhc_template.Clone()
            numu_templates['fhc_ratio']=h2_fhc_musel_template.Clone()
            swap_templates['fhc_ratio']=h2_fhc_swap_template.Clone()
        if "rhc" not in exclude:
            toClone = h_rhc_muselection_mc.Clone()
            toClone.Divide(toClone,h_rhc_selection_mc)
            mc_histsToStitch['rhc_ratio']=toClone
            toClone = h_rhc_muselection_data.Clone()
            toClone.Divide(toClone,h_rhc_selection_data)
            data_histsToStitch['rhc_ratio']=toClone
            nue_templates['rhc_ratio'] = h2_rhc_template.Clone()
            numu_templates['rhc_ratio']=h2_rhc_musel_template.Clone()
            swap_templates['rhc_ratio']=h2_rhc_swap_template.Clone()

    n_bins_new = 0
    for h in mc_histsToStitch:
        for i in range(0,mc_histsToStitch[h].GetNbinsX()+1):
            if mc_histsToStitch[h].GetBinContent(i) > minBinCont:
                n_bins_new+=1

    htemp_mc = ROOT.TH1D('mc_stitched', 'mc_stitched', n_bins_new, 0,  n_bins_new )
    htemp_mc_nue = ROOT.TH1D('mc_stitched_nue', 'mc_stitched_nue', n_bins_new, 0,  n_bins_new )
    htemp_mc_numu = ROOT.TH1D('mc_stitched_numu', 'mc_stitched_numu', n_bins_new, 0,  n_bins_new )
    htemp_mc_nutau = ROOT.TH1D('mc_stitched_nutau', 'mc_stitched_nutau', n_bins_new, 0,  n_bins_new )
    htemp_mc_nueselection = ROOT.TH1D('mc_stitched_nueselection', 'mc_stitched_nueselection', n_bins_new, 0,  n_bins_new )
    htemp_mc_ratio = ROOT.TH1D('mc_stitched_ratio', 'mc_stitched_ratio', n_bins_new, 0,  n_bins_new )
    htemp_mc_swap = ROOT.TH1D('mc_stitched_swap', 'mc_stitched_swap', n_bins_new, 0,  n_bins_new )
    htemp_data = ROOT.TH1D('data_stitched', 'data_stitched', n_bins_new, 0,  n_bins_new )

    StitchThis(htemp_mc, mc_histsToStitch,"",mc_histsToStitch)
    StitchThis(htemp_mc_nue, mc_histsToStitch,'nue',mc_histsToStitch)
    StitchThis(htemp_mc_numu, mc_histsToStitch,'numu',mc_histsToStitch)
    StitchThis(htemp_mc_nutau, mc_histsToStitch,'nutau',mc_histsToStitch)
    StitchThis(htemp_mc_nueselection, mc_histsToStitch,'nueselection',mc_histsToStitch)
    StitchThis(htemp_mc_ratio, mc_histsToStitch,'ratio',mc_histsToStitch)
    StitchThis(htemp_mc_swap, mc_histsToStitch,'swap',mc_histsToStitch)
    StitchThis(htemp_data, data_histsToStitch,"data",mc_histsToStitch,pseudodata)

    h2template_nue = ROOT.TH2D('LE_template_nue', 'LE_template_nue', n_bins_new, 0,  n_bins_new,34,0,0.495)
    h2template_numu = ROOT.TH2D('LE_template_numu', 'LE_template_numu', n_bins_new, 0,  n_bins_new,34,0,0.495)
    h2template_swap = ROOT.TH2D('LE_template_swap', 'LE_template_swap', n_bins_new, 0,  n_bins_new,34,0,0.495)

    StitchThis2D(h2template_nue, nue_templates,mc_histsToStitch,"nue")
    StitchThis2D(h2template_numu, numu_templates,mc_histsToStitch,"numu")
    StitchThis2D(h2template_swap, swap_templates,mc_histsToStitch,"swap")

    mnv_data = PlotUtils.MnvH1D( htemp_data )
    mnv_data.GetXaxis().SetTitle("Bin number")
    mnv_mc   = PlotUtils.MnvH1D( htemp_mc )
    mnv_mc_nue  = PlotUtils.MnvH1D( htemp_mc_nue )
    mnv_mc_numu = PlotUtils.MnvH1D( htemp_mc_numu )
    mnv_mc_swap = PlotUtils.MnvH1D( htemp_mc_swap )

    c1 = ROOT.TCanvas()
    mnv_mc_nue.Draw("hist")
    c1.Print("plots/nue_hist.png")
    mnv_mc_numu.Draw("hist")
    c1.Print("plots/numu_hist.png")
    mnv_mc_swap.Draw("hist")
    c1.Print("plots/nuswap_hist.png")

    htemp_mc_nutau.Draw("hist")
    c1.Print("plots/nutau_id.png")
    htemp_mc_nueselection.Draw("hist")
    c1.Print("plots/nusel_id.png")
    htemp_mc_ratio.Draw("hist")
    c1.Print("plots/nuratio_id.png")

    h2template_nue.Draw("colz")
    c1.Print("plots/nue_temp.png")
    h2template_numu.Draw("colz")
    c1.Print("plots/numu_temp.png")
    h2template_swap.Draw("colz")
    c1.Print("plots/nuswap_temp.png")

    FillErrorBandfromHist2( mnv_data, data_histsToStitch, mc_histsToStitch, pseudodata, isData=True)
    FillErrorBandfromHist2( mnv_mc, mc_histsToStitch, mc_histsToStitch, isData=False)

    FillErrorBandfromHist(mnv_mc_nue,  mc_histsToStitch, mc_histsToStitch, isData=False,sample="nue")
    FillErrorBandfromHist(mnv_mc_numu, mc_histsToStitch, mc_histsToStitch, isData=False,sample="numu")
    FillErrorBandfromHist(mnv_mc_swap, mc_histsToStitch, mc_histsToStitch, isData=False,sample="swap")

    mnv_mc.PopVertErrorBand("LowQ2Pi")
    mnv_mc_nue.PopVertErrorBand("LowQ2Pi")
    mnv_mc_numu.PopVertErrorBand("LowQ2Pi")
    mnv_mc_swap.PopVertErrorBand("LowQ2Pi")
    mnv_data.PopVertErrorBand("LowQ2Pi")

    mnv_mc.SetLineColor(ROOT.kRed)
    mc = mnv_mc.GetCVHistoWithError()
    data = mnv_data.GetCVHistoWithError()

    mc.GetXaxis().SetTitle("Bin number")
    mc.GetYaxis().SetTitle("Entries")
    MNVPLOTTER.ApplyAxisStyle(mc,True,True)

    canvas = ROOT.TCanvas()
    MNVPLOTTER.DrawErrorSummary(mnv_mc,"TR",True,True)
    canvas.Print("plots/stitched_errors_mc_Total.png")
    canvas = ROOT.TCanvas()
    MNVPLOTTER.DrawErrorSummary(mnv_mc_swap,"TR",True,True,0)
    canvas.Print("plots/stitched_errors_mc_swap_Total.png")
    canvas = ROOT.TCanvas()
    MNVPLOTTER.DrawErrorSummary(mnv_mc_numu,"TR",True,True)
    canvas.Print("plots/stitched_errors_mcnumu_Total.png")
    canvas = ROOT.TCanvas()
    MNVPLOTTER.DrawErrorSummary(mnv_mc_nue,"TR",True,True)
    canvas.Print("plots/stitched_errors_mcnue_Total.png")
    canvas = ROOT.TCanvas()
    MNVPLOTTER.DrawErrorSummary(mnv_data,"TR",True,True)
    canvas.Print("plots/stitched_errors_data_Total.png")

    dataprint = []
    mcprint = []
    for i in range(1,mnv_data.GetNbinsX()+1):
        dataprint.append(mnv_data.GetBinContent(i))
    for i in range(1,mnv_mc.GetNbinsX()+1):
        mcprint.append(mnv_mc.GetBinContent(i))
    dataprint = np.array(dataprint)
    mcprint = np.array(mcprint)
    
    np.savetxt("mc_"+ftag+".csv",mcprint,delimiter=",")
    np.savetxt("data_"+ftag+".csv",dataprint,delimiter=",")

    Nbins = mnv_mc.GetNbinsX()
    if False:
        errbands = []
        for k,v in CONSOLIDATED_ERROR_GROUPS.items():
            covMatrix = np.zeros(shape=[Nbins,Nbins],dtype='f')
            covmx = mnv_mc.GetSysErrorMatrix(v[0])
            if len(v) > 1:
                for vs in range(1,len(v)):
                    name = v[vs] 
                    covmx += mnv_mc.GetSysErrorMatrix(name)
            mnv_mc.PushCovMatrix(k,covmx)
            for i in range(0,Nbins):
                for j in range(0,Nbins):
                    covMatrix[i][j] = covmx[i+1][j+1]
            np.savetxt("covmatrix_{}_".format(k)+ftag+".csv",covMatrix,delimiter=',')
            errbands.extend(v)

        for errorband in mnv_mc.GetVertErrorBandNames():
            if errorband not in errbands:
                print("no vert error ",errorband)
        for errorband in mnv_mc.GetLatErrorBandNames():
            if errorband not in errbands:
                print("no lat error ",errorband)

        covMatrix = np.zeros(shape=[Nbins,Nbins],dtype='f')
        fluxcov = mnv_mc.GetSysErrorMatrix("Flux")
        for i in range(0,Nbins):
            for j in range(0,Nbins):
                covMatrix[i][j] = fluxcov[i+1][j+1]
        np.savetxt("covmatrix_flux_"+ftag+".csv".format(k),covMatrix,delimiter=',')

    GetCovarianceMatrix(mnv_mc,mnv_data)
    makeScales = False
    NbinsTotal = mnv_mc.GetNbinsX()
    if makeScales:
        for i in range(10):
            scale_matrix = GetScaleMatrix(covArray,NbinsTotal)
            np.savetxt('ScaleMatrix_{}.txt'.format(i),scale_matrix)
        logging.info("Made all scale matrices, exiting now.")

    if pseudodata:
        f = ROOT.TFile("/exp/minerva/app/users/rhowell/cmtuser/CCNue/FeldmanCousins/NuE_stitched_hists_pseudo.root","RECREATE")
    else:
        f = ROOT.TFile("/exp/minerva/app/users/rhowell/cmtuser/CCNue/FeldmanCousins/NuE_stitched_hists.root","RECREATE")

    mnv_data.Write()
    mnv_mc.Write()
    mnv_mc_nue.Write()
    mnv_mc_numu.Write()
    mnv_mc_swap.Write()
    #htemp_mc_numu.Write()
    #htemp_mc_nue.Write()
    htemp_mc_nutau.Write()
    htemp_mc_nueselection.Write()
    htemp_mc_ratio.Write()
    #htemp_mc_swap.Write()
    h2template_numu.Write()
    h2template_nue.Write()
    h2template_swap.Write()

    DataMCPlot(mnv_data,mnv_mc)
    print(Chi2DataMC(mnv_data,mnv_mc))
    f.Close()
