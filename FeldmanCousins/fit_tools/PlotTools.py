import os
import time
import logging, sys
import ROOT
import PlotUtils
import numpy as np
np.random.seed(0)
#np.set_printoptions(precision=1)
#np.set_printoptions(linewidth=1520)
#np.set_printoptions(threshold=sys.maxsize)
from scipy import optimize, integrate
import argparse
ccnueroot = os.environ.get('CCNUEROOT')

import math
import psutil
import multiprocessing
import threading
nthreads = 4
from array import array
from fit_tools.FitTools import *

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
#MNVPLOTTER.chi2_use_overflow_err = True
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.07
#MNVPLOTTER.extra_top_margin = -.035# go slightly closer to top of pad
#MNVPLOTTER.mc_bkgd_color = 46 
#MNVPLOTTER.mc_bkgd_line_color = 46
MNVPLOTTER.legend_n_columns = 1
#MNVPLOTTER.mc_line_width = 0
MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  
#MNVPLOTTER.axis_maximum = 500 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .055
#MNVPLOTTER.legend_offset_x           = .15
#MNVPLOTTER.mc_line_width = 0
MNVPLOTTER.axis_minimum=0.1

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

rowstodo = 100

def MakePlot(stitched_mc,stitched_data,templates,bestFit={"m":0,"ue4":0,"umu4":0,"utau4":0}):
    fitHist = GetOscillatedHistogram(stitched_mc, templates, bestFit["m"], bestFit["ue4"], bestFit["umu4"], 0.0)
    chi2_model = Chi2DataMC(stitched_data,fitHist)
    chi2_null = Chi2DataMC(stitched_data,stitched_mc)

    print(bestFit["m"], bestFit["ue4"], bestFit["umu4"],chi2_model-chi2_null)

    c1 = ROOT.TCanvas()
    margin = .12
    bottomFraction = .2
    overall = ROOT.TCanvas("Data/MC")
    top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()
    top.SetLogy()
    
    fitHist.GetXaxis().SetTitle("Bin number")
    fitHist.GetYaxis().SetTitle("Entries")
    MNVPLOTTER.ApplyAxisStyle(fitHist,True,True)

    nullRatio =  stitched_data.Clone()
    oscRatio =  fitHist.Clone()

    fitHist.SetLineColor(ROOT.kRed)
    fitHist.SetLineWidth(3)
    stitched_mc.SetLineColor(ROOT.kBlue)
    stitched_mc.SetLineWidth(3)
    stitched_mc.GetYaxis().SetTitle("Nevents")
    stitched_mc.Draw("hist")
    fitHist.Draw("hist same")
    stitched_data.Draw("same")

    osc = fitHist.GetCVHistoWithError()
    osc.SetLineColor(ROOT.kRed)
    osc.SetLineWidth(3)
    osc.SetMarkerStyle(0)
    osc.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
    osc.Draw("E2 SAME")

    null = stitched_mc.GetCVHistoWithError()
    null.SetLineColor(ROOT.kBlue)
    null.SetLineWidth(3)
    null.SetMarkerStyle(0)
    null.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
    null.Draw("E2 SAME")

    leg = ROOT.TLegend()
    leg.AddEntry(stitched_data,"Data","p")
    leg.AddEntry(fitHist,"Best Fit #chi^{2}="+"{:.1f}".format(chi2_model),"l")
    #leg.AddEntry(fitHist,"Oscillation #chi^{2}="+"{:.1f}".format(chi2_model),"l")
    #leg.AddEntry(fitHist,"RAA #chi^{2}="+"{:.1f}".format(chi2_model),"l")
    leg.AddEntry(stitched_mc,"Null Hypothesis #chi^{2}="+"{:.1f}".format(chi2_null),"l")
    leg.Draw()

    oscRatio.Divide(oscRatio, stitched_mc)
    nullRatio.Divide(nullRatio,stitched_mc)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    nullErrors = stitched_mc.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
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
    nullErrors.Draw("E2")

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
    nullRatio.SetLineColorAlpha(ROOT.kBlue+1,0.6)
    oscRatio.Draw("same")
    nullRatio.Draw("same")

    #Draw a flat line at 1 for oscRatio of MC to itself
    straightLine = nullErrors.Clone()
    straightLine.SetLineColor(ROOT.kBlue)
    straightLine.SetLineWidth(3)
    straightLine.SetFillColor(0)
    straightLine.Draw("HIST L SAME")

    leg1 = ROOT.TLegend(.5,.5,.9,.9)
    leg1.AddEntry(nullRatio,"Data/Null Hypothesis","p")
    leg1.AddEntry(oscRatio,"Model/Null Hypothesis","l")
    leg1.AddEntry(straightLine,"Null/Null Hypothesis","l")
    #leg1.Draw()

    top.cd()

    #overall.Print("fit_{}_thesis.png".format(bestFit["m"]))
    overall.Print("fit_{}_{}_{}.png".format(bestFit["m"],bestFit['ue4'],bestFit['umu4']))

def DataMCPlot(mnv_data,mnv_mc,cv_data,cv_mc):
    #mnv_data.Scale(1,"width")
    #mnv_mc.Scale(1,"width")
    #cv_data.Scale(1,"width")
    #cv_mc.Scale(1,"width")

    c1 = ROOT.TCanvas()
    margin = .12
    bottomFraction = .2
    overall = ROOT.TCanvas("Data/MC")
    top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()
    if mnv_data.GetBinContent(mnv_data.GetMaximumBin()) > 1e4:
        top.SetLogy()
    MNVPLOTTER.axis_minimum = 0.1
    cv_data.SetTitle(mnv_mc.GetName())
    cv_data.SetTitle(mnv_data.GetName())
    cv_data.GetXaxis().SetTitle("Bin number")
    cv_data.GetYaxis().SetTitle("Entries")
    MNVPLOTTER.ApplyAxisStyle(cv_data,True,True)

    ratio = cv_data.Clone()

    cv_data.SetMarkerColor(ROOT.kBlack)
    mnv_mc.SetLineColor(ROOT.kBlue)
    mnv_mc.SetLineWidth(3)
    mnv_mc.SetMarkerStyle(0)

    cv_data.Draw("P")
    mnv_mc.Draw("HIST SAME")
    cv_mc.SetLineColor(ROOT.kRed)
    cv_mc.SetLineWidth(3)
    cv_mc.SetMarkerStyle(0)
    cv_data.Draw("P SAME")
    cv_mc.Draw("HIST SAME")
    mc = mnv_mc.GetCVHistoWithError()
    mc.SetLineColor(ROOT.kBlue)
    mc.SetLineWidth(3)
    mc.SetMarkerStyle(0)
    mc.SetFillColorAlpha(ROOT.kBlue + 1, 0.4)
    mc.Draw("E2 SAME")
    cv = cv_mc.GetCVHistoWithError()
    cv.SetLineColor(ROOT.kRed)
    cv.SetLineWidth(3)
    cv.SetMarkerStyle(0)
    cv.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
    cv.Draw("E2 SAME")

    leg = ROOT.TLegend(.7,.7,1,1)
    leg.SetHeader(mnv_mc.GetName(),"C")
    leg.AddEntry(cv_data,"CV data",'p')
    leg.AddEntry(cv_mc,"CV MC",'l')
    leg.AddEntry(mnv_mc,"Universe MC","l")
    leg.Draw()

    ratio.Divide(ratio, cv_mc)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    #Now fill mcRatio with 1 for bin content and fractional error
    mcRatio = mnv_mc.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(1, mcRatio.GetXaxis().GetNbins()+1): 
        mcRatio.SetBinError(whichBin, max(mcRatio.GetBinContent(whichBin), 1e-9))
        mcRatio.SetBinContent(whichBin, 1)

    ratio.SetTitle("")
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(3)
    ratio.SetTitleSize(0)

    ratio.GetYaxis().SetTitle("#splitline{Ratio to}{Null Hypothesis}")
    ratio.GetYaxis().SetLabelSize(.13)
    ratio.GetYaxis().SetTitleSize(0.1)
    ratio.GetYaxis().SetTitleOffset(0.6)
    ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
    ratio.GetXaxis().SetTitleSize(0.16)
    ratio.GetXaxis().SetTitleOffset(0.9)
    ratio.GetXaxis().SetLabelSize(.15)
    ratio.GetXaxis().SetTitle("Bin Number")
    ratio.SetMinimum(0.5)
    ratio.SetMaximum(1.5)
    ratio.Draw()

    univratio = mnv_mc.Clone()
    univratio.Divide(univratio,cv_mc)
    univratio.SetLineColor(ROOT.kBlue)
    univratio.SetLineWidth(3)
    univratio.Draw("HIST SAME")

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
    
    MNVPLOTTER.WritePreliminary(0.4, 0.05, 5e-2, True)
    overall.Print("plots/{}.png".format(mnv_mc.GetName()))
