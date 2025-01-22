import os
import time
import logging, sys
import ROOT
import PlotUtils
import numpy as np
from root_numpy import matrix
from scipy import optimize, integrate
import argparse
ccnueroot = os.environ.get('CCNUEROOT')

import math
import psutil
import multiprocessing
import threading
nthreads = 4
from array import array
from Tools.FitTools import *
from Tools.OscHistogram import *
from Tools.PlotHelper import *

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
#MNVPLOTTER.chi2_use_overflow_err = True
MNVPLOTTER.draw_normalized_to_bin_width=False
#MNVPLOTTER.extra_top_margin = -.035# go slightly closer to top of pad
#MNVPLOTTER.mc_bkgd_color = 46 
#MNVPLOTTER.mc_bkgd_line_color = 46
MNVPLOTTER.legend_n_columns = 1
#MNVPLOTTER.mc_line_width = 0
#MNVPLOTTER.mc_line_width = 0
MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  
#MNVPLOTTER.axis_maximum = 500 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
#MNVPLOTTER.legend_text_size        = .4
#MNVPLOTTER.legend_offset_x           = .15
#MNVPLOTTER.mc_line_width = 0
MNVPLOTTER.axis_minimum=0.1

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

def PlotNorms(stitched_mc,stitched_data,fitHist,params=[]):
    chi2_model = Chi2DataMC(stitched_data,fitHist)
    chi2_null = Chi2DataMC(stitched_data,stitched_mc)

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
    nullErrors.SetMinimum(0)
    nullErrors.SetMaximum(2)

    
    #Draw the data ratios
    oscRatio.SetMinimum(0)
    oscRatio.SetMaximum(2)
    nullRatio.SetMinimum(0)
    nullRatio.SetMaximum(2)
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
    if len(params) > 0:
        leg1.SetHeader("fhc:{:.2f} rhc:{:.2f}".format(params[0],params[1]))
    leg1.AddEntry(nullRatio,"Data/Null Hypothesis","p")
    leg1.AddEntry(oscRatio,"Model/Null Hypothesis","l")
    leg1.AddEntry(straightLine,"Null/Null Hypothesis","l")
    #leg1.Draw()

    top.cd()
    if len(params) > 0:
        overall.Print("fit_norm_{:.2f}_{:.2f}.png".format(params[0],params[1]))
    else:
        overall.Print("fit_muonnorm.png")

def PlotFluxMarginalizationEffects(histogram,parameters,name="",plotSamples=False,usePseudo=False):
    histogram.SetPlottingStyle()
    
    invCov=histogram.GetInverseCovarianceMatrix()

    nullSolution,nullPen = FluxSolution(histogram,invCov=invCov,usePseudo=usePseudo)
    oscSolution,oscPen   = FluxSolution(histogram,invCov=invCov,useOsc=True,usePseudo=usePseudo)

    nullSols = ROOT.TH1D("null_solutions","Null Solutions",20,-1,1)
    oscSols = ROOT.TH1D("osc_solutions","Best Fit Solutions",20,-1,1)
    for i in range(len(nullSolution)):
        nullSols.Fill(nullSolution[i])
        oscSols.Fill(oscSolution[i])

    nullSolution = nullSolution/np.max(nullSolution)
    oscSolution = oscSolution/np.max(oscSolution)

    chi2_null_marg,null_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,usePseudo=usePseudo)
    chi2_model_marg,model_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,useOsc=True,usePseudo=usePseudo)
    chi2_null,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,usePseudo=usePseudo)
    chi2_model,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,useOsc=True,usePseudo=usePseudo)
    h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()

    h_osc = histogram.GetOscillatedHistogram()
    _,_ = Chi2DataMC(histogram,invCov=invCov,useOsc=True,setHists=True,marginalize=True,usePseudo=usePseudo)
    h_osc_marg = histogram.GetOscillatedHistogram()
    h_null = histogram.GetMCHistogram()
    _,_ = Chi2DataMC(histogram,invCov=invCov,setHists=True,marginalize=True,usePseudo=usePseudo)
    h_null_marg = histogram.GetMCHistogram()

    #c1 = ROOT.TCanvas()
    c1 = ROOT.TCanvas("CMarg", "canvas", 1024, 640)
    c1.SetFillStyle(4000)

    lMargin = 0.12
    rMargin = 0.05
    bMargin = 0.15
    tMargin = 0.05

    c1.Divide(1,3)

    # ------------------------------------ First Pad -----------------------------------------
    c1.cd(1)
    pad = ROOT.gPad
    pad.SetMargin(lMargin, rMargin, bMargin, tMargin)

    nullSols.GetYaxis().SetTitle("Frequency")
    nullSols.GetXaxis().SetTitle("Marginalization Solution Pulls")
    nullSols.GetYaxis().SetTitleOffset(1)
    nullSols.GetYaxis().SetTitleFont(43)
    nullSols.GetYaxis().SetTitleSize(16)
    nullSols.GetYaxis().SetLabelFont(43)
    nullSols.GetYaxis().SetLabelSize(16)
    nullSols.GetXaxis().SetTitleOffset(2)
    nullSols.GetXaxis().SetTitleFont(43)
    nullSols.GetXaxis().SetTitleSize(16)
    nullSols.GetXaxis().SetLabelFont(43)
    nullSols.GetXaxis().SetLabelSize(16)

    nullSols.SetLineWidth(2)
    nullSols.SetLineColor(ROOT.kRed)
    oscSols.SetLineWidth(2)
    oscSols.SetLineColor(ROOT.kBlue)

    nullSols.Draw("hist")
    #oscSols.Draw("hist l same")

    leg1 = ROOT.TLegend(.3,.2)
    leg1.AddEntry(nullSols,"Null Hypothesis Pulls","l")
    #leg1.AddEntry(oscSols,"Best Fit Pulls","l")
    #leg1.Draw()

    # ------------------------------------ Second Pad -----------------------------------------
    # plot mc - flux universe / mc

    c1.cd(2)
    pad = ROOT.gPad
    pad.SetMargin(lMargin, rMargin, bMargin, tMargin)

    straightLine = h_null.Clone()
    for whichBin in range(0, straightLine.GetXaxis().GetNbins()+1): 
        straightLine.SetBinContent(whichBin, 0)

    straightLine.SetMinimum(-.5)
    straightLine.SetMaximum(.5)
    straightLine.GetYaxis().SetTitle("#frac{Universe - CV}{CV}")
    straightLine.GetYaxis().SetTitleOffset(1)
    straightLine.GetYaxis().SetTitleFont(43)
    straightLine.GetYaxis().SetTitleSize(16)
    straightLine.GetYaxis().SetLabelFont(43)
    straightLine.GetYaxis().SetLabelSize(16)
    straightLine.GetXaxis().SetTitleOffset(2)
    straightLine.GetXaxis().SetTitleFont(43)
    straightLine.GetXaxis().SetTitleSize(16)
    straightLine.GetXaxis().SetLabelFont(43)
    straightLine.GetXaxis().SetLabelSize(16)

    straightLine.Draw("hist l")

    cv = h_null.GetCVHistoWithError()
    h_test = h_null.Clone()
    errband = h_test.GetVertErrorBand("Flux")
    for i in range(errband.GetNHists()):
        universe = errband.GetHist(i)
        universe.Add(cv,-1)
        universe.Divide(universe,cv)
        universe.SetMarkerSize(0)
        universe.DrawClone("same hist l")

    straightLine.Draw("hist l same")

    leg2 = ROOT.TLegend(.3,.2)
    leg2.AddEntry(straightLine,"Central Value","l")
    #leg2.Draw()

    # ------------------------------------ Third Pad -----------------------------------------
    # plot solution *  [mc - flux universe / mc ]

    c1.cd(3)
    pad = ROOT.gPad
    pad.SetMargin(lMargin, rMargin, bMargin, tMargin)

    line = straightLine.Clone()
    line.GetYaxis().SetTitle("#frac{solution#times(Universe - CV)}{CV}")
    line.Draw("hist l")

    cv = h_null.GetCVHistoWithError()
    errband = h_null.GetVertErrorBand("Flux")
    for i in range(errband.GetNHists()):
        universe = errband.GetHist(i)
        universe.Add(cv,-1)
        universe.Scale(nullSolution[i])
        universe.Divide(universe,cv)
        universe.SetMarkerSize(0)
        universe.DrawClone("same hist l")

    line.Draw("hist l same")

    leg3 = ROOT.TLegend(.3,.2)
    leg3.AddEntry(straightLine,"Central Value","l")
    #leg3.Draw()

    c1.Print("{}_marg_effects_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

def PlotOscillationRatios(histogram,parameters,name="",plotSamples=False,usePseudo=False):
    histogram.SetPlottingStyle()
    
    invCov=histogram.GetInverseCovarianceMatrix()

    nullSolution,nullPen = FluxSolution(histogram,invCov=invCov,usePseudo=usePseudo)
    oscSolution,oscPen   = FluxSolution(histogram,invCov=invCov,useOsc=True,usePseudo=usePseudo)

    nullSols = ROOT.TH1D("null_solutions","Null Solutions",20,-.5,.5)
    oscSols = ROOT.TH1D("osc_solutions","Best Fit Solutions",20,-.5,.5)
    for i in range(len(nullSolution)):
        nullSols.Fill(nullSolution[i])
        oscSols.Fill(oscSolution[i])

    chi2_null_marg,null_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,usePseudo=usePseudo)
    chi2_model_marg,model_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,useOsc=True,usePseudo=usePseudo)
    chi2_null,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,usePseudo=usePseudo)
    chi2_model,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,useOsc=True,usePseudo=usePseudo)
    h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()

    h_osc = histogram.GetOscillatedHistogram()
    _,_ = Chi2DataMC(histogram,invCov=invCov,useOsc=True,setHists=True,marginalize=True,usePseudo=usePseudo)
    h_osc_marg = histogram.GetOscillatedHistogram()
    h_null = histogram.GetMCHistogram()
    _,_ = Chi2DataMC(histogram,invCov=invCov,setHists=True,marginalize=True,usePseudo=usePseudo)
    h_null_marg = histogram.GetMCHistogram()

    #c1 = ROOT.TCanvas()
    n_plots = 2

    c1 = ROOT.TCanvas("C", "canvas", 1024, 640)
    c1.SetFillStyle(4000)

    Nx = 1
    Ny = 2

    lMargin = 0.12
    rMargin = 0.05
    bMargin = 0.15
    tMargin = 0.05

    CanvasPartition(c1, Nx, Ny, lMargin, rMargin, bMargin, tMargin)

    # ------------------------------------ First Pad -----------------------------------------
    pad1 = c1.FindObject("pad_{}_{}".format(0,1))
    pad1.cd()

    nullRatio = h_data.Clone()
    oscRatio  = h_osc.Clone()

    nullRatio.Divide(nullRatio,h_null)
    oscRatio.Divide(oscRatio, h_null)

    nullErrors = h_null.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
        nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
        nullErrors.SetBinContent(whichBin, 1)
    nullErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)

    straightLine = nullErrors.Clone()
    straightLine.SetLineColor(ROOT.kRed)
    straightLine.SetLineWidth(2)
    straightLine.SetFillColor(0)

    #oscErrors = h_osc.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    #for whichBin in range(0, oscErrors.GetXaxis().GetNbins()+1): 
    #    oscErrors.SetBinError(whichBin, max(oscErrors.GetBinContent(whichBin), 1e-9))
    #    oscErrors.SetBinContent(whichBin, 1)
    #oscErrors.SetFillColorAlpha(ROOT.kBlue + 1, 0.4)

    nullErrors.SetMinimum(0)
    nullErrors.SetMaximum(2)
    nullErrors.GetYaxis().SetTitle("#splitline{Ratio to}{Null Hypothesis}")
    nullErrors.GetYaxis().SetTitleOffset(1)
    nullErrors.GetYaxis().SetTitleFont(43)
    nullErrors.GetYaxis().SetTitleSize(30)
    nullErrors.GetYaxis().SetLabelFont(43)
    nullErrors.GetYaxis().SetLabelSize(30)
    nullErrors.GetXaxis().SetTitleOffset(2)
    nullErrors.GetXaxis().SetTitleFont(43)
    nullErrors.GetXaxis().SetTitleSize(30)
    nullErrors.GetXaxis().SetLabelFont(43)
    nullErrors.GetXaxis().SetLabelSize(30)

    nullRatio.SetFillColor(0)
    nullRatio.SetLineWidth(2)
    oscRatio.SetFillColor(0)
    straightLine.SetFillColor(0)

    nullErrors.Draw("E2")
    #oscErrors.Draw("E2 SAME")
    nullRatio.Draw("same")
    oscRatio.Draw("same hist l")
    straightLine.Draw("HIST L SAME")

    x1, y1, x2, y2 = XtoPad(0.05),YtoPad(0.65),XtoPad(.4),YtoPad(.95)

    leg1 = ROOT.TLegend(x1,y1,x2,y2)
    leg1.AddEntry(h_data,"Data","p")
    leg1.AddEntry(oscRatio,"Best Fit #chi^{2}="+"{:.2f}".format(chi2_model),"l")
    leg1.AddEntry(straightLine,"Null #chi^{2}="+"{:.2f}".format(chi2_null),"l")
    leg1.Draw()

    # ------------------------------------ Second Pad -----------------------------------------
    pad2 = c1.FindObject("pad_{}_{}".format(0,0))
    pad2.cd()

    nullMargRatio =  h_data.Clone()
    oscMargRatio =  h_osc_marg.Clone()

    nullMargRatio.Divide(nullMargRatio,h_null_marg)
    oscMargRatio.Divide(oscMargRatio, h_null_marg)

    h_null_marg.PopVertErrorBand("Flux")
    nullMargErrors = h_null_marg.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullMargErrors.GetXaxis().GetNbins()+1): 
        nullMargErrors.SetBinError(whichBin, max(nullMargErrors.GetBinContent(whichBin), 1e-9))
        nullMargErrors.SetBinContent(whichBin, 1)
    nullMargErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)

    #h_osc_marg.PopVertErrorBand("Flux")
    #oscMargErrors = h_osc_marg.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    #for whichBin in range(0, oscMargErrors.GetXaxis().GetNbins()+1): 
    #    oscMargErrors.SetBinError(whichBin, max(oscMargErrors.GetBinContent(whichBin), 1e-9))
    #    oscMargErrors.SetBinContent(whichBin, 1)
    #oscMargErrors.SetFillColorAlpha(ROOT.kBlue + 1, 0.4)

    nullMargErrors.SetMinimum(0)
    nullMargErrors.SetMaximum(2)
    nullMargErrors.GetYaxis().SetTitle("#splitline{Ratio to}{Marg'd Hypothesis}")
    nullMargErrors.GetYaxis().SetTitleOffset(1)
    nullMargErrors.GetYaxis().SetTitleFont(43)
    nullMargErrors.GetYaxis().SetTitleSize(30)
    nullMargErrors.GetYaxis().SetLabelFont(43)
    nullMargErrors.GetYaxis().SetLabelSize(30)
    nullMargErrors.GetXaxis().SetTitleOffset(2)
    nullMargErrors.GetXaxis().SetTitleFont(43)
    nullMargErrors.GetXaxis().SetTitleSize(30)
    nullMargErrors.GetXaxis().SetLabelFont(43)
    nullMargErrors.GetXaxis().SetLabelSize(30)
    nullMargRatio.SetFillColor(0)
    nullMargRatio.SetLineWidth(2)
    oscMargRatio.SetFillColor(0)

    nullMargErrors.Draw("E2")
    #oscMargErrors.Draw("E2 SAME")
    nullMargRatio.Draw("same")
    oscMargRatio.Draw("same hist l")
    straightLine.Draw("HIST L SAME")

    x1, y1, x2, y2 = XtoPad(0.05),YtoPad(0.65),XtoPad(.4),YtoPad(.95)

    leg2 = ROOT.TLegend(x1,y1,x2,y2)
    leg2.AddEntry(h_data,"Data","p")
    leg2.AddEntry(oscMargRatio,"Best Fit #chi^{2}="+"{:.2f}".format(chi2_model_marg),"l")
    leg2.AddEntry(straightLine,"Null #chi^{2}="+"{:.2f}".format(chi2_null_marg),"l")
    leg2.Draw()

    c1.Print("{}_ratios_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

def PlotOscillationEffects(histogram,parameters,name="",plotSamples=False,usePseudo=False):
    invCov=histogram.GetInverseCovarianceMatrix()

    h_null = histogram.GetMCHistogram()
    h_osc = histogram.GetOscillatedHistogram()
    h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()

    chi2_null,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,usePseudo=usePseudo)
    chi2_model,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,useOsc=True,usePseudo=usePseudo)

    chi2_model_marg,model_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,useOsc=True,usePseudo=usePseudo,setHists=True)
    chi2_null_marg,null_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,usePseudo=usePseudo,setHists=True)

    h_null_marg = histogram.GetMCHistogram()
    h_osc_marg = histogram.GetOscillatedHistogram()

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
    

    h_osc.SetLineColor(ROOT.kBlue)
    h_osc.SetLineWidth(2)
    h_osc.SetMarkerStyle(0)
    h_osc_marg.SetLineColor(ROOT.kBlue)
    h_osc_marg.SetLineWidth(2)
    h_osc_marg.SetMarkerStyle(0)
    h_osc_marg.SetLineStyle(2)

    h_null.SetLineColor(ROOT.kRed)
    h_null.SetLineWidth(2)
    h_null_marg.SetLineColor(ROOT.kPink)
    h_null_marg.SetLineWidth(2)
    h_null_marg.SetLineStyle(2)
    h_null.GetYaxis().SetTitle("Nevents")
    h_null.Draw("hist")
    h_null_marg.Draw("hist same")
    h_osc.Draw("hist same")
    h_osc_marg.Draw("hist_ same")
    h_data.Draw("same")

    null = h_null.GetCVHistoWithError()
    null.SetLineColor(ROOT.kRed)
    null.SetLineWidth(2)
    null.SetMarkerStyle(0)
    null.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
    null.Draw("E2 SAME")
    null_marg = h_null_marg.GetCVHistoWithError()
    null_marg.SetLineColor(ROOT.kPink)
    null_marg.SetLineWidth(2)
    null_marg.SetLineStyle(2)
    null_marg.SetMarkerStyle(0)
    null_marg.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
    null_marg.Draw("E2 SAME")

    osc = h_osc.GetCVHistoWithError()
    osc.SetLineColor(ROOT.kBlue)
    osc.SetLineWidth(2)
    osc.SetMarkerStyle(0)
    osc.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
    osc.Draw("E2 SAME")
    osc_marg = h_osc_marg.GetCVHistoWithError()
    osc_marg.SetLineColor(ROOT.kBlue)
    osc_marg.SetLineWidth(2)
    osc_marg.SetLineStyle(2)
    osc_marg.SetMarkerStyle(0)
    osc_marg.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
    osc_marg.Draw("E2 SAME")

    top_text = "#Delta m^{2}="+"{:.1f}".format(parameters["m"])+" |U_{e4}|^{2}="+"{:.4f}".format(parameters["ue4"])
    bot_text = "|U_{#mu4}|^{2}="+"{:.5f}".format(parameters["umu4"])+" |U_{#tau4}|^{2}="+"{:.1f}".format(parameters["utau4"])
    header = "#splitline{%s}{%s}" % (top_text,bot_text)
    MNVPLOTTER.AddPlotLabel(header,.35,.25)

    leg = ROOT.TLegend(.28,.4)

    #leg.SetBorderSize(0)
    leg.AddEntry(h_data,"Data","p")
    top_text = "Null Hypothesis"
    bot_text = "#chi^{2}="+"{:.2f}".format(chi2_null)
    leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    leg.AddEntry(h_null,leg_text,"l")
    top_text = "Null Hypothesis Marg."
    bot_text = "#chi^{2}="+"{:.2f} + {:.2f} penalty".format(chi2_null_marg-null_pen,null_pen)
    leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    leg.AddEntry(h_null_marg,leg_text,"l")
    top_text = "Best Osc. Fit"
    bot_text = "#chi^{2}="+"{:.2f}".format(chi2_model)
    leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    leg.AddEntry(h_osc,leg_text,"l")
    top_text = "Best Osc. Fit Marg."
    bot_text = "#chi^{2}="+"{:.2f} + {:.2f} penalty".format(chi2_model_marg-model_pen,model_pen)
    leg_text = "#splitline{%s}{%s}" % (top_text,bot_text)
    leg.AddEntry(h_osc_marg,leg_text,"l")
    leg.Draw()

    nullRatio =  h_data.Clone()
    nullMargRatio =  h_null_marg.Clone()
    oscRatio =  h_osc.Clone()
    oscMargRatio =  h_osc_marg.Clone()

    nullRatio.Divide(nullRatio,h_null)
    nullMargRatio.Divide(nullMargRatio,h_null)
    oscRatio.Divide(oscRatio, h_null)
    oscMargRatio.Divide(oscMargRatio, h_null)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    nullErrors = h_null.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
        nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
        nullErrors.SetBinContent(whichBin, 1)

    nullRatio.SetLineColor(ROOT.kBlack)
    nullRatio.SetLineWidth(3)

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
    nullErrors.SetMinimum(0)
    nullErrors.SetMaximum(2)
    
    #Draw the data ratios
    nullRatio.SetMinimum(0)
    nullRatio.SetMaximum(2)
    nullRatio.Draw("same")
    nullMargRatio.Draw("same hist l")
    oscRatio.Draw('same hist l')
    oscMargRatio.Draw('same hist l')

    #Draw a flat line at 1 for oscRatio of MC to itself
    straightLine = nullErrors.Clone()
    straightLine.SetLineColor(ROOT.kRed)
    straightLine.SetLineWidth(2)
    straightLine.SetFillColor(0)
    straightLine.Draw("HIST L SAME")

    top.cd()

    overall.Print("plots/{}_stitched_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

    plots = histogram.keys
    nullSolution,nullPen = FluxSolution(histogram,invCov=invCov,usePseudo=usePseudo)
    plots.append("fhc_ratio")
    plots.append("rhc_ratio")
    histogram.titles.append("FHC CC #nu_{#mu}/#nu_{e} Ratio")
    histogram.titles.append("RHC CC anti #nu_{#mu}/#nu_{e} Ratio")

    for i,plot in enumerate(plots):
        title = histogram.titles[i]
        margin = .12
        bottomFraction = .2
        overall = ROOT.TCanvas(plot)
        top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
        bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

        top.Draw()
        bottom.Draw()

        top.cd()

        if 'ratio' in plot:
            numu_hists,numu_total_hist = OscillateSubHistogram(histogram,plot[:3]+'_numu_selection',parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])
            nue_hists,nue_total_hist = OscillateSubHistogram(histogram,plot[:3]+'_nue_selection',parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])

            total_hist = numu_total_hist.Clone()
            total_hist.Divide(total_hist,nue_total_hist)
            swap_fraction = total_hist.Clone()
            norm_fraction = total_hist.Clone()

            hists = []
            for j,nue in enumerate(nue_hists):
                if "#mu" in nue.GetTitle():
                    swap_fraction.Divide(nue,swap_fraction)
                    norm_fraction.Add(swap_fraction,-1)

                    norm_fraction.SetTitle("#nu_{#mu}/#nu_{e}")
                    swap_fraction.SetTitle("#nu_{#mu}/"+nue.GetTitle())

                    swap_fraction.SetLineColor(ROOT.kBlue)
                    swap_fraction.SetFillColor(ROOT.kBlue)

                    hists.append(norm_fraction)
                    hists.append(swap_fraction)

            h_data = histogram.data_samples[plot[:3]+'_numu_selection'].Clone()
            h_data.Divide(h_data,histogram.data_samples[plot[:3]+'_nue_selection'])
            subSample = histogram.mc_samples[plot[:3]+'_numu_selection'].Clone()
            subSample.Divide(subSample,histogram.mc_samples[plot[:3]+'_nue_selection'])
        else:
            hists,total_hist = OscillateSubHistogram(histogram,plot,parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])
            h_data = histogram.data_samples[plot].Clone()
            subSample = histogram.mc_samples[plot].Clone()

        cv = np.array(subSample)[1:-1]
        mc = np.array(total_hist)[1:-1]
        band = subSample.GetVertErrorBand("Flux")
        nhists = band.GetNHists()
        universes = np.array([np.array(band.GetHist(l))[1:-1] for l in range(nhists)])
        cv_table = np.array([cv for l in range(len(universes))])
        A = universes - cv_table

        new_cv = mc + nullSolution @ A

        weights = total_hist.GetCVHistoWithStatError()
        for j in range(1,weights.GetNbinsX()+1):
            weight = weights.GetBinContent(j) / new_cv[j-1] if new_cv[j-1] != 0 else weights.GetBinContent(j)
            weights.SetBinContent(j,weight)
            weights.SetBinError(j,0)

        total_hist.DivideSingle(total_hist,weights)
        total_hist.PopVertErrorBand("Flux")

        TArray = ROOT.TObjArray()
        for hist in hists:
            hist.DivideSingle(hist,weights)
            if 'elastic' in plot:
                hist.Scale(2,'width')
            elif 'ratio' not in plot:
                hist.Scale(1,'width')
            TArray.Add(hist)

        if "elastic" in plot:
            if not histogram.data_scaled[plot]:
                h_data.Scale(2,'width')
            h_data.GetXaxis().SetTitle("Electron Energy [ GeV ]")
            h_data.GetYaxis().SetTitle("NEvents / 2 GeV")
            total_hist.Scale(2,'width')
        elif 'ratio' not in plot:
            if not histogram.data_scaled[plot]:
                h_data.Scale(1,'width')
            if "imd" in plot:
                h_data.GetXaxis().SetTitle("Muon Energy [ GeV ]")
            else:
                h_data.GetXaxis().SetTitle("Neutrino Energy Estimator [ GeV ]")

            h_data.GetYaxis().SetTitle("NEvents / GeV")
            total_hist.Scale(1,'width')
        else:
            h_data.GetXaxis().SetTitle("Neutrino Energy Estimator [ GeV ]")

        MNVPLOTTER.DrawDataStackedMC(h_data,TArray,1.,"TR","Data",-1,-1,-1) 
        for obj in top.GetListOfPrimitives():
            if type(obj) == ROOT.TH1D:
                obj.SetTitle(title)

        bottom.cd()
        bottom.SetTopMargin(0)
        bottom.SetBottomMargin(0.3)

        ratio = h_data.Clone()
        ratio.Divide(ratio, total_hist)

        #Now fill mcRatio with 1 for bin content and fractional error
        mcRatio = total_hist.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
        for whichBin in range(1, mcRatio.GetXaxis().GetNbins()+1): 
            mcRatio.SetBinError(whichBin, max(mcRatio.GetBinContent(whichBin), 1e-9))
            mcRatio.SetBinContent(whichBin, 1)

        ratio.SetTitle("")

        ratio.GetYaxis().SetTitle("Data/MC")
        ratio.GetYaxis().SetLabelSize(.13)
        ratio.GetYaxis().SetTitleSize(0.1)
        ratio.GetYaxis().SetTitleOffset(0.6)
        ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
        ratio.GetXaxis().SetTitleSize(0.16)
        ratio.GetXaxis().SetTitleOffset(0.9)
        ratio.GetXaxis().SetLabelSize(.15)
        ratio.SetMinimum(0)
        ratio.SetMaximum(2)
        ratio.SetLineWidth(1)
        ratio.Draw('E1 X0')

        #Error envelope for the MC
        mcRatio.SetLineColor(ROOT.kRed)
        mcRatio.SetLineWidth(3)
        mcRatio.SetMarkerStyle(0)
        mcRatio.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
        mcRatio.SetFillStyle(1001)
        mcRatio.Draw("E2 SAME")

        straightLine = mcRatio.Clone()
        straightLine.SetFillStyle(0)
        straightLine.Draw("HIST SAME")
        ROOT.gStyle.SetOptTitle(1)

        overall.Print("plots/{}_oscillated_{:.1f}_{:.3f}_{:.4f}.png".format(plot,parameters["m"],parameters['ue4'],parameters['umu4']))

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
    #top.SetTopMargin(0.1)

    if mnv_data.GetBinContent(mnv_data.GetMaximumBin()) > 1e4:
        top.SetLogy()

    MNVPLOTTER.axis_minimum = 0.1
    cv_data.GetXaxis().SetTitle("Bin number")
    cv_data.GetYaxis().SetTitle("Entries")

    ratio = cv_data.Clone()
    data = cv_data.GetCVHistoWithError()
    data.SetMarkerColor(ROOT.kBlack)
    data.SetTitle(mnv_mc.GetName())
    #data.SetTitleSize(0.15,'t')
    #data.SetTitleOffset(0.1)
    #data.SetLabelSize(0.15)
    mnv_mc.SetLineColor(ROOT.kBlue)
    mnv_mc.SetLineWidth(3)
    mnv_mc.SetMarkerStyle(0)

    data.Draw("P")
    mnv_mc.Draw("HIST SAME")
    cv_mc.SetLineColor(ROOT.kRed)
    cv_mc.SetLineWidth(3)
    cv_mc.SetMarkerStyle(0)
    #cv_data.Draw("P SAME")
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
    chi2_model = Chi2DataMC(mnv_data,mnv_mc)
    leg.SetHeader("#chi^{2}_{model}="+"{:.2f}".format(chi2_model),"C")
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
    ratio.SetMinimum(0)
    ratio.SetMaximum(2)
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
    ROOT.gStyle.SetOptTitle(1)
    c1.Update()
    #MNVPLOTTER.WritePreliminary(0.4, 0.05, 5e-2, True)
    overall.Print("plots/{}.png".format(mnv_mc.GetName()))

def DataMCCVPlot(cv_data,cv_mc,title="mc_stitched.png"):
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
    if title != "mc_stitched.png":
        cv_data.GetXaxis().SetTitle("Energy Estimator (GeV)")
        cv_data.GetYaxis().SetTitle("Muon / Electron Flavor Ratio")
    else:
        cv_data.GetXaxis().SetTitle("Bin Number")
        cv_data.GetYaxis().SetTitle("Entries")
    cv_data.GetXaxis().SetTitle("Bin Number")
    cv_data.GetYaxis().SetTitle("Entries")
    top.SetLogy()
    MNVPLOTTER.axis_minimum = 0.1

    ratio = cv_data.GetCVHistoWithError().Clone()

    cv_data.SetMarkerColor(ROOT.kBlack)
    cv_data.SetLineColor(ROOT.kBlack)
    cv_data.SetMarkerSize(1)
    cv_data.SetMarkerStyle(20)
    cv_data.SetLineWidth(3)
    cv_data.SetMinimum(0.1)

    cv_data.GetCVHistoWithError().Draw("P")
    cv_mc.SetLineColor(ROOT.kRed)
    cv_mc.SetLineWidth(3)
    cv_mc.SetMarkerStyle(0)
    cv_data.Draw("P SAME")
    cv_mc.Draw("HIST SAME")
    cv = cv_mc.GetCVHistoWithError()
    cv.SetLineColor(ROOT.kRed)
    cv.SetLineWidth(3)
    cv.SetMarkerStyle(0)
    cv.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
    cv.Draw("E2 SAME")

    leg = ROOT.TLegend(.65,.35,.8,.55)
    leg.AddEntry(cv_data,"Data",'p')
    leg.AddEntry(cv_mc,"Null Hypothesis",'l')
    leg.Draw()

    ratio.Divide(ratio, cv_mc)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    #Now fill mcRatio with 1 for bin content and fractional error
    mcRatio = cv_mc.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
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
    if title != "mc_stitched.png":
        ratio.GetXaxis().SetTitle("Energy Estimator (GeV)")
    else:
        ratio.GetXaxis().SetTitle("Bin Number")
    ratio.SetMinimum(0)
    ratio.SetMaximum(2)
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
    
    MNVPLOTTER.WritePreliminary(0.4, 0.05, 5e-2, True)
    overall.Print("plots/{}".format(title))
