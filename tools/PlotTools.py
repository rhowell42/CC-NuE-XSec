import math
import ROOT
import PlotUtils
<<<<<<< HEAD
import heapq
import ctypes
=======
import numpy as np
from scipy import optimize, integrate
import argparse
ccnueroot = os.environ.get('CCNUEROOT')
plotutils = os.environ.get('PLOTUTILSROOT')

import math
import psutil
import multiprocessing
import threading
nthreads = 4
>>>>>>> feature/sterile_neutrino
from array import array
from Tools.FitTools import *
from Tools.Histogram import *
from Tools.Helper import *

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder

<<<<<<< HEAD
from config.AnalysisConfig import AnalysisConfig
from config.SignalDef import SIGNAL_DEFINATION
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
=======
logging.basicConfig(stream=sys.stderr, level=logging.INFO)
>>>>>>> feature/sterile_neutrino

MNVPLOTTER = PlotUtils.MnvPlotter()
MNVPLOTTER.draw_normalized_to_bin_width=False
<<<<<<< HEAD
MNVPLOTTER.legend_text_size = 0.02
MNVPLOTTER.extra_top_margin = -0.035# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46 
MNVPLOTTER.mc_bkgd_line_color = 46
=======
>>>>>>> feature/sterile_neutrino

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)
legend_text_size = .025

<<<<<<< HEAD
#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .03
MNVPLOTTER.legend_n_columns = 1

CANVAS = ROOT.TCanvas("c2","c2",1200,1000)
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
MNVPLOTTER.error_summary_group_map.clear();
for k,v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v :
        vec.push_back(vs)
    MNVPLOTTER.error_summary_group_map[k]= vec
=======
def UndoBinWidthNorm(histogram):
    for i in range(0,histogram.GetNbinsX()+1):
        width = histogram.GetBinWidth(i)
        cont = histogram.GetBinContent(i)
        new_cont = cont*width
        if new_cont != 0:
            histogram.SetBinContent(i,new_cont)

        for name in histogram.GetVertErrorBandNames():
            band = histogram.GetVertErrorBand(name)
            cont = band.GetBinContent(i)
            new_cont = cont*width
            if new_cont != 0:
                band.SetBinContent(i,new_cont)
>>>>>>> feature/sterile_neutrino

            for j in range(band.GetNHists()):
                hist = band.GetHist(j)
                cont = hist.GetBinContent(i)
                new_cont = cont*width
                if new_cont != 0:
                    hist.SetBinContent(i,new_cont)


<<<<<<< HEAD
def Logz(canvas):
    canvas.SetLogz(1)
    return True

def PrepareSlicer(hist_holder):
    if hist_holder.dimension==1:
        return lambda x:[x.Clone()]
    elif hist_holder.dimension == 2:
        return Make2DSlice
    else:
        return None
        #raise KeyError("Only 1D,2D histograms are supported")
=======
def PlotSampleMarginalizationEffects(sample_histogram):
    histogram = copy.deepcopy(sample_histogram)
    samplesToFit = [["elastic","imd","numu","ratio"],["elastic"],["elastic","imd"],["numu"],["elastic","imd","numu"],["elastic","imd","ratio"]]
    names = ["All","elastic","elastic_imd","numu","elastic_imd_numu","elastic_imd_ratio"]
    titles = ["#nu+e, IMD, CC #nu_{#mu}, CC #nu_{#mu}/#nu_{e}","#nu+e","#nu+e, IMD","CC #nu_{#mu}","#nu+e, IMD, CC #nu_{#mu}","#nu+e, IMD, CC #nu_{#mu}/#nu_{e}"]
    colors = [ROOT.kRed,ROOT.kBlue,ROOT.kViolet,ROOT.kOrange,ROOT.kTeal,ROOT.kGreen]
    hists = []
    flux_hists = []
    chi2s = []
    pens = []
>>>>>>> feature/sterile_neutrino

    subSample = histogram.mc_hists["fhc_elastic"].Clone()
    band = subSample.GetVertErrorBand("Flux")
    nhists = band.GetNHists()
    universes = np.array([np.array(band.GetHist(l))[1:-1] for l in range(nhists)])
    fhc_integrals = universes.sum(axis=1)
    old_fhc = subSample.Integral()

    subSample = histogram.mc_hists["rhc_elastic"].Clone()
    band = subSample.GetVertErrorBand("Flux")
    nhists = band.GetNHists()
    universes = np.array([np.array(band.GetHist(l))[1:-1] for l in range(nhists)])
    rhc_integrals = universes.sum(axis=1)
    old_rhc = subSample.Integral()

    c0 = ROOT.TCanvas()

<<<<<<< HEAD
def PrepareSignalDecompose(data_hists,mc_hists,cates,sys_on_mc = True, sys_on_data = False,ratio=False,only_cated = False):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew

    if (data_hists.valid and mc_hists.valid ):
        hists = [data_hists.GetHist().GetCVHistoWithError() if sys_on_data else data_hists.GetHist().GetCVHistoWithStatError()]
        cate_hists,colors,titles  = mc_hists.GetCateList(cates)
        if only_cated :
            totalHist = ReSumHists(cate_hists)
=======
    mg = ROOT.TMultiGraph()
    mg.SetTitle(";#nu_{#mu}-mode Number of #nue events;#bar{#nu}_{#mu}-mode Number of #nue events")
    mg.GetHistogram().GetXaxis().SetLimits(1000,1800)
    mg.GetHistogram().GetYaxis().SetRangeUser(650,1200)

    g1 = ROOT.TGraph(len(fhc_integrals),fhc_integrals,rhc_integrals)
    ROOT.SetOwnership(g1, False)
    g1.SetTitle("Flux Universes")
    g1.SetMarkerStyle(4)
    g1.SetMarkerColorAlpha(ROOT.kBlack,1)
    g1.SetLineWidth(0)
    g2 = ROOT.TGraph()
    ROOT.SetOwnership(g2, False)
    g2.SetPoint(0,old_fhc,old_rhc)
    g2.SetTitle("Central Value")
    g2.SetMarkerStyle(41)
    g2.SetMarkerSize(3)
    g2.SetMarkerColor(ROOT.kRed)
    g2.SetLineWidth(0)
    g3 = ROOT.TGraphErrors(1,np.array([1338.0]),np.array([838.2]),np.array([99.7]),np.array([63.1]))
    ROOT.SetOwnership(g3, False)
    g3.SetTitle("Mean/RMS Before Constraint Result")
    g3.SetMarkerStyle(20)
    g3.SetMarkerColor(ROOT.kBlack)
    g4 = ROOT.TGraphErrors(1,np.array([1239.0]),np.array([778.4]),np.array([41.9]),np.array([31.2]))
    ROOT.SetOwnership(g4, False)
    g4.SetTitle("Mean/RMS After Constraint Result")
    g4.SetMarkerStyle(20)
    g4.SetMarkerColor(ROOT.kBlack)

    mg.Add(g1)
    mg.Add(g2)
    mg.Add(g3)
    mg.Add(g4)

    f_numu = ROOT.TFile.Open(plotutils+'/data/flux/flux-g4numiv6-pdg14-minervame1D1M1NWeightedAve.root')
    fhc_numu = f_numu.Get("flux_E_unweighted")
    f_numu.Close()
    f_nue = ROOT.TFile.Open(plotutils+'/data/flux/flux-g4numiv6-pdg12-minervame1D1M1NWeightedAve.root')
    fhc_nue = f_nue.Get("flux_E_unweighted")
    f_nue.Close()
    fhc_numu_univ = fhc_numu.GetVertErrorBand("Flux")
    fhc_nue_univ = fhc_nue.GetVertErrorBand("Flux")

    for i,samples in enumerate(samplesToFit):
        new_histogram = StitchedHistogram(names[i],True)

        for key in histogram.keys + ["fhc_ratio","rhc_ratio"]:
            for sample in samples:
                if sample in key:
                    if 'ratio' in sample:
                        for beam in ['fhc','rhc']:
                            h_mc_numu = histogram.mc_hists[beam+"_numu_selection"].Clone()
                            h_data_numu = histogram.data_hists[beam+"_numu_selection"].Clone()
                            h_mc_nue = histogram.mc_hists[beam+"_nue_selection"].Clone()
                            h_data_nue = histogram.data_hists[beam+"_nue_selection"].Clone()
                            h_mc_numu.Divide(h_mc_numu,h_mc_nue)
                            h_data_numu.Divide(h_data_numu,h_data_nue)
                            new_histogram.AddHistograms(beam+"_ratio",h_mc_numu.Clone(),h_data_numu.Clone())
                    else:
                        new_histogram.AddHistograms(key,histogram.mc_hists[key].Clone(),histogram.data_hists[key].Clone())

        new_histogram.Stitch()
        invCov=new_histogram.GetInverseCovarianceMatrix(sansFlux=True)

        fluxSolution,nullPen = FluxSolution(new_histogram,invCov=invCov)
        invCov = new_histogram.GetInverseCovarianceMatrix(sansFlux=True)
        plot_histogram = copy.deepcopy(histogram)
        chi2,penalty = Chi2DataMC(plot_histogram,fluxSolution=fluxSolution,invCov=invCov,setHists=True,marginalize=True)

        hist = plot_histogram.GetMCHistogram()
        hists.append(hist)
        chi2s.append(chi2)
        pens.append(penalty)

        subSample = histogram.mc_hists["fhc_elastic"].Clone()
        weights = ReweightCV(subSample,fluxSolution=fluxSolution)
        new_fhc = subSample.Integral()

        subSample = histogram.mc_hists["rhc_elastic"].Clone()
        weights = ReweightCV(subSample,fluxSolution=fluxSolution)
        new_rhc = subSample.Integral()

        new_fhc_numu = fhc_numu.Clone()
        new_fhc_nue = fhc_nue.Clone()
        weights = ReweightCV(new_fhc_numu,fluxSolution=fluxSolution)
        weights = ReweightCV(new_fhc_nue,fluxSolution=fluxSolution)

        flux_hists.append(new_fhc_nue)

        g_ = ROOT.TGraph()
        ROOT.SetOwnership(g_, False)
        g_.SetPoint(0,new_fhc,new_rhc)
        g_.SetTitle(titles[i])
        g_.SetMarkerStyle(29)
        g_.SetMarkerSize(3)
        g_.SetLineWidth(0)
        g_.SetMarkerColorAlpha(colors[i],0.4)
        mg.Add(g_)

    mg.Draw("AP")
    pad = ROOT.gPad
    pad.BuildLegend()
    c0.Print("plots/integrated_elastic_events.png")

    c1 = ROOT.TCanvas("C2", "canvas2", 1024, 640)
    c1.cd(1)
    h_null =  histogram.GetMCHistogram()
    h_data = histogram.GetDataHistogram()
    invCov = histogram.GetInverseCovarianceMatrix(sansFlux=True)
    chi2,pen = Chi2DataMC(histogram,invCov=invCov)

    nullRatio = h_data.Clone()
    nullRatio.Divide(nullRatio,h_null)

    nullErrors = h_null.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
        nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
        nullErrors.SetBinContent(whichBin, 1)
    nullErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)

    straightLine = nullErrors.Clone()
    straightLine.SetFillStyle(0)

    RatioAxis(nullErrors,MNVPLOTTER)
    nullErrors.SetMinimum(0)
    nullErrors.SetMaximum(2)
    nullErrors.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")
    nullErrors.GetXaxis().SetTitleOffset(1.5)
    nullErrors.Draw("E2")
    nullRatio.Draw("SAME")

    c1.SetTopMargin(0.35)
    c1.SetRightMargin(0.05)
    leg = ROOT.TLegend(0.05, 0.675, 0.95, .975)
    leg.SetTextSize(legend_text_size);
    leg.SetNColumns(len(titles)//2) # // median N number of rows per column
    top = "Data"
    bottom = "#chi^{2}="+"{:.2f}".format(chi2)
    legend = "#splitline{%s}{%s}" % (top,bottom)
    leg.AddEntry(nullRatio,legend,"p")

    for i,hist in enumerate(hists):
        hist.Divide(hist,h_null)
        hist.SetLineColor(colors[i])
        hist.SetFillStyle(0)
        if i == 0:
            hist.SetLineStyle(2)
        hist.Draw("same hist l")

        top = titles[i]
        bottom = "#chi^{2}="+"{:.2f}+".format(chi2s[i]-pens[i])+"{:.2f} pen.".format(pens[i])
        legend = "#splitline{%s}{%s}" % (top,bottom)
        leg.AddEntry(hist,legend,"l")

    straightLine.Draw("same hist") 

    leg.Draw()

    c1.Print("plots/stitched_flux_marg_effects.png")

    new_bins = array('d',list(range(0,21)))
    UndoBinWidthNorm(fhc_nue)
    UndoBinWidthNorm(fhc_numu)

    fhc_numu = fhc_nue.Rebin(20,"hnew",new_bins)
    fhc_numu.Scale(1,"width")
    fhc_nue = fhc_nue.Rebin(20,"hnew",new_bins)
    fhc_nue.Scale(1,"width")

    for i,hist in enumerate(flux_hists):
        UndoBinWidthNorm(flux_hists[i])
        flux_hists[i] = hist.Rebin(20,str(i),new_bins)
        flux_hists[i].Scale(1,'width')

    fhc_numu.GetXaxis().SetRangeUser(0,20)
    fhc_numu.GetXaxis().SetTitle("Neutrino Energy")
    fhc_numu.SetTitle("FHC #nu_{#mu} Flux Prediction")

    fhc_nue.GetXaxis().SetRangeUser(0,20)
    fhc_nue.GetXaxis().SetTitle("Neutrino Energy")
    fhc_nue.SetTitle("FHC #nu_{e} Flux Prediction")

    PlotWithRatio(MNVPLOTTER,"plots/FluxReweighting.png",fhc_nue,hists=flux_hists,titles=titles,colors=colors)

def PlotFluxMarginalizationEffects(sample_histogram,parameters,name="",plotSamples=False,usePseudo=False):
    histogram = copy.deepcopy(sample_histogram)
    histogram.SetPlottingStyle()
    
    invCov=histogram.GetInverseCovarianceMatrix(sansFlux=True)

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
    RatioAxis(nullSols,MNVPLOTTER)

    nullSols.SetLineWidth(2)
    nullSols.SetLineColor(ROOT.kRed)
    oscSols.SetLineWidth(2)
    oscSols.SetLineColor(ROOT.kBlue)

    nullSols.Draw("hist")
    oscSols.Draw("hist same")

    leg1 = ROOT.TLegend(.3,.2)
    leg1.AddEntry(nullSols,"Null Hypothesis Pulls","l")
    leg1.AddEntry(oscSols,"Osc. Model Pulls","l")
    leg1.Draw()

    # ------------------------------------ Second Pad -----------------------------------------
    # plot mc - flux universe / mc

    c1.cd(2)
    pad = ROOT.gPad
    pad.SetMargin(lMargin, rMargin, bMargin, tMargin)

    straightLine = h_null.Clone()
    for whichBin in range(0, straightLine.GetXaxis().GetNbins()+1): 
        straightLine.SetBinContent(whichBin, 0)

    RatioAxis(straightLine,MNVPLOTTER)
    straightLine.SetMinimum(-.5)
    straightLine.SetMaximum(.5)
    straightLine.GetYaxis().SetTitle("#frac{Universe - CV}{CV}")

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
    leg3.Draw()

    c1.Print("plots/{}_marg_effects_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

def PlotOscillationRatios(sample_histogram,parameters,name="",plotSamples=False,usePseudo=False):
    histogram = copy.deepcopy(sample_histogram)
    histogram.SetPlottingStyle()
    
    invCov=histogram.GetInverseCovarianceMatrix(sansFlux=True)

    nullSolution,nullPen = FluxSolution(histogram,invCov=invCov,usePseudo=usePseudo)
    oscSolution,oscPen   = FluxSolution(histogram,invCov=invCov,useOsc=True,usePseudo=usePseudo)

    nullSols = ROOT.TH1D("null_solutions","Null Solutions",20,-.5,.5)
    oscSols = ROOT.TH1D("osc_solutions","Best Fit Solutions",20,-.5,.5)
    for i in range(len(nullSolution)):
        nullSols.Fill(nullSolution[i])
        oscSols.Fill(oscSolution[i])

    h_null = histogram.GetMCHistogram()
    h_osc = histogram.GetOscillatedHistogram()
    h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()

    chi2_null,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,usePseudo=usePseudo)
    chi2_model,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,useOsc=True,usePseudo=usePseudo)

    chi2_model_marg,model_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,useOsc=True,usePseudo=usePseudo,setHists=True)
    chi2_null_marg,null_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,usePseudo=usePseudo,setHists=True)

    h_null_marg = histogram.GetMCHistogram()
    h_osc_marg = histogram.GetOscillatedHistogram()

    c1 = ROOT.TCanvas("C", "canvas", 1024, 640)
    c1.SetFillStyle(4000)

    Nx = 1
    Ny = 3

    lMargin = 0.12
    rMargin = 0.05
    bMargin = 0.15
    tMargin = 0.05

    CanvasPartition(c1, Nx, Ny, lMargin, rMargin, bMargin, tMargin)

    # ------------------------------------ First Pad -----------------------------------------
    pad1 = c1.FindObject("pad_{}_{}".format(0,2))
    pad1.cd()

    nullRatio = h_data.Clone()
    oscRatio  = h_osc.Clone()
    margRatio = h_null_marg.Clone()
    flux_null = h_null.Clone()

    nullRatio.Divide(nullRatio,h_null)
    oscRatio.Divide(oscRatio, h_null)
    margRatio.Divide(margRatio, h_null)

    for err in flux_null.GetVertErrorBandNames():
        if err != "Flux":
            flux_null.PopVertErrorBand(err)

    fluxErrors = flux_null.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, fluxErrors.GetXaxis().GetNbins()+1): 
        fluxErrors.SetBinError(whichBin, max(fluxErrors.GetBinContent(whichBin), 1e-9))
        fluxErrors.SetBinContent(whichBin, 1)
    fluxErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)

    fluxCV = fluxErrors.Clone()
    fluxCV.SetLineColor(ROOT.kRed)
    fluxCV.SetLineWidth(2)
    fluxCV.SetFillColor(0)

    fluxErrors.SetMinimum(.5)
    fluxErrors.SetMaximum(1.5)
    fluxErrors.SetLineWidth(0)
    fluxErrors.GetYaxis().SetTitle("#splitline{Ratio to}{Null Flux CV}")

    margRatio.SetLineStyle(2)
    oscRatio.SetFillColor(0)

    RatioAxis(fluxErrors,MNVPLOTTER)
    fluxErrors.Draw("E2")
    oscRatio.Draw("same hist l")
    margRatio.Draw("same l")
    fluxCV.Draw("same hist")

    x1, y1, x2, y2 = XtoPad(0.05),YtoPad(0.65),XtoPad(.4),YtoPad(.95)

    leg0 = ROOT.TLegend(x1,y1,x2,y2)
    leg0.AddEntry(fluxErrors,"Null Hypothesis Flux CV","f")
    leg0.AddEntry(oscRatio,"Osc. Model","l")
    leg0.AddEntry(margRatio,"Null Marg'd","l")
    leg0.Draw()

    # ------------------------------------ Second Pad -----------------------------------------
    pad1 = c1.FindObject("pad_{}_{}".format(0,1))
    pad1.cd()

    nullErrors = h_null.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
        nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
        nullErrors.SetBinContent(whichBin, 1)
    nullErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)

    straightLine = nullErrors.Clone()
    straightLine.SetLineColor(ROOT.kRed)
    straightLine.SetLineWidth(2)
    straightLine.SetFillColor(0)

    nullErrors.SetMinimum(0)
    nullErrors.SetMaximum(2)
    nullErrors.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")

    margRatio.SetLineStyle(2)
    oscRatio.SetFillColor(0)

    RatioAxis(nullErrors,MNVPLOTTER)
    nullErrors.Draw("E2")
    nullRatio.Draw("same")
    oscRatio.Draw("same hist l")
    margRatio.Draw("same l")
    straightLine.Draw("same hist")

    x1, y1, x2, y2 = XtoPad(0.05),YtoPad(0.65),XtoPad(.4),YtoPad(.95)

    leg1 = ROOT.TLegend(x1,y1,x2,y2)
    leg1.AddEntry(h_data,"Data","p")
    leg1.AddEntry(oscRatio,"Osc. Model #chi^{2}="+"{:.2f}".format(chi2_model),"l")
    leg1.AddEntry(straightLine,"Null #chi^{2}="+"{:.2f}".format(chi2_null),"l")
    leg1.AddEntry(margRatio,"Null Marg'd #chi^{2}="+"{:.2f}+{:.2f}".format(chi2_null_marg-null_pen,null_pen),"l")
    leg1.Draw()

    # ------------------------------------ Third Pad -----------------------------------------
    pad2 = c1.FindObject("pad_{}_{}".format(0,0))
    pad2.cd()

    nullMargRatio =  h_data.Clone()
    oscMargRatio =  h_osc_marg.Clone()

    nullMargRatio.Divide(nullMargRatio,h_null_marg)
    oscMargRatio.Divide(oscMargRatio, h_null_marg)

    h_null_marg.PopVertErrorBand("Flux")
    h_null_marg.AddMissingErrorBandsAndFillWithCV(h_data)
    nullMargErrors = h_null_marg.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullMargErrors.GetXaxis().GetNbins()+1): 
        nullMargErrors.SetBinError(whichBin, max(nullMargErrors.GetBinContent(whichBin), 1e-9))
        nullMargErrors.SetBinContent(whichBin, 1)
    nullMargErrors.SetFillColorAlpha(ROOT.kPink + 1, 0.4)

    RatioAxis(nullMargErrors,MNVPLOTTER)
    nullMargErrors.SetMinimum(0)
    nullMargErrors.SetMaximum(2)
    nullMargErrors.GetYaxis().SetTitle("#splitline{Ratio to Marg'd}{Null Hypothesis}")
    nullMargRatio.SetFillColor(0)
    nullMargRatio.SetLineWidth(2)
    oscMargRatio.SetFillColor(0)

    nullMargErrors.Draw("E2")
    nullMargRatio.Draw("same")
    oscMargRatio.Draw("same hist l")
    straightLine.Draw("same hist")

    x1, y1, x2, y2 = XtoPad(0.05),YtoPad(0.65),XtoPad(.4),YtoPad(.95)

    leg2 = ROOT.TLegend(x1,y1,x2,y2)
    leg2.AddEntry(h_data,"Data","p")
    leg2.AddEntry(oscMargRatio,"Osc. Model #chi^{2}="+"{:.2f}+{:.2f}".format(chi2_model_marg-model_pen,model_pen),"l")
    leg2.AddEntry(straightLine,"Null #chi^{2}="+"{:.2f}+{:.2f}".format(chi2_null_marg-null_pen,null_pen),"l")
    leg2.Draw()

    c1.Print("plots/{}_stitched_ratios_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

def PlotOscillationGrid(sample_histogram,parameters,name="",exclude="",lam=1):
    plots = histogram.keys
    nullSolution,nullPen = FluxSolution(histogram,invCov=invCov,usePseudo=usePseudo,exclude=exclude,lam=lam)
    plots.append("fhc_ratio")
    plots.append("rhc_ratio")
    histogram.titles["fhc_ratio"] = "FHC CC #nu_{#mu}/#nu_{e} Ratio"
    histogram.titles["rhc_ratio"] = "RHC CC anti #nu_{#mu}/#nu_{e} Ratio"

    mc_hists = []
    osc_hists = []
    data_hists = []
    titles = []

    for i,plot in enumerate(plots):
        title = histogram.titles[plot]
        titles.append(title)

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
                    norm_fraction.SetFillStyle(3244)
                    norm_fraction.SetLineColor(ROOT.kRed)
                    norm_fraction.SetFillColor(ROOT.kRed)

                    swap_fraction.SetTitle("#nu_{#mu}/"+nue.GetTitle())
                    swap_fraction.SetFillStyle(3244)
                    swap_fraction.SetLineColor(ROOT.kBlue)
                    swap_fraction.SetFillColor(ROOT.kBlue)

                    hists.append(norm_fraction)
                    hists.append(swap_fraction)

            h_data = histogram.data_hists[plot[:3]+'_numu_selection'].Clone()
            h_data.Divide(h_data,histogram.data_hists[plot[:3]+'_nue_selection'])
            subSample = histogram.mc_hists[plot[:3]+'_numu_selection'].Clone()
            subSample.Divide(subSample,histogram.mc_hists[plot[:3]+'_nue_selection'])
>>>>>>> feature/sterile_neutrino
        else:
            hists,total_hist = OscillateSubHistogram(histogram,plot,parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])
            h_data = histogram.data_hists[plot].Clone()
            subSample = histogram.mc_hists[plot].Clone()

        cv = np.array(subSample)[1:-1]
        mc = np.array(total_hist)[1:-1]

        weights = ReweightCV(total_hist,fluxSolution=nullSolution,cv=cv,mc=mc)
        total_hist.PopVertErrorBand("Flux")
        total_hist.AddMissingErrorBandsAndFillWithCV(h_data)

        TArray = ROOT.TObjArray()
        for hist in hists:
            hist.DivideSingle(hist,weights)
            if 'elastic' in plot:
                hist.Scale(2,'width')
            elif 'ratio' not in plot:
                hist.Scale(1,'width')
            TArray.Add(hist)

        if "elastic" in plot:
            h_data.Scale(2,'width')
            total_hist.Scale(2,'width')
        elif 'ratio' not in plot:
            h_data.Scale(1,'width')
            total_hist.Scale(1,'width')

        mc_hists.append(TArray)
        osc_hists.append(total_hist)
        data_hists.append(h_data)

    overall = plot_osc_side_by_side(mc_hists, osc_hists, data_hists, titles,MNVPLOTTER)
    overall.Print("plots/grid_oscillated_{:.1f}_{:.3f}_{:.4f}.png".format(parameters["m"],parameters['ue4'],parameters['umu4']))

def PlotOscillationEffects(sample_histogram,parameters,name="",plotSamples=False,usePseudo=False,exclude="",lam=1):
    histogram = copy.deepcopy(sample_histogram)
    invCov=histogram.GetInverseCovarianceMatrix(sansFlux=True)

    h_null = histogram.GetMCHistogram()
    h_osc = histogram.GetOscillatedHistogram()
    h_data = histogram.GetPseudoHistogram() if usePseudo else histogram.GetDataHistogram()

    chi2_null,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,usePseudo=usePseudo,exclude=exclude,lam=lam)
    chi2_model,_ = Chi2DataMC(histogram,invCov=invCov,marginalize=False,useOsc=True,usePseudo=usePseudo,exclude=exclude,lam=lam)

    chi2_model_marg,model_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,useOsc=True,usePseudo=usePseudo,setHists=True,exclude=exclude)
    chi2_null_marg,null_pen = Chi2DataMC(histogram,invCov=invCov,marginalize=True,usePseudo=usePseudo,setHists=True,lam=lam,exclude=exclude)

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
    
    h_osc_marg.SetLineStyle(2)
    h_null_marg.SetLineStyle(2)

    h_null.SetTitle(name)
    
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
    nullErrors.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")
    RatioAxis(nullErrors,MNVPLOTTER)
    nullErrors.GetXaxis().SetTitle("Bin Number")
    nullErrors.SetMinimum(0)
    nullErrors.SetMaximum(2)
    nullErrors.Draw("E2")

    #Draw the data ratios
    nullRatio.Draw("same")
    nullMargRatio.Draw("same hist l")
    oscRatio.Draw('same hist l')
    oscMargRatio.Draw('same hist l')

    #Draw a flat line at 1 for oscRatio of MC to itself
    straightLine = nullErrors.Clone()
    straightLine.SetLineColor(ROOT.kRed)
    straightLine.SetLineWidth(2)
    straightLine.SetFillColor(0)
    straightLine.Draw("HIST SAME")

    top.cd()

    overall.Print("plots/{}_stitched_{:.1f}_{:.3f}_{:.4f}.png".format(name,parameters["m"],parameters['ue4'],parameters['umu4']))

    if not plotSamples:
        return

    plots = histogram.keys
    nullSolution,nullPen = FluxSolution(histogram,invCov=invCov,usePseudo=usePseudo,exclude=exclude,lam=lam)
    plots.append("fhc_ratio")
    plots.append("rhc_ratio")
    histogram.titles["fhc_ratio"] = "FHC CC #nu_{#mu}/#nu_{e} Ratio"
    histogram.titles["rhc_ratio"] = "RHC CC anti #nu_{#mu}/#nu_{e} Ratio"

    for i,plot in enumerate(plots):
        title = histogram.titles[plot]
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
                    norm_fraction.SetFillStyle(3244)
                    norm_fraction.SetLineColor(ROOT.kRed)
                    norm_fraction.SetFillColor(ROOT.kRed)

                    swap_fraction.SetTitle("#nu_{#mu}/"+nue.GetTitle())
                    swap_fraction.SetFillStyle(3244)
                    swap_fraction.SetLineColor(ROOT.kBlue)
                    swap_fraction.SetFillColor(ROOT.kBlue)

                    hists.append(norm_fraction)
                    hists.append(swap_fraction)

            h_data = histogram.data_hists[plot[:3]+'_numu_selection'].Clone()
            h_data.Divide(h_data,histogram.data_hists[plot[:3]+'_nue_selection'])
            subSample = histogram.mc_hists[plot[:3]+'_numu_selection'].Clone()
            subSample.Divide(subSample,histogram.mc_hists[plot[:3]+'_nue_selection'])
        else:
            hists,total_hist = OscillateSubHistogram(histogram,plot,parameters["m"],parameters['ue4'],parameters['umu4'],parameters['utau4'])
            h_data = histogram.data_hists[plot].Clone()
            subSample = histogram.mc_hists[plot].Clone()

        cv = np.array(subSample)[1:-1]
        mc = np.array(total_hist)[1:-1]

        weights = ReweightCV(total_hist,fluxSolution=nullSolution,cv=cv,mc=mc)
        total_hist.PopVertErrorBand("Flux")
        total_hist.AddMissingErrorBandsAndFillWithCV(h_data)

        TArray = ROOT.TObjArray()
        for hist in hists:
            hist.DivideSingle(hist,weights)
            if 'elastic' in plot:
                hist.Scale(2,'width')
            elif 'ratio' not in plot:
                hist.Scale(1,'width')
            TArray.Add(hist)

        if "elastic" in plot:
            h_data.Scale(2,'width')
            total_hist.Scale(2,'width')
        elif 'ratio' not in plot:
            h_data.Scale(1,'width')
            total_hist.Scale(1,'width')

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

        #Error envelope for the MC
        mcRatio.SetLineColor(ROOT.kRed)
        mcRatio.SetLineWidth(2)
        mcRatio.SetMarkerStyle(0)
        mcRatio.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
        mcRatio.SetFillStyle(1001)
        mcRatio.GetYaxis().SetTitle("Data/Model")
        mcRatio.GetXaxis().SetTitle(ratio.GetXaxis().GetTitle())
        mcRatio.SetMinimum(0)
        mcRatio.SetMaximum(2)
        RatioAxis(mcRatio,MNVPLOTTER)
        mcRatio.Draw("E2")

        ratio.Draw('E1 X0 SAME')

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

    ratio.GetYaxis().SetTitle("#splitline{Ratio to Null}{Hypothesis}")
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
    straightLine.Draw("HIST SAME")
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

<<<<<<< HEAD
def PrepareSignalDecomposeRatio(data_hists,mc_hists,Grouping,only_cated = False):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew

    if not mc_hists.valid:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")

    mc_list,color,title = mc_hists.GetCateList(Grouping)
    if data_hists.valid :
        plotfunction =  lambda mnvplotter, data_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend="TR")(data_hist,mc_ints)
        hists = [data_hists.GetHist()]
    else:
        plotfunction =  lambda mnvplotter, mc_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend = "TR")(mc_hist,mc_ints)
        tmp = mc_hists.GetHist().Clone()
        tmp.Reset()
        hists =[tmp]
    hists.extend(mc_list)
    if only_cated :
        totalHist = ReSumHists(mc_list)
    else:
        totalHist = mc_hists.GetHist()
    for i in hists:
        i.Divide(i,totalHist)
    return plotfunction,hists

def PrepareComp(data_hists,mc_hists,sys_on_mc = True, sys_on_data = False,as_frac=False):
    if not (data_hists.valid or mc_hists.valid):
        raise KeyError("both data and mc is None")
    if (data_hists.valid and mc_hists.valid):
        plotfunction = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCWithErrorBand(data_hist,mc_hist,1.0,"TR")
        hists = [data_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_data else data_hists.GetHist().GetCVHistoWithStatError(),
                 mc_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_mc else mc_hist.GetHist().GetCVHistoWithStatError()]
    elif data_hists.valid:
        plotfunction = lambda mnvplotter,data_hist: data_hist.Draw("E1X0")
        hists = [data_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_data else data_hists.GetHist().GetCVHistoWithStatError()]
    else:
        plotfunction = lambda mnvplotter,mc_hist: mnvplotter.DrawMCWithErrorBand(mc_hist)
        hists = [mc_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_mc else mc_hist.GetHist().GetCVHistoWithStatError()]

    return plotfunction,hists

def PrepareRatio(data_hists,mc_hists):
    if not (data_hists.valid and mc_hists.valid):
        raise KeyError("both data and mc is Required for ratio")
    plotfunction = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCRatio(data_hist, mc_hist, 1.0 ,True,True,0,2)
    h_data = data_hists.GetHist()
    h_mc = mc_hists.GetHist()
    hists = [h_data,
             h_mc]
    return plotfunction,hists

def PrepareBkgRatio(data_hists,mc_hists):
    if not (mc_hists.valid):
        raise KeyError("mc is Required for BKG ratio")

    out_bkg = mc_hists.hists["Total"].Clone("bkgTotal")
    out_bkg.Reset()
    out_sig = mc_hists.hists["Total"].Clone("sigTotal")
    out_sig.Reset()

    for group in mc_hists.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINATION and mc_hists.hists[group]:
            out_bkg.Add(mc_hists.hists[group])
        elif group in SIGNAL_DEFINATION and mc_hists.hists[group]:
            out_sig.Add(mc_hists.hists[group])
            #print("SIGNAL ",group, out_sig.Integral())

    plotfunction = lambda mnvplotter, out_sig, out_bkg: MakeBkgRatioPlot(out_sig, out_bkg, mnvplotter)
    hists = [out_sig,
             out_bkg]

    return plotfunction,hists

def PrepareErr(data_hists,mc_hists,bkgs=None,sys_on_mc=True,sys_on_data=False,grouping=None):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew
    MNVPLOTTER.axis_maximum = 0.4
    MNVPLOTTER.legend_n_columns = 1
    plotfunction = lambda mnvplotter,data_hist: mnvplotter.DrawErrorSummary(data_hist,"TR",True,True,0.07)
    if data_hists.valid and sys_on_data :
        hist = data_hists.GetHist()
        hist.SetTitle("Background Subtracted Data")
        hist.GetXaxis().SetTitle("Energy_{estimator}")
        hists =[hist]
    elif mc_hists.valid and sys_on_mc:
        if bkgs:
            cate_hists,colors,titles  = mc_hists.GetCateList(bkgs)
            hist = ReSumHists(cate_hists)
            hist.SetTitle("Predicted Signal")
            hist.GetXaxis().SetTitle("Energy_{estimator}")
        else:
            hist = mc_hists.GetHist()

        #hist.PopVertErrorBand("SuSA_Valencia_Weight")
        hist.PopVertErrorBand("LowQ2Pi_None")
        hist.PopVertErrorBand("LowQ2Pi")
        #hist.PopVertErrorBand("MK_model")
        #hist.PopVertErrorBand("fsi_weight")
        hists =[hist]
        #hists =[hist]
    else:
        raise KeyError("Can't make error plot for systematics config mismatch")
    #updatePlotterErrorGroup(grouping)
    #AdaptivePlotterErrorGroup(hists[0],7)
    return plotfunction,hists

def PrepareErrorBand(data_hists,mc_hists, name, sys_on_mc=True,sys_on_data=False):
    plotfunction = lambda mnvplotter,hist: MakeErrorBandPlot(hist,name,mnvplotter)
    if data_hists.valid and sys_on_data:
        hists =[data_hists.GetHist()]
    elif mc_hists.valid and sys_on_mc:
        hists =[mc_hists.GetHist()]
    else:
        raise KeyError("Can't make error plot for systematics config mismatch")
    return plotfunction,hists

def Prepare2DStack(data_hists,mc_hists,Grouping = None):
    if not mc_hists.valid:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")
    mc_list,color,title = mc_hists.GetCateList(Grouping)

    if data_hists.valid:
        plotfunction =  lambda mnvplotter, data_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title,legend="TR")(data_hist,mc_ints)
        hists = [data_hists.GetHist()]
    else:
        plotfunction =  lambda mnvplotter, mc_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend = "TR")(mc_hist,mc_ints)
        tmp = mc_hists.GetHist().Clone()
        tmp.Reset()
        hists =[tmp]
    hists.extend(mc_list)
    return plotfunction,hists

def PrepareStack(data_hists,mc_hists,Grouping = None):
    if not mc_hists.valid:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")
    mc_list,color,title = mc_hists.GetCateList(Grouping)
    if data_hists.valid:
        plotfunction =  lambda mnvplotter, data_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend="TR")(data_hist,mc_ints)
        if "frontdedx" in data_hists.GetHist().GetName():
            hist = data_hists.GetHist()
            #hist.GetXaxis().SetTitle("dE/dx (MeV/cm)")
            hist.GetYaxis().SetTitle("dNEvents/d(dE/dx)")
        hists = [data_hists.GetHist()]
    else:
        plotfunction =  lambda mnvplotter, mc_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend = "TR")(mc_hist,mc_ints)
        tmp = mc_hists.GetHist().Clone()
        tmp.Reset()
        hists =[tmp]
    for hist in mc_list:
        if "frontdedx" in hist.GetName():
            #hist.GetXaxis().SetTitle("dE/dx (MeV/cm)")
            hist.GetYaxis().SetTitle("dNEvents/d(dE/dx)")
        if hist == None:
            del hist
    hists.extend(mc_list)
    return plotfunction,hists

def PrepareStackNew(data_hists,mc_hists,Grouping = None):
    if not mc_hists.valid:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")
    mc_list,color,title = mc_hists.GetCateList(Grouping)
    if data_hists.valid :
        plotfunction =  lambda mnvplotter, data_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend="TR")(data_hist,mc_ints)
        hists = [data_hists.GetHist()]
    else:
        plotfunction =  lambda mnvplotter, mc_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend = "TR")(mc_hist,mc_ints)
        tmp = mc_hists.GetHist().Clone()
        tmp.Reset()
        hists =[tmp]
    hists.extend(mc_list)
    return plotfunction,hists

def RebinPrepareStack(data_hists,mc_hists,Grouping = None):
    if not mc_hists.valid:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")
    mc_list,color,title = mc_hists.GetCateList(Grouping)
    if data_hists.valid :
        plotfunction =  lambda mnvplotter, data_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend="TR")(data_hist,mc_ints)
        data_hists.GetHist().Rebin(5)
        hists = [data_hists.GetHist()]
    else:
        plotfunction =  lambda mnvplotter, mc_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend = "TR")(mc_hist,mc_ints)
        tmp = mc_hists.GetHist().Clone()
        tmp.Reset()
        hists =[tmp]
    for hist in mc_list:
        hist.Rebin(5)
    hists.extend(mc_list)
    return plotfunction,hists

def CategoryHist(data_hists,mc_hists,category):
    if not(mc_hists.valid):
        raise KeyError("No MC histogram to plot migration")
    if mc_hists.valid and category:
        hist = mc_hists.GetHist().Clone()
        hist.Reset()
        for cate in category["cate"]:
            if cate in mc_hists.hists:
                tmp = mc_hists.hists[cate]
                hist.Add(tmp)
        #hist.SetMaximum(100)
        hists = [hist]
        #plotfunction = lambda mnvplotter,mc_hist: mnvplotter.DrawNormalizedMigrationHistogram(mc_hist,False,False,True,True)
        plotfunction = lambda mnvplotter,mc_hist: MakeCategoryHist(mc_hist)
        return plotfunction,hists
    else:
        hist = mc_hists.GetHist()
        #hist.SetMaximum(100)
        hists = [hist]
        #plotfunction = lambda mnvplotter,mc_hist: mnvplotter.DrawNormalizedMigrationHistogram(mc_hist,False,False,True,True)
        plotfunction = lambda mnvplotter,mc_hist: MakeCategoryHist(mc_hist)
        #return plotfunction,hists

def MakeCategoryHist(mc_hist):
    #tcanvas = ROOT.TCanvas()
    #SetMargin(tcanvas)
    mc_hist.DrawCopy("colz")

def PrepareDiff(data_hists,mc_hists):
    if not (data_hists.valid and mc_hists.valid):
        raise KeyError("both data and mc is Required for differece")
    hists = [data_hists.GetHist(),mc_hists.GetHist()]
    def plotDifference(mnvplotter,data_hist,mc_hist):
        #print data_hist
        #ndf = ROOT.Long()
        ndf = ctypes.c_double()
        chi2mnv=mnvplotter.Chi2DataMC(data_hist,mc_hist,ndf)
        #print(chi2mnv,ndf,data_hist.Integral("width"),mc_hist.Integral("width"))
        chi2,ndf=CalChi2(data_hist,mc_hist)
        tmp = data_hist.Clone()
        tmp.Add(mc_hist,-1)
        #print(data_hist.GetName(),tmp.Integral(0,-1,"width"))
        tmp.Draw("E1")
        size = 0.035
        #align = ROOT.Long()
        #xLabel = ROOT.Double()
        #yLabel = ROOT.Double()
        align = ctypes.c_int(1)
        xLabel = ctypes.c_double(1.)
        yLabel = ctypes.c_double(1.)
        mnvplotter.DecodePosition("TR", size, align, xLabel, yLabel )
        mnvplotter.AddPlotLabel("chi2/ndf: {:.4f}/{:d}".format(chi2,ndf),xLabel, yLabel, size, 4, 112, 22)#align)
        print(("chi2/ndf: {:.4f}/{:d}".format(chi2,ndf)))

    return plotDifference,hists

def PrepareMigration(data_hists,mc_hists,Grouping):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew
    if not(mc_hists.valid):
        raise KeyError("No MC histogram to plot migration")
    mc_list,color,title = mc_hists.GetCateList(Grouping)
    totalHist = ReSumHists(mc_list)
    totalHist.SetTitle("Combined Signal")
    totalHist.SetMaximum(100)
    hists = [totalHist]
    plotfunction = lambda mnvplotter,mc_hist: mnvplotter.DrawNormalizedMigrationHistogram(mc_hist,False,False,True,True)
    return plotfunction,hists

def PrepareHist2D(data_hists,mc_hists,Grouping):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew
    if not(mc_hists.valid):
        raise KeyError("No MC histogram to plot migration")
    mc_list,color,title = mc_hists.GetCateList(Grouping)
    totalHist = ReSumHists(mc_list)
    totalHist.SetTitle("Combined Signal")
    hists = [totalHist]
    plotfunction = lambda mnvplotter,mc_hist: mc_hist.DrawCopy("colz")
    return plotfunction,hists

def updatePlotterErrorGroup(group,mnvplotter=MNVPLOTTER):
    mnvplotter.error_summary_group_map.clear();
    for k,v in group.items():
        vec = ROOT.vector("std::string")()
        for vs in v :
            vec.push_back(vs)
        mnvplotter.error_summary_group_map[k]= vec
        
def TopNErrorBand(hist,topN):
    #find N largest errorband given histogram hist
    heap = []
    for errorband_name in hist.GetVertErrorBandNames():
        sum_error = hist.GetVertErrorBand(errorband_name).GetErrorBand(False,False).Integral()
        #print(errorband_name,hist.GetVertErrorBand(errorband_name).GetErrorBand(False,False).Integral())
        if len(heap)<topN:
            heapq.heappush(heap, (sum_error,errorband_name))
        elif sum_error>heap[0][0]:
            heapq.heappushpop(heap,(sum_error,errorband_name))
    #print(heap)
    result = []
    while heap:
        result.append(heapq.heappop(heap)[1])
    return result

#def AdaptivePlotterErrorGroup(hist,topN,mnvplotter=MNVPLOTTER):
#    result = TopNErrorBand(hist,topN)
#    rest = [i for i in hist.GetVertErrorBandNames() if i not in result]
#    updatePlotterErrorGroup({"Rest":rest},mnvplotter)

def AdaptivePlotterErrorGroup(hist,result,mnvplotter=MNVPLOTTER):
    rest = [i for i in hist.GetVertErrorBandNames() if i not in result]
    updatePlotterErrorGroup({"Rest":rest},mnvplotter)

def CalMXN(N_plots,max_horizontal_plot = 4):
    height = math.ceil(1.0*N_plots/max_horizontal_plot)  #hard code max 3 plots per line
    width = math.ceil(N_plots/height) #number of plots per line
    return int(width),int(height)

def Make2DSlice(hist2D_o, X_slice=True, bin_start = 1 , bin_end = 0,interval = 1):
    # x slice true produce 1d histograms of Y variable for each x variable bins.
    # bin_end = 0 : all except overflow, = -1: all including overflow, other: bins before(i.e. *not* including ) bin_end

    slicing_hists = []
    hist2D = hist2D_o.Clone()
    axis = hist2D.GetXaxis() if X_slice else hist2D.GetYaxis()
    Nbins = axis.GetNbins()
    start = max(0,bin_start)
    end = Nbins - bin_end + 1 if bin_end <=0 else bin_end
    for i in range(start,end,interval):
        slicing_hists.append(hist2D.ProjectionX(hist2D.GetName()+str(i),i,i+interval-1,"o") if not X_slice else hist2D.ProjectionY(hist2D.GetName()+str(i),i,i+interval-1,"o"))
        slicing_hists[-1].SetTitle("%.2f<%s<%.2f"%(axis.GetBinLowEdge(i),axis.GetTitle(),axis.GetBinUpEdge(i+interval-1)))
        slicing_hists[-1].GetYaxis().SetTitle(hist2D.GetZaxis().GetTitle())
        if "frontdedx" in hist2D.GetName():
            slicing_hists[-1].Rebin(2)

    del hist2D
    return slicing_hists

def MakeDataMCPlot(data_hist, mc_hist, pot_scale=1, sys_on_mc = True, sys_on_data = False, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    local_mc = (mc_hist.GetCVHistoWithError() if sys_on_mc else mc_hist.GetCVHistoWithStatError()) if mc_hist else None
    local_data = (data_hist.GetCVHistoWithError() if sys_on_data else data_hist.GetCVHistoWithStatError()) if data_hist else None
    if not (local_mc or local_data):
        raise KeyError("both data and mc is None")
    if local_mc and local_data:
        mnvplotter.DrawDataMCWithErrorBand(local_data,local_mc,pot_scale,"TR")
    elif local_mc:
        mnvplotter.DrawMCWithErrorBand(local_mc,pot_scale)
    else:
        local_data.Draw("E1X0")

def MakeModelVariantPlot(data_hist, mc_hists, color=None, title=None,legend ="TR",pot_scale=1.0,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    if not mc_hists:
        raise KeyError("Doesn't make sense to plot model variation without MC")
    TArray = ROOT.TObjArray()
    for i in range(len(mc_hists)):
        if color:
            mc_hists[i].SetLineColor(color[i])
        if title:
            mc_hists[i].SetTitle(title[i])
        TArray.Add(mc_hists[i])
    #mnvplotter.DrawDataMCVariations(data_hist,TArray,pot_scale,legend,True,True,False,False,False)
    mnvplotter.DrawDataMCVariations(data_hist,TArray,pot_scale,legend,True,True,False,False)

def MakeDataMCStackedPlot(data_hist, mc_hists, legend = "TR", pot_scale=1, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    if not mc_hists:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")
    TArray = ROOT.TObjArray()
    for i in range(len(mc_hists)):
        if color:
            mc_hists[i].SetFillColor(color[i])
        if title:
            mc_hists[i].SetTitle(title[i])

        if mc_hists[i]:
            TArray.Add(mc_hists[i])
    if data_hist is not None:
        mnvplotter.DrawDataStackedMC(data_hist,TArray,pot_scale,legend,"Data",0,0,1001)
    #else:
    elif TArray is not None:
        mnvplotter.DrawStackedMC(TArray,pot_scale,legend,0,0,1001)

def MakeSignalDecomposePlot(data_hist, mc_hist, mc_hists, title, color, pot_scale = 1.0, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    #if data_hist is not None:
    #    mnvplotter.DrawDataMCWithErrorBand(data_hist,mc_hist,pot_scale,"TR")
    #else:
    #    mnvplotter.DrawMCWithErrorBand(mc_hist,pot_scale)
    mnvplotter.axis_minimum = 0
    tcanvas = canvas.GetPad(1)
    #tcanvas.SetLogy()
    mc_hists[0].SetLineWidth(4)
    mc_hists[0].SetLineColor(color[0])
    #mc_hists[0].Scale(pot_scale)
    mc_hists[0].SetMaximum(1)
    mc_hists[0].DrawCopy("HIST")
    for i in range(1,len(mc_hists)):
        mc_hists[i].SetLineWidth(4)
        mc_hists[i].SetLineColor(color[i])
        #mc_hists[i].Scale(pot_scale)
        mc_hists[i].DrawCopy("HIST SAME")

def MakeRatioPlot(data_hist,mc_hist,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    if not (data_hist and mc_hist):
        raise KeyError("both data and mc is Required for ratio")
    mnvplotter.DrawDataMCRatio(data_hist, mc_hist, 1.0 ,True,True,0,2)

def MakeErrPlot(hist,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    mnvplotter.DrawErrorSummary(hist)

def MakeErrorBandPlot(hist,name,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    #hist.PopVertErrorBand("SuSA_Valencia_Weight")
    #hist.PopVertErrorBand("LowQ2Pi_None")
    errorband = hist.GetVertErrorBand(name)
    errorband.DrawAll("",True)

def MakeBkgRatioPlot(out_sig, out_bkg, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    signal = out_sig.GetCVHistoWithError().Clone()
    background = out_bkg.GetCVHistoWithError().Clone()
    tune = -1
    best_cut = 0
    s = 0
    b = 0
    S = 0
    B = 0
    x = array('d',[])
    y1 = array('d',[])
    width = array('d',[])
    for i in range(signal.GetNbinsX()):
        s+=signal.GetBinContent(i)
        b+=background.GetBinContent(i)
        ratio = (s/math.sqrt(s+b)) if (s + b > 0) else 0
        x.append(signal.GetBinLowEdge(i)+signal.GetBinWidth(i))
        y1.append(ratio)
        width.append(signal.GetBinWidth(i))
        if ratio > tune:
            tune = ratio
            best_cut = i
            S=s
            B=b

    x[0]=0
    x.append(x[-1]+signal.GetBinWidth(signal.GetNbinsX()))
    h1 = ROOT.TH1F("h1", "Signal Significance", signal.GetNbinsX(), x)
    for i in range(signal.GetNbinsX()):
        h1.Fill(x[i]+width[i]/2,y1[i])
    h1.GetXaxis().SetTitle("E_{available}")
    h1.GetYaxis().SetTitle("s/sqrt(s+b)")

    #mnvplotter.mc_error_color = 0
    #total = background/(signal + background)
    total = out_bkg.Clone()
    #total.PopVertErrorBand("SuSA_Valencia_Weight")
    #total.PopVertErrorBand("LowQ2Pi_None")
    #total.PopVertErrorBand("MK_model")
    total.Scale(1/total.Integral())
    total.GetXaxis().SetTitle("E_{available} + E_{lepton}")
    #total.GetYaxis().SetTitle("Background Fraction per Bin")
    mnvplotter.DrawMCWithErrorBand(total)
    #mnvplotter.DrawMCWithErrorBand(h1)
    #mnvplotter.DrawDataMCRatio(out_sig, out_bkg, 1.0 ,True,True,0,4,"Signal/Background")
    #mnvplotter.AddPlotLabel("Max S/sqrt(S+B) at E_avail = {:.2f}".format(x[best_cut]+width[best_cut]),0.65,0.8,0.15)
    #mnvplotter.AddPlotLabel("{:.3f}/sqrt({:.3f}+{:.3f}) = {:.4f}".format(S,S,B,tune),0.65,0.625,0.15)
    #mnvplotter.AddPlotLabel("Max S/sqrt(S+B) at E_lepton = {:.3f}".format(signal.GetBinLowEdge(best_cut)+signal.GetBinWidth(best_cut)),0.65,0.8,0.08)
    #mnvplotter.AddPlotLabel("{:.3f}/sqrt({:.3f}+{:.3f} = {:.4f}".format(S,S,B,tune),0.65,0.625,0.08)

#def MakeGridPlot(MakingSlice,MakingEachPlot,input_hists,CanvasConfig=lambda canvas:True, draw_seperate_legend = False, mnvplotter=MNVPLOTTER,canvas=CANVAS):
#    if canvas is CANVAS:
#        canvas.Clear()
#    slices = list(map(lambda *args: args, *list(map(MakingSlice,input_hists))))
#    N_plots = len(slices)
#    canvas.Divide(*CalMXN(N_plots+int(draw_seperate_legend)))
#    for i in range(N_plots):
#        canvas.cd(i+1)
#        if not CanvasConfig(canvas.GetPad(i+1)):
#            print("Warning: failed setting canvas.")
#        if N_plots>1:
#            SetMargin(canvas.GetPad(i+1))
#        MakingEachPlot(mnvplotter,*slices[i])
#        mnvplotter.AddHistoTitle(slices[i][0].GetTitle())
#    if draw_seperate_legend:
#        for i in range(N_plots):
#            Tleg = GetTLegend(canvas.GetPad(i+1))
#        if Tleg:
#            canvas.cd(N_plots+1)
#            Tleg.SetX1(0)
#            Tleg.SetX2(1)
#            Tleg.SetY1(0)
#            Tleg.SetY2(1)
#            Tleg.SetTextSize(2*MNVPLOTTER.legend_text_size);
#            Tleg.Draw()

def SumGridPlots(MakingSlice,MakingEachPlot,input_hists,CanvasConfig=lambda canvas:True, draw_seperate_legend = False, mnvplotter=MNVPLOTTER,canvas=CANVAS,outname=None,title=None):
    slices = list(map(lambda *args: args, *list(map(MakingSlice,input_hists)))) # MakingSlice is Make2DSlice, takes argument to slice along x or y axis 
    new_slices = slices[0]
    for i in range(1,len(slices)):
        for j in range(len(slices[0])):
            new_slices[j].Add(slices[i][j])
    new_slices = [new_slices]
    del slices
    N_plots = 1
    canvas=CANVAS
    canvas.Divide(*CalMXN(N_plots+int(draw_seperate_legend)))
    for i in range(N_plots):
        canvas.cd(i+1)
        tcanvas = canvas.GetPad(i+1)
        if not CanvasConfig(tcanvas):
            print("Warning: failed setting canvas.")
        if N_plots <= 1:
            SetMargin(tcanvas)
            #tcanvas.SetLogz()
        else:
            SetMarginMulti(tcanvas)
            #i+=34
        MakingEachPlot(mnvplotter,*new_slices[i])
        #if title:
        #    mnvplotter.AddHistoTitle(title)
        #else:
        #    mnvplotter.AddHistoTitle(new_slices[i][0].GetTitle())

        #code.interact(local=locals())
    if draw_seperate_legend:
        Tleg = None
        for i in range(N_plots):
            Tleg = GetTLegend(canvas.GetPad(i+1)) or Tleg
        if Tleg:
            canvas.cd(N_plots+1)
            Tleg.SetX1(0.16)
            Tleg.SetX2(0.98)
            Tleg.SetY1(0.08)
            Tleg.SetY2(0.92)
            Tleg.SetTextSize(2*MNVPLOTTER.legend_text_size);
            Tleg.Draw()
    if outname:
        mnvplotter.MultiPrint(canvas,outname,"png")


def MakeGridPlot(MakingSlice,MakingEachPlot,input_hists,CanvasConfig=lambda canvas:True, draw_seperate_legend = False, mnvplotter=MNVPLOTTER,canvas=CANVAS,outname=None,title=None):
    if canvas is CANVAS:
        canvas.Clear()
        #canvas.SetLeftMargin(0)
    slices = list(map(lambda *args: args, *list(map(MakingSlice,input_hists))))
    N_plots = len(slices)
    #if N_plots > 1:
    #    N_plots = 9
    canvas.Divide(*CalMXN(N_plots+int(draw_seperate_legend)))
    for i in range(N_plots):
        canvas.cd(i+1)
        tcanvas = canvas.GetPad(i+1)
        if not CanvasConfig(tcanvas):
            print("Warning: failed setting canvas.")
        if N_plots <= 1:
            SetMargin(tcanvas)
            #tcanvas.SetLogz()
        else:
            SetMarginMulti(tcanvas)
            #i+=34
        MakingEachPlot(mnvplotter,*slices[i])
        if title:
            mnvplotter.AddHistoTitle(title)
        else:
            mnvplotter.AddHistoTitle(slices[i][0].GetTitle())

        #code.interact(local=locals())
    if draw_seperate_legend:
        Tleg = None
        for i in range(N_plots):
            Tleg = GetTLegend(canvas.GetPad(i+1)) or Tleg
        if Tleg:
            canvas.cd(N_plots+1)
            Tleg.SetX1(0.16)
            Tleg.SetX2(0.98)
            Tleg.SetY1(0.08)
            Tleg.SetY2(0.92)
            Tleg.SetTextSize(2*MNVPLOTTER.legend_text_size);
            Tleg.Draw()
    if outname:
        mnvplotter.MultiPrint(canvas,outname,"png")


def Print(outname,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    mnvplotter.MultiPrint(canvas,outname)

def SetMargin(pad):
    pad.SetRightMargin(0.15)
    pad.SetLeftMargin(.15)
    pad.SetTopMargin(0.08)
    pad.SetBottomMargin(0.2)
    ROOT.TGaxis.SetExponentOffset(0.04,-0.1,"y")

def SetMarginMulti(pad):
    pad.SetRightMargin(0.1)
    pad.SetLeftMargin(0.1)
    pad.SetTopMargin(0.08)
    pad.SetBottomMargin(0.05)
    ROOT.TGaxis.SetExponentOffset(0.04,-0.1,"y")

def GetTLegend(pad):
    tlist = pad.GetListOfPrimitives()
    for i in tlist:
        if (isinstance(i,ROOT.TLegend)) :
            pad.RecursiveRemove(i)
            pad.Update()
            return i
    return None

def MakeDataMCStackedPlot(data_hist, mc_hists, color=None, title=None, legend = "TR", pot_scale=1, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    TArray = ROOT.TObjArray()
    for i in range(len(mc_hists)):
        if not mc_hists[i] or mc_hists[i].Integral() <= 0:
            continue
        if color is not None:
            mc_hists[i].SetFillColor(color[i])
        if title is not None:
            mc_hists[i].SetTitle(title[i])
        TArray.Add(mc_hists[i])
    if data_hist.Integral() > 0:
        try:
            mnvplotter.DrawDataStackedMC(data_hist,TArray,pot_scale,legend,"Data",0,0,1001)
        except:
            pass
    else:
        mnvplotter.DrawStackedMC(TArray,pot_scale,legend,0,0,1001)


def MakeMigrationPlots(hist, output, no_text = False, fix_width = True, mnvplotter= MNVPLOTTER, canvas = CANVAS):
    if canvas is CANVAS:
        canvas.Clear()
    #make sure migration matrix has fix bin width
    if fix_width:
        hist.GetXaxis().Set(hist.GetNbinsX(),0,hist.GetNbinsX())
        hist.GetXaxis().SetTitle("reco bin number")
        hist.GetYaxis().Set(hist.GetNbinsY(),0,hist.GetNbinsY())
        hist.GetYaxis().SetTitle("truth bin number")
    mnvplotter.DrawNormalizedMigrationHistogram(hist,False,False,True,no_text)
    mnvplotter.MultiPrint(canvas,output)

def Make2DPlot(hist,output,mnvplotter= MNVPLOTTER, canvas = CANVAS):
    if canvas is CANVAS:
        canvas.Clear()

    plotfunction = lambda mnvplotter,mc_hist: mc_hist.DrawCopy("colz")
    return plotfunction,hists

def CalChi2(data_hist,mc_hist,pot_scale=1.0):
    chi2 = 0
    ndf = 0
    for i in range (0,data_hist.GetSize()):
        data = data_hist.GetBinContent(i)
        mc = mc_hist.GetBinContent(i)*pot_scale
        sig = data_hist.GetBinError(i)
        #print (i,data,mc,chi2)
        chi2 += (data-mc)**2/sig/sig if data>0 else 0
        ndf += 1 if data>0 else 0
    print("data/mc integral: {}/{}".format(data_hist.Integral(),mc_hist.Integral()))
    print("chi2/ndf of histogram {} is {}/{}.".format(data_hist.GetName(),chi2,ndf))
    return chi2,ndf
=======
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
>>>>>>> feature/sterile_neutrino

    #Draw a flat line at 1 for ratio of MC to itself
    straightLine = mcRatio.Clone()
    straightLine.SetFillStyle(0)
    straightLine.Draw("HIST SAME")
    top.cd()
    
    MNVPLOTTER.WritePreliminary(0.4, 0.05, 5e-2, True)
    overall.Print("plots/{}".format(title))
