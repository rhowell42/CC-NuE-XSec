"""
Get Crossection Plots.
Inputs are unfolded data histogram, efficiency histogram(or seperated true signal and selected signal),
"""

import ROOT
import PlotUtils
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities, PlotTools
from tools.PlotLibrary import PLOT_SETTINGS
from config.SystematicsConfig import USE_NUE_CONSTRAINT,CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS
from config import DrawingConfig
from config.CutConfig import NEUTRINO_ENERGY_RANGE,FIDUCIAL_Z_RANGE
from functools import partial

ROOT.TH1.AddDirectory(False)

XSEC_TO_MAKE = [
    #"q0 vs q3",
    "Visible Energy vs q3",
    "Visible Energy vs Lepton Pt"
]

TARGET_UTILS = PlotUtils.TargetUtils.Get()

def GetXSectionHistogram(unfolded,efficiency,is_mc):
    #divide by efficiency
    unfolded.Divide(unfolded,efficiency)
    #divide by flux
    DivideFlux(unfolded,is_mc)
    #divide by N nucleaon
    Nnucleon = TARGET_UTILS.GetTrackerNNucleons(FIDUCIAL_Z_RANGE[0],FIDUCIAL_Z_RANGE[1],is_mc)
    unfolded.Scale(1.0/Nnucleon)
    return unfolded

def DivideFlux(unfolded,is_mc):
    frw= PlotUtils.flux_reweighter("minervame1A",12,USE_NUE_CONSTRAINT) #playlist is dummy for now
    flux = frw.GetIntegratedFluxReweighted(12,unfolded,NEUTRINO_ENERGY_RANGE[0],NEUTRINO_ENERGY_RANGE[1],False)
    flux.PopVertErrorBand("Flux_BeamFocus")
    flux.PopVertErrorBand("ppfx1_Total")
    flux.Scale(1e-4*(mc_pot if is_mc else data_pot)) #change unit to nu/cm^2
    unfolded.Divide(unfolded,flux)

def GetEfficiency(ifile, ifile_truth, plot):
    if ifile_truth:
        print ("using separate efficiency file")
        ifile = ifile_truth
    num = Utilities.GetHistogram(ifile,PLOT_SETTINGS[plot+" Migration"]["name"]+"_truth")
    den = Utilities.GetHistogram(ifile,PLOT_SETTINGS["True Signal "+plot]["name"])
    hist_out = num.Clone()
    hist_out.Divide(num,den,1.0,1.0,"B")
    #plot efficiency
    hist_out.GetYaxis().SetTitle(den.GetYaxis().GetTitle())
    hist_out.GetXaxis().SetTitle("Eavail (GeV)") #hard coding for now
    hist_out.GetZaxis().SetTitle("Efficiency")

    plotter = lambda mnvplotter, hist: mnvplotter.DrawMCWithErrorBand(hist.GetCVHistoWithError())
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[hist_out],lambda x: True,False)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"eff"))
    #PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[hist_out],AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"eff"))
    plotter = lambda mnvplotter, hist: mnvplotter.DrawErrorSummary(hist)
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[hist_out],lambda x: True,False)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"eff_err"))
    return hist_out

def GetMCXSectionHistogram(mc_file,plot):
    # mc unfolded is truth_signal:
    unfolded = Utilities.GetHistogram(mc_file,PLOT_SETTINGS["True Signal "+plot]["name"])
    DivideFlux(unfolded,True)
    Nnucleon = TARGET_UTILS.GetTrackerNNucleons(FIDUCIAL_Z_RANGE[0],FIDUCIAL_Z_RANGE[1],True)
    unfolded.Scale(1.0/Nnucleon)
    sig_dep = []
    colors = []
    titles = []
    for chan in ["CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"]:
        tmp =Utilities.GetHistogram(mc_file,"{}_{}".format(PLOT_SETTINGS["True Signal "+plot]["name"],chan))
        DivideFlux(tmp,True)
        tmp.Scale(1.0/Nnucleon)
        sig_dep.append(tmp)
        colors.append(DrawingConfig.SignalDecomposition[chan]["color"])
        titles.append(DrawingConfig.SignalDecomposition[chan]["title"])
    return unfolded,sig_dep,colors,titles

if __name__ == "__main__":
    playlist= AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_reco_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    unfolded_file = ROOT.TFile.Open(AnalysisConfig.UnfoldedHistoPath(playlist,AnalysisConfig.bkgTune_tag,False))
    #mc_reco_file = ROOT.TFile.Open(AnalysisCmc_reco_file = ROOT.TFile.Open(AnalysisConfig.SelectionHistoPath(playlist,False,False))
    mc_truth_file = ROOT.TFile.Open(AnalysisConfig.SelectionHistoPath(playlist+"-BigNuE",False))
    xsec_file = ROOT.TFile.Open(AnalysisConfig.XSecHistoPath(playlist),"RECREATE")
    for plot in XSEC_TO_MAKE:
        PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
        unfolded = Utilities.GetHistogram(unfolded_file,PLOT_SETTINGS[plot]["name"]+"_bkg_unfolding")
        efficiency = GetEfficiency(mc_reco_file,mc_truth_file,plot)
        xsec = GetXSectionHistogram(unfolded,efficiency,False)
        mc_xsec,sig_dep,colors,titles = GetMCXSectionHistogram(mc_reco_file,plot)
        ylabel = xsec.GetZaxis().GetTitle().replace("NEvents","#sigma")
        #ylabel += "(cm^2/Gev^2/c^2/Nucleon)"
        xsec.GetZaxis().SetTitle(ylabel)
        mc_xsec.GetZaxis().SetTitle(ylabel)
        mc_xsec.GetXaxis().SetTitle("Eavail (GeV)") #hard coding for now

        for h in [xsec,mc_xsec,*sig_dep]:
            h.Scale(1.0,"width")

        plotter = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCWithErrorBand(data_hist.GetCVHistoWithError(),mc_hist.GetCVHistoWithStatError(),1.0,"TR")
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec,mc_xsec],draw_seperate_legend=True)
        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec"))
        plotter = lambda mnvplotter, data_hist,mc_hist : mnvplotter.DrawDataMCRatio(data_hist.GetCVHistoWithError(),mc_hist.GetCVHistoWithStatError(),1.0,True,0,2)
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec,mc_xsec],draw_seperate_legend=True)
        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_ratio"))
        plotter = lambda mnvplotter, hist: mnvplotter.DrawErrorSummary(hist)
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec],draw_seperate_legend=True)
        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_err"))
        PlotTools.AdaptivePlotterErrorGroup(xsec,18)
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec],draw_seperate_legend=True)
        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_err_top7"))

        plotter = lambda mnvplotter,data_hist, mc_hist, *mc_ints : partial(PlotTools.MakeSignalDecomposePlot,color=colors,title=titles)(data_hist,mc_hist,mc_ints)
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec,mc_xsec,*sig_dep],draw_seperate_legend=True)
        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_sigdep"))
        xsec_file.cd()
        xsec.Write(PLOT_SETTINGS[plot]["name"]+"_dataxsec")
        mc_xsec.Write(PLOT_SETTINGS[plot]["name"]+"_mcxsec")
        print("done")
    xsec_file.Close()
