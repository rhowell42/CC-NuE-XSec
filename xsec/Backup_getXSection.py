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
from config.CutConfig import NEUTRINO_ENERGY_RANGE,FIDUCIAL_Z_RANGE
from config import DrawingConfig

from functools import partial
ROOT.TH1.AddDirectory(False)

XSEC_TO_MAKE = [
    #"q0 vs q3",
    "Visible Energy vs q3",
    "Visible Energy vs Lepton Pt"
]

threshold = 100
TARGET_UTILS = PlotUtils.TargetUtils.Get()

def GetXSectionHistogram(unfolded,efficiency,is_mc):
    #divide by efficiency
    efficiency.AddMissingErrorBandsAndFillWithCV(unfolded)
    unfolded.Divide(unfolded,efficiency)
    #divide by flux
    DivideFlux(unfolded,is_mc)
    #divide by N nucleaon
    Nnucleon = TARGET_UTILS.GetTrackerNNucleons(FIDUCIAL_Z_RANGE[0],FIDUCIAL_Z_RANGE[1],is_mc)
    h_nucleon = GetNnucleonError(unfolded,Nnucleon)
    
    return unfolded

def GetNnucleonError(hist,ntargets):
    hist_target = hist.Clone("number_of_targets")
    hist_target.ClearAllErrorBands()
    hist_target.Reset()
    errband_name = "Target_Mass_CH"
    band_err = 0.014
    hist_target.AddVertErrorBand(errband_name,2)

    for i in range(hist_target.GetSize()):
        hist_target.SetBinContent(i,ntargets)
        hist_target.SetBinError(i,0)
        hist_target.GetVertErrorBand(errband_name).SetBinContent(i,ntargets)
        hist_target.GetVertErrorBand(errband_name).GetHist(0).SetBinContent(i,ntargets*(1-band_err))
        hist_target.GetVertErrorBand(errband_name).GetHist(1).SetBinContent(i,ntargets*(1+band_err))
    hist_target.AddMissingErrorBandsAndFillWithCV(hist)
    print ("Target Normalization: {:.4E},{:.4E}".format(ntargets,ntargets*0.014))
    return hist_target

def ProjectBinContent(num, den):
    xNBins = num.GetNbinsX()
    yNBins = num.GetNbinsY() 
    for j in range(1,yNBins+1):
        for i in range(1,xNBins+1): 
            if 0 <= den.GetBinContent(i,j) < 100: #treshold 
                num, den = ProjectUniverseContent(num, den, i, j) #need to do this for every universe
                denNewContent = den.GetBinContent(i-1,j)
                numNewContent = num.GetBinContent(i-1,j) 
                den.SetBinContent(i,j,denNewContent)
                num.SetBinContent(i,j,numNewContent)
    return num, den


def ProjectUniverseContent(num, den, i, j):
    for bandname in den.GetVertErrorBandNames():
        nHists = den.GetVertErrorBand(bandname).GetNHists()
        for k in range(0,nHists):
            denNewContent = den.GetVertErrorBand(bandname).GetHist(k).GetBinContent(i-1,j)
            numNewContent = num.GetVertErrorBand(bandname).GetHist(k).GetBinContent(i-1,j)
            den.GetVertErrorBand(bandname).GetHist(k).SetBinContent(i,j,denNewContent)
            num.GetVertErrorBand(bandname).GetHist(k).SetBinContent(i,j,numNewContent) 
    return  num, den
   
def DivideFlux(unfolded,is_mc):
    print(unfolded)
    frw= PlotUtils.flux_reweighter("minervame6A",-12,USE_NUE_CONSTRAINT) #playlist is dummy for now 
    flux = frw.GetIntegratedFluxReweighted(-12,unfolded,NEUTRINO_ENERGY_RANGE[0],NEUTRINO_ENERGY_RANGE[1],False)
    
    flux.PopVertErrorBand("Flux_BeamFocus")
    flux.PopVertErrorBand("ppfx1_Total")
    
    flux.PopVertErrorBand("SuSA_Valencia_Weight")
    unfolded.PopVertErrorBand("SuSA_Valencia_Weight")

    flux.PopVertErrorBand("MK_model")
    unfolded.PopVertErrorBand("MK_model")

    flux.PopVertErrorBand("fsi_weight")
    unfolded.PopVertErrorBand("fsi_weight")

    flux.Scale(1e-4*(mc_pot if is_mc else data_pot)) #change unit to nu/cm^2 
    unfolded.Divide(unfolded,flux)

def GetEfficiency(ifile, ifile_truth, plot):
    if ifile_truth:
        print ("using separate efficiency file")
        ifile = ifile_truth
    num = Utilities.GetHistogram(ifile,PLOT_SETTINGS[plot+" Migration"]["name"]+"_truth")
    den = Utilities.GetHistogram(ifile,PLOT_SETTINGS["True Signal "+plot]["name"])
    hist_out,den_new = ProjectBinContent(num,den)
    hist_out.Divide(hist_out,den_new,1.0,1.0,"B")
 
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
    mc_truth_file = ROOT.TFile.Open(AnalysisConfig.SelectionHistoPath("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcEePtCollabRedo_nx-BigNuE_collab1_fspline.root",False))
    #mc_truth_file = ROOT.TFile.Open(AnalysisConfig.SelectionHistoPath("NueOnly_nx",False))
    xsec_file = ROOT.TFile.Open(AnalysisConfig.XSecHistoPath(playlist),"RECREATE")
    for plot in XSEC_TO_MAKE:
        PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)
        unfolded = Utilities.GetHistogram(unfolded_file,PLOT_SETTINGS[plot]["name"]+"_bkg_unfolding") 
        efficiency = GetEfficiency(mc_reco_file,mc_truth_file,plot) 
        xsec = GetXSectionHistogram(unfolded,efficiency,False)
     
        mc_xsec,sig_dep,colors,titles = GetMCXSectionHistogram(mc_reco_file,plot)
        ylabel = xsec.GetZaxis().GetTitle().replace("NEvents","d^{{2}}#sigma/d{X}d{Y}".format(X="E_{avail}",Y=plot.split()[-1]))
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
        #PlotTools.MNVPLOTTER.axis_maximum = 0.1
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec],draw_seperate_legend=False)

        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_err"))
        #PlotTools.AdaptivePlotterErrorGroup(xsec,18)
        #PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec],draw_seperate_legend=True)
        #PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_err_top7"))
        PlotTools.MNVPLOTTER.axis_maximum = -1111


        plotter = lambda mnvplotter,data_hist, mc_hist, *mc_ints : partial(PlotTools.MakeSignalDecomposePlot,color=colors,title=titles)(data_hist.GetCVHistoWithError(),mc_hist.GetCVHistoWithStatError(),mc_ints)
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[xsec,mc_xsec,*sig_dep],draw_seperate_legend=True)
        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_sigdep"))
        xsec_file.cd()
        xsec.Write(PLOT_SETTINGS[plot]["name"]+"_dataxsec")
        mc_xsec.Write(PLOT_SETTINGS[plot]["name"]+"_mcxsec")
        print("done")
    xsec_file.Close()
