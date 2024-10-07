"""
Get Crossection Plots.
Inputs are unfolded data histogram, efficiency histogram(or seperated true signal and selected signal),
"""

import ROOT
import PlotUtils
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities, PlotTools
from tools.PlotLibrary import PLOT_SETTINGS
from tools.PlotLibrary import HistHolder
from config.SystematicsConfig import USE_NUE_CONSTRAINT,CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS
from config import DrawingConfig
from config.DrawingConfig import PLOTS_TO_MAKE,Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer
from config.CutConfig import NEUTRINO_ENERGY_RANGE,FIDUCIAL_Z_RANGE
from functools import partial
mnvplotter = PlotUtils.MnvPlotter()

ROOT.TH1.AddDirectory(False)

XSEC_TO_MAKE = [
    #"q0 vs q3",
    "Visible Energy vs q3",
    "Visible Energy vs Lepton Pt"
]
USE_BIGNUE = True
threshold = 100 if USE_BIGNUE else 1
TARGET_UTILS = PlotUtils.TargetUtils.Get()
warping_errorband = ["fsi_weight","SuSA_Valencia_Weight","MK_model"]#,"LowQ2Pi_Joint","LowQ2Pi_NUPI0"]
#warping_errorband = ["SuSA_Valencia_Weight"]

FLUX="minervame1A"
#FLUX="minervame1d1m1nweightedave"

COLORS=ROOT.MnvColors.GetColors()
MODELS = {
    "MnvTune v1": {
        "errorband":(None,None),
        "color":COLORS[0]
    },
    "2p2h Tune (QE)": {
        "errorband":("Low_Recoil_2p2h_Tune",2),
        "color":COLORS[1]
    },
    "SuSA 2p2h" : {
        "errorband":("SuSA_Valencia_Weight",0),
        "color":COLORS[2]
    },
    "MK model": {
        "errorband": ("MK_model",0),
        "color":COLORS[7]
    },
#    "Low Q2 Pion Joint": {
#        "errorband" : ("LowQ2Pi_Joint",0),
#        "color": COLORS[5]
#    },
#    "Low Q2 Pion NuPi0": {
#        "errorband" : ("LowQ2Pi_NUPI0",0),
#        "color": COLORS[6]
#    }
}


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


def DivideFlux(unfolded,is_mc):
    frw= PlotUtils.flux_reweighter(FLUX,14,USE_NUE_CONSTRAINT) #playlist is dummy for now 
    flux = frw.GetIntegratedFluxReweighted(14,unfolded,NEUTRINO_ENERGY_RANGE[0],NEUTRINO_ENERGY_RANGE[1],False)
    flux.PopVertErrorBand("Flux_BeamFocus")
    flux.PopVertErrorBand("ppfx1_Total")
    flux.Scale(1e-4*(mc_pot if is_mc else data_pot)) #change unit to nu/cm^2
    print ("Flux Normalization: {:.4E},{:.4E}".format(flux.GetBinContent(1,1),flux.GetTotalError(False).GetBinContent(1,1)))
    unfolded.Divide(unfolded,flux)

def ProjectUniverseContent(num, den, i, j,num_o,den_o):
    l = 0
    denNewContent = 0
    numNewContent = 0
    while i-l>=0:
        denNewContent += den_o.GetBinContent(i-l,j)
        numNewContent += num_o.GetBinContent(i-l,j)
        l+=1
        if denNewContent>threshold:
            break

    den.SetBinContent(i,j,denNewContent)
    num.SetBinContent(i,j,numNewContent)
    #need to do this for every universe
    for bandname in den.GetVertErrorBandNames():
        nHists = den.GetVertErrorBand(bandname).GetNHists()
        for k in range(0,nHists):
            denNewContent = sum (den_o.GetVertErrorBand(bandname).GetHist(k).GetBinContent(i-_,j) for _ in range(l))
            numNewContent = sum (num_o.GetVertErrorBand(bandname).GetHist(k).GetBinContent(i-_,j) for _ in range(l))
            den.GetVertErrorBand(bandname).GetHist(k).SetBinContent(i,j,denNewContent)
            num.GetVertErrorBand(bandname).GetHist(k).SetBinContent(i,j,numNewContent)
    return  i-l>=0

def ProjectBinContent(num, den):
    num_new = num.Clone("{}_smooth".format(num.GetName()))
    den_new = den.Clone("{}_smooth".format(den.GetName()))
    xNBins = num.GetNbinsX()
    yNBins = num.GetNbinsY()
    for j in range(0,yNBins+2):
        for i in range(0,xNBins+2):
            if 0 < den.GetBinContent(i,j) < threshold: #treshold
                if not ProjectUniverseContent(num_new, den_new, i, j,num,den):
                    print ("not enough events in {}-th y bin".format(j))
    return num_new, den_new
    

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
    hist_out.GetXaxis().SetTitle("E_{avail} (GeV)") #hard coding for now
    hist_out.GetZaxis().SetTitle("Efficiency")

    plotter = lambda mnvplotter, hist: mnvplotter.DrawMCWithErrorBand(hist.GetCVHistoWithError())
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[hist_out],lambda x: True,False)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"eff"))
    #PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[hist_out],AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"eff"))
    plotter = lambda mnvplotter, hist: mnvplotter.DrawErrorSummary(hist)
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[hist_out],lambda x: True,False)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"eff_err"))
    return hist_out

def GetXSectionHistogram(unfolded,is_mc):
    #divide by flux
    DivideFlux(unfolded,is_mc)
    #divide by N nucleaon
    Nnucleon = TARGET_UTILS.GetTrackerNNucleons(FIDUCIAL_Z_RANGE[0],FIDUCIAL_Z_RANGE[1],is_mc)
    print("nucleon: ",Nnucleon)
    unfolded.Scale(1.0/Nnucleon)
    #h_nucleon = GetNnucleonError(unfolded,Nnucleon)  
    #unfolded.Divide(unfolded,h_nucleon) 
    return unfolded

def GetMCXSectionHistogram(mc_file,plot):
    # mc unfolded is truth_signal:
    
    unfolded = Utilities.GetHistogram(mc_file,PLOT_SETTINGS["True Signal "+plot]["name"])  
    DivideFlux(unfolded,True)
    Nnucleon = TARGET_UTILS.GetTrackerNNucleons(FIDUCIAL_Z_RANGE[0],FIDUCIAL_Z_RANGE[1],True)
    unfolded.Scale(1.0/Nnucleon)
    sig_dep = []
    colors = []
    titles = []
    for chan in DrawingConfig.SignalDecomposition.keys():#["CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"]: 
        if chan == "Background":
            continue

        tmp =Utilities.GetHistogram(mc_file,"{}_{}".format(PLOT_SETTINGS["True Signal "+plot]["name"],chan))
        
        DivideFlux(tmp,True)
        tmp.Scale(1.0/Nnucleon)
        sig_dep.append(tmp)
        colors.append(DrawingConfig.SignalDecomposition[chan]["color"])
        titles.append(DrawingConfig.SignalDecomposition[chan]["title"])
    return unfolded,sig_dep,colors,titles

def DrawModelComparison(data_hist,mc_hist,models=MODELS,band_on_mc=True):
    _cate = []
    _mc_models = []
    _colors = []
    for k,v in models.items():
        _cate.append(k)
        _colors.append(v["color"])
        htmp = mc_hist if band_on_mc else data_hist
        if v["errorband"][0]:
            try:
                _mc_models.append(PlotUtils.MnvH2D(htmp.GetVertErrorBand(v["errorband"][0]).GetHist(v["errorband"][1])))
            except ReferenceError:
                continue
        else:
            _mc_models.append(PlotUtils.MnvH2D(htmp.GetCVHistoWithStatError()))
    plotter = lambda mnvplotter,data_hist, *mc_ints : partial(PlotTools.MakeModelVariantPlot,color=_colors,title=_cate)(data_hist, mc_ints)
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[data_hist,*_mc_models],draw_seperate_legend=False)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_models"))
    # data_hist.GetZaxis().SetTitle("Data/MC")
    # mc_hist.GetZaxis().SetTitle("Data/MC")
    h0 = PlotUtils.MnvH2D(mc_hist.GetCVHistoWithStatError())
    for i in _mc_models:
        i.Divide(i,h0)
        i.GetZaxis().SetTitle("Data/MC")
    h0.AddMissingErrorBandsAndFillWithCV(data_hist)
    data_hist.Divide(data_hist,h0)
    data_hist.GetZaxis().SetTitle("Data/MC")
    PlotTools.MakeGridPlot(PlotTools.Make2DSlice,plotter,[data_hist,*_mc_models],draw_seperate_legend=True)
    PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"xsec_models_ratio"))

def setAxisTitles(h,plot):
    xlabel = "E_{avail}"
    ylabel = "q_{3}" if plot.split()[-1]=="q3" else "P^{t}_{lep}"
    zlabel = "d^{{2}}#sigma/d{X}d{Y}".format(X=xlabel,Y=ylabel)
    zlabel += " (#times 10^{39} cm^{2}/GeV^{2})"
    h.GetXaxis().SetTitle("{} (GeV)".format(xlabel))
    h.GetYaxis().SetTitle("{} (GeV)".format(ylabel))
    h.GetZaxis().SetTitle(zlabel)

if __name__ == "__main__":
    playlist= AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    standPOT = data_pot if data_pot is not None else mc_pot 

    AaronHist = HistHolder("True Signal Lepton Pt Aaron",mc_file,"Signal",True,mc_pot,standPOT)
    CCQEHist = HistHolder("True Signal Lepton Pt CCQE",mc_file,"Signal",True,mc_pot,standPOT)
    Default_Scale(AaronHist)
    Default_Scale(CCQEHist)

    aHist = AaronHist.hists["SinglePion"]
    ccHist = CCQEHist.hists["CCQE"]

    GetXSectionHistogram(aHist,True)
    GetXSectionHistogram(ccHist,True)
    h1_aHist = aHist.GetCVHistoWithStatError()
    h1_ccHist = ccHist.GetCVHistoWithStatError()
    mnv_aHist = aHist.Clone()
    mnv_ccHist = ccHist.Clone()

    mnv_aHist.DivideSingle(mnv_aHist,h1_aHist)
    c1 = ROOT.TCanvas()
    mnv_aHist.GetVertErrorBand("GENIE_D2_MaRES").DrawAll("HIST e",True)
    c1.Print("singlepion_ratio.png")
    c1.Clear()
    aHist.GetVertErrorBand("GENIE_D2_MaRES").DrawAll("HIST e",True)
    c1.Print("singlepion.png")
    
    mnv_ccHist.DivideSingle(mnv_ccHist,h1_ccHist)
    c1 = ROOT.TCanvas()
    mnv_ccHist.GetVertErrorBand("GENIE_MaZExpCCQE").DrawAll("HIST e",True)
    c1.Print("ccqe_ratio.png")
    c1.Clear()
    ccHist.GetVertErrorBand("GENIE_MaZExpCCQE").DrawAll("HIST e",True)
    c1.Print("ccqe.png")

    mnvplotter.DrawErrorSummary(ccHist,"TR",True,True,0.075)
    c1.Print("ccqe_errorband.png")
    mnvplotter.DrawErrorSummary(aHist,"TR",True,True,0.075)
    c1.Print("singlepion_errorband.png")

    mnvplotter.DrawMCWithErrorBand(ccHist.GetCVHistoWithError())
    c1.Print("ccqe_hist.png")
    mnvplotter.DrawMCWithErrorBand(aHist.GetCVHistoWithError())
    c1.Print("singlepion_hist.png")

