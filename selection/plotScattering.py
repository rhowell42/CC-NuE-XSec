import os
import sys
import ROOT
import PlotUtils
import numpy as np
from array import array
#insert path for modules of this package.
from config import PlotConfig
from config.AnalysisConfig import AnalysisConfig
from config.DrawingConfig import PLOTS_TO_MAKE,Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer
from tools import Utilities,PlotTools
from tools.PlotLibrary import HistHolder

MNVPLOTTER = PlotUtils.MnvPlotter()
# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
SELECTED_SIDEBANDS = AnalysisConfig.sidebands

def MakePlot(sideband_map,config):
    #scale
    muonHolder = sideband_map["True_Muon"]
    elecHolder = sideband_map["True_Electron"]

    Default_Scale(muonHolder)
    Default_Scale(elecHolder)
    CanvasConfig = config.setdefault("canvasconfig",lambda x:True)
    muonHist = muonHolder.GetHist().Clone()
    elecHist = elecHolder.GetHist().Clone()

    cm = ROOT.TCanvas("cm")
    muonHist.DrawCopy("colz")
    cm.Print("muon_scattering_template.png")
    ce = ROOT.TCanvas("ce")
    elecHist.DrawCopy("colz")
    ce.Print("elec_scattering_template.png")

    NUE_SCATTERING_BINNING = array('d',[0.8,2.,3.,5.,7.,9.,20.])
    cd = ROOT.TCanvas("cd")
    muonHist.Add(elecHist)
    profile = muonHist.ProfileX()
    h1 = ROOT.TH1D('h1','nue scattering',6,NUE_SCATTERING_BINNING)
    for i in range(profile.GetNbinsX()):
        h1.Fill(profile.GetBinContent(i),np.random.rand())

    h1.DrawCopy()
    cd.Print("nue_scattering_distr.png")

    return True

if __name__ == "__main__":
    #input knobs
    playlist=AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    sideband_map = {}
    c = 0
    for config in PLOTS_TO_MAKE:
        sideband_group =  config.setdefault("sideband_group",["Signal"]+SELECTED_SIDEBANDS)
        for sideband in sideband_group:
            mc_hists = HistHolder(config["name"] if "name" in config else config,mc_file,sideband,True,mc_pot,standPOT)
            if sideband != "Signal":
                sideband_map[sideband] = mc_hists
                c = config
        MakePlot(sideband_map,c)
