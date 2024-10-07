import ROOT
import math
from functools import partial
import PlotUtils

from tools.PlotLibrary import HistHolder
from config.AnalysisConfig import AnalysisConfig
from config.BackgroundFitConfig import CATEGORY_FACTORS
from config.SignalDef import SIGNAL_DEFINATION
from config.DrawingConfig import NuEElasticCategory
from tools import Utilities,PlotTools
ROOT.TH1.AddDirectory(False)

def SubtractPoissonHistograms(h,h1):
    h.AddMissingErrorBandsAndFillWithCV (h1)
    errors = []
    for i in range(h.GetSize()):
        errors.append(math.sqrt(h.GetBinError(i)**2 + h1.GetBinError(i)**2))
    h.Add(h1,-1)
    for i in range(h.GetSize()):
        h.SetBinError(i,errors[i])

    return h

def FitLeakage(mnvplotter,data_hist,mc_hist):
    mnvplotter.DrawDataMCWithErrorBand(data_hist,mc_hist,1.0,"TR")
    size = 0.05
    align = ROOT.Long()
    xLabel = ROOT.Double()
    yLabel = ROOT.Double()
    mnvplotter.DecodePosition("TR", size, align, xLabel, yLabel )
    data_hist.Fit("gaus","M0+","",0,100)
    fit = data_hist.GetFunction("gaus")
    if not fit:
        return None
    fit.Draw("SAME")
    u =fit.GetParameter(1)
    sig = fit.GetParameter(2)
    s = "data: {:.2f}#pm{:.2f}".format(u,sig)
    mnvplotter.AddPlotLabel(s,xLabel,yLabel,size,4,112,align)
    mc_hist.Fit("gaus","M0+","",0,100)
    fit = mc_hist.GetFunction("gaus")
    fit.Draw("SAME")
    u =fit.GetParameter(1)
    sig = fit.GetParameter(2)
    s = "mc: {:.2f}#pm{:.2f}".format(u,sig)
    mnvplotter.AddPlotLabel(s,xLabel,yLabel-0.05,size,4,112,align)
    
    


if __name__ == "__main__":
    #input knobs
    playlist=AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    standPOT = data_pot
    config = {"variables":["Visible Zoomin","Lepton Energy"]}
    data_hists = HistHolder(config["name"] if "name" in config else config,data_file,"Signal",False,data_pot,standPOT)
    mc_hists = HistHolder(config["name"] if "name" in config else config,mc_file,"Signal",True,mc_pot,standPOT)
    data_hists.POTScale(False)
    mc_hists.POTScale(False)
    mc_ints, color,title = mc_hists.GetCateList(NuEElasticCategory)
    for i,v in enumerate(title):
        if v !="nu+e":
            SubtractPoissonHistograms(data_hists.GetHist(),mc_ints[i])
        else:
            MC_NuE = mc_ints[i]
    slicer =  lambda hist: PlotTools.Make2DSlice(hist,interval=5)
    hists = [data_hists.GetHist(),MC_NuE]
    PlotTools.MakeGridPlot(slicer,FitLeakage,hists,lambda x:True ,True)
    PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","leakage_fit_testing"))


