import ROOT
import math
from functools import partial
import PlotUtils
import UnfoldUtils
import re

from tools.PlotLibrary import PLOT_SETTINGS,HistHolder
from config.AnalysisConfig import AnalysisConfig
from config.BackgroundFitConfig import CATEGORY_FACTORS
from config.SignalDef import SIGNAL_DEFINATION
from config import DrawingConfig
from tools import Utilities,PlotTools
from config.UnfoldingConfig import HISTOGRAMS_TO_UNFOLD,REGULATION_PARAMETERS
from config.DrawingConfig import Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer
from tools.unfoldingwrapper import DeOverflowWrapper,WithBackgroundWrapper,IdentityWrapper
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS

MNVUNFOLD = UnfoldUtils.MnvUnfold()
ROOT.TH1.AddDirectory(False)
USE_BIGNUE=True

def SubtractPoissonHistograms(h,h1,pseudo_data=False):
    if pseudo_data:
        h.ClearAllErrorBands ()
    h.AddMissingErrorBandsAndFillWithCV(h1)
    errors = []
    for i in range(h.GetSize()):
        errors.append(math.sqrt(h.GetBinError(i)**2 + h1.GetBinError(i)**2))
    h.Add(h1,-1)
    for i in range(h.GetSize()):
        if errors[i]==0:
            continue
        h.SetBinError(i,errors[i])
    return h


def MakePlot(data_hists,mc_hists,config):
    #scale
    if "scale" in config:
        config["scale"](data_hists)
        config["scale"](mc_hists)
    else:
        Default_Scale(data_hists)
        Default_Scale(mc_hists)
    CanvasConfig = config.setdefault("canvasconfig",lambda x:True)
    PlotType = config.setdefault("plot_type",Default_Plot_Type)
    slicer = config.setdefault("slicer", DefaultSlicer(data_hists))
    draw_seperate_legend = config.setdefault("draw_seperate_legend",data_hists.dimension!=1 and PlotType != "migration")
    try:
        custom_tag = config["tag"]+PlotType if "tag" in config else PlotType+"unfolded"
        if PlotType == "custom":
            plotfunction,hists=config["getplotters"](data_hists,mc_hists)
        else:
            if "args" in config:
                args = config["args"]
            elif "args" in DefaultPlotters[PlotType]:
                args = DefaultPlotters[PlotType]["args"]
            else:
                args = None
            if args is None:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists)
            else:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists,*args)
            PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
            PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,sideband,custom_tag))
            print("plot {} made.".format(data_hists.plot_name))
    except KeyError as e:
        print("plot {} not made.".format(data_hists.plot_name))
        print(e)
        return False
    return True


class unfolder(object):
    def __init__(**kwargs):
        pass

    def unfold():
        pass

    def backgroundFit():
        pass

def GetMigrationHistograms(mc_file, name):
    migration_hist = Utilities.GetHistogram(mc_file,PLOT_SETTINGS[name]["name"])
    migration_reco = Utilities.GetHistogram(mc_file,PLOT_SETTINGS[name]["name"]+"_reco")
    migration_truth = Utilities.GetHistogram(mc_file,PLOT_SETTINGS[name]["name"]+"_truth")
    return (migration_hist,migration_reco,migration_truth)

def unfolding(hist_to_unfold,migration,mc_reco,true_signal,mc_truth):
    cov = ROOT.TMatrixD(1,1)
    hist_to_unfold.AddMissingErrorBandsAndFillWithCV(mc_reco)
    migration.AddMissingErrorBandsAndFillWithCV(hist_to_unfold)
    mc_reco.AddMissingErrorBandsAndFillWithCV(hist_to_unfold)
    mc_truth.AddMissingErrorBandsAndFillWithCV(hist_to_unfold)
    data = hist_to_unfold.Clone()
    unfolded3 = true_signal.Clone()
    MNVUNFOLD.UnfoldHistoWithFakes(unfolded3,cov,migration,data,mc_reco,ROOT.nullptr,ROOT.nullptr,REGULATION_PARAMETERS[plot],True,True)
    return unfolded3

def DrawPostUnfolding(data_hist1,data_hist2):
    CanvasConfig = lambda x:True
    slicer = DefaultSlicer(data_hists)
    draw_seperate_legend = data_hists.dimension!=1
    plotfunction = lambda mnvplotter, data_hist, mc_hist: mnvplotter.DrawDataMCWithErrorBand(data_hist,mc_hist,1.0,"TR")
    hists = [data_hist2.GetCVHistoWithError(),data_hist1.GetCVHistoWithError()]
    PlotTools.updatePlotterErrorGroup(CONSOLIDATED_ERROR_GROUPS)

    PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
    PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","comp_unfolded_"+background_scale_tag))
    plotfunction = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCRatio(data_hist, mc_hist, 1.0 ,True,0,2)
    PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
    PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","ratio_unfolded_"+background_scale_tag))

    PlotTools.MNVPLOTTER.axis_maximum = 0.5
    plotfunction = lambda mnvplotter,data_hist: mnvplotter.DrawErrorSummary(data_hist)
    hists = [data_hist1]
    PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
    PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","err1_unfolded_"+background_scale_tag))
    hists = [data_hist2]
    PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
    PlotTools.MNVPLOTTER.axis_maximum = -1111
    PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","err2_unfolded_"+background_scale_tag))
    data_hist1.GetTotalErrorMatrix().Draw("COLZ")
    PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","err1matrix_unfolded_"+background_scale_tag))
    data_hist2.GetTotalErrorMatrix().Draw("COLZ")
    PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","err2matrix_unfolded_"+background_scale_tag))


if __name__ == "__main__":
    #input knobs
    playlist= AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_file,pot_scale = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag)
    background_scale_tag = AnalysisConfig.bkgTune_tag
    scale_file=ROOT.TFile.Open(AnalysisConfig.BackgroundFitPath(playlist,background_scale_tag))
    if not scale_file:
        scale_file = mc_file
        print ("Didn't find scale file. using mc file instead")
    signal_rich_file = ROOT.TFile.Open(AnalysisConfig.SelectionHistoPath(re.sub("-tune[0-9]","",playlist)+"-BigNuE",False)) if USE_BIGNUE else None
    if not signal_rich_file:
        signal_rich_file = mc_file
        print ("Didn't find signal rich file, using mc file instead")
    unfolded_file = ROOT.TFile.Open(AnalysisConfig.UnfoldedHistoPath(playlist,background_scale_tag),"RECREATE")


    for plot in HISTOGRAMS_TO_UNFOLD:
        data_hists= HistHolder(plot,data_file,"Signal",False,1.0)
        #mc_bkg = Utilities.GetHistogram(scale_file,data_hists.plot_name+"_mcbkg")
        mc_hists = HistHolder(plot,scale_file,"Signal",True,pot_scale)

        migration_name = plot + " Migration"
        migration_hists = GetMigrationHistograms(signal_rich_file,migration_name)
        true_signal = HistHolder("True Signal "+plot,mc_file,"Signal",True,pot_scale)
        signal_background,_,titles = mc_hists.GetCateList(DrawingConfig.SignalBackground)

        # at this point, three different normalization histogram could exist:
        # 1. data normalization, 2. standard MC (4x data) normalization, 3. signal rich normalization (back ground prediction is wrong)
        # if background including in unfolding: scale signal rich normaliztion to standard MC POT
        # if background excluded in unfolding: scale standard MC to data.



        #Nevent1 = migration_hists[1].Integral(0,-1,0,-1)
        #Nevent2 = signal_background[1].Integral(0,-1,0,-1)
        #print("Nevent",Nevent1,Nevent2)
        #migration_hists[1].Add(signal_background[0],Nevent1/Nevent2)
        signal_background[0].Scale(pot_scale)
        data = data_hists.GetHist()
        SubtractPoissonHistograms(data, signal_background[0],True)
        unfolded = unfolding(data,migration_hists[0],migration_hists[1],true_signal.GetHist(),migration_hists[2])
        #unfolded = unfold(plot,data)
        unfolded_file.cd()
        unfolded.Write(data_hists.plot_name+"_bkg_unfolding")
        #unfolded.Scale(1,"width")
        DrawPostUnfolding(unfolded,migration_hists[2])

    unfolded_file.Close()
