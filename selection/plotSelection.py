import os
import sys
import ROOT
import PlotUtils

#insert path for modules of this package.
from config import PlotConfig
from config.AnalysisConfig import AnalysisConfig
from config.DrawingConfig import PLOTS_TO_MAKE,Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer
from tools import Utilities,PlotTools
from tools.PlotLibrary import HistHolder

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
SELECTED_SIDEBANDS = AnalysisConfig.sidebands

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
    slicer = config.setdefault("slicer", DefaultSlicer(data_hists)) if PlotType!="migration" else PlotTools.IdentitySlicer 
    draw_seperate_legend = config.setdefault("draw_seperate_legend",data_hists.dimension!=1 and PlotType != "migration")
    try:
        custom_tag = config["tag"]+PlotType if "tag" in config else PlotType
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

if __name__ == "__main__":
    #input knobs
    playlist=AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    standPOT = data_pot if data_pot is not None else mc_pot

    sideband_map = {}
    for config in PLOTS_TO_MAKE:
        sideband_group =  config.setdefault("sideband_group",["Signal"]+SELECTED_SIDEBANDS)
        if isinstance(sideband_group,list):
            for sideband in sideband_group:
                data_hists = HistHolder(config["name"] if "name" in config else config,data_file,sideband,False,data_pot,standPOT)
                mc_hists = HistHolder(config["name"] if "name" in config else config,mc_file,sideband,True,mc_pot,standPOT)
                MakePlot(data_hists,mc_hists,config)
        else:
            #assuing sideband_group is a tuple of name, and list of sidebands
            sideband = sideband_group[0]
            sidebands = sideband_group[1]
            data_hists = HistHolder(config["name"] if "name" in config else config,data_file,sidebands[0],False,data_pot,standPOT)
            mc_hists = HistHolder(config["name"] if "name" in config else config,mc_file,sidebands[0],True,mc_pot,standPOT)
            for _ in range(1,len(sidebands)):
                data_hists.Add(HistHolder(config["name"] if "name" in config else config,data_file,sidebands[_],False,data_pot,standPOT))
                mc_hists.Add(HistHolder(config["name"] if "name" in config else config,mc_file,sidebands[_],True,mc_pot,standPOT))
            MakePlot(data_hists,mc_hists,config)
