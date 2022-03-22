import os
import sys
import ROOT
import PlotUtils
from collections import OrderedDict
ROOT.TH1.AddDirectory(False)

#insert path for modules of this package.
from config import PlotConfig
from config.AnalysisConfig import AnalysisConfig
from config.DrawingConfig import PLOTS_TO_MAKE,Default_Plot_Type,Default_Scale,DefaultPlotters
from tools import Utilities,PlotTools
from tools.PlotLibrary import HistHolder
from config.CutConfig import SAMPLE_CUTS

config = {
    "variables": ["W"]
}
COLORS=ROOT.MnvColors.GetColors()
Categories = OrderedDict()
Categories["Signal"] = {
    "title":"Signal",
    "cate":{"CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE","CCNuE2p2h"},
    "color": COLORS[0]
}
Categories["Background"] = {
    "title":"Background",
    "color": COLORS[1],
    "cate" :{"Other","NonFiducial","NuEElastic","NCCOH","NCPi0","NCOther","CCPi0","CCOther","CCNuEAntiNu"}
}

if __name__ == "__main__":
    playlist=AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    # standPOT=data_pot
    # sideband = "Signal"
    
    # data_hists = HistHolder(config,data_file,sideband,False,data_pot,standPOT)
    # mc_hists = HistHolder(config,mc_file,sideband,True,mc_pot,standPOT)

    # mc_ints,_,titles = mc_hists.GetCateList(Categories)
    # data_hist = data_hists.GetHist() if data_hists else None
    # total_data = data_hist.Integral(0,-1) if data_hist else 1
    # total_signal = mc_ints[1].Integral(0,-1)
    # total_bkg = mc_ints[0].Integral(0,-1)

    # print(("Total data: {}, MC signal: {}, MC background: {}".format(total_data,total_signal,total_bkg)))
    # for i in range(1,26):
    #     cutted_data = data_hist.Integral(1,i) if data_hist else 1
    #     cutted_signal = mc_ints[1].Integral(1,i)
    #     cutted_bkg = mc_ints[0].Integral(1,i)
    #     #print("W cut at 0-{}GeV, remaining data: {}, MC signal: {}, MC background: {}".format(0.1*i,cutted_data,cutted_signal,cutted_bkg))
    #     #print("fraction of data remaining: {}, relative efficiency: {}, fractional backround remaining: {}".format(cutted_data/total_data,cutted_signal/total_signal,cutted_bkg/total_bkg))
    #     #print("metric: eff*pur: {}".format(cutted_signal/total_signal*cutted_signal/(cutted_bkg+cutted_signal)))
    #     print(("{} & {:.2f} & {:.2f} & {:.2f} & {:.4f}\\\\".format(0.1*i,cutted_signal/total_signal,cutted_bkg/total_bkg,cutted_signal/(cutted_bkg+cutted_signal),cutted_signal/total_signal*cutted_signal/(cutted_bkg+cutted_signal))))

    data_tree =ROOT.RDataFrame("cut_stat",data_file)
    mc_tree =ROOT.RDataFrame("cut_stat",mc_file)
    cut_stat = {}
    for cut in SAMPLE_CUTS["Signal"]:
        cut_stat[cut] = (data_tree.Sum("{}_{}".format(cut,"selected")).GetValue(),
                              mc_tree.Sum("{}_{}".format(cut,"selected")).GetValue(),
                              mc_tree.Sum("{}_{}".format(cut,"signal")).GetValue())
        #print (cut_stat[cut])
    for cut,v in sorted(cut_stat.items(),key=lambda x: x[1][0],reverse=True):
        print ("{} & {:.0f} & {:.0f} & {:.0f} \\".format(cut,*v))
