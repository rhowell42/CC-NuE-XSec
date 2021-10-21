import ROOT
import math
from functools import partial
import PlotUtils
import UnfoldUtils

from tools.PlotLibrary import PLOT_SETTINGS,HistHolder
from config.AnalysisConfig import AnalysisConfig
from config.BackgroundFitConfig import CATEGORY_FACTORS
from config.SignalDef import SIGNAL_DEFINATION
from config import PlotConfig
from tools import Utilities,PlotTools
from config.UnfoldingConfig import HISTOGRAMS_TO_UNFOLD,REGULATION_PARAMETER
from config.DrawingConfig import Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer
from tools.unfoldingwrapper import DeOverflowWrapper,WithBackgroundWrapper,IdentityWrapper

MNVUNFOLD = UnfoldUtils.MnvUnfold()
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
            print "plot {} made.".format(data_hists.plot_name)
    except KeyError as e:
        print "plot {} not made.".format(data_hists.plot_name)
        print e
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


# def unfold(name,hist_to_unfold = None):
#     print "unfolding"
#     migration_name = name + " Migration"
#     migration_hists = GetMigrationHistograms(mc_file,migration_name)
#     unfolded = PlotUtils.MnvH2D()
#     if hist_to_unfold is None:
#         hist_to_unfold = data_hists[PLOT_SETTINGS[name]["name"]]

#     MNVUNFOLD.UnfoldHisto2D(unfolded,migration_hists[0],migration_hists[1],migration_hists[2],hist_to_unfold,REGULATION_PARAMETER,True,True)
#     # MnvH2D method don't return covariance matrix for some reason, do it again by TH2D and push it back.

#     cov_matrix = ROOT.TMatrixD()
#     dummy_unfolded = ROOT.TH2D()
#     MNVUNFOLD.UnfoldHisto2D(dummy_unfolded,cov_matrix,migration_hists[0].GetCVHistoWithStatError(),migration_hists[1].GetCVHistoWithStatError(),migration_hists[2].GetCVHistoWithStatError(),hist_to_unfold.GetCVHistoWithStatError(),REGULATION_PARAMETER)
#     unfolded.SetTitle(hist_to_unfold.GetTitle())
#     # RooUnfold bug adds two extra bins to cov matrix. kill them.
#     cov_matrix.ResizeTo(dummy_unfolded.fN,dummy_unfolded.fN)
#     unfolded.PushCovMatrix("unfolding_cov_matrix_{}".format(name),cov_matrix)
#     return unfolded

# def unfolding(hist_to_unfold,migration,mc_reco,mc_truth,mc_bkg,true_signal):
#     cov = ROOT.TMatrixD(1,1)
#     #hist_to_unfold.AddMissingErrorBandsAndFillWithCV(mc_bkg)
#     data = hist_to_unfold.Clone()
#     mc_reco_new = mc_bkg.Clone()
#     mc_reco_new.Add(mc_reco)
#     unfolded3 = mc_truth.Clone()
#     MNVUNFOLD.UnfoldHistoWithFakes(unfolded3,cov,migration,data,mc_reco_new,ROOT.nullptr,ROOT.nullptr,REGULATION_PARAMETER,False,False)
#     del data
#     del mc_reco_new
#     print "done with fakes"
#     data = hist_to_unfold.Clone()
#     mc_bkg.Scale(pot_scale)
#     SubtractPoissonHistograms(data,mc_bkg)
#     unfolded2 = mc_truth.Clone()
#     MNVUNFOLD.UnfoldHistoWithFakes(unfolded2,cov,migration,data,mc_reco,ROOT.nullptr,ROOT.nullptr,REGULATION_PARAMETER,False,False)
#     del data
#     print unfolded2.GetBinContent(19),unfolded3.GetBinContent(19)
#     return unfolded2,unfolded3

def unfolding(hist_to_unfold,migration,mc_reco,true_signal,mc_truth):
    cov = ROOT.TMatrixD(1,1)
    hist_to_unfold.AddMissingErrorBandsAndFillWithCV(mc_reco)
    data = hist_to_unfold.Clone()
    unfolded3 = data.Clone()
    MNVUNFOLD.UnfoldHistoWithFakes(unfolded3,cov,migration,data,mc_reco,true_signal,ROOT.nullptr,REGULATION_PARAMETER,True,True)
    print true_signal.Integral(),unfolded3.Integral()
    return unfolded3

def DrawPreUnfolding(data_hists,mc_hists):
    data_hists.POTScale()
    mc_hists.POTScale()
    CanvasConfig = lambda x:True
    slicer = DefaultSlicer(data_hists)
    draw_seperate_legend = data_hists.dimension!=1
    plotfunction = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCWithErrorBand(data_hist,mc_hist,1.0,"TR")
    hists = [data_hists.GetHist().GetCVHistoWithError(),mc_hists.GetSumCate(["CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"]).GetCVHistoWithError()]
    PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
    PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","bkg_subbed_data_"+background_scale_tag))

    plotfunction = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCRatio(data_hist, mc_hist, 1.0 ,True,0,2)
    PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
    PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","bkg_subbed_data_ratio_"+background_scale_tag))
     
    plotfunction,hists = PlotTools.PrepareErr(data_hists,mc_hists,False,True)
    PlotTools.MNVPLOTTER.axis_maximum = 0.5
    PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
    PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","err_"+background_scale_tag))
    PlotTools.MNVPLOTTER.axis_maximum = -1111
    # plotfunction,hists = PlotTools.PrepareErrorBand(data_hists,mc_hists,"bkgfit",False,True)
    # PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
    # PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","errband_"+background_scale_tag))

def DrawPostUnfolding(data_hist1,data_hist2):
    CanvasConfig = lambda x:True
    slicer = DefaultSlicer(data_hists)
    draw_seperate_legend = data_hists.dimension!=1
    plotfunction = lambda mnvplotter, data_hist, mc_hist: mnvplotter.DrawDataMCWithErrorBand(data_hist,mc_hist,1.0,"TR")
    hists = [data_hist2.GetCVHistoWithError(),data_hist1.GetCVHistoWithError()]
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
    unfolded_file = ROOT.TFile.Open(AnalysisConfig.UnfoldedHistoPath(playlist,background_scale_tag),"RECREATE")


    for plot in HISTOGRAMS_TO_UNFOLD:
        data_hists= HistHolder(plot,data_file,"Signal",False,1.0)
        #mc_bkg = Utilities.GetHistogram(scale_file,data_hists.plot_name+"_mcbkg")
        mc_hists = HistHolder(plot,scale_file,"Signal",True,pot_scale)
        #mc_fixedBkgtmp.Scale(pot_scale)
        #mc_fixedBkg =mc_bkg.Clone()
        #mc_fixedBkg.Add(mc_fixedBkgtmp,-1)
        #data = SubtractPoissonHistograms(data_hists.hists["Total"],mc_bkg)
        #mc_hists = HistHolder(plot,mc_file,"Signal",True,pot_scale)
        #DrawPreUnfolding(data_hists,mc_hists)

        # plotfunction,hists = PlotTools.PrepareErrorBand(data_hists,mc_hists,"bkgfit",False,True)
        # PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
        # PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,"Signal","errband_"+background_scale_tag))
        #mc_unfolded.Write(data_hists.plot_name+"_{}closure".format("" if closure else "pseudo"))
        #data_unfolded = unfold(plot, data)
        # data_unfolded.Write(data_hists.plot_name+"_unfolded")
        # mc_truth = Utilities.GetHistogram(mc_file,PLOT_SETTINGS[plot+" Migration"]["name"]+"_truth")
        # mc_truth.GetXaxis().SetTitle(data_unfolded.GetXaxis().GetTitle())
        # mc_truth.GetYaxis().SetTitbble(data_unfolded.GetYaxis().GetTitle())
        # mc_truth.GetZaxis().SetTitle(data_unfolded.GetZaxis().GetTitle())
        # data_unfolded.Scale(1.0,"width")
        # mc_truth.Scale(pot_scale,"width")
        #PlotDataMC(data_unfolded,mc_truth,playlist,"unfolding" if not closure else "closure")

        migration_name = plot + " Migration"
        migration_hists = GetMigrationHistograms(mc_file,migration_name)
        true_signal = HistHolder("True Signal "+plot,mc_file,"Signal",True,pot_scale)
        unfolded = unfolding(data_hists.GetHist(),migration_hists[0],mc_hists.GetHist(),true_signal.GetHist(),migration_hists[2])
        #unfolded = unfold(plot,data)
        unfolded_file.cd()
        unfolded.Write(data_hists.plot_name+"_bkg_unfolding")
        #unfolded.Scale(1,"width")
        #DrawPostUnfolding(unfolded)

    unfolded_file.Close()
