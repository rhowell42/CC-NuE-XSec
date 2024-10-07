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


def GetMigrationHistograms(mc_file, name):
    migration_hist = Utilities.GetHistogram(mc_file,PLOT_SETTINGS[name]["name"])
    migration_reco = Utilities.GetHistogram(mc_file,PLOT_SETTINGS[name]["name"]+"_reco")
    migration_truth = Utilities.GetHistogram(mc_file,PLOT_SETTINGS[name]["name"]+"_truth")
    return (migration_hist,migration_reco,migration_truth)

def unfolding(hist_to_unfold,migration,mc_reco,true_signal,mc_truth,iteration_i):
    cov = ROOT.TMatrixD(1,1)
    hist_to_unfold.AddMissingErrorBandsAndFillWithCV(mc_reco)
    migration.AddMissingErrorBandsAndFillWithCV(hist_to_unfold)
    mc_reco.AddMissingErrorBandsAndFillWithCV(hist_to_unfold)
    mc_truth.AddMissingErrorBandsAndFillWithCV(hist_to_unfold)
    data = hist_to_unfold.Clone()
    unfolded3 = true_signal.Clone()
    MNVUNFOLD.UnfoldHistoWithFakes(unfolded3,cov,migration,data,mc_reco,ROOT.nullptr,ROOT.nullptr,iteration_i,True,True)
    del data
    return unfolded3


def IterationDrawer(mnvplotter,mc_ref,*iters,titles):
    colors = ROOT.MnvColors.GetColors()
    leg = ROOT.TLegend(0.6, 0.1, 0.9, 0.9)
    ROOT.SetOwnership(leg,False)
    for i,v in enumerate(iters):
        v.Divide(v,mc_ref)
        v.SetLineColor(colors[i])
        v.SetMaximum(2)
        #i.SetMinimum(0)
        v.Draw("HIST SAME")
        leg.AddEntry(v,titles[i])
    leg.Draw()

    # size = len(iters)+1
    # titles = ["mc"]+[h.GetTitle() for h in iters]
    # PlotTools.MakeModelVariantPlot(mc_ref,iters,color=colors[:size],title=titles)


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
        _,_,mc_ref = GetMigrationHistograms(mc_file,migration_name)
        mc_ref.GetYaxis().SetTitle("q_{3}" if plot.split()[-1]=="q3" else "P^{t}_{lep}")
        mc_ref.Scale(pot_scale)
        unfolded_list = []
        titles = []
        for i in range(1,16,2):
            unfolded_list.append( unfolding(data,migration_hists[0],migration_hists[1],true_signal.GetHist(),migration_hists[2],i))
            titles.append("iter{}".format(i))

        drawer = partial(IterationDrawer,titles=titles)

        print("making plot")
        PlotTools.MakeGridPlot(PlotTools.Make2DSlice,drawer,[mc_ref,*unfolded_list],draw_seperate_legend=True)
        PlotTools.Print(AnalysisConfig.PlotPath(PLOT_SETTINGS[plot]["name"],playlist,"unfolding_iters"))
        #unfolded = unfold(plot,data)
        # unfolded_file.cd()
        # unfolded.Write(data_hists.plot_name+"_bkg_unfolding")
        #unfolded.Scale(1,"width")
        del unfolded_list
    print("done")
    unfolded_file.Close()
