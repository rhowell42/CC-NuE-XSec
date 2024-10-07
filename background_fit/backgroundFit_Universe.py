import os
import sys
import ROOT
import PlotUtils
import math
import copy
import gc
from array import array
from collections import OrderedDict
import psutil

from tools.PlotLibrary import HistHolder
from config.AnalysisConfig import AnalysisConfig
from config import BackgroundFitConfig
from tools import Utilities,PlotTools
from config.UnfoldingConfig import HISTOGRAMS_TO_UNFOLD
from config.DrawingConfig import SignalOnly,Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer,PLOTS_TO_MAKE,SignalChargedBackground
from config.SignalDef import SIGNAL_DEFINATION
mnvplotter = PlotUtils.MnvPlotter()

from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
mnvplotter.error_summary_group_map.clear();
for k,v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v :
        vec.push_back(vs)
    mnvplotter.error_summary_group_map[k]= vec
# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

def Get1DScaleFactor(variable_hists,scale_hists):
    scale_dict = {}
    comparable_scale = MakeComparableMnvHXD(variable_hists.GetHist(), scale_hists, False)
    for cate in variable_hists.hists:
        
        if variable_hists.hists[cate] is None:
            continue
        scaled =variable_hists.hists[cate].Clone()
        try:
            scale = comparable_scale[BackgroundFitConfig.CATEGORY_FACTORS[cate]]
            scaled.Multiply(scaled,scale)  
            scale_dict[BackgroundFitConfig.CATEGORY_FACTORS[cate]] = scaled.Integral()/variable_hists.hists[cate].Integral() if variable_hists.hists[cate].Integral() != 0 else 0
        except KeyError:
            pass
        del scaled
    del comparable_scale
    return scale_dict

def MakeComparableMnvHXD(hist, scale_hist, y_axis=False):
    new_scale = {} 
    for cate in scale_hist: 
        new_scale[cate] = hist.Clone()
    xbins = hist.GetNbinsX()+2 #including under/overflows.
    for i in range(0,hist.GetSize()):
        nx = i%xbins
        ny = i//xbins
        scale_bin_entry = hist.GetYaxis().GetBinCenter(ny) if y_axis else hist.GetXaxis().GetBinCenter(nx)
        for cate in scale_hist:
            k = scale_hist[cate].FindBin(scale_bin_entry)
            new_scale[cate].SetBinContent(i,scale_hist[cate].GetBinContent(k))
            new_scale[cate].SetBinError(i,scale_hist[cate].GetBinError(k))
            
            for bandname in new_scale[cate].GetErrorBandNames():
                errorband = new_scale[cate].GetVertErrorBand(bandname)
                errorband.SetBinContent(i,scale_hists[cate].GetCVHistoWithStatError().GetBinContent(i))
                for ith in range(errorband.GetNHists()):
                    errorband.GetHist(ith).SetBinContent(i,scale_hist[cate].GetVertErrorBand(bandname).GetHist(ith).GetBinContent(k))
                    errorband.GetHist(ith).SetBinError(i,scale_hist[cate].GetVertErrorBand(bandname).GetHist(ith).GetBinError(k))

    return new_scale

def WriteScaleToMnvH1D(hist, scale, scale_err = None,  errorband=None,i=None):
    for group in scale:	
        if errorband is None:
            universe_hist = hist[group]
        elif i is None:
            universe_hist = hist[group].GetVertErrorBand(errorband)
        else:
            universe_hist= hist[group].GetVertErrorBand(errorband).GetHist(i)
        
        for q in range(0,universe_hist.GetNbinsX()+1):
            universe_hist.SetBinContent(q,scale[group].GetBinContent(q))
            universe_hist.SetBinError(q,scale[group].GetBinError(q))

def RunUniverseMinimizer(datasideband_histholders, datasignal_histholders, mcsideband_histholders, mcsignal_histholders, error_band = None, i = None):
    index = 0
    data_sideband = datasideband_histholders[index].GetHist()
    data_signal = datasignal_histholders[index].GetHist()
    mc_sidebandBKG = mcsideband_histholders[index].GetHist().Clone()
    mc_sidebandBKG.Reset()
    mc_signalBKG = mcsignal_histholders[index].GetHist().Clone()
    mc_signalBKG.Reset()
    mc_sidebandSIG = mcsideband_histholders[index].GetHist().Clone()
    mc_sidebandSIG.Reset()
    mc_signalSIG = mcsignal_histholders[index].GetHist().Clone()
    mc_signalSIG.Reset()

    for cate in mcsignal_histholders[index].hists:
        if cate == "Total":
            continue
        if cate not in SIGNAL_DEFINATION and cate != "NuEElastic":
            mc_sidebandBKG.Add(mcsideband_histholders[index].hists[cate])
            mc_signalBKG.Add(mcsignal_histholders[index].hists[cate])
        if cate in SIGNAL_DEFINATION:
            mc_sidebandSIG.Add(mcsideband_histholders[index].hists[cate])
            mc_signalSIG.Add(mcsignal_histholders[index].hists[cate])

    mc_signalNUEEL = mcsignal_histholders[index].hists["NuEElastic"].Clone()
    mc_sidebandNUEEL = mcsideband_histholders[index].hists["NuEElastic"].Clone()

    ### want a MC-like pseudodata signal region to avoid preliminary unblinding
    if AnalysisConfig.pseudodata:
        for q in range(0,data_signal.GetNbinsX()+1):
            data_signal.SetBinContent(q,mcsignal_histholders[index].hists["Total"].GetBinContent(q))
        for q in range(0,data_sideband.GetNbinsX()+1):
            data_sideband.SetBinContent(q,mcsideband_histholders[index].hists["Total"].GetBinContent(q))

    if error_band is not None and i is not None:
        mc_sidebandBKG = mc_sidebandBKG.GetVertErrorBand(error_band).GetHist(i).Clone()
        mc_signalBKG = mc_signalBKG.GetVertErrorBand(error_band).GetHist(i).Clone()
        mc_sidebandSIG = mc_sidebandSIG.GetVertErrorBand(error_band).GetHist(i).Clone()
        mc_signalSIG = mc_signalSIG.GetVertErrorBand(error_band).GetHist(i).Clone()
        mc_sidebandNUEEL = mc_sidebandNUEEL.GetVertErrorBand(error_band).GetHist(i).Clone()
        mc_signalNUEEL= mc_signalNUEEL.GetVertErrorBand(error_band).GetHist(i).Clone()
    elif error_band is not None:
        mc_sidebandBKG = mc_sidebandBKG.GetVertErrorBand(error_band).Clone()
        mc_signalBKG = mc_signalBKG.GetVertErrorBand(error_band).Clone()
        mc_sidebandSIG = mc_sidebandSIG.GetVertErrorBand(error_band).Clone()
        mc_signalSIG = mc_signalSIG.GetVertErrorBand(error_band).Clone()
        mc_sidebandNUEEL = mc_sidebandNUEEL.GetVertErrorBand(error_band).Clone()
        mc_signalNUEEL= mc_signalNUEEL.GetVertErrorBand(error_band).Clone()

    bkgscale = (mc_sidebandSIG * (mc_signalNUEEL - data_signal) + mc_signalSIG * (data_sideband - mc_sidebandNUEEL))/(mc_sidebandBKG * mc_signalSIG - mc_sidebandSIG * mc_signalBKG)
    sigscale = (mc_sidebandBKG * (data_signal - mc_signalNUEEL) + mc_signalBKG * (mc_sidebandNUEEL - data_sideband)) / (mc_sidebandBKG * mc_signalSIG - mc_sidebandSIG * mc_signalBKG)
    predscale = (data_sideband - mc_sidebandSIG - mc_sidebandNUEEL) / mc_sidebandBKG
    scales = {"signal":sigscale,"background":bkgscale,"prediction":predscale}
    del mc_sidebandBKG
    del mc_sidebandSIG
    del mc_sidebandNUEEL
    del mc_signalBKG
    del mc_signalSIG
    del mc_signalNUEEL

    return scales

def RunMinimizer(datasideband_histholders,datasignal_histholders, mcsideband_histholders, mcsignal_histholders,scale_hists):
    scales = RunUniverseMinimizer(datasideband_histholders,datasignal_histholders,mcsideband_histholders,mcsignal_histholders) 
    WriteScaleToMnvH1D(scale_hists,scales,None)

    #errorbands:
    for error_band in (mcsideband_histholders[0].GetHist().GetErrorBandNames()):
        #do errorband hist
        scales = RunUniverseMinimizer(datasideband_histholders,datasignal_histholders,mcsideband_histholders,mcsignal_histholders,error_band) 
        WriteScaleToMnvH1D(scale_hists,scales,None,error_band)

        for i in range(mcsideband_histholders[0].GetHist().GetVertErrorBand(error_band).GetNHists()):
            #do errorband universes 
            scales = RunUniverseMinimizer(datasideband_histholders,datasignal_histholders,mcsideband_histholders,mcsignal_histholders,error_band,i) 
            WriteScaleToMnvH1D(scale_hists,scales,None,error_band,i)

def TuneMC(hist_holder, scale_hists, x_axis=False, y_axis=False, prediction=False):
    if (x_axis and y_axis):
        return None # shouldnt happend
    elif not (x_axis or y_axis):
        try:
            ScaleCategories1D(hist_holder,scale_hists,prediction) #scale_dict is a global variable
        except AttributeError:
            return False
    else:
        comparable_scale = MakeComparableMnvHXD(hist_holder.GetHist(),scale_hists,y_axis)
        try:
            #ScaleCategories(hist_holder,comparable_scale)
            ScaleCategories(hist_holder,comparable_scale,prediction)
        except AttributeError:
            del comparable_scale
            return False
        del comparable_scale
    hist_holder.ResumTotal()
    return True

def ScaleCategories(hist_holder,scale_hists,prediction=False):
    for cate in hist_holder.hists:
        if cate == "Total":
            continue
        try:
            if not prediction:
                if cate not in SIGNAL_DEFINATION and cate != "NuEElastic":
                    scale = scale_hists["background"]
                    hist_holder.hists[cate].Multiply(hist_holder.hists[cate],scale)
                if cate in SIGNAL_DEFINATION:
                    scale = scale_hists["signal"]
                    hist_holder.hists[cate].Multiply(hist_holder.hists[cate],scale)
            else:
                if cate not in SIGNAL_DEFINATION and cate != "NuEElastic":
                    scale = scale_hists["background"]
                    hist_holder.hists[cate].Multiply(hist_holder.hists[cate],scale)

        except KeyError:
            print("KeyError with {} in {}".format(cate,hist_holder.sideband))
            continue

def ScaleCategories1D(hist_holder,scale_dict):
    for cate in hist_holder.hists:
        try:
            scale = scale_dict[BackgroundFitConfig.CATEGORY_FACTORS[cate]]
            hist_holder.hists[cate].Scale(scale,bin_width_normalize=True) 
        except KeyError:
            continue

def BackgroundSubtraction(data_hists, mc_hists, pred_hists, errs = None):
    data_hists.POTScale(False)
    mc_hists.POTScale(False)
    pred_hists.POTScale(False)
    out_data = data_hists.GetHist().Clone()
    out_mc = pred_hists.hists["Total"].Clone()
    out_data.AddMissingErrorBandsAndFillWithCV(out_mc)

    for group in mc_hists.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINATION:
            SubtractPoissonHistograms(out_data,mc_hists.hists[group]) #data tuned signal
            SubtractPoissonHistograms(out_mc,pred_hists.hists[group]) #no oscillation predicted signal

    return out_data,out_mc

def GetBackground(mc_hists):
    out_bkg = mc_hists.hists["Total"].Clone("bkgTotal")
    out_bkg.Reset()

    for group in mc_hists.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINATION:
            out_bkg.Add(mc_hists.hists[group])
    return out_bkg

def SubtractPoissonHistograms(h,h1):
    errors = []
    for i in range(h.GetSize()):
        errors.append(math.sqrt(h.GetBinError(i)**2 + h1.GetBinError(i)**2))
    h.Add(h1,-1)
    for i in range(h.GetSize()):
        h.SetBinError(i,errors[i])
    return h

def GetScaledDataMC(hist,datafile,mcfile,region):
    data_hist = HistHolder(hist,datafile,region,False,pot_scale)
    mc_hist = HistHolder(hist,mcfile,region,True,pot_scale)
    pred_hist = HistHolder(hist,mcfile,region,True,pot_scale)
    fit_on_axis = scaled_hist_name.upper() in data_hist.plot_name.upper() or "estimator" in data_hist.plot_name.lower()
    if fit_on_axis: # fit_on_axis = True
        fit_on_yaxis = ("_"+scaled_hist_name).upper() in data_hist.plot_name.upper() or "_estimator" in data_hist.plot_name.lower()
        #print(("fit {} on {} axis".format(data_hist.plot_name, "y" if fit_on_yaxis else "x")))
        TuneMC(mc_hist, scale_hists, not fit_on_yaxis , fit_on_yaxis )
        TuneMC(pred_hist, scale_hists, not fit_on_yaxis , fit_on_yaxis, True)
    else:
        #print(("not fitting {} on any axis".format(data_hist.plot_name)))
        variable_hist = HistHolder(BackgroundFitConfig.HIST_TO_FIT,mcfile,region,True,pot_scale) 
        scale_dict = Get1DScaleFactor(variable_hist,scale_hists)
        TuneMC(mc_hist, scale_dict, False , False )
        TuneMC(pred_hist, scale_hists, False, False, True)
        del scale_dict
        del variable_hist
    return data_hist,mc_hist,pred_hist

if __name__ == "__main__":
    #input knobs
    playlist=AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_path = type_path_map["data"]
    mc_path = type_path_map["mc"]

    dfile,mfile,pot_scale = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag)
    #output knobs:
    background_fit_tag = AnalysisConfig.bkgTune_tag
    BackgroundFitConfig.SetGlobalParameter(background_fit_tag)

    sel_histholder = HistHolder(BackgroundFitConfig.HIST_TO_FIT,mfile,"Signal",True,pot_scale)
    sid_histholder = HistHolder(BackgroundFitConfig.HIST_TO_FIT,mfile,"dEdX",True,pot_scale)
    scaled_hist_name = sel_histholder.plot_name

    datafile,mcfile,pot_scale = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag)
    datafile.Close()
    datafile.DeleteAll()
    mcfile.Close()
    mcfile.DeleteAll()
    dfile.Close()
    dfile.DeleteAll()
    mfile.Close()
    mfile.DeleteAll()
    scalefile=ROOT.TFile(AnalysisConfig.BackgroundFitPath(playlist,background_fit_tag),"RECREATE")
    scalefile.Close()
    scalefile.DeleteAll()

    for name in sel_histholder.GetHist().GetVertErrorBandNames():
        print("Running over {}".format(name))
        n_univ = sel_histholder.GetHist().GetVertErrorBand(name).GetNHists()
        if name == "Flux":
            n_univ = 100
        for universe in range(n_univ):
            datafile = ROOT.TFile(data_path)
            mcfile = ROOT.TFile(mc_path)
            scalefile=ROOT.TFile(AnalysisConfig.BackgroundFitPath(playlist,background_fit_tag),"UPDATE")
            print("Doing universe {}".format(universe))

            datasideband_histholders = []
            mcsideband_histholders = []
            datasignal_histholders = []
            mcsignal_histholders = []

            #fit scale histograms
            scale_hists = {"signal":None,"background":None}
            scale_hists["signal"] = sel_histholder.GetHist().Clone()
            scale_hists["signal"].Reset()
            scale_hists["signal"].GetYaxis().SetTitle("Scale Factor")
            scale_hists["signal"].SetTitle("Signal Scale Factor")
            scale_hists["background"] = sid_histholder.GetHist().Clone()
            scale_hists["background"].Reset()
            scale_hists["background"].GetYaxis().SetTitle("Scale Factor")
            scale_hists["background"].SetTitle("Background Scale Factor")
            scale_hists["prediction"] = sid_histholder.GetHist().Clone()
            scale_hists["prediction"].Reset()
            scale_hists["prediction"].GetYaxis().SetTitle("Scale Factor")
            scale_hists["prediction"].SetTitle("Background Scale Factor")

            for region in AnalysisConfig.sidebands:
                datasideband_histholders.append(HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,datafile,region,False))
                mcsideband_histholders.append(HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,mcfile,region,True,pot_scale))
                mcsideband_histholders[-1].POTScale(False)

            datasignal_histholders.append(HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,datafile,"Signal",False))

            mcsignal_histholders.append(HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,mcfile,"Signal",True,pot_scale))
            mcsignal_histholders[-1].POTScale(False)

            signalHist = HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,mcfile,"Signal",True,pot_scale)
            signalHist.POTScale(False)

            mc_prediction = signalHist.GetHist().Clone()
            mc_prediction.Reset()
            for cate in signalHist.hists:
                if cate in SIGNAL_DEFINATION:
                    mc_prediction.Add(signalHist.hists[cate])

            for h in mcsignal_histholders:
                ratio = h.GetHist().GetVertErrorBand(name).GetHist(universe).Clone()
                CV = h.GetHist().GetCVHistoWithStatError().Clone()
                ratio.Divide(ratio,CV)
                ratio = PlotUtils.MnvH1D(ratio)
                ratio.AddMissingErrorBandsAndFillWithCV(h.GetHist())
                for i in range(0,ratio.GetNbinsX()+1):
                    ratio.SetBinError(i,0)
                h.GetHist().Multiply(h.GetHist(),ratio)
                h.GetHist().RenameHistosAndErrorBands(ratio.GetName()+"_{}_{}".format(name,universe))

            ratio = mc_prediction.GetVertErrorBand(name).GetHist(universe).Clone()
            CV = mc_prediction.GetCVHistoWithStatError().Clone()
            ratio.Divide(ratio,CV)
            ratio = PlotUtils.MnvH1D(ratio)
            ratio.AddMissingErrorBandsAndFillWithCV(mc_prediction)
            for i in range(0,ratio.GetNbinsX()+1):
                ratio.SetBinError(i,0)
            mc_prediction.Multiply(mc_prediction,ratio)
            mc_prediction.RenameHistosAndErrorBands(ratio.GetName()+"_{}_{}".format(name,universe))


            RunMinimizer(datasideband_histholders,datasignal_histholders,mcsideband_histholders,mcsignal_histholders,scale_hists)

            region = "Signal"
            for factor in scale_hists: 
                hist = scale_hists[factor]
                hist.SetXTitle("E_{estimator}")
                hist.Write("{}_Scale_Factor_{}_{}".format(factor,name,universe))

            for hist in HISTOGRAMS_TO_UNFOLD:
                data_hist,mc_hist,pred_hist = GetScaledDataMC(hist,datafile,mcfile,region)
                mc_hist.GetHist().Write(data_hist.plot_name+"_{}_{}".format(name,universe))
                subbedData, subbedMC = BackgroundSubtraction(data_hist,mc_hist,pred_hist)
                subbedData.Write(data_hist.plot_name+"_data_bkgSubbed_{}_{}".format(name,universe)) #added here
                mc_prediction.Write(data_hist.plot_name+"_predicted_Signal_{}_{}".format(name,universe)) #added here

            scalefile.Close()
            datafile.Close()
            mcfile.Close()
            scalefile.DeleteAll()
            datafile.DeleteAll()
            mcfile.DeleteAll()
