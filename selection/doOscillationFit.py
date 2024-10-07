import os
import sys
import ROOT
import PlotUtils
import numpy as np
import math
from array import array

#insert path for modules of this package.
from config import PlotConfig
from config.AnalysisConfig import AnalysisConfig
from config.DrawingConfig import PLOTS_TO_MAKE,Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer
from config.SignalDef import SIGNAL_DEFINATION
from tools import Utilities,PlotTools
from tools.PlotLibrary import HistHolder

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.04
#MNVPLOTTER.extra_top_margin = -.035# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46 
MNVPLOTTER.mc_bkgd_line_color = 46

MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .04

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
SELECTED_SIDEBANDS = AnalysisConfig.sidebands

def GetOscillatedHistograms(m, theta_nue, theta_numu):
    #nue_hist  = unoscillated_normal.Clone().GetCVHistoWithError()
    #numu_hist = unoscillated_swapped.Clone().GetCVHistoWithError()
    nue_hist  = unoscillated_normal.Clone()
    numu_hist = unoscillated_swapped.Clone()

    for q in range(nue_hist.GetNbinsX()):
        P_ee = 1 - theta_nue*np.sin(1.27 * m * dis_template.GetBinContent(q))**2
        P_mue = theta_numu*np.sin(1.27 * m * app_template.GetBinContent(q))**2
        e_before = nue_hist.GetBinContent(q)
        mu_before = numu_hist.GetBinContent(q)
        left_nue = P_ee * e_before
        new_nue = P_mue * mu_before
        try:
            newerr = ((nue_hist.GetBinError(q)*left_nue/e_before)**2 + (left_nue**.5)**2)**.5
            appearednewerr = ((numu_hist.GetBinError(q)*new_nue/mu_before)**2 + (new_nue**.5)**2)**.5
        except:
            newerr = ((nue_hist.GetBinError(q))**2 + (left_nue**.5)**2)**.5 
            appearednewerr = ((numu_hist.GetBinError(q))**2 + (new_nue**.5)**2)**.5 
        nue_hist.SetBinContent(q,left_nue)
        numu_hist.SetBinContent(q,new_nue)

        nue_hist.SetBinError(q,newerr)
        numu_hist.SetBinError(q,appearednewerr)

    numu_hist.SetBinContent(numu_hist.GetNbinsX(),0)
    numu_hist.SetBinError(numu_hist.GetNbinsX(),0)
    return(nue_hist, numu_hist)

def FeldmanCousins():
    #scale
    theta_ee = np.linspace(0,1,20)
    theta_mue = np.linspace(0,0.1,9)
    deltam   = np.linspace(0,50,1)
    unoscillated = unoscillated_normal.Clone().GetCVHistoWithError()

    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(errors,True,True,False)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_fulldata_ErrorMatrix_Corr.pdf"))
    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(errors,True,False,False)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_fulldata_ErrorMatrix_Cov.pdf"))
    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(errors,True,True,True)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_fulldata_ErrorMatrix_Corr_asFrac.pdf"))

    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(unoscillated_normal,True,True,False)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_normal_nue_ErrorMatrix_Corr.pdf"))
    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(unoscillated_normal,True,False,False)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_normal_nue_ErrorMatrix_Cov.pdf"))
    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(unoscillated_normal,True,True,True)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_normal_nue_ErrorMatrix_Corr_asFrac.pdf"))

    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(unoscillated_swapped,True,True,False)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_swapped_nue_ErrorMatrix_Corr.pdf"))
    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(unoscillated_swapped,True,False,False)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_swapped_nue_ErrorMatrix_Cov.pdf"))
    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(unoscillated_swapped,True,True,True)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_swapped_nue_ErrorMatrix_Corr_asFrac.pdf"))

    #dis, app = GetOscillatedHistograms(m, ee, mue)
    #total_oscillation = dis.Clone()
    #total_oscillation.Add(app)
    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(total_oscillation,True,True,False)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_full_oscillation_ErrorMatrix_Corr.pdf"))
    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(total_oscillation,True,False,False)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_full_oscillation_ErrorMatrix_Cov.pdf"))
    #errorsC = ROOT.TCanvas()
    #MNVPLOTTER.DrawErrorMatrices(total_oscillation,True,True,True)
    #errorsC.Print(AnalysisConfig.PlotPath(datasignal.plot_name,sideband,"_full_oscillation_ErrorMatrix_Corr_asFrac.pdf"))
    fake_data = errors.Clone()
    for i in range(fake_data.GetNbinsX()):
        newerr = ((errors.GetCVHistoWithError(False).GetBinError(i))**2 + (fake_data.GetBinContent(i)**.5)**2)**.5
        fake_data.SetBinError(i,newerr)

    chi2minData = 1e10
    bestm = 0
    bestee = 0
    bestmue = 0
    for datafit_m in deltam:
        for datafit_ee in theta_ee:
            for datafit_mue in theta_mue:
                datadisTest, dataappTest = GetOscillatedHistograms(datafit_m,datafit_ee,datafit_mue)
                data_oscillatedFit = datadisTest+dataappTest
                chi2 = MNVPLOTTER.Chi2DataMC(data_oscillatedFit,fake_data)/fake_data.GetNbinsX()
                print(chi2)
                if chi2 < chi2minData:
                    chi2minData = chi2
                    bestm = datafit_m
                    bestee = datafit_ee
                    bestmue = datafit_mue

    print("best oscillation fit for m:{},s_ee:{},s_mue:{}".format(m,ee,mue))
    print("->"+" chi2:{},m:{},s_ee:{},s_mue:{}".format(chi2minData,bestm,bestee,bestmue))

    test_chi2s = []
    data_chi2s = []
    N_universes = 100
    alpha = 0.35
    exclusion_points = {"m":[],"s_ee":[],"s_mue":[]}

    for test_m in deltam:
        for test_ee in theta_ee:
            for test_mue in theta_mue:
                disMC, appMC = GetOscillatedHistograms(test_m, test_ee, test_mue)
                oscillatedMC = disMC + appMC
                for N in range(N_universes):
                    chi2minMC = 1e10
                    universe = oscillatedMC.Clone()

                    for i in range(universe.GetNbinsX()):
                        fluctuation = np.random.normal(oscillatedMC.GetBinContent(i),oscillatedMC.GetBinError(i))
                        if fluctuation <= 0:
                            fluctuation = 0
                        universe.SetBinContent(i,fluctuation)

                    for fit_m in deltam:
                        for fit_ee in theta_ee:
                            for fit_mue in theta_mue:
                                disTest, appTest = GetOscillatedHistograms(fit_m,fit_ee,fit_mue)
                                oscillatedFit = disTest+appTest
                                chi2 = MNVPLOTTER.Chi2DataMC(oscillatedFit,universe)/universe.GetNbinsX()
                                if chi2 < chi2minMC:
                                    chi2minMC = chi2

                    deltaChi2_test = MNVPLOTTER.Chi2DataMC(universe,unoscillated)/universe.GetNbinsX() - chi2minMC
                    test_chi2s.append(deltaChi2_test)

                deltaChi2_data = MNVPLOTTER.Chi2DataMC(oscillatedMC,fake_data)/universe.GetNbinsX() - chi2minData
                count = 0
                for entry in deltaChi2_test:
                    if deltaChi2_data < entry:
                        count += 1
                if count/N_universes < alpha:
                    exclusion_points["m"].append(test_m)
                    exclusion_points["s_ee"].append(test_ee)
                    exclusion_points["s_mue"].append(test_mue)
    print(exclusion_points)

def MakePlot():
    theta_ee = [0.05]
    theta_mue = [0.1]
    deltam   = [7.34]
    chi2     = 0
    unoscillated = unoscillated_normal.Clone().GetCVHistoWithError()

    for m in deltam:
        for ee in theta_ee:
            for mue in theta_mue:
                dis, app = GetOscillatedHistograms(m, ee, mue)
                total_oscillation = dis + app
                for i in range(total_oscillation.GetNbinsX()):
                    newerr = ((errors.GetCVHistoWithError(False).GetBinError(i))**2 + (total_oscillation.GetBinContent(i)**.5)**2)**.5
                    total_oscillation.SetBinError(i,newerr)

                margin = .12
                bottomFraction = .2
                overall = ROOT.TCanvas("Data/MC")
                top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
                bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

                top.Draw()
                bottom.Draw()

                top.cd()

                MNVPLOTTER.DrawDataMCWithErrorBand(total_oscillation,unoscillated,1,"TR")
                MNVPLOTTER.AddChi2Label(total_oscillation,unoscillated,1,"BC",.04,.45)

                bottom.cd()
                bottom.SetTopMargin(0)
                bottom.SetBottomMargin(0.3)

                ratio =  total_oscillation.Clone()
                ratio.Divide(ratio, unoscillated)

                #Now fill mcRatio with 1 for bin content and fractional error
                mcRatio = unoscillated_normal.GetTotalError(False, True, False) #The second "true" makes this fractional error
                for whichBin in range(1, mcRatio.GetXaxis().GetNbins()+1): 
                    mcRatio.SetBinError(whichBin, max(mcRatio.GetBinContent(whichBin), 1e-9))
                    mcRatio.SetBinContent(whichBin, 1)

                ratio.SetTitle("")
                ratio.SetLineColor(ROOT.kBlack)
                ratio.SetLineWidth(3)
                ratio.SetTitleSize(0)

                ratio.GetYaxis().SetTitle("Data / MC")
                ratio.GetYaxis().SetLabelSize(.15)
                ratio.GetYaxis().SetTitleSize(0.16)
                ratio.GetYaxis().SetTitleOffset(0.4)
                ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.

                ratio.GetXaxis().SetTitleSize(0.16)
                ratio.GetXaxis().SetTitleOffset(0.9)
                ratio.GetXaxis().SetLabelSize(.15)

                ratio.SetMinimum(.5)
                ratio.SetMaximum(1.5)
                ratio.Draw()

                #Error envelope for the MC
                mcRatio.SetLineColor(ROOT.kRed)
                mcRatio.SetLineWidth(3)
                mcRatio.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
                mcRatio.Draw("E2 SAME")

                #Draw a flat line at 1 for ratio of MC to itself
                straightLine = mcRatio.Clone()
                straightLine.SetFillStyle(0)
                straightLine.Draw("HIST SAME")
                top.cd()
                title = ROOT.TPaveText(0.3, 0.91, 0.7, 1.0, "nbNDC") #no border and use Normalized Device Coordinates to place it
                title.SetFillStyle(0)
                title.SetLineColor(ROOT.kWhite)
                title.AddText("Pseudo Data Oscillation")
                title.Draw()

                MNVPLOTTER.WritePreliminary(0.3, 0.82, 5e-2, True)
                overall.Print(AnalysisConfig.PlotPath(mcsignal.plot_name,sideband,"fit_oscillated_fullsample_fullsample.png"))            

                oscillatedC = ROOT.TCanvas("oscillated")
                dis.SetFillColor(ROOT.kBlue)
                #unoscillated.SetMarkerStyle(21)
                dis.SetMarkerColor(ROOT.kBlue)
                dis.SetLineColor(ROOT.kBlue)
                app.SetFillColor(ROOT.kRed)
                #app.SetMarkerStyle(21)
                app.SetMarkerColor(ROOT.kRed)
                app.SetLineColor(ROOT.kRed)
                dis.SetTitle("Remaining #nu_{e}")
                app.SetTitle("Oscillated #nu_{#mu}")
                stacked = ROOT.THStack("hs1","stacked")
                stacked.Add(dis)
                stacked.Add(app)
                unoscillated.SetMarkerStyle(MNVPLOTTER.data_marker);
                unoscillated.SetMarkerSize(MNVPLOTTER.data_marker_size);
                unoscillated.SetLineWidth(MNVPLOTTER.data_line_width);
                unoscillated.SetLineStyle(MNVPLOTTER.data_line_style);
                unoscillated.SetLineColor(MNVPLOTTER.data_color);
                unoscillated.SetMinimum(0)
                unoscillated.SetTitle("Neutrino-4 Hypothesis")
                unoscillated.SetXTitle("E_{e}+E_{available}")
                unoscillated.Draw("E1 X0")
                stacked.Draw("HIST SAME")
                unoscillated.Draw("SAME E1 X0")

                TLegend = ROOT.TLegend(0.68,0.72,0.98,0.98)
                TLegend.AddEntry(unoscillated,"no oscillation")
                TLegend.AddEntry(dis,"remaining #nu_{e}")
                TLegend.AddEntry(app,"oscillated #nu_{#mu}")
                TLegend.Draw("SAME")

                oscillatedC.Print(AnalysisConfig.PlotPath(mcsignal.plot_name,sideband,"full_stacked_oscillated_fullsample.png"))

                normalC = ROOT.TCanvas("Normal")
                out_sig.Draw("colz")
                normalC.Print(AnalysisConfig.PlotPath(mcsignal.plot_name,sideband,"normal_2dHist_fullsample_fullsample.png"))
                appearedC = ROOT.TCanvas("Swapped")
                appeared_out_sig.Draw("colz")
                appearedC.Print(AnalysisConfig.PlotPath(mcsignal.plot_name,sideband,"appeared_2dHist_fullsample_fullsample.png")) 

                aC = ROOT.TCanvas("aC")
                appDiv = app/unoscillated
                appDiv.SetMinimum(-0.1)
                appDiv.SetMaximum(1.2)
                appDiv.SetXTitle("E_{e}+E_{available}")
                appDiv.SetYTitle("ratio with no oscillation sample")
                appDiv.SetTitle("Neutrino-4 Hypothesis")
                appDiv.SetLineWidth(3)
                disDiv = dis/unoscillated
                disDiv.SetLineWidth(3)
                appDiv.Draw()
                disDiv.Draw("SAME")
                TLegend = ROOT.TLegend(0.40,0.40,0.70,0.550)
                TLegend.AddEntry(appDiv,"appearance ratio")
                TLegend.AddEntry(disDiv,"survival ratio")
                TLegend.Draw("SAME")
                aC.Print(AnalysisConfig.PlotPath(mcsignal.plot_name,sideband,"over_unoscillated_ratios_fullsample.png")) 

                profC = ROOT.TCanvas("profC")
                appeared_L_OVER_E.SetLineWidth(3)
                appeared_L_OVER_E.SetLineColor(ROOT.kRed)
                L_OVER_E.SetLineWidth(3)
                L_OVER_E.SetLineColor(ROOT.kBlue)
                appeared_L_OVER_E.SetXTitle("E_{e}+E_{available}")
                appeared_L_OVER_E.SetYTitle("True L/E Profiles")
                appeared_L_OVER_E.Draw()
                L_OVER_E.Draw("SAME")
                TLegend = ROOT.TLegend(0.5,0.60,0.8,0.750)
                TLegend.AddEntry(L_OVER_E,"Normal L/E Profile")
                TLegend.AddEntry(appeared_L_OVER_E,"Flavor Swapped L/E Profile")
                TLegend.Draw("SAME")
                profC.Print(AnalysisConfig.PlotPath(mcsignal.plot_name,sideband,"L_OVER_E_profiles_fullsample.png")) 


                diffC = ROOT.TCanvas("diffC")
                difference = disDiv - appDiv
                difference.SetMinimum(0)
                difference.SetMaximum(1.25)
                difference.SetXTitle("E_{e}+E_{available}")
                difference.SetYTitle("disappearence/normal - appearance/normal")
                difference.SetTitle("Neutrino-4 Hypothesis")
                difference.SetLineWidth(3)
                difference.Draw()
                diffC.Print(AnalysisConfig.PlotPath(mcsignal.plot_name,sideband,"disappearance_minus_appearance_fullsample.png"))

    return True

if __name__ == "__main__":
    #input knobs
    playlist=AnalysisConfig.playlist
    appeared_playlist = "me5A_swap"

    bkg_file_path = AnalysisConfig.BackgroundFitPath(playlist,AnalysisConfig.bkgTune_tag,False)

    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    appeared_type_path_map = { t:"/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcme5A_swap_flavorswapped_fspline.root" for t in AnalysisConfig.data_types}

    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    appeared_data_file,appeared_mc_file,appeared_pot_scale,appeared_data_pot,appeared_mc_pot = Utilities.getFilesAndPOTScale(appeared_playlist,appeared_type_path_map,AnalysisConfig.ntuple_tag,True)

    bkg = ROOT.TFile.Open(bkg_file_path)

    standPOT = data_pot if data_pot is not None else mc_pot 
    appeared_standPOT = appeared_data_pot if appeared_data_pot is not None else appeared_mc_pot 
    #signalstandPOT = signaldata_pot if signaldata_pot is not None else signalmc_pot 
    sideband_map = {}
    sideband = "Signal"
    L_OVER_E = HistHolder("Reco Energy vs L/E",mc_file,"Signal",True,mc_pot,standPOT)
    appeared_L_OVER_E = HistHolder("Reco Energy vs L/E",appeared_mc_file,"Signal",True,appeared_mc_pot,appeared_standPOT)
    mcsignal = HistHolder("Background Subbed MC",bkg,"Signal",True,mc_pot,standPOT)
    appeared_mcsignal = HistHolder("Biased Neutrino Energy",appeared_mc_file,"Signal",True,appeared_mc_pot,appeared_standPOT)
    datasignal = HistHolder("Background Subbed Data",bkg,"Signal",False,mc_pot,standPOT)

    Default_Scale(appeared_L_OVER_E)
    Default_Scale(L_OVER_E)
    Default_Scale(appeared_mcsignal)
    Default_Scale(mcsignal)
    Default_Scale(datasignal)

    out_sig = L_OVER_E.hists["Total"].Clone("sigTotal")
    out_sig.Reset()
    appeared_out_sig = appeared_L_OVER_E.hists["Total"].Clone("sigTotal")
    appeared_out_sig.Reset()
    appeared_signal = appeared_mcsignal.hists["Total"].Clone("sigTotal")
    appeared_signal.Reset()

    for group in L_OVER_E.hists:
        if group == "Total":
                continue
        elif group in SIGNAL_DEFINATION:
            if L_OVER_E.hists[group]:
                out_sig.Add(L_OVER_E.hists[group])
            if appeared_L_OVER_E.hists[group]:
                appeared_out_sig.Add(appeared_L_OVER_E.hists[group])
            if appeared_mcsignal.hists[group]:
                appeared_signal.Add(appeared_mcsignal.hists[group])

    L_OVER_E = out_sig.Clone().ProfileY()
    appeared_L_OVER_E = appeared_out_sig.Clone().ProfileY()
    mc = mcsignal.GetHist().Clone()
    data = datasignal.GetHist().Clone()

    unoscillated_normal = mc.Clone()
    unoscillated_swapped = appeared_signal.Clone()
    errors   = data.Clone()
    errors.PopVertErrorBand("SuSA_Valencia_Weight")
    errors.PopVertErrorBand("MK_model")
    errors.PopVertErrorBand("LowQ2Pi_None")
    unoscillated_normal.PopVertErrorBand("SuSA_Valencia_Weight")
    unoscillated_normal.PopVertErrorBand("MK_model")
    unoscillated_normal.PopVertErrorBand("LowQ2Pi_None")
    unoscillated_swapped.PopVertErrorBand("SuSA_Valencia_Weight")
    unoscillated_swapped.PopVertErrorBand("MK_model")
    unoscillated_swapped.PopVertErrorBand("LowQ2Pi_None")
    dis_template = L_OVER_E.Clone()
    app_template = appeared_L_OVER_E.Clone()
        
    MakePlot()
