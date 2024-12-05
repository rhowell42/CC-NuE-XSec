import os
from collections import OrderedDict
import argparse
import logging, sys
import ROOT
import PlotUtils
from Tools.FitTools import *
from Tools.PlotTools import *
from Tools.StitchTools import *
from Tools.OscHistogram import *

import numpy as np
#np.random.seed(0)
#np.set_printoptions(precision=3)
#np.set_printoptions(linewidth=1520)
#np.set_printoptions(threshold=sys.maxsize)

import math
from array import array

#insert path for modules of this package.
#from config import PlotConfig
#from config.AnalysisConfig import AnalysisConfig
#from config.DrawingConfig import PLOTS_TO_MAKE,Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer
from config.SignalDef import SWAP_SIGNAL_DEFINITION, SIGNAL_DEFINITION
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
from tools import Utilities
from tools.PlotLibrary import HistHolder

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
MNVPLOTTER.draw_normalized_to_bin_width=False
#MNVPLOTTER.extra_top_margin = -.033# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46 
MNVPLOTTER.mc_bkgd_line_color = 46
#MNVPLOTTER.mc_line_width = 0

MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.1
MNVPLOTTER.width_xspace_per_letter = .3
MNVPLOTTER.legend_text_size        = .07
MNVPLOTTER.legend_n_columns        = 1
MNVPLOTTER.axis_draw_grid_x = True
MNVPLOTTER.axis_draw_grid_y = True

MNVPLOTTER.error_summary_group_map.clear();
for k,v in CONSOLIDATED_ERROR_GROUPS.items():
    vec = ROOT.vector("std::string")()
    for vs in v :
        vec.push_back(vs)
    MNVPLOTTER.error_summary_group_map[k]= vec

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

muonList =  ["Muon_Energy_Resolution","Muon_Energy_MINOS","Muon_Energy_MINERvA"]
errsToRemove = ["LowQ2Pi"]

def SwapUniverseCV(h_in,err,univ):
    if err in h_in.GetVertErrorBandNames():
        ratio = h_in.GetVertErrorBand(err).GetHist(univ).Clone()
    else:
        print("Histogram doesn't have {} universe in errorband {}, skipping...".format(univ,err))
        return(h_in)

    CV = h_in.GetCVHistoWithStatError().Clone()
    ratio.Divide(ratio,CV)
    ratio = PlotUtils.MnvH1D(ratio)
    ratio.AddMissingErrorBandsAndFillWithCV(h_in)
    for i in range(0,ratio.GetNbinsX()+1):
        ratio.SetBinError(i,0)

    h_in.Multiply(h_in,ratio)
    h_in.RenameHistosAndErrorBands(ratio.GetName()+"_{}_{}".format(err,univ))
    return(h_in)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--exclude",
                        dest ="exclude",
                        default = "None",
    )
    parser.add_argument("--ratio",
                        dest ="ratio",
                        default = False,
                        action = "store_true",
    )
    parser.add_argument("--fit_muons",
                        dest ="fit_muons",
                        default = False,
                        action = "store_true",
    )

    args = parser.parse_args()
    exclude = str(args.exclude).lower()
    doratio = args.ratio
    fit_muons = args.fit_muons
    if doratio:
        ftag = "ratio"
    else:
        ftag = "stitched"

    # ---------------------- NuEEl Selection Block -----------------------------
    h_fhc_nueel_mc = ROOT.TFile('/exp/minerva/data/users/rhowell/nueel/DeepikaMc_LauraData_Mnv.root').Get('mnv_eff_cor_mc')
    h_fhc_nueel_data = ROOT.TFile('/exp/minerva/data/users/rhowell/nueel/DeepikaMc_LauraData_Mnv.root').Get('mnv_eff_cor_data')

    h_rhc_nueel_mc = ROOT.TFile("/exp/minerva/data/users/rhowell/nueel/electronE_bkgsub_FullSetFinalv2_nobinnorm.root").Get('mc_fluxonly')   #.Get('mc_bkgsub_effcor')
    h_rhc_nueel_data =  ROOT.TFile("/exp/minerva/data/users/rhowell/nueel/electronE_bkgsub_FullSetFinalv2_nobinnorm.root").Get('data_bkgsub_effcor')

    # ---------------------- IMD Selection Block -----------------------------
    h_fhc_imd_mc = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("FHC_MC")
    h_fhc_imd_data = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("FHC_Data")

    h_rhc_imd_mc = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("RHC_MC")
    h_rhc_imd_data = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("RHC_Data")

    # ---------------------- NuMu Selection Block -----------------------------
    type_path_map = {'data':'/exp/minerva/data/users/rhowell/nu_mu/kin_dist_dataFHC_Selection_100Univ_thesis_muon_MAD.root','mc':'/exp/minerva/data/users/rhowell/nu_mu/kin_dist_mcFHC_Selection_100Univ_thesis_muon_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    fhc_muselection_mcHold = HistHolder("Biased Neutrino Energy",mc_file,"Signal",True,mc_pot,standPOT)
    fhc_muselection_dataHold = HistHolder("Biased Neutrino Energy",data_file,"Signal",False,data_pot,standPOT)

    type_path_map = {'data':'/exp/minerva/data/users/rhowell/antinu_mu/kin_dist_dataRHC_Selection_1000Univ_thesis_muon_MAD.root','mc':'/exp/minerva/data/users/rhowell/antinu_mu/kin_dist_mcRHC_Selection_1000Univ_thesis_muon_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_100Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    rhc_muselection_mcHold = HistHolder("Biased Neutrino Energy",mc_file,"Signal",True,mc_pot,standPOT)
    rhc_muselection_dataHold = HistHolder("Biased Neutrino Energy",data_file,"Signal",False,data_pot,standPOT)

    h_fhc_muselection_mc = fhc_muselection_mcHold.GetHist()
    h_rhc_muselection_mc = rhc_muselection_mcHold.GetHist()
    h_fhc_muselection_data = fhc_muselection_dataHold.GetHist()
    h_rhc_muselection_data = rhc_muselection_dataHold.GetHist()

    fhc_muselection_mcHold.POTScale(False)
    rhc_muselection_mcHold.POTScale(False)
    fhc_muselection_dataHold.POTScale(False)
    rhc_muselection_dataHold.POTScale(False)

    h_fhc_muselection_mc.Reset()
    h_rhc_muselection_mc.Reset()

    cates = ["CCQE","CCDelta","CCDIS","CC2p2h","CCOther","CCWrongSign"]
    for group in fhc_muselection_mcHold.hists:
        if group == "Total":
            continue
        elif group in cates:
            h_fhc_muselection_mc.Add(fhc_muselection_mcHold.hists[group])

    for group in rhc_muselection_mcHold.hists:
        if group == "Total":
            continue
        elif group in cates:
            h_rhc_muselection_mc.Add(rhc_muselection_mcHold.hists[group])

    fhc_universe = "/exp/minerva/data/users/rhowell/nu_e/bkgfit_FHC_Selection_100Univ_Sigmas_N4_tune_thesis_MAD.root"
    rhc_universe = "/exp/minerva/data/users/rhowell/antinu_e/bkgfit_RHC_Selection_100Univ_Sigmas_N4_tune_thesis_MAD.root"

    chi2filename = "chi2_values_{}_{}.txt".format(ftag,exclude.replace(" ","_"))
    with open(chi2filename,'w') as file:
        file.write("errorband")
        for i in range(1,101):
            file.write(", Universe {}".format(i))
        file.write("\n")

        # ---------------------- NuE CV Selection Block -----------------------------
        fhc_cv = "/exp/minerva/data/users/rhowell/nu_e/bkgfit_FHC_Selection_100Univ_N4_tune_thesis_MAD.root"
        h_fhc_selection_mc   = ROOT.TFile.Open(fhc_cv).Get("EN4_predicted_Signal")
        h_fhc_selection_data = ROOT.TFile.Open(fhc_cv).Get("EN4_data_bkgSubbed")

        rhc_cv = "/exp/minerva/data/users/rhowell/antinu_e/bkgfit_RHC_Selection_100Univ_N4_tune_thesis_MAD.root"
        h_rhc_selection_mc   = ROOT.TFile.Open(rhc_cv).Get("EN4_predicted_Signal")
        h_rhc_selection_data = ROOT.TFile.Open(rhc_cv).Get("EN4_data_bkgSubbed")

        # ---------------------- Add Electron Energy Scale Systematic to Selection -----------------------------
        fhc_scale_CV = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Electron_Scale_electron_scale_MAD.root").Get("EN4")
        fhc_scale_p1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_p1sig")
        fhc_scale_m1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_m1sig")
        fhc_scale_p1sig = fhc_scale_p1sig/fhc_scale_CV
        fhc_scale_m1sig = fhc_scale_m1sig/fhc_scale_CV
        rhc_scale_CV = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Electron_Scale_electron_scale_MAD.root").Get("EN4")
        rhc_scale_p1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_p1sig")
        rhc_scale_m1sig = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Electron_Scale_electron_scale_MAD.root").Get("EN4_m1sig")
        rhc_scale_p1sig = rhc_scale_p1sig/rhc_scale_CV
        rhc_scale_m1sig = rhc_scale_m1sig/rhc_scale_CV

        h_fhc_selection_mc.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
        for i in range(h_fhc_selection_mc.GetNbinsX()+1):
            h_fhc_selection_mc.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,h_fhc_selection_mc.GetBinContent(i) * fhc_scale_m1sig.GetBinContent(i))
            h_fhc_selection_mc.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,h_fhc_selection_mc.GetBinContent(i) * fhc_scale_p1sig.GetBinContent(i))

        h_rhc_selection_mc.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
        for i in range(h_rhc_selection_mc.GetNbinsX()+1):
            h_rhc_selection_mc.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,h_rhc_selection_mc.GetBinContent(i) * rhc_scale_m1sig.GetBinContent(i))
            h_rhc_selection_mc.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,h_rhc_selection_mc.GetBinContent(i) * rhc_scale_p1sig.GetBinContent(i))

        # ---------------------- Create Stitched CV Histograms -----------------------------
        cv_histogram = StitchedHistogram("cv_histogram")

        # ----- Initialize histogram objects with all samples ----- #
        cv_histogram.AddHistogram('fhc_nueel',h_fhc_nueel_mc,h_fhc_nueel_data)
        cv_histogram.AddHistogram('fhc_imd',h_fhc_imd_mc,h_fhc_imd_data)
        cv_histogram.AddHistogram('fhc_muselection',h_fhc_muselection_mc,h_fhc_muselection_data)
        cv_histogram.AddHistogram('fhc_selection',h_fhc_selection_mc,h_fhc_selection_data)
        cv_histogram.AddHistogram('rhc_nueel',h_rhc_nueel_mc,h_rhc_nueel_data)
        cv_histogram.AddHistogram('rhc_imd',h_rhc_imd_mc,h_rhc_imd_data)
        cv_histogram.AddHistogram('rhc_muselection',h_rhc_muselection_mc,h_rhc_muselection_data)
        cv_histogram.AddHistogram('rhc_selection',h_rhc_selection_mc,h_rhc_selection_data)

        # ----- Remove samples that we want to exclude from analysis ----- #
        cv_histogram.ApplyExclusion(exclude)

        # ----- Process Systematics and Synchronize across histograms ----- #
        cv_histogram.CleanErrorBands(errsToRemove)

        if doratio: # do we want to replace selection samples with flavor ratios
            if "fhc" not in exclude: # do we care about the fhc component
                cv_histogram.MakeRatio('fhc')
            if "rhc" not in exclude: # do we care about the rhc component
                cv_histogram.MakeRatio('rhc')
            if not fit_muons: # do we want to keep the muon selections in addition to flavor ratios
                cv_histogram.RemoveHistogram('fhc_muselection')
                cv_histogram.RemoveHistogram('rhc_muselection')

        # ----- Stitch histograms together ----- #
        cv_histogram.Stitch()

        # ---------------------- Start Swapping Universes -----------------------------
        err_names = h_rhc_selection_mc.GetVertErrorBandNames()
        for err in err_names:
            err = str(err)
            n_universes = cv_histogram.mc_hist.GetVertErrorBand(err).GetNHists()
            print("loading {} errorband that has {} n_universes".format(err,n_universes))
            chi2s = []
            for univ in range(n_universes):
                if err != "ElectronScale" and err in rhc_scale_CV.GetVertErrorBandNames():
                    # ---------------------- NuE Universe Selection Block -----------------------------
                    h_fhc_selection_mc_universe   = ROOT.TFile.Open(fhc_universe).Get("EN4_predicted_Signal_{}_{}".format(err,univ))
                    h_fhc_selection_data_universe = ROOT.TFile.Open(fhc_universe).Get("EN4_data_bkgSubbed_{}_{}".format(err,univ))

                    h_rhc_selection_mc_universe   = ROOT.TFile.Open(rhc_universe).Get("EN4_predicted_Signal_{}_{}".format(err,univ))
                    h_rhc_selection_data_universe = ROOT.TFile.Open(rhc_universe).Get("EN4_data_bkgSubbed_{}_{}".format(err,univ))

                    h_fhc_selection_mc_universe.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
                    for i in range(h_fhc_selection_mc_universe.GetNbinsX()+1):
                        h_fhc_selection_mc_universe.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,h_fhc_selection_mc_universe.GetBinContent(i) * fhc_scale_m1sig.GetBinContent(i))
                        h_fhc_selection_mc_universe.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,h_fhc_selection_mc_universe.GetBinContent(i) * fhc_scale_p1sig.GetBinContent(i))

                    h_rhc_selection_mc_universe.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
                    for i in range(h_rhc_selection_mc_universe.GetNbinsX()+1):
                        h_rhc_selection_mc_universe.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,h_rhc_selection_mc_universe.GetBinContent(i) * rhc_scale_m1sig.GetBinContent(i))
                        h_rhc_selection_mc_universe.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,h_rhc_selection_mc_universe.GetBinContent(i) * rhc_scale_p1sig.GetBinContent(i))

                elif err in muonList:
                    h_fhc_selection_mc_universe = h_fhc_selection_mc.Clone()
                    h_rhc_selection_mc_universe = h_rhc_selection_mc.Clone()
                    h_fhc_selection_data_universe = h_fhc_selection_data.Clone()
                    h_rhc_selection_data_universe = h_rhc_selection_data.Clone()
                else:
                    h_fhc_selection_mc_universe = h_fhc_selection_mc.Clone()
                    h_rhc_selection_mc_universe = h_rhc_selection_mc.Clone()
                    h_fhc_selection_data_universe = h_fhc_selection_data.Clone()
                    h_rhc_selection_data_universe = h_rhc_selection_data.Clone()
                    h_fhc_selection_mc_universe = SwapUniverseCV(h_fhc_selection_mc_universe,err,univ)
                    h_rhc_selection_mc_universe = SwapUniverseCV(h_rhc_selection_mc_universe,err,univ)

                # ---------------------- NuMu Universe Selection Block -----------------------------
                h_fhc_muselection_mc_universe = h_fhc_muselection_mc.Clone()
                h_rhc_muselection_mc_universe = h_rhc_muselection_mc.Clone()
                h_fhc_muselection_data_universe = h_fhc_muselection_data.Clone()
                h_rhc_muselection_data_universe = h_rhc_muselection_data.Clone()

                # ---------------------- Scattering Universe Selection Block -----------------------------
                h_fhc_nueel_mc_universe = h_fhc_nueel_mc.Clone()
                h_rhc_nueel_mc_universe = h_rhc_nueel_mc.Clone()
                h_fhc_imd_mc_universe = h_fhc_imd_mc.Clone()
                h_rhc_imd_mc_universe = h_rhc_imd_mc.Clone()

                h_fhc_nueel_data_universe = h_fhc_nueel_data.Clone()
                h_rhc_nueel_data_universe = h_rhc_nueel_data.Clone()
                h_fhc_imd_data_universe = h_fhc_imd_data.Clone()
                h_rhc_imd_data_universe = h_rhc_imd_data.Clone()
               
                # ---------------------- Swap Each Systematic Universe with the CV -----------------------------
                # Do monte carlo histograms
                h_fhc_muselection_mc_universe = SwapUniverseCV(h_fhc_muselection_mc_universe,err,univ)
                h_rhc_muselection_mc_universe = SwapUniverseCV(h_rhc_muselection_mc_universe,err,univ)
                h_fhc_nueel_mc_universe = SwapUniverseCV(h_fhc_nueel_mc_universe,err,univ)
                h_rhc_nueel_mc_universe = SwapUniverseCV(h_rhc_nueel_mc_universe,err,univ)
                h_fhc_imd_mc_universe = SwapUniverseCV(h_fhc_imd_mc_universe,err,univ)
                h_rhc_imd_mc_universe = SwapUniverseCV(h_rhc_imd_mc_universe,err,univ)

                # Do data histograms to mimic background subtraction effect, muon selections have no data systematics
                h_fhc_nueel_data_universe = SwapUniverseCV(h_fhc_nueel_data_universe,err,univ)
                h_rhc_nueel_data_universe = SwapUniverseCV(h_rhc_nueel_data_universe,err,univ)
                h_fhc_imd_data_universe = SwapUniverseCV(h_fhc_imd_data_universe,err,univ)
                h_rhc_imd_data_universe = SwapUniverseCV(h_rhc_imd_data_universe,err,univ)

                # ---------------------- Create Stitched Universe Histograms -----------------------------
                universe_histogram = StitchedHistogram(err+"_{}".format(univ))

                # ----- Initialize histogram objects with all samples ----- #
                universe_histogram.AddHistogram('fhc_nueel',h_fhc_nueel_mc_universe,h_fhc_nueel_data_universe)
                universe_histogram.AddHistogram('fhc_imd',h_fhc_imd_mc_universe,h_fhc_imd_data_universe)
                universe_histogram.AddHistogram('fhc_muselection',h_fhc_muselection_mc_universe,h_fhc_muselection_data_universe)
                universe_histogram.AddHistogram('fhc_selection',h_fhc_selection_mc_universe,h_fhc_selection_data_universe)
                universe_histogram.AddHistogram('rhc_nueel',h_rhc_nueel_mc_universe,h_rhc_nueel_data_universe)
                universe_histogram.AddHistogram('rhc_imd',h_rhc_imd_mc_universe,h_rhc_imd_data_universe)
                universe_histogram.AddHistogram('rhc_muselection',h_rhc_muselection_mc_universe,h_rhc_muselection_data_universe)
                universe_histogram.AddHistogram('rhc_selection',h_rhc_selection_mc_universe,h_rhc_selection_data_universe)
                
                # ----- Remove samples that we want to exclude from analysis ----- #
                universe_histogram.ApplyExclusion(exclude)

                # ----- Process Systematics and Synchronize across histograms ----- #
                universe_histogram.CleanErrorBands(errsToRemove)

                if doratio: # do we want to replace selection samples with flavor ratios
                    if "fhc" not in exclude: # do we care about the fhc component
                        universe_histogram.MakeRatio('fhc')
                    if "rhc" not in exclude: # do we care about the rhc component
                        universe_histogram.MakeRatio('rhc')

                    if not fit_muons: # do we want to keep the muon selections in addition to flavor ratios
                        universe_histogram.RemoveHistogram('fhc_muselection')
                        universe_histogram.RemoveHistogram('rhc_muselection')

                # ----- Stitch histograms together ----- #
                universe_histogram.Stitch()
                universe_histogram.SyncHistograms(cv_histogram)

                chi2 = Chi2DataMC(universe_histogram.data_hist,universe_histogram.mc_hist)
                chi2s.append(chi2)
                print("universe chi2: ",chi2)

                DataMCPlot(universe_histogram.data_hist,universe_histogram.mc_hist,cv_histogram.data_hist,cv_histogram.mc_hist)
                c1 = ROOT.TCanvas()
                MNVPLOTTER.DrawErrorSummary(universe_histogram.mc_hist,"TR",True,True,0)
                c1.Print('plots/'+universe_histogram.name+"_mc_summary.png")
                MNVPLOTTER.DrawErrorSummary(universe_histogram.data_hist,"TR",True,True,0)
                c1.Print('plots/'+universe_histogram.name+"_data_summary.png")

            print("writing {} universe for {} errorband".format(n_universes,err))
            file.write(err)
            for i in range(0,n_universes):
                file.write(", {}".format(chi2s[i]))

            file.write("\n")
    print("dnof: ",cv_histogram.data_hist.GetNbinsX()," chi2: ",Chi2DataMC(cv_histogram.data_hist,cv_histogram.mc_hist))

