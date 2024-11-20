import os
import copy
from collections import OrderedDict
import argparse
import logging, sys
import ROOT
import PlotUtils
from Tools.FitTools import *
from Tools.PlotTools import *
from Tools.StitchTools import *
import numpy as np

import math
from array import array

from config.SignalDef import SWAP_SIGNAL_DEFINATION, SIGNAL_DEFINATION
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
from tools import Utilities
from tools.PlotLibrary import HistHolder
ccnueroot = os.environ.get('CCNUEROOT')
MNVPLOTTER = PlotUtils.MnvPlotter()
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bin_width",
                        dest = "binwidth",
                        default = False,
                        action="store_true",
    )
    parser.add_argument("--thesis_plots",
                        dest = "thesis_plots",
                        default = False,
                        action="store_true",
    )
    parser.add_argument("--pseudodata",
                        dest = "pseudodata",
                        default = False,
                        action="store_true",
    )
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
    binwidthScale = args.binwidth
    pseudodata = args.pseudodata
    exclude = str(args.exclude).lower()
    doratio = args.ratio
    thesis_plots = args.thesis_plots
    fit_muons = args.fit_muons
    if thesis_plots:
        binwidthScale = True
    if doratio:
        ftag = "ratio"
        binwidthScale = False
    else:
        ftag = "stitched"

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

    type_path_map = {'data':'/exp/minerva/data/users/rhowell/nu_e/kin_dist_dataFHC_Selection_100Univ_thesis_MAD.root','mc':'/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Selection_100Univ_thesis_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    fhc_selection_mc = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/bkgfit_FHC_Selection_100Univ_N4_tune_thesis_MAD.root")
    fhc_selection_mcHold = HistHolder("Predicted MC",fhc_selection_mc,"Signal",True,mc_pot,standPOT)
    h_fhc_selection_data = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/bkgfit_FHC_Selection_100Univ_N4_tune_thesis_MAD.root").Get("EN4_data_bkgSubbed")
    fhc_sel = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_Selection_100Univ_Genie_thesis_MAD.root")
    fhc_sel_template = HistHolder("Reco Energy vs L/E",fhc_sel,"Signal",True,mc_pot,standPOT)


    type_path_map = {'mc':'/exp/minerva/data/users/rhowell/nu_e_swap/kin_dist_mcFHC_Selection_100Univ_thesis_swap_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    print("fhc_swap: ",mc_pot,data_pot)
    fhc_sel_swap = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/nu_e_swap/kin_dist_mcFHC_Selection_100Univ_thesis_swap_MAD.root")
    fhc_swap_sel_template = HistHolder("Reco Energy vs L/E",fhc_sel_swap,"Signal",True,mc_pot,standPOT)
    fhc_swap_selection_mc = HistHolder("Biased Neutrino Energy",fhc_sel_swap,"Signal",True,mc_pot,standPOT)
    
    type_path_map = {'data':'/exp/minerva/data/users/rhowell/antinu_e/kin_dist_dataRHC_Selection_1000Univ_thesis_MAD.root','mc':'/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Selection_1000Univ_thesis_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_1000Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    rhc_selection_mc = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/bkgfit_RHC_Selection_1000Univ_N4_tune_thesis_MAD.root")
    rhc_selection_mcHold = HistHolder("Predicted MC",rhc_selection_mc,"Signal",True,mc_pot,standPOT)
    h_rhc_selection_data = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/bkgfit_RHC_Selection_1000Univ_N4_tune_thesis_MAD.root").Get("EN4_data_bkgSubbed")
    rhc_sel = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Selection_Genie_thesis_MAD.root")
    rhc_sel_template = HistHolder("Reco Energy vs L/E",rhc_sel,"Signal",True,mc_pot,standPOT)

    type_path_map = {'mc':'/exp/minerva/data/users/rhowell/antinu_e_swap/kin_dist_mcRHC_Selection_1000Univ_thesis_swap_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_1000Univ",type_path_map,"MAD",True)
    print("rhc_swap: ",mc_pot,data_pot)
    rhc_sel_swap = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e_swap/kin_dist_mcRHC_Selection_1000Univ_thesis_swap_MAD.root")
    rhc_swap_sel_template = HistHolder("Reco Energy vs L/E",rhc_sel_swap,"Signal",True,mc_pot,standPOT)
    rhc_swap_selection_mc = HistHolder("Biased Neutrino Energy",rhc_sel_swap,"Signal",True,mc_pot,standPOT)

    ### NuMu Selection ###
    cates = ["CCQE","CCDelta","CCDIS","CC2p2h","CCOther","CCWrongSign"]

    type_path_map = {'data':'/exp/minerva/data/users/rhowell/nu_mu/kin_dist_dataFHC_Selection_100Univ_thesis_muon_MAD.root','mc':'/exp/minerva/data/users/rhowell/nu_mu/kin_dist_mcFHC_Selection_100Univ_thesis_muon_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("FHC_Selection_100Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    fhc_muselection_mcHold = HistHolder("Biased Neutrino Energy",mc_file,"Signal",True,mc_pot,standPOT)
    fhc_muselection_dataHold = HistHolder("Biased Neutrino Energy",data_file,"Signal",False,data_pot,standPOT)
    fhc_musel_template = HistHolder("Reco Energy vs L/E",mc_file,"Signal",True,mc_pot,standPOT)

    type_path_map = {'data':'/exp/minerva/data/users/rhowell/antinu_mu/kin_dist_dataRHC_Selection_1000Univ_thesis_muon_MAD.root','mc':'/exp/minerva/data/users/rhowell/antinu_mu/kin_dist_mcRHC_Selection_1000Univ_thesis_muon_MAD.root'}
    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale("RHC_Selection_1000Univ",type_path_map,"MAD",True)
    standPOT = data_pot if data_pot is not None else mc_pot 
    rhc_muselection_mcHold = HistHolder("Biased Neutrino Energy",mc_file,"Signal",True,mc_pot,standPOT)
    rhc_muselection_dataHold = HistHolder("Biased Neutrino Energy",data_file,"Signal",False,data_pot,standPOT)
    rhc_musel_template = HistHolder("Reco Energy vs L/E",mc_file,"Signal",True,mc_pot,standPOT)

    h2_fhc_nueel_template_nue = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_nue')
    h2_fhc_nueel_template_numu = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_numu')
    h2_fhc_nueel_template_anue = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinue')
    h2_fhc_nueel_template_anumu = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_e/kin_dist_mcFHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinumu')

    h_fhc_imd_template = ROOT.TFile('/exp/minerva/data/users/rhowell/nu_mu/FHC_imd_scattering_template.root').Get('imd_scattering_template')
    h_rhc_imd_template = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_mu/RHC_imd_scattering_template.root').Get('imd_scattering_template')

    h2_rhc_nueel_template_nue = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_nue')
    h2_rhc_nueel_template_numu = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_numu')
    h2_rhc_nueel_template_anue = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinue')
    h2_rhc_nueel_template_anumu = ROOT.TFile('/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_scattering_scatter_MAD.root').Get('nue_scattering_template_NoCuts_antinumu')

    h_fhc_swap_selection = fhc_swap_selection_mc.hists["Total"].Clone("sigTotal") 
    h_rhc_swap_selection = rhc_swap_selection_mc.hists["Total"].Clone("sigTotal") 

    h_fhc_muselection_mc = fhc_muselection_mcHold.GetHist()
    h_rhc_muselection_mc = rhc_muselection_mcHold.GetHist()
    h_fhc_muselection_data = fhc_muselection_dataHold.GetHist()
    h_rhc_muselection_data = rhc_muselection_dataHold.GetHist()
    h2_fhc_musel_template = fhc_musel_template.GetHist()
    h2_rhc_musel_template = rhc_musel_template.GetHist()

    h2_fhc_template = fhc_sel_template.hists["Total"].Clone("sigTotal")
    h2_rhc_template = rhc_sel_template.hists["Total"].Clone("sigTotal")
    h2_fhc_swap_template = fhc_swap_sel_template.hists["Total"].Clone("sigTotal")
    h2_rhc_swap_template = rhc_swap_sel_template.hists["Total"].Clone("sigTotal")

    fhc_muselection_mcHold.POTScale(binwidthScale)
    rhc_muselection_mcHold.POTScale(binwidthScale)
    fhc_muselection_dataHold.POTScale(binwidthScale)
    rhc_muselection_dataHold.POTScale(binwidthScale)
    fhc_swap_selection_mc.POTScale(binwidthScale)
    rhc_swap_selection_mc.POTScale(binwidthScale)

    h_fhc_swap_selection.Reset()
    h_rhc_swap_selection.Reset()

    h2_fhc_template.Reset()
    h2_rhc_template.Reset()
    h2_fhc_swap_template.Reset()
    h2_rhc_swap_template.Reset()

    h_fhc_muselection_mc.Reset()
    h_rhc_muselection_mc.Reset()
    h2_fhc_musel_template.Reset()
    h2_rhc_musel_template.Reset()

    for group in fhc_sel_template.hists:
        if group == "Total":
            continue
        elif group in SIGNAL_DEFINATION:
            if fhc_sel_template.hists[group]:
                h2_fhc_template.Add(fhc_sel_template.hists[group])

    for group in rhc_sel_template.hists:
        if group == "Total":
            continue
        elif group in SIGNAL_DEFINATION:
            if rhc_sel_template.hists[group]:
                h2_rhc_template.Add(rhc_sel_template.hists[group])

    for group in fhc_swap_sel_template.hists:
        if group == "Total":
            continue
        elif group in SWAP_SIGNAL_DEFINATION:
            if fhc_swap_sel_template.hists[group]:
                h_fhc_swap_selection.Add(fhc_swap_selection_mc.hists[group])
                h2_fhc_swap_template.Add(fhc_swap_sel_template.hists[group])

    for group in rhc_swap_sel_template.hists:
        if group == "Total":
            continue
        elif group in SWAP_SIGNAL_DEFINATION:
            if rhc_swap_sel_template.hists[group]:
                h_rhc_swap_selection.Add(rhc_swap_selection_mc.hists[group])
                h2_rhc_swap_template.Add(rhc_swap_sel_template.hists[group])

    for group in fhc_muselection_mcHold.hists:
        if group == "Total":
            continue
        elif group in cates:
            h_fhc_muselection_mc.Add(fhc_muselection_mcHold.hists[group])
            h2_fhc_musel_template.Add(fhc_musel_template.hists[group])

    for group in fhc_muselection_mcHold.hists:
        if group == "Total":
            continue
        elif group in cates:
            h_rhc_muselection_mc.Add(rhc_muselection_mcHold.hists[group])
            h2_rhc_musel_template.Add(rhc_musel_template.hists[group])

    h_fhc_selection_mc = fhc_selection_mcHold.GetHist()
    h_rhc_selection_mc = rhc_selection_mcHold.GetHist()
    
    h_fhc_selection_mc.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    h_fhc_swap_selection.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    sys_p = h_fhc_selection_mc.GetCVHistoWithError() * fhc_scale_p1sig
    sys_m = h_fhc_selection_mc.GetCVHistoWithError() * fhc_scale_m1sig
    swap_p = h_fhc_swap_selection.GetCVHistoWithError() * fhc_scale_p1sig
    swap_m = h_fhc_swap_selection.GetCVHistoWithError() * fhc_scale_m1sig
    for i in range(h_fhc_selection_mc.GetNbinsX()+1):
        h_fhc_selection_mc.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,sys_p.GetBinContent(i))
        h_fhc_selection_mc.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,sys_m.GetBinContent(i))

        h_fhc_swap_selection.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,swap_p.GetBinContent(i))
        h_fhc_swap_selection.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,swap_m.GetBinContent(i))

    h_rhc_selection_mc.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    h_rhc_swap_selection.AddVertErrorBandAndFillWithCV("ElectronScale", 2)
    sys_p = h_rhc_selection_mc.GetCVHistoWithError() * rhc_scale_p1sig
    sys_m = h_rhc_selection_mc.GetCVHistoWithError() * rhc_scale_m1sig
    swap_p = h_rhc_swap_selection.GetCVHistoWithError() * rhc_scale_p1sig
    swap_m = h_rhc_swap_selection.GetCVHistoWithError() * rhc_scale_m1sig
    for i in range(h_rhc_selection_mc.GetNbinsX()+1):
        h_rhc_selection_mc.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,sys_p.GetBinContent(i))
        h_rhc_selection_mc.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,sys_m.GetBinContent(i))

        h_rhc_swap_selection.GetVertErrorBand("ElectronScale").GetHist(0).SetBinContent(i,swap_p.GetBinContent(i))
        h_rhc_swap_selection.GetVertErrorBand("ElectronScale").GetHist(1).SetBinContent(i,swap_m.GetBinContent(i))

    #get FHC electron energy hist
    h_fhc_nueel_mc = ROOT.TFile('/exp/minerva/data/users/rhowell/nueel/DeepikaMc_LauraData_Mnv.root').Get('mnv_eff_cor_mc')
    h_fhc_nueel_data = ROOT.TFile('/exp/minerva/data/users/rhowell/nueel/DeepikaMc_LauraData_Mnv.root').Get('mnv_eff_cor_data')

    #get RHC electron energy hist
    h_rhc_nueel_mc = ROOT.TFile("/exp/minerva/data/users/rhowell/nueel/electronE_bkgsub_FullSetFinalv2_nobinnorm.root").Get('mc_fluxonly')   #.Get('mc_bkgsub_effcor')
    h_rhc_nueel_data =  ROOT.TFile("/exp/minerva/data/users/rhowell/nueel/electronE_bkgsub_FullSetFinalv2_nobinnorm.root").Get('data_bkgsub_effcor')

    #get IMD histograms
    h_fhc_imd_mc = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("FHC_MC")
    h_fhc_imd_data = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("FHC_Data")

    h_rhc_imd_mc = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("RHC_MC")
    h_rhc_imd_data = ROOT.TFile.Open("IMD/PubResult_v18_Combined.root").Get("RHC_Data")

    fhcnueelnue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_nue')
    fhcnueelnumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_numu')
    fhcnueelanue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_anue')
    fhcnueelanumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_FHC_withoutconstraint.root").Get('electron_energy_anumu')
    rhcnueelnue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_nue')
    rhcnueelnumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_numu')
    rhcnueelanue = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_anue')
    rhcnueelanumu = ROOT.TFile.Open("NuEnumberEvents/ElectronEnergySpectrum_RHC_withoutconstraint.root").Get('electron_energy_anumu')

    #Lateral errorbands are deprecated, convert to vertical in all histograms that have them
    LateralToVertical([h_fhc_nueel_data,h_rhc_nueel_data,h_fhc_imd_mc,h_rhc_imd_mc,h_fhc_imd_data,h_rhc_imd_data])
    #Replace beam angle systematic with 4 shifts for x and y, +-1sigma, with two errorbands for each position
    SeparateBeamAngle([h_fhc_selection_mc,h_fhc_muselection_mc,h_fhc_selection_data,
        h_rhc_selection_mc,h_rhc_muselection_mc,h_rhc_selection_data])

    f1 = fhcnueelnue.Clone()
    f2 = fhcnueelanue.Clone()
    f3 = fhcnueelnumu.Clone()
    f4 = fhcnueelanumu.Clone()
    f1.Add(f2)
    f1.Add(f3)
    f1.Add(f4)
    r1 = rhcnueelnue.Clone()
    r2 = rhcnueelanue.Clone()
    r3 = rhcnueelnumu.Clone()
    r4 = rhcnueelanumu.Clone()
    r1.Add(r2)
    r1.Add(r3)
    r1.Add(r4)
    fhcweight = f1.GetCVHistoWithError()
    rhcweight = r1.GetCVHistoWithError()
    for i in range(0,f1.GetNbinsX()+1):
        f_ratio = f1.GetBinContent(i)/h_fhc_nueel_mc.GetBinContent(i) if h_fhc_nueel_mc.GetBinContent(i) != 0 else 1
        r_ratio = r1.GetBinContent(i)/h_rhc_nueel_mc.GetBinContent(i) if h_rhc_nueel_mc.GetBinContent(i) != 0 else 1
        fhcweight.SetBinContent(i,f_ratio)
        fhcweight.SetBinError(i,0)
        rhcweight.SetBinContent(i,r_ratio)
        rhcweight.SetBinError(i,0)
    fhcnueelnumu.DivideSingle(fhcnueelnumu,fhcweight)
    fhcnueelanumu.DivideSingle(fhcnueelanumu,fhcweight)
    fhcnueelanue.DivideSingle(fhcnueelanue,fhcweight)
    fhcnueelnue.DivideSingle(fhcnueelnue,fhcweight)
    rhcnueelnumu.DivideSingle(rhcnueelnumu,rhcweight)
    rhcnueelanumu.DivideSingle(rhcnueelanumu,rhcweight)
    rhcnueelanue.DivideSingle(rhcnueelanue,rhcweight)
    rhcnueelnue.DivideSingle(rhcnueelnue,rhcweight)
    r1 = rhcnueelnue.Clone()
    r2 = rhcnueelanue.Clone()
    r3 = rhcnueelnumu.Clone()
    r4 = rhcnueelanumu.Clone()
    r1.Add(r2)
    r1.Add(r3)
    r1.Add(r4)
    
    #Make all hist have the same errorbands. Filled with CV if they don't have it.
    SyncErrorBandsv2([h_fhc_imd_mc,h_rhc_imd_mc,h_fhc_nueel_mc,h_rhc_nueel_mc,h_fhc_selection_mc,
        h_rhc_selection_mc,h_fhc_muselection_mc,h_rhc_muselection_mc,h_fhc_imd_data,h_rhc_imd_data,
        h_fhc_nueel_data,h_rhc_nueel_data,h_fhc_selection_data,h_rhc_selection_data,
        h_fhc_muselection_data,h_rhc_muselection_data,fhcnueelanumu,fhcnueelanue,fhcnueelnumu,
        fhcnueelnue,rhcnueelanumu,rhcnueelanue,rhcnueelnumu,rhcnueelnue,h_fhc_swap_selection,h_rhc_swap_selection])

    mc_histsToStitch = OrderedDict()
    mc_histsToStitch['fhc_nueel'] = h_fhc_nueel_mc
    mc_histsToStitch['fhc_imd'] = h_fhc_imd_mc
    mc_histsToStitch['fhc_muselection']=h_fhc_muselection_mc
    mc_histsToStitch['fhc_selection']=h_fhc_selection_mc
    mc_histsToStitch['rhc_nueel'] = h_rhc_nueel_mc
    mc_histsToStitch['rhc_imd'] = h_rhc_imd_mc
    mc_histsToStitch['rhc_muselection']=h_rhc_muselection_mc
    mc_histsToStitch['rhc_selection']=h_rhc_selection_mc

    nue_templates = OrderedDict()
    h2_fhc_nueel_template_nue.Add(h2_fhc_nueel_template_anue)
    nue_templates['fhc_nueel'] = h2_fhc_nueel_template_nue.Clone()
    nue_templates['fhc_imd'] = h_fhc_imd_template.Clone() # placeholder
    nue_templates['fhc_imd'].Reset()
    nue_templates['fhc_muselection']=h2_fhc_musel_template.Clone() # placeholder
    nue_templates['fhc_muselection'].Reset()
    nue_templates['fhc_selection']=h2_fhc_template.Clone()
    h2_rhc_nueel_template_anue.Add(h2_rhc_nueel_template_nue)
    nue_templates['rhc_nueel'] = h2_rhc_nueel_template_anue.Clone()
    nue_templates['rhc_imd'] = h_rhc_imd_template.Clone() #  placeholder
    nue_templates['rhc_imd'].Reset()
    nue_templates['rhc_muselection']=h2_rhc_musel_template.Clone()# placeholder
    nue_templates['rhc_muselection'].Reset()
    nue_templates['rhc_selection']=h2_rhc_template.Clone()
    
    numu_templates = OrderedDict()
    h2_fhc_nueel_template_numu.Add(h2_fhc_nueel_template_anumu)
    numu_templates['fhc_nueel'] = h2_fhc_nueel_template_numu.Clone()
    numu_templates['fhc_imd'] = h_fhc_imd_template.Clone() # placeholder
    numu_templates['fhc_muselection']=h2_fhc_musel_template.Clone()
    numu_templates['fhc_selection']=h2_fhc_swap_template.Clone()
    numu_templates['fhc_selection'].Reset()
    h2_rhc_nueel_template_anumu.Add(h2_rhc_nueel_template_numu)
    numu_templates['rhc_nueel'] = h2_rhc_nueel_template_anumu.Clone()
    numu_templates['rhc_imd'] = h_rhc_imd_template.Clone() # placeholder
    numu_templates['rhc_muselection']=h2_rhc_musel_template.Clone()
    numu_templates['rhc_selection']=h2_rhc_swap_template.Clone()
    numu_templates['rhc_selection'].Reset()

    swap_templates = OrderedDict()
    h2_fhc_nueel_template_numu.Add(h2_fhc_nueel_template_anumu)
    swap_templates['fhc_nueel'] = h2_fhc_nueel_template_numu.Clone()
    swap_templates['fhc_imd'] = h_fhc_imd_template.Clone() # placeholder
    swap_templates['fhc_muselection']=h2_fhc_musel_template.Clone()
    swap_templates['fhc_selection']=h2_fhc_swap_template.Clone()
    h2_rhc_nueel_template_anumu.Add(h2_rhc_nueel_template_numu)
    swap_templates['rhc_nueel'] = h2_rhc_nueel_template_anumu.Clone()
    swap_templates['rhc_imd'] = h_rhc_imd_template.Clone() # placeholder
    swap_templates['rhc_muselection']=h2_rhc_musel_template.Clone()
    swap_templates['rhc_selection']=h2_rhc_swap_template.Clone()
    for h in swap_templates.keys():
        if "_selection" not in h and "nueel" not in h:
            swap_templates[h].Reset()

    data_histsToStitch = OrderedDict()
    data_histsToStitch['fhc_nueel'] = h_fhc_nueel_data
    data_histsToStitch['fhc_imd'] = h_fhc_imd_data
    data_histsToStitch['fhc_muselection']=h_fhc_muselection_data
    data_histsToStitch['fhc_selection']=h_fhc_selection_data
    data_histsToStitch['rhc_nueel'] = h_rhc_nueel_data
    data_histsToStitch['rhc_imd'] = h_rhc_imd_data
    data_histsToStitch['rhc_muselection']=h_rhc_muselection_data
    data_histsToStitch['rhc_selection']=h_rhc_selection_data
    
    for h in list(mc_histsToStitch.keys()):
        try:
            if "fhc" in exclude and "selection" in h and "fhc" in h:
                print('deleting '+h)
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
            if "rhc" in exclude and "selection" in h and "rhc" in h:
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
            if "numu" in exclude and "muselection" in h:
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
            if "nue" in exclude and "_selection" in h:
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
            if "elastic" in exclude and "nueel" in h:
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
            if "imd" in exclude and "imd" in h:
                del(mc_histsToStitch[h])
                del(nue_templates[h])
                del(numu_templates[h])
                del(swap_templates[h])
                del(data_histsToStitch[h])
        except:
            continue

    if doratio:
        for h in list(mc_histsToStitch.keys()):
            if not fit_muons:
                if "selection" in h:
                    del(mc_histsToStitch[h])
                    del(nue_templates[h])
                    del(numu_templates[h])
                    del(swap_templates[h])
                    del(data_histsToStitch[h])
            else:
                if "_selection" in h:
                    del(mc_histsToStitch[h])
                    del(nue_templates[h])
                    del(numu_templates[h])
                    del(swap_templates[h])
                    del(data_histsToStitch[h])

        if "fhc" not in exclude:
            toClone = h_fhc_muselection_mc.Clone()
            toClone.Divide(toClone,h_fhc_selection_mc)
            mc_histsToStitch['fhc_ratio']=toClone
            toClone = h_fhc_muselection_data.Clone()
            toClone.Divide(toClone,h_fhc_selection_data)
            data_histsToStitch['fhc_ratio']=toClone
            nue_templates['fhc_ratio'] = h2_fhc_template.Clone()
            numu_templates['fhc_ratio']=h2_fhc_musel_template.Clone()
            swap_templates['fhc_ratio']=h2_fhc_swap_template.Clone()
        if "rhc" not in exclude:
            toClone = h_rhc_muselection_mc.Clone()
            toClone.Divide(toClone,h_rhc_selection_mc)
            mc_histsToStitch['rhc_ratio']=toClone
            toClone = h_rhc_muselection_data.Clone()
            toClone.Divide(toClone,h_rhc_selection_data)
            data_histsToStitch['rhc_ratio']=toClone
            nue_templates['rhc_ratio'] = h2_rhc_template.Clone()
            numu_templates['rhc_ratio']=h2_rhc_musel_template.Clone()
            swap_templates['rhc_ratio']=h2_rhc_swap_template.Clone()

    nue_histsToStitch           = copy.deepcopy(mc_histsToStitch)
    numu_histsToStitch          = copy.deepcopy(mc_histsToStitch)
    nutau_histsToStitch         = copy.deepcopy(mc_histsToStitch)
    fhc_histsToStitch           = copy.deepcopy(mc_histsToStitch)
    nue_selection_histsToStitch = copy.deepcopy(mc_histsToStitch)
    ratio_histsToStitch         = copy.deepcopy(mc_histsToStitch)
    swap_histsToStitch          = copy.deepcopy(mc_histsToStitch)

    for h in mc_histsToStitch.keys():
        if 'fhc' in h:
            UnityHist(fhc_histsToStitch[h])
        elif 'rhc' in h:
            EmptyHist(fhc_histsToStitch[h])

        if 'nueel' in h:
            EmptyHist(nue_selection_histsToStitch[h])
            EmptyHist(ratio_histsToStitch[h])
            UnityHist(nutau_histsToStitch[h])
            if 'fhc' in h:
                nue_histsToStitch[h] = fhcnueelnue.Clone()
                nue_histsToStitch[h].Add(fhcnueelanue)
                numu_histsToStitch[h] = fhcnueelnumu.Clone()
                numu_histsToStitch[h].Add(fhcnueelanumu)
                swap_histsToStitch[h] = fhcnueelnumu.Clone()
                swap_histsToStitch[h].Add(fhcnueelanumu)
            elif 'rhc' in h:
                nue_histsToStitch[h] = rhcnueelnue.Clone()
                nue_histsToStitch[h].Add(rhcnueelanue)
                numu_histsToStitch[h] = rhcnueelnumu.Clone()
                numu_histsToStitch[h].Add(rhcnueelanumu)
                swap_histsToStitch[h] = rhcnueelnumu.Clone()
                swap_histsToStitch[h].Add(rhcnueelanumu)
        elif 'imd' in h:
            EmptyHist(nue_selection_histsToStitch[h])
            EmptyHist(nue_histsToStitch[h])
            EmptyHist(swap_histsToStitch[h])
            EmptyHist(ratio_histsToStitch[h])
            EmptyHist(nutau_histsToStitch[h])
        elif 'muselection' in h:
            EmptyHist(nue_selection_histsToStitch[h])
            EmptyHist(nue_histsToStitch[h])
            EmptyHist(swap_histsToStitch[h])
            EmptyHist(ratio_histsToStitch[h])
            EmptyHist(nutau_histsToStitch[h])
        elif '_selection' in h:
            UnityHist(nue_selection_histsToStitch[h])
            EmptyHist(numu_histsToStitch[h])
            EmptyHist(nutau_histsToStitch[h])
            EmptyHist(ratio_histsToStitch[h])
            if 'fhc' in h:
                swap_histsToStitch[h] = h_fhc_swap_selection
            elif 'rhc' in h:
                swap_histsToStitch[h] = h_rhc_swap_selection
        elif 'ratio' in h:
            EmptyHist(nue_selection_histsToStitch[h])
            EmptyHist(nutau_histsToStitch[h])
            UnityHist(ratio_histsToStitch[h])
            if 'fhc' in h:
                swap_histsToStitch[h] = h_fhc_swap_selection
                numu_histsToStitch[h] = h_fhc_muselection_mc
                nue_histsToStitch[h] = h_fhc_selection_mc
            elif 'rhc' in h:
                swap_histsToStitch[h] = h_rhc_swap_selection
                numu_histsToStitch[h] = h_rhc_muselection_mc
                nue_histsToStitch[h] = h_rhc_selection_mc

    n_bins_new = 0
    for h in mc_histsToStitch:
        for i in range(0,mc_histsToStitch[h].GetNbinsX()+1):
            if mc_histsToStitch[h].GetBinContent(i) > minBinCont:
                n_bins_new+=1

    htemp_mc = ROOT.TH1D('mc_stitched', 'mc_stitched', n_bins_new, 0,  n_bins_new )
    htemp_mc_nue = ROOT.TH1D('mc_stitched_nue', 'mc_stitched_nue', n_bins_new, 0,  n_bins_new )
    htemp_mc_numu = ROOT.TH1D('mc_stitched_numu', 'mc_stitched_numu', n_bins_new, 0,  n_bins_new )
    htemp_mc_nutau = ROOT.TH1D('mc_stitched_nutau', 'mc_stitched_nutau', n_bins_new, 0,  n_bins_new )
    htemp_mc_fhc = ROOT.TH1D('mc_stitched_fhc', 'mc_stitched_fhc', n_bins_new, 0,  n_bins_new )
    htemp_mc_nueselection = ROOT.TH1D('mc_stitched_nueselection', 'mc_stitched_nueselection', n_bins_new, 0,  n_bins_new )
    htemp_mc_ratio = ROOT.TH1D('mc_stitched_ratio', 'mc_stitched_ratio', n_bins_new, 0,  n_bins_new )
    htemp_mc_swap = ROOT.TH1D('mc_stitched_swap', 'mc_stitched_swap', n_bins_new, 0,  n_bins_new )
    htemp_data = ROOT.TH1D('data_stitched', 'data_stitched', n_bins_new, 0,  n_bins_new )

    StitchThis(htemp_mc, mc_histsToStitch,mc_histsToStitch)
    StitchThis(htemp_mc_nue, nue_histsToStitch,mc_histsToStitch)
    StitchThis(htemp_mc_numu, numu_histsToStitch,mc_histsToStitch)
    StitchThis(htemp_mc_nutau, nutau_histsToStitch,mc_histsToStitch)
    StitchThis(htemp_mc_fhc, fhc_histsToStitch,mc_histsToStitch)
    StitchThis(htemp_mc_nueselection, nue_selection_histsToStitch,mc_histsToStitch)
    StitchThis(htemp_mc_ratio, ratio_histsToStitch,mc_histsToStitch)
    StitchThis(htemp_mc_swap, swap_histsToStitch,mc_histsToStitch)
    StitchThis(htemp_data, data_histsToStitch,mc_histsToStitch,True,pseudodata)

    h2template_nue = ROOT.TH2D('LE_template_nue', 'LE_template_nue', n_bins_new, 0,  n_bins_new,34,0,0.495)
    h2template_numu = ROOT.TH2D('LE_template_numu', 'LE_template_numu', n_bins_new, 0,  n_bins_new,34,0,0.495)
    h2template_swap = ROOT.TH2D('LE_template_swap', 'LE_template_swap', n_bins_new, 0,  n_bins_new,34,0,0.495)

    StitchThis2D(h2template_nue, nue_templates,mc_histsToStitch,"nue")
    StitchThis2D(h2template_numu, numu_templates,mc_histsToStitch,"numu")
    StitchThis2D(h2template_swap, swap_templates,mc_histsToStitch,"swap")

    mnv_data = PlotUtils.MnvH1D( htemp_data )
    mnv_data.GetXaxis().SetTitle("Bin number")
    mnv_mc   = PlotUtils.MnvH1D( htemp_mc )
    mnv_mc_nue  = PlotUtils.MnvH1D( htemp_mc_nue )
    mnv_mc_numu = PlotUtils.MnvH1D( htemp_mc_numu )
    mnv_mc_swap = PlotUtils.MnvH1D( htemp_mc_swap )

    c1 = ROOT.TCanvas()
    mnv_mc_nue.Draw("hist")
    c1.Print("plots/nue_hist.png")
    mnv_mc_numu.Draw("hist")
    c1.Print("plots/numu_hist.png")
    mnv_mc_swap.Draw("hist")
    c1.Print("plots/nuswap_hist.png")

    htemp_mc_nutau.Draw("hist")
    c1.Print("plots/nutau_id.png")
    htemp_mc_nueselection.Draw("hist")
    c1.Print("plots/nusel_id.png")
    htemp_mc_ratio.Draw("hist")
    c1.Print("plots/nuratio_id.png")
    htemp_mc_fhc.Draw("hist")
    c1.Print("plots/fhc_id.png")

    h2template_nue.Draw("colz")
    c1.Print("plots/nue_temp.png")
    h2template_numu.Draw("colz")
    c1.Print("plots/numu_temp.png")
    h2template_swap.Draw("colz")
    c1.Print("plots/nuswap_temp.png")

    FillErrorBandfromHist2( mnv_data, data_histsToStitch, mc_histsToStitch, True, pseudodata)
    FillErrorBandfromHist2( mnv_mc, mc_histsToStitch, mc_histsToStitch, False)

    FillErrorBandfromHist2(mnv_mc_nue,  nue_histsToStitch, mc_histsToStitch)
    FillErrorBandfromHist2(mnv_mc_numu, numu_histsToStitch, mc_histsToStitch)
    FillErrorBandfromHist2(mnv_mc_swap, swap_histsToStitch, mc_histsToStitch)

    # ----- remove errorbands from MnvH1Ds and convert systematics to errorbands to save time in chi2 calculation ----- #
    #ConsolidateErrorMatrices(mnv_data)
    #ConsolidateErrorMatrices(mnv_mc)
    #ConsolidateErrorMatrices(mnv_mc_nue)
    #ConsolidateErrorMatrices(mnv_mc_numu)
    #ConsolidateErrorMatrices(mnv_mc_swap)

    mnv_mc.SetLineColor(ROOT.kRed)
    mc = mnv_mc.GetCVHistoWithError()
    data = mnv_data.GetCVHistoWithError()

    mc.GetXaxis().SetTitle("Bin number")
    mc.GetYaxis().SetTitle("Entries")
    MNVPLOTTER.ApplyAxisStyle(mc,True,True)

    canvas = ROOT.TCanvas()
    MNVPLOTTER.DrawErrorSummary(mnv_mc,"TR",True,True)
    canvas.Print("plots/stitched_errors_mc_Total.png")
    MNVPLOTTER.DrawErrorSummary(mnv_data,"TR",True,True)
    canvas.Print("plots/stitched_errors_data_Total.png")
    MNVPLOTTER.DrawErrorSummary(mnv_mc_nue,"TR",True,True)
    canvas.Print("plots/stitched_errors_nue_Total.png")
    MNVPLOTTER.DrawErrorSummary(mnv_mc_numu,"TR",True,True)
    canvas.Print("plots/stitched_errors_numu_Total.png")
    MNVPLOTTER.DrawErrorSummary(mnv_mc_swap,"TR",True,True)
    canvas.Print("plots/stitched_errors_swap_Total.png")

    dataprint = []
    mcprint = []
    for i in range(1,mnv_data.GetNbinsX()+1):
        dataprint.append(mnv_data.GetBinContent(i))
    for i in range(1,mnv_mc.GetNbinsX()+1):
        mcprint.append(mnv_mc.GetBinContent(i))
    dataprint = np.array(dataprint)
    mcprint = np.array(mcprint)
    
    np.savetxt("mc_"+ftag+".csv",mcprint,delimiter=",")
    np.savetxt("data_"+ftag+".csv",dataprint,delimiter=",")
                
    for name in mnv_mc.GetVertErrorBandNames():
        n_univ = mnv_mc.GetVertErrorBand(name).GetNHists()
        err_hists = []
        for univ in range(n_univ):
            hist_vals = []
            h_univ = mnv_mc.GetVertErrorBand(name).GetHist(univ)
            for i in range(h_univ.GetNbinsX()+1):
                hist_vals.append(h_univ.GetBinContent(i))
            err_hists.append(hist_vals)
        np.savetxt("errorband_hists_{}.csv".format(name),np.array(err_hists),delimiter=',')
            

    GetCovarianceMatrix(mnv_mc,mnv_data,ftag)

    if pseudodata:
        f = ROOT.TFile("{}/FeldmanCousins/NuE_stitched_hists_pseudo.root".format(ccnueroot),"RECREATE")
    else:
        f = ROOT.TFile("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot),"RECREATE")

    mnv_data.Write()
    mnv_mc.Write()
    mnv_mc_nue.Write()
    mnv_mc_numu.Write()
    mnv_mc_swap.Write()
    htemp_mc_nutau.Write()
    htemp_mc_fhc.Write()
    htemp_mc_nueselection.Write()
    htemp_mc_ratio.Write()
    h2template_numu.Write()
    h2template_nue.Write()
    h2template_swap.Write()

    DataMCCVPlot(mnv_data,mnv_mc)
    print(Chi2DataMC(mnv_data,mnv_mc))
    f.Close()
