import os
from collections import OrderedDict
import argparse
import logging, sys
import ROOT
import PlotUtils
from fit_tools.FitTools import *
from fit_tools.PlotTools import *
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
from config.SignalDef import SWAP_SIGNAL_DEFINATION, SIGNAL_DEFINATION
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

#MNVPLOTTER.error_summary_group_map.clear();
#for k,v in CONSOLIDATED_ERROR_GROUPS.items():
#    vec = ROOT.vector("std::string")()
#    for vs in v :
#        vec.push_back(vs)
#    MNVPLOTTER.error_summary_group_map[k]= vec

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
minBinCont = 1

def SyncErrorBandsv2(hists):
    for _h1 in range(len(hists)):
        for _h2 in range(len(hists)):
            h1 = hists[_h1]
            h2 = hists[_h2]
            if h1 == h2:
                continue
            for name in h1.GetVertErrorBandNames():
                if name == "Flux":
                    n_universes = 100
                elif name in ["Muon_Energy_Resolution","Muon_Energy_MINOS","Muon_Energy_MINERvA"]:
                    continue
                elif(name not in h2.GetVertErrorBandNames()):
                    n_universes = h1.GetVertErrorBand(name).GetNHists()
                    h2.AddVertErrorBandAndFillWithCV(name, n_universes)

            for name in h1.GetLatErrorBandNames():
                if name == "Flux":
                    n_universes = 100
                elif name in ["Muon_Energy_Resolution","Muon_Energy_MINOS","Muon_Energy_MINERvA"]:
                    continue
                elif(name not in h2.GetLatErrorBandNames()):
                    n_universes = h1.GetLatErrorBand(name).GetNHists()
                    h2.AddLatErrorBandAndFillWithCV(name, n_universes)

            for name in h2.GetVertErrorBandNames():
                if name == "Flux":
                    n_universes = 100
                elif name in ["Muon_Energy_Resolution","Muon_Energy_MINOS","Muon_Energy_MINERvA"]:
                    continue
                elif(name not in h1.GetVertErrorBandNames()):
                    n_universes = h2.GetVertErrorBand(name).GetNHists()
                    h1.AddVertErrorBandAndFillWithCV(name, n_universes)

            for name in h2.GetLatErrorBandNames():
                if name == "Flux":
                    n_universes = 100
                elif name in ["Muon_Energy_Resolution","Muon_Energy_MINOS","Muon_Energy_MINERvA"]:
                    continue
                elif(name not in h1.GetLatErrorBandNames()):
                    n_universes = h2.GetLatErrorBand(name).GetNHists()
                    h1.AddLatErrorBandAndFillWithCV(name, n_universes)

def FillErrorBandfromHist2(h_new,h_olds,mchists = None, isData=False):
    offset = 1
    for h in h_olds:
        h_old = h_olds[h]
        Nbins = h_old.GetNbinsX()+1 if not mchists else mchists[h].GetNbinsX()+1

        errorband_names_vert = h_old.GetVertErrorBandNames()
        errorband_names_lat = h_old.GetLatErrorBandNames()
        n_univ = 0
        sys_bc = 0.0

        for error_band in errorband_names_vert:
            #print("looping over vert error ", error_band)
            if error_band == 'Flux' and isData:# and h_old.GetVertErrorBand(error_band).GetNHists() < 1000:            #hacked for the flux since fhc files have 100 and rhc have 1000 n_universes
                n_univ = 100
            elif error_band == 'Flux':
                n_univ = 100
            else:
                n_univ = h_old.GetVertErrorBand(error_band).GetNHists()
                if error_band == "GENIE_MaZExpCCQE":
                    n_univ = 100
                elif error_band == "GENIE_MaCCQE":
                    n_univ = 2

            if error_band == "Muon_Energy_MINERvA" or error_band == "Muon_Energy_MINOS" or error_band == "Muon_Energy_Resolution":
                new_band = "v_{}".format(error_band)
            else:
                new_band = error_band

            if not h_new.HasVertErrorBand( new_band ) and h_old.HasVertErrorBand( error_band ):
                h_new.AddVertErrorBandAndFillWithCV( new_band, n_univ )
                h_new.GetVertErrorBand(new_band).SetUseSpreadError( h_old.GetVertErrorBand(error_band).GetUseSpreadError())

            for universe in range(0, n_univ):
                bin_offset = offset
                for b in range(1, Nbins):
                    bin_c = h_old.GetBinContent(b)
                    mcbin_c = mchists[h].GetBinContent(b)
                    sys_bc = h_old.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)
                    ratio = sys_bc/bin_c if bin_c != 0 else 0
                    if mcbin_c <= minBinCont:
                        bin_offset += -1
                        continue
                    bin_new = b + bin_offset - 1
                    h_new.GetVertErrorBand(new_band).GetHist(universe).SetBinContent( bin_new, sys_bc )

        for error_band in errorband_names_lat:
            if error_band == "Muon_Energy_MINERvA" or error_band == "Muon_Energy_MINOS" or error_band == "Muon_Energy_Resolution":
                new_band = "h_{}".format(error_band)
            else:
                new_band = error_band
            #print("looping over lat error ", error_band)
            if not h_new.HasLatErrorBand(new_band) and h_old.HasLatErrorBand(error_band):
                n_univ = h_old.GetLatErrorBand(error_band).GetNHists()
            else:
                continue
            if not h_new.HasLatErrorBand( new_band ):
                h_new.AddLatErrorBandAndFillWithCV( new_band , n_univ )
                h_new.GetLatErrorBand(new_band).SetUseSpreadError( h_old.GetLatErrorBand(error_band).GetUseSpreadError())
            for universe in range(0, n_univ):
                bin_offset = offset
                for b in range(1, Nbins):  #bin 0 and 1 have no entries. Bin seven is the last bin
                    bin_c = h_old.GetBinContent(b)
                    mcbin_c = mchists[h].GetBinContent(b)
                    sys_bc = h_old.GetLatErrorBand(error_band).GetHist(universe).GetBinContent(b)
                    ratio = sys_bc/bin_c if bin_c !=0 else 0
                    if mcbin_c <= minBinCont:
                        bin_offset += -1
                        continue
                    bin_new = b  + bin_offset - 1
                    h_new.GetLatErrorBand(new_band).GetHist(universe).SetBinContent( bin_new, sys_bc )
        for i in range(1,Nbins):
            if mchists[h].GetBinContent(i) <= minBinCont:
                continue
            offset+=1

def StitchThis(h_new,h_olds,sample="",mchists=None):
    i_new = 0
    for h in h_olds:
        h_old = h_olds[h]
        for i in range(1,h_old.GetNbinsX()+1):
            bin_c = h_old.GetBinContent(i)
            
            if mchists[h].GetBinContent(i) <= minBinCont:
                continue
            i_new += 1

            h_new.SetBinContent(i_new,bin_c)
            if sample=="data":
                continue

            # no statistical error on elastic scattering special production
            if 'nueel' in h:
                h_new.SetBinError(i_new,0)

def SwapUniverseCV(h_in,err,univ):
    try:
        ratio = h_in.GetVertErrorBand(err).GetHist(univ).Clone()
    except:
        print("Histogram doesn't have systematic errorband {}, skipping.".format(err))
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

    args = parser.parse_args()
    exclude = str(args.exclude).lower()
    doratio = args.ratio
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
        
        err_names = h_rhc_selection_mc.GetVertErrorBandNames()

        #Make all hist have the same errorbands. Filled with CV if they don't have it.
        SyncErrorBandsv2([h_fhc_imd_mc,h_rhc_imd_mc,h_fhc_nueel_mc,h_rhc_nueel_mc,h_fhc_selection_mc,h_rhc_selection_mc,h_fhc_muselection_mc,h_rhc_muselection_mc])
        SyncErrorBandsv2([h_fhc_imd_data,h_rhc_imd_data,h_fhc_nueel_data,h_rhc_nueel_data,h_fhc_selection_data,h_rhc_selection_data,h_fhc_muselection_data,h_rhc_muselection_data])
        
        mc_cvToStitch = OrderedDict()
        mc_cvToStitch['fhc_nueel'] = h_fhc_nueel_mc
        mc_cvToStitch['fhc_imd'] = h_fhc_imd_mc
        mc_cvToStitch['fhc_muselection']=h_fhc_muselection_mc
        mc_cvToStitch['fhc_selection']=h_fhc_selection_mc
        mc_cvToStitch['rhc_nueel'] = h_rhc_nueel_mc
        mc_cvToStitch['rhc_imd'] = h_rhc_imd_mc
        mc_cvToStitch['rhc_muselection']=h_rhc_muselection_mc
        mc_cvToStitch['rhc_selection']=h_rhc_selection_mc

        data_cvToStitch = OrderedDict()
        data_cvToStitch['fhc_nueel'] = h_fhc_nueel_data
        data_cvToStitch['fhc_imd'] = h_fhc_imd_data
        data_cvToStitch['fhc_muselection']=h_fhc_muselection_data
        data_cvToStitch['fhc_selection']=h_fhc_selection_data
        data_cvToStitch['rhc_nueel'] = h_rhc_nueel_data
        data_cvToStitch['rhc_imd'] = h_rhc_imd_data
        data_cvToStitch['rhc_muselection']=h_rhc_muselection_data
        data_cvToStitch['rhc_selection']=h_rhc_selection_data

        for h in list(mc_cvToStitch.keys()):
            if "fhc" in exclude and "selection" in h and "fhc" in h:
                del(mc_cvToStitch[h])
                del(data_cvToStitch[h])
                continue
            if "rhc" in exclude and "selection" in h and "rhc" in h:
                del(mc_cvToStitch[h])
                del(data_cvToStitch[h])
                continue
            if "numu" in exclude and "muselection" in h:
                del(mc_cvToStitch[h])
                del(data_cvToStitch[h])
                continue
            if "nue" in exclude and "_selection" in h:
                del(mc_cvToStitch[h])
                del(data_cvToStitch[h])
                continue
            if "elastic" in exclude and "nueel" in h:
                del(mc_cvToStitch[h])
                del(data_cvToStitch[h])
                continue
            if "imd" in exclude and "imd" in h:
                del(mc_cvToStitch[h])
                del(data_cvToStitch[h])
                continue

        if doratio:
            for h in list(mc_cvToStitch.keys()):
                if "selection" in h:
                    del(mc_cvToStitch[h])
                    del(data_cvToStitch[h])
            if "fhc" not in exclude:
                toClone = h_fhc_muselection_mc.Clone()
                toClone.Divide(toClone,h_fhc_selection_mc)
                mc_cvToStitch['fhc_ratio']=toClone
                toClone = h_fhc_muselection_data.Clone()
                toClone.Divide(toClone,h_fhc_selection_data)
                data_cvToStitch['fhc_ratio']=toClone
            if "rhc" not in exclude:
                toClone = h_rhc_muselection_mc.Clone()
                toClone.Divide(toClone,h_rhc_selection_mc)
                mc_cvToStitch['rhc_ratio']=toClone
                toClone = h_rhc_muselection_data.Clone()
                toClone.Divide(toClone,h_rhc_selection_data)
                data_cvToStitch['rhc_ratio']=toClone

        n_bins_new = 0
        for h in mc_cvToStitch:
            for i in range(0,mc_cvToStitch[h].GetNbinsX()+1):
                if mc_cvToStitch[h].GetBinContent(i) > minBinCont:
                    n_bins_new+=1

        htemp_mc = ROOT.TH1D('mc_stitched', 'mc_stitched', n_bins_new, 0,  n_bins_new )
        htemp_data = ROOT.TH1D('data_stitched', 'data_stitched', n_bins_new, 0,  n_bins_new )

        StitchThis(htemp_mc, mc_cvToStitch,"",mc_cvToStitch)
        StitchThis(htemp_data, data_cvToStitch,"data",mc_cvToStitch)

        mnv_cv_data = PlotUtils.MnvH1D( htemp_data )
        mnv_cv_mc   = PlotUtils.MnvH1D( htemp_mc )

        FillErrorBandfromHist2( mnv_cv_data, data_cvToStitch, mc_cvToStitch, isData=True)
        FillErrorBandfromHist2( mnv_cv_mc, mc_cvToStitch, mc_cvToStitch, isData=False)
        
        for err in err_names:
            err = str(err)
            n_universes = mnv_cv_mc.GetVertErrorBand(err).GetNHists()
            print("loading {} errorband that has {} n_universes".format(err,n_universes))
            chi2s = []
            for univ in range(n_universes):
                # ---------------------- NuE Universe Selection Block -----------------------------
                h_fhc_universe_mc   = ROOT.TFile.Open(fhc_universe).Get("EN4_predicted_Signal_{}_{}".format(err,univ))
                h_fhc_universe_data = ROOT.TFile.Open(fhc_universe).Get("EN4_data_bkgSubbed_{}_{}".format(err,univ))

                h_rhc_universe_mc   = ROOT.TFile.Open(rhc_universe).Get("EN4_predicted_Signal_{}_{}".format(err,univ))
                h_rhc_universe_data = ROOT.TFile.Open(rhc_universe).Get("EN4_data_bkgSubbed_{}_{}".format(err,univ))
                
                # ---------------------- NuMu Universe Selection Block -----------------------------
                h_fhc_mu_universe_mc = h_fhc_muselection_mc.Clone()
                h_rhc_mu_universe_mc = h_rhc_muselection_mc.Clone()

                # ---------------------- Scattering Universe Selection Block -----------------------------
                h_fhc_nueel_universe_mc = h_fhc_nueel_mc.Clone()
                h_rhc_nueel_universe_mc = h_rhc_nueel_mc.Clone()
                h_fhc_imd_universe_mc = h_fhc_imd_mc.Clone()
                h_rhc_imd_universe_mc = h_rhc_imd_mc.Clone()

                h_fhc_nueel_universe_data = h_fhc_nueel_data.Clone()
                h_rhc_nueel_universe_data = h_rhc_nueel_data.Clone()
                h_fhc_imd_universe_data = h_fhc_imd_data.Clone()
                h_rhc_imd_universe_data = h_rhc_imd_data.Clone()
                
                #Make all hist have the same errorbands. Filled with CV if they don't have it.
                SyncErrorBandsv2([h_fhc_imd_universe_mc,h_rhc_imd_universe_mc,h_fhc_nueel_universe_mc,h_rhc_nueel_universe_mc,h_fhc_universe_mc,h_rhc_universe_mc,h_fhc_mu_universe_mc,h_rhc_mu_universe_mc])
                SyncErrorBandsv2([h_fhc_imd_universe_data,h_rhc_imd_universe_data,h_fhc_nueel_universe_data,h_rhc_nueel_universe_data,h_fhc_universe_data,h_rhc_universe_data,h_fhc_muselection_data,h_rhc_muselection_data])
               
                # ---------------------- Swap Each Systematic Universe with the CV -----------------------------
                # Do monte carlo histograms
                h_fhc_mu_universe_mc = SwapUniverseCV(h_fhc_mu_universe_mc,err,univ)
                h_rhc_mu_universe_mc = SwapUniverseCV(h_rhc_mu_universe_mc,err,univ)
                h_fhc_nueel_universe_mc = SwapUniverseCV(h_fhc_nueel_universe_mc,err,univ)
                h_rhc_nueel_universe_mc = SwapUniverseCV(h_rhc_nueel_universe_mc,err,univ)
                h_fhc_imd_universe_mc = SwapUniverseCV(h_fhc_imd_universe_mc,err,univ)
                h_rhc_imd_universe_mc = SwapUniverseCV(h_rhc_imd_universe_mc,err,univ)

                # Do data histograms to mimic background subtraction effect
                h_fhc_nueel_universe_data = SwapUniverseCV(h_fhc_nueel_universe_data,err,univ)
                h_rhc_nueel_universe_data = SwapUniverseCV(h_rhc_nueel_universe_data,err,univ)
                h_fhc_imd_universe_data = SwapUniverseCV(h_fhc_imd_universe_data,err,univ)
                h_rhc_imd_universe_data = SwapUniverseCV(h_rhc_imd_universe_data,err,univ)

                mc_histsToStitch = OrderedDict()
                mc_histsToStitch['fhc_nueel'] = h_fhc_nueel_universe_mc
                mc_histsToStitch['fhc_imd'] = h_fhc_imd_universe_mc
                mc_histsToStitch['fhc_muselection']=h_fhc_mu_universe_mc
                mc_histsToStitch['fhc_selection']=h_fhc_universe_mc
                mc_histsToStitch['rhc_nueel'] = h_rhc_nueel_universe_mc
                mc_histsToStitch['rhc_imd'] = h_rhc_imd_universe_mc
                mc_histsToStitch['rhc_muselection']=h_rhc_mu_universe_mc
                mc_histsToStitch['rhc_selection']=h_rhc_universe_mc

                data_histsToStitch = OrderedDict()
                data_histsToStitch['fhc_nueel'] = h_fhc_nueel_universe_data
                data_histsToStitch['fhc_imd'] = h_fhc_imd_universe_data
                data_histsToStitch['fhc_muselection']=h_fhc_muselection_data
                data_histsToStitch['fhc_selection']=h_fhc_universe_data
                data_histsToStitch['rhc_nueel'] = h_rhc_nueel_universe_data
                data_histsToStitch['rhc_imd'] = h_rhc_imd_universe_data
                data_histsToStitch['rhc_muselection']=h_rhc_muselection_data
                data_histsToStitch['rhc_selection']=h_rhc_universe_data

                for h in list(mc_histsToStitch.keys()):
                    if "fhc" in exclude and "selection" in h and "fhc" in h:
                        del(mc_histsToStitch[h])
                        del(data_histsToStitch[h])
                        continue
                    if "rhc" in exclude and "selection" in h and "rhc" in h:
                        del(mc_histsToStitch[h])
                        del(data_histsToStitch[h])
                        continue
                    if "numu" in exclude and "muselection" in h:
                        del(mc_histsToStitch[h])
                        del(data_histsToStitch[h])
                        continue
                    if "nue" in exclude and "_selection" in h:
                        del(mc_histsToStitch[h])
                        del(data_histsToStitch[h])
                        continue
                    if "elastic" in exclude and "nueel" in h:
                        del(mc_histsToStitch[h])
                        del(data_histsToStitch[h])
                        continue
                    if "imd" in exclude and "imd" in h:
                        del(mc_histsToStitch[h])
                        del(data_histsToStitch[h])
                        continue

                if doratio:
                    for h in list(mc_histsToStitch.keys()):
                        if "selection" in h:
                            del(mc_histsToStitch[h])
                            del(data_histsToStitch[h])
                    if "fhc" not in exclude:
                        toClone = h_fhc_mu_universe_mc.Clone()
                        toClone.Divide(toClone,h_fhc_universe_mc)
                        mc_histsToStitch['fhc_ratio']=toClone
                        toClone = h_fhc_muselection_data.Clone()
                        toClone.Divide(toClone,h_fhc_universe_data)
                        data_histsToStitch['fhc_ratio']=toClone
                    if "rhc" not in exclude:
                        toClone = h_rhc_mu_universe_mc.Clone()
                        toClone.Divide(toClone,h_rhc_universe_mc)
                        mc_histsToStitch['rhc_ratio']=toClone
                        toClone = h_rhc_muselection_data.Clone()
                        toClone.Divide(toClone,h_rhc_universe_data)
                        data_histsToStitch['rhc_ratio']=toClone

                huniverse_mc = ROOT.TH1D('mc_stitched_{}_{}'.format(err,univ), 'mc_stitched_{}_{}'.format(err,univ), n_bins_new, 0,  n_bins_new )
                huniverse_data = ROOT.TH1D('data_stitched_{}_{}'.format(err,univ), 'data_stitched_{}_{}'.format(err,univ), n_bins_new, 0,  n_bins_new )

                StitchThis(huniverse_mc, mc_histsToStitch,"",mc_cvToStitch)
                StitchThis(huniverse_data, data_histsToStitch,"data",mc_cvToStitch)

                mnv_universe_data = PlotUtils.MnvH1D( huniverse_data )
                mnv_universe_mc   = PlotUtils.MnvH1D( huniverse_mc )

                FillErrorBandfromHist2( mnv_universe_data, data_histsToStitch, mc_cvToStitch, isData=True)
                FillErrorBandfromHist2( mnv_universe_mc, mc_histsToStitch, mc_cvToStitch, isData=False)
                
                chi2 = Chi2DataMC(mnv_universe_data,mnv_universe_mc)
                chi2s.append(chi2)
                DataMCPlot(mnv_universe_data,mnv_universe_mc,mnv_cv_data.Clone(),mnv_cv_mc.Clone())

            file.write(err)
            for i in range(0,n_universes):
                print("writing {}th universe for {} errorband".format(i,err))
                file.write(", {}".format(chi2s[i]))
            file.write("\n")
    print("dnof: ",mnv_cv_data.GetNbinsX()," chi2: ",Chi2DataMC(mnv_cv_data,mnv_cv_mc))

