#python 2/3 compatible lines
from __future__ import print_function
if hasattr(__builtins__, 'raw_input'):
    input = raw_input

import os
import sys
import ROOT
import PlotUtils
import math
import copy
from array import array
from collections import OrderedDict
from tools.PlotLibrary import PLOT_SETTINGS,VariantPlotsNamingScheme, HistHolder
from config.AnalysisConfig import AnalysisConfig
from config import BackgroundFitConfig
from tools import Utilities,PlotTools


# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

fitPath = '/exp/minerva/data/users/shenry/antinu_e'
nue_mcfile = 'kin_dist_mcFHCScaled_nx_collab1_MAD.root'
nue_datafile = 'kin_dist_dataFHCScaled_nx_collab1_MAD.root'
#PLOTPATH = "/exp/minerva/data/users/shenry/WrongSign/"
PLOTPATH = "/exp/minerva/app/users/shenry/cmtuser/gitMAT/CCNue/background_fit/" 

#nue_file = 'kin_dist_mcFHCFinalScaled_nx_new_fspline.root'.format(fitPath)
#nue_datafile =  'kin_dist_dataFHCFinalScaled_nx_new_fspline.root'.format(fitPath)
#nue_mcfile = 'kin_dist_mcFHCFinalScaled_nx_new_fspline.root'.format(fitPath)

# nue_background = ROOT.TFile.Open("{}/{}".format(fitPath,nue_file))
# data_background = ROOT.TFile.Open("{}/{}".format(fitPath,nue_datafile))

# nuedata_map_path = "{}/{}".format(fitPath,nue_datafile)
# nuemc_map_path = "{}/{}".format(fitPath,nue_mcfile)

HISTOGRAMS_TO_NUE_ADD = [
  "Visible Energy vs q3",
  "Visible Energy vs Lepton Pt",
#  "Lepton Energy",
#  "Lepton Theta",
#  "Q3",
#  "Visible Energy",
#  "Lepton Pt",
]



def DataSubtraction(data_hist,nue_hist,background_hist,pot_scale):
    out_data = data_hist.Clone()
    out_nue = nue_hist.Clone()
    out_background = background_hist.Clone()  
     
    out_nue.AddMissingErrorBandsAndFillWithCV(nue_hist)
    out_background.AddMissingErrorBandsAndFillWithCV(background_hist)
    
    out_nue.Scale(pot_scale) #FHC data pot normalize
    out_background.Scale(pot_scale) #FHC data pot normalize
    out_data.AddMissingErrorBandsAndFillWithCV(background_hist) 
    out_data.Add(out_nue,-1)
    out_data.Add(out_background,-1)
    print("Data Subtraction:") 
    return out_data

def DataSubtractionFHC(data_hist,mc_hist,wrong_sign_hist,ref): 
    data_hist.AddMissingErrorBandsAndFillWithCV(ref)
    mc_hist.AddMissingErrorBandsAndFillWithCV(ref)
    wrong_sign_hist.AddMissingErrorBandsAndFillWithCV(ref)
    data_hist.Add(mc_hist,-1)
    data_hist.Add(wrong_sign_hist) 
    return data_hist

def DrawWrongSign(mnvplotter,data_hist,mc_hist):
    data_hist.SetLineColor(ROOT.kRed-2)
    data_hist.Draw("HIST")
    mc_hist.SetLineColor(ROOT.kBlue-2)
    mc_hist.Draw("HIST SAME")
    leg = ROOT.TLegend(0.6,0.6,0.9,0.9)
    leg.AddEntry(data_hist,"Constrained wrong sign")
    leg.AddEntry(mc_hist,"GENIE wrong sign")
    ROOT.SetOwnership(leg,False)
    leg.Draw()
    

if __name__ == "__main__":
    
    #input knobs
    playlist= AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    data_file,mc_reco_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    
    nue_files = []
    nue_pots = []
    for i in [nue_datafile,nue_mcfile]:
        path = "{}/{}".format(fitPath,i)
        nue_files.append(ROOT.TFile.Open(path))
        nue_pots.append(Utilities.getPOTFromFile(path))

   
    regions = ["Signal"] + AnalysisConfig.sidebands
    newplaylist=ROOT.TFile.Open('/exp/minerva/data/users/shenry/antinu_e/WrongSign.root',"UPDATE")
    
    for region in regions:
        for config in HISTOGRAMS_TO_NUE_ADD: 
            data_nue = HistHolder(config["name"] if "name" in config else config,nue_files[0],region,False,nue_pots[0],nue_pots[0])
           
            mc_nue = HistHolder(config["name"] if "name" in config else config,nue_files[1],region,True,nue_pots[1],nue_pots[0])
            ref = HistHolder(config["name"] if "name" in config else config,mc_reco_file,region,True,mc_pot,mc_pot)
            data_nue.POTScale(False)
            mc_nue.POTScale(False)
           
            pred = DataSubtractionFHC(data_nue.GetHist(),mc_nue.GetHist(),mc_nue.hists["CCNu"],ref.GetHist())
            newplaylist.cd()
            pred.Scale(mc_pot/nue_pots[0])
            pred.Write("{}_{}".format(pred.GetName(),"CCNu"))
  
            #drawing
            Slicer = PlotTools.Make2DSlice 
            PlotTools.MakeGridPlot(Slicer,DrawWrongSign,[pred,ref.hists["CCNu"]],draw_seperate_legend = True)
            PlotTools.CANVAS.Print("{}{}.png".format(PLOTPATH,pred.GetName()))
               
    if input("Warning: Will update mc file, continue? (y/N)").upper()!="Y":
        newplaylist.Close()
        exit(1)
    mc_reco_file.Close()
    mc_reco_file=ROOT.TFile.Open(type_path_map["mc"],"UPDATE")
    #mc_reco_file=ROOT.TFile.Open("/exp/minerva/data/users/shenry/antinu_e/kin_dist_mcNoScaleSys_nx_collab1_fspline.root","UPDATE")
    newplaylist.ReOpen("READ")
    region="Signal" 
    for config in HISTOGRAMS_TO_NUE_ADD:
        ref = HistHolder(config["name"] if "name" in config else config,mc_reco_file,region,True,1)
        pred = newplaylist.Get(ref.plot_name+"_CCNu")
        ref.hists["CCNu"]=pred.Clone(ref.hists["CCNu"].GetName())
        print (ref.hists["CCNu"].GetName(),config,HISTOGRAMS_TO_NUE_ADD)
        ref.ResumTotal()
        mc_reco_file.cd()
        
       
        for i in ["Total","CCNu"]:
            h = pred.hists[i]
            print(h.GetName(),pred.hists[i]) 
            pred.Write(pred.GetName(),ROOT.TObject.kOverwrite)
        print(ref.hists["CCNu"].GetBinContent(1,1))




    newplaylist.Close()
    mc_reco_file.Write()
    mc_reco_file.Close()
    data_file.Close()



