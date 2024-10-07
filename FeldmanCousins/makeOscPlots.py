import os
import logging, sys
import ROOT
import PlotUtils
import numpy as np
np.random.seed(0)
np.set_printoptions(precision=1)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)
from scipy import optimize, integrate
import argparse
ccnueroot = os.environ.get('CCNUEROOT')

import math
import psutil
import multiprocessing
import threading
nthreads = 4
from array import array

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from fit_tools.FitTools import *

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
#MNVPLOTTER.chi2_use_overflow_err = True
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.04
#MNVPLOTTER.extra_top_margin = -.035# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46 
MNVPLOTTER.mc_bkgd_line_color = 46
MNVPLOTTER.legend_n_columns = 3
#MNVPLOTTER.mc_line_width = 0
MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  
MNVPLOTTER.axis_maximum = 500 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .03
MNVPLOTTER.legend_offset_x           = .15

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

if __name__ == "__main__":
    BLUEARC = "/exp/minerva/data/users/{}/surfaces".format(os.environ["USER"])
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--grid",
                        action="store_true",
                        default = False,
                        help = "Run macro on grid, Input/Output path must be updated to avoid direct access to BlueArc"
    )
    parser.add_argument("-o", "--output",
                        dest = "output_dir",
                        help="Use alternate location for output file.",
                        default=BLUEARC
    )
    parser.add_argument("-m", "--delta_m",
                        dest = "delta_m",
                        help="Delta m^2 value to probe.",
                        type=float,
                        default=0
    )
    parser.add_argument("-U", "--U_tau4",
                        dest = "U_tau4",
                        help="U_tau4 parameter to probe.",
                        type=float,
                        default=0
    )
    parser.add_argument("--pseudodata",
                        dest = "pseudodata",
                        default = False,
                        action="store_true",
    )

    args = parser.parse_args()
    outdir_surface = args.output_dir
    delta_m = args.delta_m
    runongrid = args.grid
    U_tau4 = args.U_tau4
    pseudodata = args.pseudodata

    filename = ""
    if pseudodata:
        filename = "NuE_stitched_hists_pseudo.root"
    else:
        filename = "NuE_stitched_hists.root"

    stitched_data = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('data_stitched')
    stitched_mc = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched')
    
    stitched_nueTemp = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('LE_template_nue')
    stitched_numuTemp = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('LE_template_numu')
    stitched_swapTemp = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('LE_template_swap')

    stitched_nue_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_nue')
    stitched_numu_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_numu')
    stitched_nutau_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_nutau')
    stitched_nueselection_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_nueselection')
    stitched_ratio_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_ratio')
    stitched_swap_energy = ROOT.TFile.Open("{}/FeldmanCousins/{}".format(ccnueroot,filename)).Get('mc_stitched_swap')

    templates = {
            "nue":stitched_nueTemp,
            "numu":stitched_numuTemp,
            "swap":stitched_swapTemp,
            "nue_energy":stitched_nue_energy,
            "numu_energy":stitched_numu_energy,
            "nutau_energy":stitched_nutau_energy,
            "swap_energy":stitched_swap_energy,
            "ratio_energy":stitched_ratio_energy,
            "nueselection_energy":stitched_nueselection_energy,
    }

    #bestFit = {"m":[7,7,7],"ue4":[0.12,.12,.12],"umu4":[0,0.005,0.01]}  # RAA
    #bestFit = {"m":[3,2.5,2.5],"ue4":[0.1,.054,0.015],"umu4":[.013,.031,.073]} # best fits
    #bestFit = {"m":[4,8,10,80,10,10,10],"ue4":[.03,.03,0.03,.03,.02,.1,.12],"umu4":[.031,.031,0.031,0.031,.005,.005,0]} # m effect
    #bestFit = {"m":[10,10],"ue4":[.02,.1],"umu4":[.005,0.005]} # Ue4 effects
    #bestFit = {"m":[10,10],"ue4":[.03,.03],"umu4":[.073,.173]} # tau/Umu4 effects
    #bestFit = {"m":[10,10],"ue4":[.03,.03],"umu4":[.073,.173]} # tau/Umu4 effects
    bestFit = {"m":[3,3],"ue4":[.1,.02],"umu4":[.013,.073]} # tau/Umu4 effects

    if True:
        for fit in range(len(bestFit["m"])):
            fitHist = GetOscillatedHistogram(stitched_mc, templates, bestFit["m"][fit], bestFit["ue4"][fit], bestFit["umu4"][fit], 0.0)
            chi2_model = Chi2DataMC(stitched_data,fitHist)
            chi2_null = Chi2DataMC(stitched_data,stitched_mc)

            print(bestFit["m"][fit], bestFit["ue4"][fit], bestFit["umu4"][fit],chi2_model-chi2_null)

            c1 = ROOT.TCanvas()
            margin = .12
            bottomFraction = .2
            overall = ROOT.TCanvas("Data/MC")
            top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
            bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

            top.Draw()
            bottom.Draw()

            top.cd()
            top.SetLogy()
            
            fitHist.GetXaxis().SetTitle("Bin number")
            fitHist.GetYaxis().SetTitle("Entries")
            MNVPLOTTER.ApplyAxisStyle(fitHist,True,True)

            nullRatio =  stitched_data.Clone()
            oscRatio =  fitHist.Clone()

            fitHist.SetLineColor(ROOT.kRed)
            fitHist.SetLineWidth(3)
            stitched_mc.SetLineColor(ROOT.kBlue)
            stitched_mc.SetLineWidth(3)
            stitched_mc.GetYaxis().SetTitle("Nevents")
            stitched_mc.Draw("hist")
            fitHist.Draw("hist same")
            stitched_data.Draw("same")

            osc = fitHist.GetCVHistoWithError()
            osc.SetLineColor(ROOT.kRed)
            osc.SetLineWidth(3)
            osc.SetMarkerStyle(0)
            osc.SetFillColorAlpha(ROOT.kPink + 1, 0.3)
            osc.Draw("E2 SAME")

            null = stitched_mc.GetCVHistoWithError()
            null.SetLineColor(ROOT.kBlue)
            null.SetLineWidth(3)
            null.SetMarkerStyle(0)
            null.SetFillColorAlpha(ROOT.kBlue + 1, 0.3)
            null.Draw("E2 SAME")

            leg = ROOT.TLegend()
            leg.AddEntry(stitched_data,"Data","p")
            leg.AddEntry(fitHist,"Best Fit #chi^{2}="+"{:.1f}".format(chi2_model),"l")
            #leg.AddEntry(fitHist,"Oscillation #chi^{2}="+"{:.1f}".format(chi2_model),"l")
            #leg.AddEntry(fitHist,"RAA #chi^{2}="+"{:.1f}".format(chi2_model),"l")
            leg.AddEntry(stitched_mc,"Null Hypothesis #chi^{2}="+"{:.1f}".format(chi2_null),"l")
            leg.Draw()

            oscRatio.Divide(oscRatio, stitched_mc)
            nullRatio.Divide(nullRatio,stitched_mc)

            bottom.cd()
            bottom.SetTopMargin(0)
            bottom.SetBottomMargin(0.3)

            nullErrors = stitched_mc.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
            for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
                nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
                nullErrors.SetBinContent(whichBin, 1)

            oscRatio.SetTitle("")
            oscRatio.SetLineColor(ROOT.kRed)
            nullRatio.SetLineColor(ROOT.kBlue)
            oscRatio.SetLineWidth(3)
            nullRatio.SetLineWidth(3)
            oscRatio.SetTitleSize(0)

            #Error envelope for the MC
            nullErrors.SetLineWidth(0)
            nullErrors.SetMarkerStyle(0)
            nullErrors.SetFillColorAlpha(ROOT.kBlue + 1, 0.4)
            nullErrors.Draw("E2")

            nullErrors.GetYaxis().SetTitle("#splitline{Ratio to}{Null Hypothesis}")
            nullErrors.GetYaxis().SetLabelSize(.13)
            nullErrors.GetYaxis().SetTitleSize(0.1)
            nullErrors.GetYaxis().SetTitleOffset(0.6)
            nullErrors.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
            nullErrors.GetXaxis().SetTitleSize(0.16)
            nullErrors.GetXaxis().SetTitleOffset(0.9)
            nullErrors.GetXaxis().SetLabelSize(.15)
            nullErrors.GetXaxis().SetTitle("Bin Number")
            nullErrors.SetMinimum(0.5)
            nullErrors.SetMaximum(1.5)

            
            #Draw the data ratios
            oscRatio.SetMinimum(0.5)
            oscRatio.SetMaximum(1.5)
            nullRatio.SetMinimum(0.5)
            nullRatio.SetMaximum(1.5)
            nullRatio.SetLineColorAlpha(ROOT.kBlue+1,0.6)
            oscRatio.Draw("same")
            nullRatio.Draw("same")

            #Draw a flat line at 1 for oscRatio of MC to itself
            straightLine = nullErrors.Clone()
            straightLine.SetLineColor(ROOT.kBlue)
            straightLine.SetLineWidth(3)
            straightLine.SetFillColor(0)
            straightLine.Draw("HIST L SAME")

            leg1 = ROOT.TLegend(.5,.5,.9,.9)
            leg1.AddEntry(nullRatio,"Data/Null Hypothesis","p")
            leg1.AddEntry(oscRatio,"Model/Null Hypothesis","l")
            leg1.AddEntry(straightLine,"Null/Null Hypothesis","l")
            #leg1.Draw()

            top.cd()

            #overall.Print("fit_{}_thesis.png".format(bestFit["m"][fit]))
            overall.Print("fit_{}_{}_{}_thesis.png".format(bestFit["m"][fit],bestFit['ue4'][fit],bestFit['umu4'][fit]))
