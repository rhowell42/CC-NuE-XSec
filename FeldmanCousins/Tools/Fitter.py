import os
import time
import logging, sys
import ROOT
import PlotUtils
import numpy as np
np.random.seed(0)
#np.set_printoptions(precision=1)
#np.set_printoptions(linewidth=1520)
#np.set_printoptions(threshold=sys.maxsize)
from scipy import optimize, integrate
import argparse
ccnueroot = os.environ.get('CCNUEROOT')

import math
import psutil
from array import array

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 

#logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
#MNVPLOTTER.chi2_use_overflow_err = True
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.04
#MNVPLOTTER.extra_top_margin = -.035# go slightly closer to top of pad
#MNVPLOTTER.mc_bkgd_color = 46 
#MNVPLOTTER.mc_bkgd_line_color = 46
MNVPLOTTER.legend_n_columns = 1
#MNVPLOTTER.mc_line_width = 0
MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  
#MNVPLOTTER.axis_maximum = 500 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .03
MNVPLOTTER.legend_offset_x           = .15
#MNVPLOTTER.mc_line_width = 0
#MNVPLOTTER.axis_minimum=0.1
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
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

def constraint(t):
    U_e4  = t[1]*0.15
    U_mu4 = t[2]*0.41
    U_tau4= t[3]*0.66
    return 1-(U_e4 + U_mu4 + U_tau4)

class MyTakeStep(object):
    def __init__(self, stepsize=0.5):
        self.stepsize = stepsize
    def __call__(self, x):
        s = self.stepsize
        rand = np.random.rand()
        direction = -1 if rand < 0.5 else 1
        old = x[0]
        x[0] += np.random.uniform(0.1*s, 0.2*s)*direction
        print("stepping from {:.1f} eV2 to {:.1f} eV2".format(old*100,x[0]*100))
        #x[1] += np.random.uniform(-.5*s, .5*s)
        #x[2] += np.random.uniform(-0.1*s, .1*s)
        return(x)

def doFit(stitched_data, templates, stitched_mc, makePlot = False, plotArgs = []):
    x0 = [0.0,0.0,0.0,0.0]
    bounds = np.array([[0.0,1],[0.0,1],[0.0,1],[0,1]], dtype = float)
    cons = [{"type":"ineq","fun":constraint}]
    mc = stitched_mc.Clone()
    data = stitched_data.Clone()
    mystep = MyTakeStep()
    minimizer_kwargs = {"bounds":bounds, "tol":1e-5, "options":{"maxiter":1000},"constraints":cons, "method":"SLSQP","args":(mc,templates,data)}
    res = optimize.basinhopping(CalChi2, x0,interval=10, take_step=mystep, stepsize=1, niter = 100,disp=True,minimizer_kwargs=minimizer_kwargs)
    chi2 = float(res.fun)
    
    if makePlot:
        fitms = res.x[0]*100
        fitUe4 = res.x[1]*0.15
        fitUmu4 = res.x[2]*0.41
        fitUtau4 = res.x[3]*0.66
        nit = res.nit
        GetOscillatedHistogramCategories(stitched_mc,stitched_data,templates,fitms,fitUe4,fitUmu4,fitUtau4,plotArgs)

    return(chi2,{"m":res.x[0]*100,"ue4":res.x[1],"umu4":res.x[2],"utau4":res.x[3]})

def fitNorm(stitched_data, templates, stitched_mc):
    x0 = [1.5,0.5]
    bounds = np.array([[0,10],[0,10]])
    mc = stitched_mc.Clone()
    data = stitched_data.Clone()

    res = optimize.minimize(NormCalc, x0,tol=1,bounds=bounds,method="Nelder-Mead",args=(mc,templates,data),options={"disp":True})
    chi2 = float(res.fun)
    bestfit = MuonNorm(stitched_mc,templates,res.x[0],res.x[1])
    return(chi2,bestfit,{"fhc":res.x[0],"rhc":res.x[1]})

def NormCalc(x,stitched_mc,templates,stitched_data):
    fhc_norm = x[0]
    rhc_norm = x[1]
    normalized = MuonNorm(stitched_mc,templates,fhc_norm,rhc_norm)
    PlotNorms(stitched_mc,stitched_data,normalized,[fhc_norm,rhc_norm])
    chi2 = Chi2DataMC(stitched_data,normalized)
    return(chi2)

def PlotNorms(stitched_mc,stitched_data,fitHist,params=[]):
    chi2_model = Chi2DataMC(stitched_data,fitHist)
    chi2_null = Chi2DataMC(stitched_data,stitched_mc)

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
    if len(params) > 0:
        leg.SetHeader("fhc:{:.2f} rhc:{:.2f}".format(params[0],params[1]))
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
    if len(params) > 0:
        overall.Print("fit_norm_{:.2f}.png".format(chi2_model))
    else:
        overall.Print("fit_muonnorm.png")

def MuonNorm(stitched_mc,templates,fhc_norm,rhc_norm):
    stitched_nueTemp = templates["nue"]
    stitched_numuTemp = templates["numu"]
    stitched_swapTemp = templates["swap"]

    stitched_nue_energy = templates["nue_energy"].Clone()
    stitched_numu_energy = templates["numu_energy"].Clone()
    stitched_nutau_id = templates["nutau_energy"].Clone()
    stitched_ratio_id = templates["ratio_energy"].Clone()
    stitched_fhc_id = templates["fhc_energy"].Clone()
    stitched_rhc_id = stitched_fhc_id.Clone()
    stitched_inv_ratio_id = stitched_ratio_id.Clone()
    stitched_swap_energy = templates["swap_energy"].Clone()
    stitched_nueselection_id = templates["nueselection_energy"].Clone()

    InvertID(stitched_inv_ratio_id)
    InvertID(stitched_rhc_id)

    fhc_numu_energy = stitched_numu_energy.Clone()
    rhc_numu_energy = stitched_numu_energy.Clone()

    fhc_numu_energy.MultiplySingle(fhc_numu_energy,stitched_fhc_id)
    rhc_numu_energy.MultiplySingle(rhc_numu_energy,stitched_rhc_id)

    fhc_numu_energy.Scale(1/fhc_norm)
    rhc_numu_energy.Scale(1/rhc_norm)

    numu = fhc_numu_energy.Clone()
    numu.Add(rhc_numu_energy)
    nue = stitched_nue_energy.Clone()

    ratio = numu.Clone()
    ratio.Divide(ratio,nue)
    ratio.MultiplySingle(ratio,stitched_ratio_id)

    nue.MultiplySingle(nue,stitched_inv_ratio_id)
    numu.MultiplySingle(numu,stitched_inv_ratio_id)

    hist = numu.Clone()
    hist.Add(nue)
    hist.Add(ratio)
    return(hist)

def Chi2DataMC(dataHist,mcHist):
    #We get the number of bins and make sure it's compatible with the NxN matrix given
    if dataHist.GetNbinsX() != mcHist.GetNbinsX():
        print("breaking error in Chi2DataMC")
        #logging.error("MnvPlotter::Chi2DataMC", "The number of bins from Data and MC histograms differ. Returning -1.")
        return(-1)

    Nbins = dataHist.GetNbinsX()

    #get the covariance matrix
    covMatrix = np.zeros(shape=[Nbins,Nbins],dtype='f')
    useOnlyShapeErrors = False
    includeStatError   = True
    errorAsFraction    = False

    covMatrixTmp  =   mcHist.GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)
    covMatrixTmp += dataHist.GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)

    for i in range(0,Nbins):
        for j in range(0,Nbins):
            if errorAsFraction:
                covMatrix[i][j] = covMatrixTmp[i+1][j+1] * (dataHist.GetBinContent(i+1)*dataHist.GetBinContent(j+1))
            else:
                covMatrix[i][j] = covMatrixTmp[i+1][j+1]

    errorMatrix = np.linalg.inv(covMatrix)

    chi2 = 0.
    # under/overflow bins not taken into account in the chi2 calculation
    for i in range(0,errorMatrix.shape[0]):
        x_data_i = dataHist.GetBinContent(i+1);
        x_mc_i   = mcHist.GetBinContent(i+1);

        for j in range(0,errorMatrix.shape[0]):
            x_data_j = dataHist.GetBinContent(j+1)
            x_mc_j   = mcHist.GetBinContent(j+1)
            chi2_ij = (x_data_i - x_mc_i) * errorMatrix[i][j] * (x_data_j - x_mc_j)
            chi2 += chi2_ij 

    return(chi2)

def CalChi2(x,stitched_mc,templates,stitched_data):
    ms = x[0]*100
    U_e4 = x[1]*0.15
    U_mu4 = x[2]*0.41
    U_tau4 = x[3]*0.66
    oscillated = GetOscillatedHistogram(stitched_mc,templates,ms,U_e4,U_mu4,U_tau4)
    chi2 = Chi2DataMC(stitched_data,oscillated)
    return(chi2)

def InvertID(hist):
    for i in range(hist.GetNbinsX()+1):
        if hist.GetBinContent(i) == 0:
            hist.SetBinContent(i,1.0)
            hist.SetBinError(i,0.0)
        elif hist.GetBinContent(i) == 1:
            hist.SetBinContent(i,0.0)
            hist.SetBinError(i,0.0)

def GetOscillatedHistogram(histogram, templates, m, U_e4, U_mu4, U_tau4,doSys = False):
    # stitched_nueselection_energy is 1 only for the nueselection sample
    # stitched_nutau_energy is 1 only for the nueel selection
    # stitched_ratio_energy is 1 only for the ratio sample
    # if 0 is False, if 1 is True

    stitched_nueTemp = templates["nue"]
    stitched_numuTemp = templates["numu"]
    stitched_swapTemp = templates["swap"]

    stitched_nue_energy = templates["nue_energy"].Clone()
    stitched_numu_energy = templates["numu_energy"].Clone()
    stitched_nutau_id = templates["nutau_energy"].Clone()
    stitched_ratio_id = templates["ratio_energy"].Clone()
    stitched_inv_ratio_id = stitched_ratio_id.Clone()
    stitched_swap_energy = templates["swap_energy"].Clone()
    stitched_nueselection_id = templates["nueselection_energy"].Clone()

    InvertID(stitched_inv_ratio_id)

    hist = histogram.Clone()
    nue_weights = hist.GetCVHistoWithStatError().Clone()
    numu_weights = hist.GetCVHistoWithStatError().Clone()
    nuenutau_weights = hist.GetCVHistoWithStatError().Clone()
    numunue_weights = hist.GetCVHistoWithStatError().Clone()
    numunutau_weights = hist.GetCVHistoWithStatError().Clone()

    for i in range(0,hist.GetNbinsX() + 1):
        nue_sin = sin_average(i,m,stitched_nueTemp)
        numu_sin = sin_average(i,m,stitched_numuTemp)
        swap_sin = sin_average(i,m,stitched_swapTemp)

        P_ee = float(1 - 4*U_e4*(1-U_e4)*nue_sin)

        P_mue = float(4*(U_e4)*(U_mu4)*swap_sin)

        P_mumu = float(1 - 4*U_mu4*(1-U_mu4)*numu_sin)

        P_mutau = float(4*U_tau4*U_mu4*numu_sin)

        P_etau = float(4*U_e4*U_tau4*nue_sin)

        nue_weights.SetBinContent(i,P_ee)
        numu_weights.SetBinContent(i,P_mumu)
        nuenutau_weights.SetBinContent(i,P_etau)
        numunue_weights.SetBinContent(i,P_mue)
        numunutau_weights.SetBinContent(i,P_mutau)

        nue_weights.SetBinError(i,0)
        numu_weights.SetBinError(i,0)
        nuenutau_weights.SetBinError(i,0)
        numunue_weights.SetBinError(i,0)
        numunutau_weights.SetBinError(i,0)
        stitched_ratio_id.SetBinError(i,0)
        stitched_nutau_id.SetBinError(i,0)
        stitched_nueselection_id.SetBinError(i,0)

    c1 = ROOT.TCanvas()
    numu = stitched_numu_energy.Clone()
    numu.MultiplySingle(numu,numu_weights)

    nue = stitched_nue_energy.Clone()
    nue.MultiplySingle(nue,nue_weights)

    numunue = stitched_numu_energy.Clone()
    numunue.MultiplySingle(numunue,numunue_weights)
    nue.Add(numunue)

    numunutau = stitched_numu_energy.Clone()
    numunutau.MultiplySingle(numunutau,numunutau_weights)

    nuenutau  = stitched_nue_energy.Clone()
    nuenutau.MultiplySingle(nuenutau,nuenutau_weights)
    
    nutau = numunutau.Clone()
    nutau.Add(nuenutau)
    nutau.MultiplySingle(nutau,stitched_nutau_id)

    ratio = numu.Clone()
    ratio.Divide(ratio,nue)
    ratio.MultiplySingle(ratio,stitched_ratio_id)

    nue.MultiplySingle(nue,stitched_inv_ratio_id)
    numu.MultiplySingle(numu,stitched_inv_ratio_id)

    hist = numu.Clone()
    hist.Add(nue)
    hist.Add(ratio)
    if nutau.Integral() > 0:
        hist.Add(nutau)

    return(hist)

def sin_average(q = 0,dm2 = 0,template = None):
    avgsin = 0
    total_N = 0
    for x in range(template.GetNbinsY()+1):
        lowEdge = template.GetYaxis().GetBinLowEdge(x)
        upEdge = template.GetYaxis().GetBinUpEdge(x)
        bin_width = upEdge-lowEdge
        bin_center = (upEdge+lowEdge)/2
        N_bin = template.GetBinContent(q,x)
        total_N+=N_bin
        if N_bin == 0:
            continue
        nue_sin = np.sin(1.27*dm2*bin_center)**2
        avgsin += nue_sin * N_bin
    if total_N != 0:
        avgsin = avgsin / total_N
    return(avgsin)
