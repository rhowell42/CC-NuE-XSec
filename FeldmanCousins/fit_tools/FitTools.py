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
import multiprocessing
import threading
nthreads = 4
from array import array

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
#MNVPLOTTER.chi2_use_overflow_err = True
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.04
#MNVPLOTTER.extra_top_margin = -.035# go slightly closer to top of pad
#MNVPLOTTER.mc_bkgd_color = 46 
#MNVPLOTTER.mc_bkgd_line_color = 46
MNVPLOTTER.legend_n_columns = 3
#MNVPLOTTER.mc_line_width = 0
MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  
#MNVPLOTTER.axis_maximum = 500 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .03
MNVPLOTTER.legend_offset_x           = .15
MNVPLOTTER.mc_line_width = 0
MNVPLOTTER.axis_minimum=0.1

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

rowstodo = 100

def constraint(t):
    U_e4  = t[1]
    U_mu4 = t[2]
    U_tau4= t[3]
    return 1-(U_e4 + U_mu4 + U_tau4)

def getFits(matrix,stitched_data,name,templates,stitched_mc):
    c = 0
    dchi2s = np.zeros(rowstodo)
    logging.debug("\nI want to loop through array of size: {}".format(matrix.shape[0]))
    for row in matrix:
        universe = stitched_data.Clone()
        for i in range(len(row)):
            before = universe.GetBinContent(i+1)
            scale = row[i]
            newCont = before*(1+scale)
            if newCont < 0:
                newCont = 0
            universe.SetBinContent(i+1,newCont)
        minchi2,_ = doFit(universe,templates,stitched_data)
        chi2 = Chi2DataMC(universe,stitched_mc)
        dchi2 = chi2 - minchi2
        dchi2s[c] = dchi2
        c+=1
        #logging.debug("done with row {}".format(c))
    #logging.info("writing chi2_{}.txt to disk".format(name))
    np.savetxt('chi2s/chi2_{}.txt'.format(name),dchi2s)

class MyTakeStep(object):
    def __init__(self, stepsize=0.5):
        self.stepsize = stepsize
    def __call__(self, x):
        s = self.stepsize
        x[0] += np.random.uniform(-2.0*s, 2.0*s)
        x[1] += np.random.uniform(-.5*s, .5*s)
        x[2] += np.random.uniform(-0.1*s, .1*s)
        return(x)

def doFit(stitched_data, templates, stitched_mc, makePlot = False, plotArgs = []):
    x0 = [0.0,0.0,0.0,0.0]
    bounds = np.array([[0.0,1.0],[0.0,0.15],[0.0,0.41],[0,0.66]], dtype = float)
    cons = [{"type":"ineq","fun":constraint}]
    mc = stitched_mc.Clone()
    data = stitched_data.Clone()

    minimizer_kwargs = {"bounds":bounds, "tol":1e-12, "options":{"maxiter":1000},"constraints":cons, "method":"SLSQP","args":(mc,templates,data)}
    res = optimize.basinhopping(CalChi2, x0,interval=5, niter = 60,disp=True,minimizer_kwargs=minimizer_kwargs)
    chi2 = float(res.fun)
    
    if makePlot:
        fitms = res.x[0]*100
        fitUe4 = res.x[1]
        fitUmu4 = res.x[2]
        fitUtau4 = res.x[3]
        nit = res.nit
        GetOscillatedHistogramCategories(stitched_mc,stitched_data,templates,fitms,fitUe4,fitUmu4,fitUtau4,plotArgs)

    return(chi2,{"m":res.x[0]*100,"ue4":res.x[1],"umu4":res.x[2],"utau4":res.x[3]})

def Chi2DataMC(dataHist,mcHist):
    #We get the number of bins and make sure it's compatible with the NxN matrix given
    if dataHist.GetNbinsX() != mcHist.GetNbinsX():
        print("breaking error in Chi2DataMC")
        logging.error("MnvPlotter::Chi2DataMC", "The number of bins from Data and MC histograms differ. Returning -1.")
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
    U_e4 = x[1]
    U_mu4 = x[2]
    U_tau4 = x[3]
    oscillated = GetOscillatedHistogram(stitched_mc,templates,ms,U_e4,U_mu4,U_tau4)
    chi2 = Chi2DataMC(stitched_data,oscillated)
    #print("chi2: {:.2f}".format(chi2))
    print("{:.12f},{:.4f},{:.4f},{:.4f},{:.4f}".format(chi2,ms,U_e4,U_mu4,U_tau4))
    return(chi2)

def GetOscillatedHistogramCategories(histogram,data_histogram,templates, m, U_e4, U_mu4,U_tau4,plotArgs = []):
    stitched_nueTemp = templates["nue"]
    stitched_numuTemp = templates["numu"]

    stitched_nue_energy = templates["nue_energy"]
    stitched_numu_energy = templates["numu_energy"]
    stitched_nutau_energy = templates["nutau_energy"]
    stitched_nueselection_energy = templates["nueselection_energy"]

    hist = histogram.Clone()

    nueCat = hist.Clone()
    numuCat = hist.Clone()
    nueswapCat = hist.Clone()
    nutauswapCat = hist.Clone()

    # stitched_nueselection_energy is 1 for every sample except the nueselection
    # stitched_nutau_energy is 1 only for the nueel selection
    # if 0 is False, if 1 is True

    for q in range(0,hist.GetNbinsX()+1):
        nue_sin = sin_average(q,m,stitched_nueTemp)
        numu_sin = sin_average(q,m,stitched_numuTemp)
        
        P_ee = 1 - 4*U_e4*(1-U_e4)*nue_sin
        eeCont = P_ee * stitched_nue_energy.GetBinContent(q)

        P_mue = 4*U_e4*U_mu4*numu_sin
        mueCont = P_mue * stitched_numu_energy.GetBinContent(q)
        if stitched_nueselection_energy.GetBinContent(q) == 1 and stitched_nutau_energy.GetBinContent(q) != 1:
            mueCont = 0

        P_mumu = 1 - 4*U_mu4*(1-U_mu4)*numu_sin
        mumuCont = P_mumu * stitched_numu_energy.GetBinContent(q)
        mumuCont*=stitched_nueselection_energy.GetBinContent(q)

        P_mutau = 4*U_tau4*U_mu4*numu_sin
        mutauCont = P_mutau * stitched_numu_energy.GetBinContent(q)

        P_etau = 4*U_e4*U_tau4*nue_sin
        etauCont = P_etau * stitched_nue_energy.GetBinContent(q)

        nueCat.SetBinContent(q,eeCont)
        numuCat.SetBinContent(q,mumuCont)
        nueswapCat.SetBinContent(q,mueCont)
        tauCont = (etauCont + mutauCont) * stitched_nutau_energy.GetBinContent(q)
        nutauswapCat.SetBinContent(q,tauCont)

        totCont = eeCont+mumuCont+mueCont+tauCont
        weight = totCont/hist.GetBinContent(q) if hist.GetBinContent(q) > 0 else 0

        hist.SetBinContent(q,totCont)


        ### oscillate CV bin content along with systematic universe contents ###
        errorband_names_vert = hist.GetVertErrorBandNames()
        errorband_names_lat = hist.GetLatErrorBandNames()
        for error_band in errorband_names_vert:
            errorband = hist.GetVertErrorBand(error_band)
            errorband.SetBinContent(q,errorband.GetBinContent(q)*weight)
            n_univ = hist.GetVertErrorBand(error_band).GetNHists()
            for universe in range(0, n_univ):
                sys_bc = hist.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(q)
                hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(q, sys_bc*weight)

        for error_band in errorband_names_lat:
            errorband = hist.GetLatErrorBand(error_band)
            errorband.SetBinContent(q,errorband.GetBinContent(q)*weight)
            n_univ = hist.GetLatErrorBand(error_band).GetNHists()
            for universe in range(0,  n_univ):
                sys_bc = hist.GetLatErrorBand(error_band).GetHist(universe).GetBinContent(q)
                hist.GetLatErrorBand(error_band).GetHist(universe).SetBinContent(q, sys_bc*weight)

    mc_hists = [nueCat,numuCat,nueswapCat,nutauswapCat]
    title = ["remaining #nu_{e}","remaining #nu_{#mu}","appeared #nu_{e}","appeared #nu_{#tau}"]
    color = [ROOT.kRed,ROOT.kBlue,ROOT.kTeal,ROOT.kGray]
    TArray = ROOT.TObjArray()
    for i in range(len(mc_hists)):
        if color is not None:
            mc_hists[i].SetFillColor(color[i])
        if title is not None:
            mc_hists[i].SetTitle(title[i])
        TArray.Add(mc_hists[i])

    margin = .12
    bottomFraction = .2
    overall = ROOT.TCanvas("Data/MC")
    top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()

    ratio =  hist.Clone()
    if len(plotArgs) > 0:
        oscdata = GetOscillatedHistogram(data_histogram,templates,plotArgs[0],plotArgs[1],plotArgs[2],plotArgs[3])
        MNVPLOTTER.DrawDataStackedMC(oscdata,TArray,1,"TR","Data",0,0,1001)
        MNVPLOTTER.AddChi2Label(oscdata,hist,1,"BC",.04,.15)
        ratio.Divide(oscdata,ratio)
    else:
        MNVPLOTTER.DrawDataStackedMC(histogram,TArray,1,"TR","Data",0,0,1001)
        MNVPLOTTER.AddChi2Label(histogram,hist,1,"BC",.04,.15)
        ratio.Divide(histogram,ratio)
    top.SetLogy()

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    #Now fill mcRatio with 1 for bin content and fractional error
    mcRatio = hist.GetTotalError(False, True, False) #The second "true" makes this fractional error
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

    ratio.SetMinimum(0.5)
    ratio.SetMaximum(1.5)
    ratio.Draw()

    #Error envelope for the MC
    mcRatio.SetLineColor(ROOT.kRed)
    mcRatio.SetLineWidth(3)
    mcRatio.SetMarkerStyle(0)
    mcRatio.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
    mcRatio.Draw("E2 SAME")

    #Draw a flat line at 1 for ratio of MC to itself
    straightLine = mcRatio.Clone()
    straightLine.SetFillStyle(0)
    straightLine.Draw("HIST L SAME")
    top.cd()
    title = ROOT.TPaveText(0.3, 0.91, 0.7, 1.0, "nbNDC") #no border and use Normalized Device Coordinates to place it
    title.SetFillStyle(0)
    title.SetLineColor(ROOT.kWhite)
    title.AddText("Pseudo Data Oscillation")
    title.Draw("SAME")
    
    title = ROOT.TPaveText(0.35, 0.54, 0.65, 0.72, "nbNDC") #no border and use Normalized Device Coordinates to place it
    title.SetFillStyle(0)
    title.SetLineColor(ROOT.kWhite)
    if len(plotArgs) > 0:
        title.AddText("#Delta m^{2}: "+"{:.3f}".format(plotArgs[0])+" -> "+"{:.3f}".format(m))
        title.AddText("U_{e4}: "+"{:.3f}".format(plotArgs[1])+" -> "+"{:.3f}".format(U_e4))
        title.AddText("U_{#mu 4}: "+"{:.3f}".format(plotArgs[2])+" -> "+"{:.3f}".format(U_mu4))
        title.AddText("U_{#tau 4}: "+"{:.3f}".format(plotArgs[3])+" -> "+"{:.3f}".format(U_tau4))
    else:
        title.AddText("#Delta m^{2}: "+"{:.3f}".format(m))
        title.AddText("U_{e4}: "+"{:.3f}".format(U_e4))
        title.AddText("U_{#mu 4}: "+"{:.3f}".format(U_mu4))
        title.AddText("U_{#tau 4}: "+"{:.3f}".format(U_tau4))
    title.Draw("SAME")

    MNVPLOTTER.WritePreliminary(0.4, 0.05, 5e-2, True)
    if len(plotArgs) > 0:
        overall.Print("categories_dm_{:.1f}_Ue4_{:.3f}_Umu4_{:.3f}_Utau4_{:.3f}.png".format(plotArgs[0],plotArgs[1],plotArgs[2],plotArgs[3]))
    else:
        overall.Print("categories_dm_{}_Ue4_{}_Umu4_{}_Utau4_{}.png".format('result','result','result','result'))

    return

def GetOscillatedHistogram(histogram, templates, m, U_e4, U_mu4, U_tau4,doSys = False):
    # stitched_nueselection_energy is 1 for every sample except the nueselection
    # stitched_nutau_energy is 1 only for the nueel selection
    # if 0 is False, if 1 is True

    stitched_nueTemp = templates["nue"]
    stitched_numuTemp = templates["numu"]
    stitched_swapTemp = templates["swap"]

    stitched_nue_energy = templates["nue_energy"]
    stitched_numu_energy = templates["numu_energy"]
    stitched_nutau_energy = templates["nutau_energy"]
    stitched_ratio_energy = templates["ratio_energy"]
    stitched_swap_energy = templates["swap_energy"]
    stitched_nueselection_energy = templates["nueselection_energy"]

    hist = histogram.Clone()
    weights = hist.GetCVHistoWithStatError().Clone()

    for q in range(0,hist.GetNbinsX() + 1):
        nue_sin = sin_average(q,m,stitched_nueTemp)
        numu_sin = sin_average(q,m,stitched_numuTemp)
        swap_sin = sin_average(q,m,stitched_swapTemp)

        P_ee = 1 - 4*U_e4*(1-U_e4)*nue_sin
        eeCont = P_ee * stitched_nue_energy.GetBinContent(q)

        P_mue = 4*(U_e4)*(U_mu4)*swap_sin
        mueCont = P_mue * stitched_swap_energy.GetBinContent(q) #TODO P_mue is 0 when it shouldn't be

        P_mumu = 1 - 4*U_mu4*(1-U_mu4)*numu_sin
        mumuCont = P_mumu * stitched_numu_energy.GetBinContent(q)

        P_mutau = 4*U_tau4*U_mu4*numu_sin
        mutauCont = P_mutau * stitched_numu_energy.GetBinContent(q)

        P_etau = 4*U_e4*U_tau4*nue_sin
        etauCont = P_etau * stitched_nue_energy.GetBinContent(q)

        tauCont = (etauCont + mutauCont) * stitched_nutau_energy.GetBinContent(q)

        if stitched_ratio_energy.GetBinContent(q) == 1:
            if eeCont + mueCont == 0:
                new_cont = 0
            else:
                new_cont = mumuCont/(mueCont+eeCont)

        else:
            new_cont = mumuCont*stitched_nueselection_energy.GetBinContent(q) + eeCont + mueCont + tauCont

        new_cont = new_cont if new_cont >= 0 else 0
        old_cont = histogram.GetBinContent(q)
        weight = new_cont/old_cont if old_cont > 0 else 1
        if weight == 0:
            weight = 1
        weights.SetBinContent(q,1/weight)
        weights.SetBinError(q,0)
        #hist.SetBinContent(q,new_cont)
        
        ### oscillate CV bin content along with systematic universe contents ###
        if doSys:
            errorband_names_vert = hist.GetVertErrorBandNames()
            errorband_names_lat = hist.GetLatErrorBandNames()
            for error_band in errorband_names_vert:
                errorband = hist.GetVertErrorBand(error_band)
                errorband.SetBinContent(q,errorband.GetBinContent(q)*weight)
                n_univ = hist.GetVertErrorBand(error_band).GetNHists()
                for universe in range(0, n_univ):
                    sys_bc = hist.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(q)
                    hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(q, sys_bc*weight)

            for error_band in errorband_names_lat:
                n_univ = hist.GetLatErrorBand(error_band).GetNHists()
                errorband = hist.GetLatErrorBand(error_band)
                errorband.SetBinContent(q,errorband.GetBinContent(q)*weight)
                for universe in range(0,  n_univ):
                    sys_bc = hist.GetLatErrorBand(error_band).GetHist(universe).GetBinContent(q)
                    hist.GetLatErrorBand(error_band).GetHist(universe).SetBinContent(q, sys_bc*weight)
            hist.SetBinError(q,hist.GetBinError(q)*np.sqrt(weight))

    hist.DivideSingle(hist,weights)
    return(hist)

def makeChi2Surface(templates,stitched_data,stitched_mc,outdir_surface,dodeltachi2 = False,deltam = 1,U_tau4 = 0):
    U_e4s = np.linspace(0,0.15,40)
    U_mu4s = 0.41*np.logspace(-3,0,9)
    U_mu4s[0] = 0
    surface = np.zeros(shape=[np.shape(U_e4s)[0],np.shape(U_mu4s)[0]],dtype='f')
    count = 0
    for j in range(U_e4s.shape[0]):
        U_e4 = U_e4s[j]
        for k in range(U_mu4s.shape[0]):
            count+=1
            U_mu4 = U_mu4s[k]
            if dodeltachi2:
                testHist = GetOscillatedHistogram(stitched_data, templates, deltam, U_e4, U_mu4, U_tau4)
                minchi2,_ = doFit(testHist,templates,stitched_mc)
                chi2_null = Chi2DataMC(stitched_data,stitched_mc)
                surface[j,k] = minchi2 - chi2_null
                logging.info("{:.2f}% done with deltachi2s. Current dchi2 = {:.4f}".format(100*count/(U_e4s.shape[0]*U_mu4s.shape[0]),surface[j,k]))
            else:
                testHist = GetOscillatedHistogram(stitched_mc, templates, deltam, U_e4, U_mu4, U_tau4)
                chi2_model = Chi2DataMC(stitched_data,testHist)
                chi2_null = Chi2DataMC(stitched_data,stitched_mc)
                surface[j,k] = chi2_model - chi2_null
        
    if dodeltachi2:
        np.save('{}/deltachi2_surface_m_{}_Utau4_{}.dat'.format(outdir_surface,deltam,U_tau4),surface)
    else:
        np.save('{}/chi2_surface_m_{}_Utau4_{}.dat'.format(outdir_surface,deltam,U_tau4),surface)

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
