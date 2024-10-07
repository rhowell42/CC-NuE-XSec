import os
import logging, sys
import ROOT
import PlotUtils
import numpy as np
np.random.seed(0)
np.set_printoptions(precision=1)
np.set_printoptions(linewidth=1520)
np.set_printoptions(threshold=sys.maxsize)
from scipy import optimize
from scipy import integrate

import math
import psutil
import multiprocessing
import threading
nthreads = 4
from array import array

#insert path for modules of this package.
from config import PlotConfig
from config.AnalysisConfig import AnalysisConfig
from config.DrawingConfig import PLOTS_TO_MAKE,Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer
from config.SignalDef import SIGNAL_DEFINATION
from tools import Utilities,PlotTools
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
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

rowstodo = 100

def getFinalFit(m,theta_mue,j):
    k = 0
    for U_e4 in U_e4s:
        temp_normal = unoscillated_normal.Clone()
        temp_swap = unoscillated_swapped.Clone()
        testHist = GetOscillatedHistogram(temp_normal, temp_swap, m, U_e4, theta_mue)

        if (m == 0 or m == 50 or 3.3 < m < 3.35 or 7.6 < m < 7.7) and (U_e4 == 0 or theta_ee == 1 or 0.25 < theta_ee < 0.256 or 0.44 < theta_ee < 0.447):
            minchi2 = doFit(testHist,makePlot=True, plotArgs = [m,U_e4,theta_mue])
        else:
            minchi2 = doFit(testHist,makePlot=False, plotArgs = [m,U_e4,theta_mue])
        chi2 = MNVPLOTTER.Chi2DataMC(unoscillated_normal,testHist)
        deltachi2 = chi2 - minchi2
        passed = np.count_nonzero(deltachi2 < test_dchi2s)
        fraction = passed/test_dchi2s.size
        results[k,j] = fraction
        k+=1
        logging.debug("Done with {:.2f}% U_e4s".format(100*k/U_e4s.shape[0]))

def constraint(t):
    ee  = t[0]
    mue = t[1]
    return(-(4 * mue / ee) + 1)

def getFits(matrix,hist,name):
    c = 0
    dchi2s = np.zeros(rowstodo)
    logging.debug("\nI want to loop through array of size: {}".format(matrix.shape[0]))
    for row in matrix:
        universe = hist.Clone()
        for i in range(len(row)):
            before = universe.GetBinContent(i)
            scale = row[i]
            universe.SetBinContent(i,before*(1+scale))
        minchi2 = doFit(universe)
        chi2 = MNVPLOTTER.Chi2DataMC(unoscillated_normal.Clone(),universe.Clone())
        dchi2 = chi2 - minchi2
        dchi2s[c] = dchi2
        c+=1
        logging.debug("done with row {}".format(c))
    logging.info("writing chi2_{}.txt to disk".format(name))
    np.savetxt('chi2s/chi2_{}.txt'.format(name),dchi2s)

class MyTakeStep(object):
    def __init__(self, stepsize=0.5):
        self.stepsize = stepsize
    def __call__(self, x):
        s = self.stepsize
        x[0] += np.random.uniform(-2.0*s, 2.0*s)
        x[1] += np.random.uniform(-0.01*s, 0.01*s)
        x[2] += np.random.uniform(-2.0*s, 2.0*s)
        return(x)

def doFit(universe, makePlot = False, plotArgs = [0,0,0]):
    logging.debug("    initializing list")
    x0 = [0.0,0.0,0.0]
    logging.debug("    initializing array")
    bnds = np.array([[0.0,1.0],[0.0,1.0],[0.0,1.0]], dtype = float)
    cons = {"type":"ineq","fun":constraint}
    logging.debug("    minimizing...")
    minimizer_kwargs = {"bounds":bnds,"tol":0.1, "options":{"maxiter":100}, "method":"SLSQP","args":(universe,)}

    before = psutil.virtual_memory()[3]/1000000000
    mytakestep = MyTakeStep()
    rng = np.random.default_rng()
    res = optimize.basinhopping(CalChi2, x0, take_step = mytakestep, interval=25, niter = 200, seed=rng, minimizer_kwargs=minimizer_kwargs)
    after = psutil.virtual_memory()[3]/1000000000
    logging.debug("RAM Added by CalChi2 (GB): {}".format(after-before))
    logging.debug('RAM Used (GB):', psutil.virtual_memory()[3]/1000000000)

    logging.debug("    returning results")
    chi2 = float(res.fun)
    if makePlot:
        fitee = res.x[0]
        fitmue = res.x[1]
        fitms = res.x[2]*100
        nit = res.nit
        temp_normal = unoscillated_normal.Clone()
        temp_swap = unoscillated_swapped.Clone()
        ROOT.SetOwnership(temp_normal, True)
        ROOT.SetOwnership(temp_swap, True)
        fitHist = GetOscillatedHistogram(temp_normal,temp_swap,fitms,fitee,fitmue)
        ROOT.SetOwnership(fitHist, True)
        oscillatedC = ROOT.TCanvas("oscillated")
        MNVPLOTTER.DrawDataMC(fitHist,universe)
        logging.info(res.message)
        oscillatedC.Print("fitm_{:.3f}_testm_{:.3f}__fitee_{:.3f}_testee_{:.3f}__fitmue_{:.3f}_testmue_{:.3f}_iter_{}.png".format(fitms,plotArgs[0],fitee,plotArgs[1],fitmue,plotArgs[2],nit))
    return(chi2)

def Chi2DataMC(dataHist,mcHist,useDataErrorMatrix = False, useOnlyShapeErrors = False, useModelStat=True):
    #We get the number of bins and make sure it's compatible with the NxN matrix given
    if dataHist.GetNbinsX() != mcHist.GetNbinsX():
        logging.error("MnvPlotter::Chi2DataMC", "The number of bins from Data and MC histograms differ. Returning -1.")
        return(-1)

    #only consider the plotted range
    lowBin  = mcHist.GetXaxis().GetFirst()
    highBin = mcHist.GetXaxis().GetLast()

    #Scaling MC to Data
    tmpMCHist = mcHist.Clone("tmpMCHist")
    tmpMCHist.Scale(1.0)

    # Defining Error Matrix dimensions
    # Either use the requested range or the full error matrix with under/overflow
    Nbins = highBin - lowBin + 1;
    if MNVPLOTTER.chi2_use_overflow_err:
        Nbins = tmpMCHist.GetNbinsX()+2

    #get the covariance matrix
    covMatrix = np.zeros(shape=[Nbins,Nbins],dtype='f')
    covMatrixTmp = None
    NbinsTotal = tmpMCHist.GetNbinsX() + 2 #+2 for under/overflow
    includeStatError = True
    errorAsFraction  = False

    # Use Error Matrix from Data or MC?
    if useDataErrorMatrix:
        covMatrixTmp = dataHist.GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)
        if useModelStat:
            covMatrixTmp+=tmpMCHist.GetStatErrorMatrix()

    else:
        covMatrixTmp = tmpMCHist.GetTotalErrorMatrix(includeStatError, errorAsFraction, useOnlyShapeErrors)
        if useModelStat:
            covMatrixTmp+=dataHist.GetStatErrorMatrix()

    if MNVPLOTTER.chi2_use_overflow_err:
        covMatrix = covMatrixTmp
    else:
        for i in range(0,Nbins):
            for j in range(0,Nbins):
                covMatrix[i][j] = covMatrixTmp[i+lowBin][j+lowBin];

    covMatrix*=1e20
    #errors,S,Vh = np.linalg.svd(covMatrix)
    #errorMatrix = np.linalg.inv(errors)
    errorMatrix = np.linalg.pinv(covMatrix)
    errorMatrix*=1e20

    chi2 = 0.
    # under/overflow bins not taken into account in the chi2 calculation
    for i in range(lowBin,highBin):
        iErrBin = i - lowBin;
        x_data_i = dataHist.GetBinContent(i);
        x_mc_i   = tmpMCHist.GetBinContent(i);

        for j in range(lowBin,highBin):
            jErrBin = j - lowBin
            x_data_j = dataHist.GetBinContent(j)
            x_mc_j   = tmpMCHist.GetBinContent(j)
            chi2_ij = (x_data_i - x_mc_i) * errorMatrix[iErrBin][jErrBin] * (x_data_j - x_mc_j)
            chi2 += chi2_ij 

    return(chi2)

def CalChi2(x,args):
    U_e4 = x[0]
    theta_numu = x[1]
    ms = x[2]*100
    target = args
    temp_normal = unoscillated_normal.Clone()
    temp_swap = unoscillated_swapped.Clone()
    oscillated = GetOscillatedHistogram(temp_normal,temp_swap,ms,U_e4,theta_numu)
    chi2 = MNVPLOTTER.Chi2DataMC(oscillated,target)
    return(chi2)

def GetOscillatedHistogram(nue_hist, numu_hist, m, U_e4, U_mu4):
    for q in range(nue_hist.GetNbinsX()):
        avgsin = sin_average(q,m,rhc_swap_sel_template)
        avgsinSwapped = sin_average(q,m,rhc_swap_sel_template)
        P_ee = 1 - 4*U_e4*(1-U_e4)*avgsin
        P_mue = 4*U_e4*U_mu4*avgsinSwapped
        e_before = nue_hist.GetBinContent(q)
        mu_before = numu_hist.GetBinContent(q)
        left_nue = P_ee * e_before
        new_nue = P_mue * mu_before
        combined = left_nue+new_nue
        if e_before > 0 and mu_before > 0:
            newerr = ((nue_hist.GetBinError(q)*combined/e_before)**2 + (combined**.5)**2)**.5
            appearednewerr = ((numu_hist.GetBinError(q)*new_nue/mu_before)**2 + (new_nue**.5)**2)**.5
        else:
            newerr = ((nue_hist.GetBinError(q))**2 + (combined**.5)**2)**.5 
            appearednewerr = ((numu_hist.GetBinError(q))**2 + (new_nue**.5)**2)**.5 

        nue_hist.SetBinContent(q,combined)
        numu_hist.SetBinContent(q,new_nue)

        nue_hist.SetBinError(q,newerr)
        numu_hist.SetBinError(q,appearednewerr)
    return(nue_hist)

def GetScaleMatrix(hist):
    lowBin = hist.GetXaxis().GetFirst()
    highBin = hist.GetXaxis().GetLast()
    NbinsTotal = highBin - lowBin + 1
    covMatrix = np.zeros(shape=[NbinsTotal,NbinsTotal],dtype='f')

    for i in range(0,NbinsTotal):
        for j in range(0,NbinsTotal):
            covMatrix[i][j] = covMatrixTmp[1+i+lowBin][1+j+lowBin]
    chol = np.linalg.cholesky(covMatrix)

    exit()
    
    N = 1000 # number of pseudo experiments
    scales = np.random.normal(loc=0,scale=1,size=[N,NbinsTotal]).astype('f')
    means = np.mean(scales,axis=0,dtype=np.float32)
    for i in range(NbinsTotal):
        scales[:,i] = scales[:,i] - means[i]

    scaleCovars = np.zeros(shape=[NbinsTotal,NbinsTotal])
    for i in range(NbinsTotal):
        for j in range(NbinsTotal):
            covar = 0
            for n in range(N):
                covar += scales[n][i] * scales[n][j]
            scaleCovars[i][j] = covar/N

    toInvert = np.linalg.cholesky(scaleCovars)
    inverse = np.linalg.inv(toInvert)
    scales = np.dot(scales,inverse)
    scales = np.dot(scales,chol)
    logging.debug("   ***   Returning {} universe scales   ***   ".format(N))
    return(scales) 

def makeChi2Surface(dodeltachi2 = False,deltams = 10):
    surface = np.zeros(shape=[100,9],dtype='f')
    U_e4s = np.linspace(0,1,100)
    theta_mues = np.linspace(0,0.025,9)
    ms = deltams
    for j in range(100):
        U_e4 = U_e4s[j]
        for k in range(9):
            theta_mue = theta_mues[k]
            temp_normal = unoscillated_normal.Clone()
            temp_swap = unoscillated_swapped.Clone()
            testHist = GetOscillatedHistogram(temp_normal, temp_swap, m, U_e4, theta_mue)
            if dodeltachi2:
                minchi2 = doFit(testHist)
                chi2 = MNVPLOTTER.Chi2DataMC(unoscillated_normal,testHist)
                surface[j,k] = chi2 - minchi2
            else:
                chi2 = MNVPLOTTER.Chi2DataMC(unoscillated_normal,testHist)
                surface[j,k] = chi2
        print("{}% done with the U_e4s".format(100*(1+j)/100))
        
    if dodelta:
        logging.info("writing deltachi2 surface file for deltam2 = {}".format(m))
        np.save("deltachi2_surface_{}.dat".format(m),surface)
    else:
        logging.info("writing chi2 surface file for deltam2 = {}".format(m))
        np.save("chi2_surface_{}.dat".format(f),surface)

def integrand(x,beta,N_bin,bin_width):
    return(np.sin(1.27*beta*x)**2 * N_bin)

def makeOscTemplates():
    nueTemplate = rhc_sel_template.Clone().ProjectionY()
    numuTemplate = rhc_swap_sel_template.Clone().ProjectionY()

    avgSwap0 = nueTemplate.Clone()
    avgHist0 = nueTemplate.Clone()

    avgSwap1 = nueTemplate.Clone()
    avgHist1 = nueTemplate.Clone()

    avgSwap2 = nueTemplate.Clone()
    avgHist2 = nueTemplate.Clone()

    avgSwap3 = nueTemplate.Clone()
    avgHist3 = nueTemplate.Clone()

    avgSwap4 = nueTemplate.Clone()
    avgHist4 = nueTemplate.Clone()

    avgSwap5 = nueTemplate.Clone()
    avgHist5 = nueTemplate.Clone()

    avgSwap6 = nueTemplate.Clone()
    avgHist6 = nueTemplate.Clone()

    avgSwap7 = nueTemplate.Clone()
    avgHist7 = nueTemplate.Clone()

    avgSwap8 = nueTemplate.Clone()
    avgHist8 = nueTemplate.Clone()

    avgSwap9 = nueTemplate.Clone()
    avgHist9 = nueTemplate.Clone()

    avgSwap10 = nueTemplate.Clone()
    avgHist10 = nueTemplate.Clone()

    c1 = ROOT.TCanvas()
    rhc_sel_template.Draw("colz")
    c1.Print("sel_template.png")
    rhc_swap_sel_template.Draw("colz")
    c1.Print("sel_swap_template.png")

    for q in range(rhc_sel_template.GetNbinsY()+1):
        m = 100
        avgsin = sin_average(q,m,rhc_sel_template)    
        avgHist1.SetBinContent(q,avgsin)
        avgHist1.SetBinError(q,0)
        avgsin = sin_average(q,m,rhc_swap_sel_template)    
        avgSwap1.SetBinContent(q,avgsin)
        avgSwap1.SetBinError(q,0)

        m = 20
        avgsin = sin_average(q,m,rhc_sel_template)    
        avgHist2.SetBinContent(q,avgsin)
        avgHist2.SetBinError(q,0)
        avgsin = sin_average(q,m,rhc_swap_sel_template)    
        avgSwap2.SetBinContent(q,avgsin)
        avgSwap2.SetBinError(q,0)

        m = 25
        avgsin = sin_average(q,m,rhc_sel_template)    
        avgHist3.SetBinContent(q,avgsin)
        avgHist3.SetBinError(q,0)
        avgsin = sin_average(q,m,rhc_swap_sel_template)    
        avgSwap3.SetBinContent(q,avgsin)
        avgSwap3.SetBinError(q,0)

        m = 50
        avgsin = sin_average(q,m,rhc_sel_template)    
        avgHist4.SetBinContent(q,avgsin)
        avgHist4.SetBinError(q,0)
        avgsin = sin_average(q,m,rhc_swap_sel_template)    
        avgSwap4.SetBinContent(q,avgsin)
        avgSwap4.SetBinError(q,0)

        m = 100
        avgsin = sin_average(q,m,rhc_sel_template)    
        avgHist5.SetBinContent(q,avgsin)
        avgHist5.SetBinError(q,0)
        avgsin = sin_average(q,m,rhc_swap_sel_template)    
        avgSwap5.SetBinContent(q,avgsin)
        avgSwap5.SetBinError(q,0)

    #ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)

    avgSwap1.SetLineWidth(3)
    avgSwap1.SetLineColor(ROOT.kBlack)
    avgSwap1.GetXaxis().SetTitle("Energy Estimator [ GeV ]")
    avgSwap1.GetYaxis().SetTitle("sin^{2}(#Delta m^{2} L_{#nu_{e}}/E)")
    avgSwap1.GetYaxis().SetRangeUser(0,1)

    avgSwap2.SetLineWidth(3)
    avgSwap3.SetLineWidth(3)
    avgSwap4.SetLineWidth(3)
    avgSwap5.SetLineWidth(3)
    avgSwap2.SetLineColor(ROOT.kRed)
    avgSwap3.SetLineColor(ROOT.kGreen)
    avgSwap4.SetLineColor(ROOT.kBlue)
    avgSwap5.SetLineColor(6)
    avgHist1.SetLineWidth(3)
    avgHist2.SetLineWidth(3)
    avgHist3.SetLineWidth(3)
    avgHist4.SetLineWidth(3)
    avgHist5.SetLineWidth(3)
    avgHist1.SetLineColor(ROOT.kBlack)
    avgHist1.GetXaxis().SetTitle("Energy Estimator [ GeV ]")
    avgHist1.GetYaxis().SetTitle("sin^{2}(1.27 #Delta m^{2} L/E)")
    avgHist1.GetYaxis().SetRangeUser(0,1)
    avgHist2.SetLineColor(ROOT.kRed)
    avgHist3.SetLineColor(ROOT.kGreen)
    avgHist4.SetLineColor(ROOT.kBlue)
    avgHist5.SetLineColor(6)
    avgSwap1.SetLineStyle(2)
    avgSwap2.SetLineStyle(2)
    avgSwap3.SetLineStyle(2)
    avgSwap4.SetLineStyle(2)
    avgSwap5.SetLineStyle(2)

    TLegend = ROOT.TLegend(0.5,0.75,0.75,0.9)
    #TLegend.SetNColumns(2)
    TLegend.SetHeader("#Delta m^{2} = 100 eV^{2}","C")
    TLegend.AddEntry(avgSwap1,"#nu_{#mu}")
    #TLegend.AddEntry(avgSwap2,"m = 20 eV^{2}")
    #TLegend.AddEntry(avgSwap3,"m = 25 eV^{2}")
    #TLegend.AddEntry(avgSwap4,"m = 50 eV^{2}")
    #TLegend.AddEntry(avgSwap5,"m = 1000 eV^{2}")
    TLegend.AddEntry(avgHist1,"#nu_{e}")
    #TLegend.AddEntry(avgHist2,"m = 20 eV^{2}")
    #TLegend.AddEntry(avgHist3,"m = 25 eV^{2}")
    #TLegend.AddEntry(avgHist4,"m = 50 eV^{2}")
    #TLegend.AddEntry(avgHist5,"m = 1000 eV^{2}")

    canv = ROOT.TCanvas("avgSwap")
    avgHist1.DrawClone("hist L")
    #avgHist2.DrawClone("SAME L hist ")
    #avgHist3.DrawClone("SAME L hist ")
    #avgHist4.DrawClone("SAME L hist ")
    #avgHist5.DrawClone("SAME L hist ")
    avgSwap1.Draw("SAME L hist")
    #avgSwap2.Draw("SAME L hist")
    #avgSwap3.Draw("SAME L hist ")
    #avgSwap4.Draw("SAME L hist ")
    #avgSwap5.Draw("SAME L hist ")
    TLegend.Draw("SAME")
    canv.Print("avg_sin_comparison_100.png")

if __name__ == "__main__":
    #input knobs
    appeared_playlist = "me5A_swap"
    AnalysisConfig.bkgTune_tag = "N4_tune"
    AnalysisConfig.selection_tag = "IMD"
    AnalysisConfig.playlist = "RHC_Selection"
    AnalysisConfig.ntuple_tag = "MAD"
    AnalysisConfig.truth = True
    AnalysisConfig.mc_only = True
    #playlist=AnalysisConfig.playlist
    playlist = "RHC_Selection"
    bkg_file_path = AnalysisConfig.BackgroundFitPath(playlist,AnalysisConfig.bkgTune_tag,False)

    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    appeared_type_path_map = { t:"/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcme5A_swap_flavorswapped_fspline.root" for t in AnalysisConfig.data_types}

    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    appeared_data_file,appeared_mc_file,appeared_pot_scale,appeared_data_pot,appeared_mc_pot = Utilities.getFilesAndPOTScale(appeared_playlist,appeared_type_path_map,AnalysisConfig.ntuple_tag,True)

    bkg = ROOT.TFile(bkg_file_path)
    ROOT.SetOwnership(bkg, True)

    standPOT = data_pot if data_pot is not None else mc_pot 
    appeared_standPOT = appeared_data_pot if appeared_data_pot is not None else appeared_mc_pot 

    rhc_sel = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e/kin_dist_mcRHC_Selection_Genie_thesis_MAD.root")
    rhc_sel_template = HistHolder("Reco Energy vs L/E",rhc_sel,"Signal",True,mc_pot,standPOT)
    rhc_sel_swap = ROOT.TFile.Open("/exp/minerva/data/users/rhowell/antinu_e_swap/kin_dist_mcRHC_Selection_1000Univ_thesis_swap_MAD.root")
    rhc_swap_sel_template = HistHolder("Reco Energy vs L/E",rhc_sel_swap,"Signal",True,appeared_mc_pot,appeared_standPOT)

    Default_Scale(rhc_sel_template)
    Default_Scale(rhc_swap_sel_template)

    rhc_sel_template = rhc_sel_template.hists["Total"]
    rhc_swap_sel_template = rhc_swap_sel_template.hists["Total"]

    if False:
        L_OVER_E = HistHolder("Reco Energy vs L/E",mc_file,"Signal",True,mc_pot,standPOT)
        appeared_L_OVER_E = HistHolder("Reco Energy vs L/E",appeared_mc_file,"Signal",True,appeared_mc_pot,appeared_standPOT)
        mcsignal = HistHolder("Background Subbed MC",bkg,"Signal",True,mc_pot,standPOT)
        appeared_mcsignal = HistHolder("Biased Neutrino Energy",appeared_mc_file,"Signal",True,appeared_mc_pot,appeared_standPOT)
        datasignal = HistHolder("Background Subbed Data",bkg,"Signal",False,data_pot,standPOT)

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

        mc = mcsignal.GetHist().Clone()
        data = datasignal.GetHist().Clone()

        unoscillated_normal = mc.Clone()
        unoscillated_swapped = appeared_signal.Clone()
        debugs   = data.Clone()
        debugs.PopVertErrorBand("SuSA_Valencia_Weight")
        debugs.PopVertErrorBand("MK_model")
        debugs.PopVertErrorBand("LowQ2Pi_None")
        unoscillated_normal.PopVertErrorBand("SuSA_Valencia_Weight")
        unoscillated_normal.PopVertErrorBand("MK_model")
        unoscillated_normal.PopVertErrorBand("LowQ2Pi_None")
        unoscillated_swapped.PopVertErrorBand("SuSA_Valencia_Weight")
        unoscillated_swapped.PopVertErrorBand("MK_model")
        unoscillated_swapped.PopVertErrorBand("LowQ2Pi_None")
        
        pseudodata = unoscillated_normal.Clone()
    
    plotRatios = True
    if plotRatios:
        makeOscTemplates()
        exit()

    makeSurfaces = False
    if makeSurfaces: # surface plot
        dodelta = True
        for m in range(1,100,3):
            logging.info("initializing multiprocess threads")
            t1 = multiprocessing.Process(target=makeChi2Surface, args=(dodelta,m), name='t1')
            t2 = multiprocessing.Process(target=makeChi2Surface, args=(dodelta,m+1), name='t2')
            t3 = multiprocessing.Process(target=makeChi2Surface, args=(dodelta,m+2), name='t3')

            logging.info("starting multiprocess threads")
            t1.start()
            t2.start()
            t3.start()

            logging.info("joining multiprocess threads")
            t1.join()
            t2.join()
            t3.join()
        exit()

    if False: # test fitter
        temp_normal = unoscillated_normal.Clone()
        temp_swap = unoscillated_swapped.Clone()
        for m in [5,7,50]:
            for U_e4 in [0.1,0.4,0.8]:
                for theta_mue in [0.005]:
                    testHist = GetOscillatedHistogram(temp_normal, temp_swap, m, U_e4, theta_mue)
                    minchi2 = doFit(testHist,makePlot=True, plotArgs = [m,U_e4,theta_mue])
        exit()

    try:
        logging.info("loading scale matrices...")
        scale_matrix = np.loadtxt("ScaleMatrix_RHC_0.txt",dtype=np.float32)
        for i in range(1,10):
            matrix = np.loadtxt("ScaleMatrix_RHC_{}.txt".format(i),dtype=np.float32)
            scale_matrix = np.concatenate((scale_matrix,matrix))
    except:
        logging.info("Couldn't find scale matrices, making them now.")
        for i in range(10):
            scale_matrix = GetScaleMatrix(pseudodata)
            np.savetxt('ScaleMatrix_RHC_{}.txt'.format(i),scale_matrix)
        logging.info("Made all scale matrices, exiting now.")
        exit()

    num_rows, num_cols = scale_matrix.shape
    logging.info("Got scale matrix with shape {} x {}".format(num_rows,num_cols))

    makeChi2s = True
    if makeChi2s:
        logging.info("Making chi2 text files...")
        for p in range(0,num_rows,rowstodo*nthreads):
            upto = p+rowstodo
            
            print("going to loop through: {} -> {}".format(p,upto))
            print("going to loop through: {} -> {}".format(upto,upto+rowstodo))
            print("going to loop through: {} -> {}".format(upto+rowstodo,upto+2*rowstodo))
            print("going to loop through: {} -> {}".format(upto+2*rowstodo,upto+3*rowstodo))
            t1 = multiprocessing.Process(target=getFits, args=(scale_matrix[p:upto,:],pseudodata,p), name='t1')
            t2 = multiprocessing.Process(target=getFits, args=(scale_matrix[upto:upto+rowstodo,:],pseudodata,upto), name='t2')
            t3 = multiprocessing.Process(target=getFits, args=(scale_matrix[upto+rowstodo:upto+2*rowstodo,:],pseudodata,upto+rowstodo), name='t3')
            t4 = multiprocessing.Process(target=getFits, args=(scale_matrix[upto+2*rowstodo:upto+3*rowstodo,:],pseudodata,upto+rowstodo*2), name='t4')

            t1.start()
            t2.start()
            t3.start()
            t4.start()

            t1.join()
            t2.join()
            t3.join()
            t4.join()

        logging.info("Done making chi2 text files, exiting...")
        exit()

    deltams = np.linspace(0,10,100)
    deltams = np.append(deltams,[15,20,50,100])
    
    U_e4s = np.linspace(0,1,100)
    theta_mues = np.linspace(0,0.025,9)

    test_dchi2s = np.loadtxt("chi2s/chi2_{}.txt".format(rowstodo))
    for f in range(rowstodo,num_rows,rowstodo):
        try:
            f_chi2 = np.loadtxt("chi2s/chi2_{}.txt".format(f))
            test_dchi2s = np.concatenate((test_dchi2s,f_chi2))
        except:
            logging.info("Couldn't find chi2s/chi2_{}.txt for some reason... skipping.".format(f))

    i = 0
    for theta_mue in theta_mues:
        j = 0
        results = np.zeros((U_e4s.shape[0],deltams.shape[0]))
        for m in deltams:
            getFinalFit(m,theta_mue,j)
            j+=1
        np.save("feldman_cousins_thetamue_{}.dat".format(theta_mue),results)
        i+=1
        logging.info("{}% done calculating Feldman Cousin Fractions".format(100*i/theta_mues.size))


    if False: # plot things
        pseudodata.SetLineWidth(2)
        pseudodata.SetLineColor(ROOT.kBlack)

        hs = ROOT.THStack("hs","")
        TLegend = ROOT.TLegend(0.5,0.60,0.8,0.750)
        TLegend.AddEntry(pseudodata,"No Oscillation Pseudodata")

        histN4 = GetOscillatedHistogram(unoscillated_normal, unoscillated_swapped, 7, .44, .002)
        count = 0
        for universe in scale_matrix:
            randPseudo = pseudodata.Clone()
            for i in range(len(universe)):
                before = randPseudo.GetBinContent(i)
                scale = universe[i]
                randPseudo.SetBinContent(i,before*(1+scale))
            randPseudo.SetLineWidth(2)
            hs.Add(randPseudo)
            count += 1
            if count == 1000:
                break

        hs.Draw("nostack")
        histN4.SetLineWidth(2)
        histN4.SetLineColor(ROOT.kRed)
        unoscillated_normal.SetLineWidth(2)
        unoscillated_normal.SetLineColor(ROOT.kBlue)
        unoscillated_normal.Draw("SAME")
        histN4.Draw("SAME")
        TLegend.AddEntry(histN4,"Neutrino-4 Hypothesis")
        TLegend.Draw("SAME")

        canv.Print("test.png")

        canv = ROOT.TCanvas("err")
        MNVPLOTTER.DrawDataMC(data,debugs)
        canv.Print("test_debugs.png")
