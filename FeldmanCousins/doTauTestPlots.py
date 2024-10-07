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
MNVPLOTTER.mc_line_width = 0
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

rowstodo = 100

def constraint(t):
    U_e4  = t[1]
    U_mu4 = t[2]
    U_tau4= t[3]
    
    return(-1*(U_e4**2 + U_mu4**2 + U_tau4**2 - 1))

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
        chi2 = MNVPLOTTER.Chi2DataMC(stitched_data,universe)
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
        x[1] += np.random.uniform(-0.1*s, 0.1*s)
        x[2] += np.random.uniform(-0.1*s, 0.1*s)
        x[3] += 0
        return(x)

def doFit(universe, makePlot = False, plotArgs = [0,0,0,0]):
    logging.debug("    initializing list")
    x0 = [0.0,0.0,0.0,0.0]
    logging.debug("    initializing array")
    bnds = np.array([[0.0,1.0],[0.0,0.394],[0.0,0.489],[0,0.718]], dtype = float)
    cons = ({"type":"ineq","fun":constraint})
    logging.debug("    minimizing...")
    minimizer_kwargs = {"bounds":bnds,"tol":0.001, "options":{"maxiter":10},"constraints":cons, "method":"SLSQP","args":(universe,)}

    before = psutil.virtual_memory()[3]/1000000000
    mytakestep = MyTakeStep()
    rng = np.random.default_rng()
    #res = optimize.basinhopping(CalChi2, x0, take_step = mytakestep, interval=25, niter = 100, seed=rng, minimizer_kwargs=minimizer_kwargs)
    res = optimize.basinhopping(CalChi2, x0, interval=20, niter = 80, seed=rng, minimizer_kwargs=minimizer_kwargs)
    after = psutil.virtual_memory()[3]/1000000000
    logging.debug("RAM Added by CalChi2 (GB): %s",str(after-before))

    logging.debug("    returning results")
    chi2 = float(res.fun)
    if makePlot:
        fitms = res.x[0]*100
        fitUe4 = res.x[1]
        fitUmu4 = res.x[2]
        fitUtau4 = res.x[3]
        nit = res.nit
        GetOscillatedHistogramCategories(stitched_mc,fitms,fitUe4,fitUmu4,fitUtau4,plotArgs)
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
    ms = x[0]*100
    U_e4 = x[1]
    U_mu4 = x[2]
    U_tau4 = x[3]
    target = args
    oscillated = GetOscillatedHistogram(stitched_mc,ms,U_e4,U_mu4,U_tau4)
    chi2 = MNVPLOTTER.Chi2DataMC(oscillated,target)
    return(chi2)

def GetOscillatedHistogramCategories(histogram, m, U_e4, U_mu4,U_tau4,plotArgs = []):
    mc_histogram = histogram.Clone()
 
    nueCat = mc_histogram.Clone()
    anueCat = mc_histogram.Clone()
    numuCat = mc_histogram.Clone()
    anumuCat = mc_histogram.Clone()
    nueswapCat = mc_histogram.Clone()
    anueswapCat = mc_histogram.Clone()
    numuswapCat = mc_histogram.Clone()
    anumuswapCat = mc_histogram.Clone() 
    nutauswapCat = mc_histogram.Clone()
    anutauswapCat = mc_histogram.Clone() 
    selCat = mc_histogram.Clone()
    selswapCat = mc_histogram.Clone()

    for q in range(1,mc_histogram.GetNbinsX()+1):
        if q <= 12:
            nue_sin = sin_average(q,m,stitched_nueTemp)
            anue_sin = sin_average(q,m,stitched_anueTemp)
            numu_sin = sin_average(q,m,stitched_numuTemp)
            anumu_sin = sin_average(q,m,stitched_anumuTemp)

            P_ee = 1 - 4*U_e4**2*(1-U_e4**2)*nue_sin
            aP_ee = 1 - 4*U_e4**2*(1-U_e4**2)*anue_sin
            eeCont = P_ee * stitched_nue_energy.GetBinContent(q)
            a_eeCont = aP_ee * stitched_anue_energy.GetBinContent(q)

            P_mue = 4*U_e4**2*U_mu4**2*numu_sin
            aP_mue = 4*U_e4**2*U_mu4**2*anumu_sin
            mueCont = P_mue * stitched_numu_energy.GetBinContent(q)
            a_mueCont = aP_mue * stitched_anumu_energy.GetBinContent(q)

            P_mumu = 1 - 4*U_mu4**2*(1-U_mu4**2)*numu_sin
            aP_mumu = 1 - 4*U_mu4**2*(1-U_mu4**2)*anumu_sin
            mumuCont = P_mumu * stitched_numu_energy.GetBinContent(q)
            a_mumuCont = aP_mumu * stitched_anumu_energy.GetBinContent(q)

            P_mutau = 4*U_tau4**2*U_mu4**2*numu_sin
            aP_mutau = 4*U_tau4**2*U_mu4**2*anumu_sin
            mutauCont = P_mutau * stitched_numu_energy.GetBinContent(q)
            a_mutauCont = aP_mutau * stitched_anumu_energy.GetBinContent(q)

            P_etau = 4*U_e4**2*U_tau4**2*nue_sin
            aP_etau = 4*U_e4**2*U_tau4**2*anue_sin
            etauCont = P_etau * stitched_nue_energy.GetBinContent(q)
            a_etauCont = aP_etau * stitched_anue_energy.GetBinContent(q)

            nueCat.SetBinContent(q,eeCont)
            anueCat.SetBinContent(q,a_eeCont)
            numuCat.SetBinContent(q,mumuCont)
            anumuCat.SetBinContent(q,a_mumuCont)

            nueswapCat.SetBinContent(q,mueCont)
            anueswapCat.SetBinContent(q,a_mueCont)
            nutauswapCat.SetBinContent(q,mutauCont+etauCont)
            anutauswapCat.SetBinContent(q,a_mutauCont+a_etauCont)

            selCat.SetBinContent(q,0)
            selswapCat.SetBinContent(q,0)
            mc_histogram.SetBinContent(q,eeCont+a_eeCont+mumuCont+a_mumuCont+mueCont+a_mueCont+mutauCont+etauCont+a_mutauCont+a_etauCont)
        else:
            selection_sin = sin_average(q,m,stitched_nueTemp)
            swapped_sin = sin_average(q,m,stitched_numuTemp)

            P_ee = 1 - 4*U_e4**2*(1-U_e4**2)*selection_sin
            P_mue = 4*U_e4**2*U_mu4**2*swapped_sin

            nueCont = P_ee * mc_histogram.GetBinContent(q)
            numuCont = P_mue * stitched_numu_energy.GetBinContent(q)

            selCat.SetBinContent(q,nueCont)
            selswapCat.SetBinContent(q,numuCont)

            nueCat.SetBinContent(q,0)
            anueCat.SetBinContent(q,0)
            numuCat.SetBinContent(q,0)
            anumuCat.SetBinContent(q,0)
            nueswapCat.SetBinContent(q,0)
            anueswapCat.SetBinContent(q,0)
            numuswapCat.SetBinContent(q,0)
            anumuswapCat.SetBinContent(q,0) 
            nutauswapCat.SetBinContent(q,0)
            anutauswapCat.SetBinContent(q,0) 

            mc_histogram.SetBinContent(q,nueCont+numuCont)


    mc_hists = [nueCat,anueCat,numuCat,anumuCat,nueswapCat,anueswapCat,nutauswapCat,anutauswapCat,selCat,selswapCat]
    title = ["#nu_{e}","anti #nu_{e}","#nu_{#mu}","anti #nu_{#mu}","appeared #nu_{e}","appeared anti #nu_{e}","appeared #nu_{#tau}","appeared anti #nu_{#tau}","CC #nu_{e}","appeared CC #nu_{e}"]
    color = [ROOT.kRed,ROOT.kOrange,ROOT.kBlue,ROOT.kTeal,ROOT.kMagenta,ROOT.kPink,ROOT.kGray,ROOT.kBlack,ROOT.kRed,ROOT.kViolet]
    TArray = ROOT.TObjArray()
    for i in range(len(mc_hists)):
        if color is not None:
            mc_hists[i].SetFillColor(color[i])
        if title is not None:
            mc_hists[i].SetTitle(title[i])
        TArray.Add(mc_hists[i])
    #c1 = ROOT.TCanvas()
    #MNVPLOTTER.DrawDataStackedMC(stitched_data,TArray,1,"TR","Data",0,0,1001)

    margin = .12
    bottomFraction = .2
    overall = ROOT.TCanvas("Data/MC")
    top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()

    #MNVPLOTTER.DrawDataMCWithErrorBand(mc_histogram,stitched_data,1,"TR")
    ratio =  mc_histogram.Clone()
    ratio2 = stitched_mc.Clone()
    if len(plotArgs) > 0:
        oscdata = GetOscillatedHistogram(stitched_mc,plotArgs[0],plotArgs[1],plotArgs[2],plotArgs[3])
        MNVPLOTTER.DrawDataStackedMC(oscdata,TArray,1,"TR","Data",0,0,1001)
        MNVPLOTTER.AddChi2Label(oscdata,mc_histogram,1,"BC",.04,.15)
        ratio.Divide(ratio, oscdata)
        ratio2.Divide(ratio2,oscdata)
    else:
        MNVPLOTTER.DrawDataStackedMC(stitched_data,TArray,1,"TR","Data",0,0,1001)
        MNVPLOTTER.AddChi2Label(stitched_data,mc_histogram,1,"BC",.04,.15)
        ratio.Divide(ratio, stitched_data)



    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    #Now fill mcRatio with 1 for bin content and fractional error
    mcRatio = mc_histogram.GetTotalError(False, True, False) #The second "true" makes this fractional error
    for whichBin in range(1, mcRatio.GetXaxis().GetNbins()+1): 
        mcRatio.SetBinError(whichBin, max(mcRatio.GetBinContent(whichBin), 1e-9))
        mcRatio.SetBinContent(whichBin, 1)

    ratio.SetTitle("")
    ratio.SetLineColor(ROOT.kBlack)
    ratio2.SetLineColor(ROOT.kBlue)
    ratio.SetLineWidth(3)
    ratio2.SetLineWidth(3)
    ratio.SetTitleSize(0)
    ratio2.SetTitleSize(0)

    ratio.GetYaxis().SetTitle("Data / MC")
    ratio.GetYaxis().SetLabelSize(.15)
    ratio.GetYaxis().SetTitleSize(0.16)
    ratio.GetYaxis().SetTitleOffset(0.4)
    ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.

    ratio.GetXaxis().SetTitleSize(0.16)
    ratio.GetXaxis().SetTitleOffset(0.9)
    ratio.GetXaxis().SetLabelSize(.15)

    ratio.SetMinimum(0)
    ratio.SetMaximum(2)
    ratio.Draw()
    ratio2.Draw("SAME HIST L")

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
    
    title = ROOT.TPaveText(0.35, 0.5, 0.65, 0.7, "nbNDC") #no border and use Normalized Device Coordinates to place it
    title.SetFillStyle(0)
    title.SetLineColor(ROOT.kWhite)
    title.AddText("#Delta m^{2}: "+"{:.3f}".format(plotArgs[0])+" -> "+"{:.3f}".format(m))
    title.AddText("U_{e4}: "+"{:.3f}".format(plotArgs[1])+" -> "+"{:.3f}".format(U_e4))
    title.AddText("U_{#mu 4}: "+"{:.3f}".format(plotArgs[2])+" -> "+"{:.3f}".format(U_mu4))
    title.AddText("U_{#tau 4}: "+"{:.3f}".format(plotArgs[3])+" -> "+"{:.3f}".format(U_tau4))
    title.Draw("SAME")

    MNVPLOTTER.WritePreliminary(0.4, 0.05, 5e-2, True)
    overall.Print("categories_dm_{:.2f}_Ue4_{:.3f}_Umu4_{:.3f}_Utau4_{:.3f}.png".format(m,U_e4,U_mu4,U_tau4))

    return

def GetOscillatedHistogram(histogram, m, U_e4, U_mu4, U_tau4):
    mc_histogram = histogram.Clone()
    for q in range(1,mc_histogram.GetNbinsX()):
        if q <= 12:
            nue_sin = sin_average(q,m,stitched_nueTemp)
            anue_sin = sin_average(q,m,stitched_anueTemp)
            numu_sin = sin_average(q,m,stitched_numuTemp)
            anumu_sin = sin_average(q,m,stitched_anumuTemp)

            P_ee = 1 - 4*U_e4**2*(1-U_e4**2)*nue_sin
            aP_ee = 1 - 4*U_e4**2*(1-U_e4**2)*anue_sin
            eeCont = P_ee * stitched_nue_energy.GetBinContent(q)
            a_eeCont = aP_ee * stitched_anue_energy.GetBinContent(q)

            P_mue = 4*(U_e4**2)*(U_mu4**2)*numu_sin
            aP_mue = 4*(U_e4**2)*(U_mu4**2)*anumu_sin
            mueCont = P_mue * stitched_numu_energy.GetBinContent(q)
            a_mueCont = aP_mue * stitched_anumu_energy.GetBinContent(q)

            P_mumu = 1 - 4*U_mu4**2*(1-U_mu4**2)*numu_sin
            aP_mumu = 1 - 4*U_mu4**2*(1-U_mu4**2)*anumu_sin
            mumuCont = P_mumu * stitched_numu_energy.GetBinContent(q)
            a_mumuCont = aP_mumu * stitched_anumu_energy.GetBinContent(q)

            P_mutau = 4*(U_tau4**2)*(U_mu4**2)*numu_sin
            aP_mutau = 4*(U_tau4**2)*(U_mu4**2)*anumu_sin
            mutauCont = P_mutau * stitched_numu_energy.GetBinContent(q)
            a_mutauCont = aP_mutau * stitched_anumu_energy.GetBinContent(q)

            P_etau = 4*(U_e4**2)*(U_tau4**2)*nue_sin
            aP_etau = 4*(U_e4**2)*(U_tau4**2)*anue_sin
            etauCont = P_etau * stitched_nue_energy.GetBinContent(q)
            a_etauCont = aP_etau * stitched_anue_energy.GetBinContent(q)

            nueel_bincontent = eeCont+a_eeCont+mueCont+a_mueCont+mumuCont+a_mumuCont+etauCont+a_etauCont+mutauCont+a_mutauCont
            mc_histogram.SetBinContent(q,nueel_bincontent)
        else:
            selection_sin = sin_average(q,m,stitched_nueTemp)
            swapped_sin = sin_average(q,m,stitched_numuTemp)

            P_ee = 1 - 4*U_e4**2*(1-U_e4**2)*selection_sin
            P_mue = 4*(U_e4**2)*(U_mu4**2)*swapped_sin

            nueCont = P_ee * mc_histogram.GetBinContent(q)
            numuCont = P_mue * stitched_numu_energy.GetBinContent(q)
            selection_bincontent = nueCont+numuCont

            mc_histogram.SetBinContent(q,selection_bincontent)

    return(mc_histogram)

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

def makeChi2Surface(dodeltachi2 = False,deltam = 1,U_tau4 = 0):
    U_e4s = np.linspace(0,0.394,40)
    U_mu4s = np.linspace(0,0.489,9)
    surface = np.zeros(shape=[np.shape(U_e4s)[0],np.shape(U_mu4s)[0]],dtype='f')
    count = 0
    for j in range(U_e4s.shape[0]):
        U_e4 = U_e4s[j]
        for k in range(U_mu4s.shape[0]):
            count+=1
            U_mu4 = U_mu4s[k]
            testHist = GetOscillatedHistogram(stitched_mc, deltam, U_e4, U_mu4, U_tau4)
            if dodeltachi2:
                minchi2 = doFit(testHist)
                chi2 = MNVPLOTTER.Chi2DataMC(stitched_data,testHist)
                surface[j,k] = chi2 - minchi2
                logging.info("{:.2f}% done with deltachi2s. Current dchi2 = {:.4f}".format(100*count/(U_e4s.shape[0]*U_mu4s.shape[0]),surface[j,k]))
            else:
                chi2 = MNVPLOTTER.Chi2DataMC(stitched_data,testHist)
                surface[j,k] = chi2
        
    if dodeltachi2:
        logging.info("writing deltachi2 surface file for deltam2 = {}".format(m))
        #np.save("deltachi2_surface_{}.dat".format(m),surface)
        np.save('{}/deltachi2_surface_m_{}_Utau4_{}.dat'.format(outdir_surface,deltam,U_tau4),surface)
    else:
        logging.info("writing chi2 surface file for deltam2 = {}".format(m))
        np.save('{}/deltachi2_surface_m_{}_Utau4_{}.dat'.format(outdir_surface,deltam,U_tau4),surface)

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
                        default=0
    )
    parser.add_argument("-U", "--U_tau4",
                        dest = "U_tau4",
                        help="U_tau4 parameter to probe.",
                        default=0
    )
    args = parser.parse_args()
    outdir_surface = args.output_dir
    delta_m = int(args.delta_m)
    runongrid = args.grid
    U_tau4 = float(args.U_tau4)

    stitched_data = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('data_stitched')
    stitched_mc = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched')
    
    stitched_nueTemp = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('LE_template_nue')
    stitched_anueTemp = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('LE_template_anue')
    stitched_numuTemp = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('LE_template_numu')
    stitched_anumuTemp = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('LE_template_anumu')

    stitched_nue_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_nue')
    stitched_anue_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_anue')
    stitched_numu_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_numu')
    stitched_anumu_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_anumu')

    for i in range(1,stitched_mc.GetNbinsX()):
        summed = stitched_nue_energy.GetBinContent(i) + stitched_anue_energy.GetBinContent(i) + stitched_numu_energy.GetBinContent(i) + stitched_anumu_energy.GetBinContent(i)
        ratio = stitched_nue_energy.GetBinContent(i) / summed
        new_nue = stitched_mc.GetBinContent(i) * ratio
        ratio = stitched_anue_energy.GetBinContent(i) / summed
        new_anue = stitched_mc.GetBinContent(i) * ratio
        ratio = stitched_numu_energy.GetBinContent(i) / summed
        new_numu = stitched_mc.GetBinContent(i) * ratio
        ratio = stitched_anumu_energy.GetBinContent(i) / summed
        new_anumu = stitched_mc.GetBinContent(i) * ratio
        if i <= 12:
            stitched_nue_energy.SetBinContent(i,new_nue)
            stitched_anue_energy.SetBinContent(i,new_anue)
            stitched_numu_energy.SetBinContent(i,new_numu)
            stitched_anumu_energy.SetBinContent(i,new_anumu)

    for m in [10,50]:
        for U_e4 in [0.02,.15]:
            for U_mu4 in [.122,.306,.489]:
                for U_tau4 in [0,.3,.7]:
                    universe = GetOscillatedHistogram(stitched_mc,m,U_e4,U_mu4,U_tau4)
                    doFit(universe, makePlot = True, plotArgs = [m,U_e4,U_mu4,U_tau4])
    exit()

    dodelta = False
    if not runongrid: # surface plot
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
    else:
        m = delta_m / 2
        makeChi2Surface(dodelta,m,U_tau4)
