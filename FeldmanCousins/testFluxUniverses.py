import os
import time
import logging, sys
import ROOT
import PlotUtils
import numpy as np
from scipy import optimize, integrate
import argparse
ccnueroot = os.environ.get('CCNUEROOT')
plotutils = os.environ.get('PLOTUTILSROOT')

import math
import psutil
import multiprocessing
import threading
nthreads = 4
from array import array
from Tools.FitTools import *
from Tools.Histogram import *
from Tools.Helper import *

#insert path for modules of this package.
from tools.PlotLibrary import HistHolder

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()
MNVPLOTTER.draw_normalized_to_bin_width=False

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)
legend_text_size = .035

def UndoBinWidthNorm(histogram):
    for i in range(0,histogram.GetNbinsX()+1):
        width = histogram.GetBinWidth(i)
        cont = histogram.GetBinContent(i)
        new_cont = cont*width
        if new_cont != 0:
            histogram.SetBinContent(i,new_cont)

        for name in histogram.GetVertErrorBandNames():
            band = histogram.GetVertErrorBand(name)
            cont = band.GetBinContent(i)
            new_cont = cont*width
            if new_cont != 0:
                band.SetBinContent(i,new_cont)

            for j in range(band.GetNHists()):
                hist = band.GetHist(j)
                cont = hist.GetBinContent(i)
                new_cont = cont*width
                if new_cont != 0:
                    hist.SetBinContent(i,new_cont)

if __name__ == "__main__":
    f_numu = ROOT.TFile.Open(plotutils+'/data/flux/flux-g4numiv6-pdg14-minervame1D1M1NWeightedAve.root')
    fhc_numu = f_numu.Get("flux_E_unweighted")
    f_numu.Close()
    f_nue = ROOT.TFile.Open(plotutils+'/data/flux/flux-g4numiv6-pdg12-minervame1D1M1NWeightedAve.root')
    fhc_nue = f_nue.Get("flux_E_unweighted")
    f_nue.Close()

    UndoBinWidthNorm(fhc_nue)
    UndoBinWidthNorm(fhc_numu)

    new_bins = array('d',list(range(0,21)))
    #new_bins = array('d',[0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9,10,12.5,15,17.5,20])

    fhc_numu = fhc_numu.Rebin(len(new_bins)-1,"hnew",new_bins)
    fhc_nue = fhc_nue.Rebin(len(new_bins)-1,"hnew",new_bins)

    fhc_numu.Scale(1e7,"width")
    fhc_nue.Scale(1e7,"width")

    fhc_numu_univ = fhc_numu.GetVertErrorBand("Flux")
    fhc_nue_univ = fhc_nue.GetVertErrorBand("Flux")

    nue_arr = np.array(fhc_nue)[1:-1]
    nhists = fhc_nue_univ.GetNHists()
    flux_universes = np.array([np.array(fhc_nue_univ.GetHist(i))[1:-1] for i in range(nhists)])
    A = flux_universes - np.array([nue_arr for i in range(nhists)])

    fhc_nue.PopVertErrorBand("Flux_BeamFocus")
    fhc_nue.PopVertErrorBand("ppfx1_Total")
    fhc_nue.PopVertErrorBand("Flux")

    h_data = fhc_nue.Clone()
    for i in range(1,h_data.GetNbinsX()+1):
        excess = 2.718281**-h_data.GetXaxis().GetBinCenter(i)
        fake_data = nue_arr[i-1] + 2*nue_arr[i-1] * excess
        h_data.SetBinContent(i,fake_data)
        h_data.SetBinError(i,np.sqrt(fake_data))
        fhc_nue.SetBinContent(i,fhc_nue.GetBinContent(i))
        fhc_nue.SetBinError(i,np.sqrt(fhc_nue.GetBinContent(i)))

    data = np.array(h_data)[1:-1]
    mc = nue_arr

    C = data - mc
    ctest = ROOT.TCanvas()
    h_data.GetXaxis().SetTitle("Neutrino Energy [ GeV ]")
    h_data.Draw("hist p")

    leg = ROOT.TLegend(.4,.3)
    leg.AddEntry(h_data,"fake data","p")
    hists = []
    pens = []
    dataHist = h_data.Clone()
    dataHist.Add(fhc_nue,-1)
    V = TMatrix_to_Numpy(dataHist.GetTotalErrorMatrix(True,False,False))[1:-1,1:-1]
    for testUniverses in [10,50,100,500,1000]:
        cloned = fhc_nue.Clone()

        I = np.identity(testUniverses)
        Atest = A[:testUniverses,:]
        L = 2 * Atest @ V @ C
        lam = 1
        Q = Atest @ V @ Atest.T + I * lam
        solution = np.linalg.inv(Q) @ L/2
        penalty = solution @ solution * lam
        new_cv = mc + solution @ Atest
        for i in range(1,cloned.GetNbinsX()+1):
            cloned.SetBinContent(i,Atest[testUniverses-1][i-1])
            cloned.SetBinContent(i,new_cv[i-1])
        hists.append(cloned)
        pens.append(penalty)
        print("{} universes with penalty chi2 = {:.2f}".format(testUniverses,penalty))

    fhc_nue.SetLineColor(ROOT.kBlack)
    fhc_nue.Draw("hist same")
    leg.AddEntry(fhc_nue,"CV prediction","l")
    for i,cloned in enumerate(hists):
        cloned.SetLineColor([ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kTeal,ROOT.kGray][i])
        cloned.Draw("hist same")
        leg.AddEntry(cloned,"{} universe profile with penalty = {:.2f}".format([10,50,100,500,1000][i],pens[i]),"l")

    leg.Draw()
    ctest.SetLogy()
    ctest.Print("plots/fakedata.png")
