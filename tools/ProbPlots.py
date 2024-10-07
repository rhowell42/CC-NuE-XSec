import os
import sys
import ROOT
import PlotUtils
import numpy as np

E_BINNING = PlotConfig.NEUTRINO4_EE_BINNING
E_BINNING_INV = PlotConfig.NEUTRINO4_EE_BINNING_INV
L_BINNING = PlotConfig.NEUTRINO4_LENGTH_BINNING

def TrueSignalPlot():
    m = [0,2,7.34]
    sines = [0,.15,.45]
    sin_mue = 0.002
    k = 0

    canvas = ROOT.TCanvas("truesig")
    Nx = len(m);
    Ny = len(sines);
    canvas.Divide(Nx,Ny) 
    canvas.SetLeftMargin(0)

    for i in range(Nx):
        for j in range(Ny):
            d_msqr = m[i]
            sinee = sines[j]
            print("**** working on loop "+str(k)+" out of "+str(Nx*Ny-1)+" ****")
            signalHist = reco_energy_mc_hists_truesignal.Clone()
            nueHist = true_energy_mc_hists_nue.Clone()
            numuHist = true_energy_mc_hists_numu.Clone()
            for i_E in range(len(E_BINNING)):
                E = E_BINNING[i_E]
                signal_projection = length_vs_energy_mc_hists_truesignal.ProjectionX(str(E),i_E,i_E+1)
                nue_projection = length_vs_energy_mc_hists_nue.ProjectionX(str(E),i_E,i_E+1)
                numu_projection = length_vs_energy_mc_hists_numu.ProjectionX(str(E),i_E,i_E+1)

                NE_Signal = round(reco_energy_mc_hists_truesignal.GetBinContent(i_E))
                NE_NuE = round(true_energy_mc_hists_nue.GetBinContent(i_E))
                NE_NuMu = round(true_energy_mc_hists_numu.GetBinContent(i_E))

                if NE_Signal == 0:
                  cut_effics = 0
                else:
                  cut_effics = NE_Signal/NE_NuE

                NE_NuE_cuts = round(NE_NuE*cut_effics)
                NuE_Lengths = np.array([nue_projection.GetRandom() for _ in range(NE_NuE_cuts)])
                NuMu_Lengths = np.array([numu_projection.GetRandom() for _ in range(NE_NuMu)])
                P_ee = 1 - sinee*np.sin(1.27*d_msqr*NuE_Lengths/E)**2
                P_mue = sin_mue*np.sin(1.27*d_msqr*NuMu_Lengths/E)**2
              
                remaining_nue = np.sum(np.random.uniform(0,1,NE_NuE_cuts) < P_ee)
                new_nue = np.sum(np.random.uniform(0,1,NE_NuMu) < P_mue)*cut_effics

                signalHist.SetBinContent(i_E,NE_Signal)
