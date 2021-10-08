import ROOT
import os, time, sys, math
import datetime as dt
import PlotUtils
import UnfoldUtils
from multiprocessing import Process

ROOT.TH1.AddDirectory(False)
numufile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme1A_nx_muon_incmuon.root"
nuefile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme1A_nx_electron_incmuon.root"

hist_name = "Enu"

MU_SIGNALS = ["CCQE","CCDelta","CC2p2h","CCDIS","CCOther"]
E_SIGNALS = ["CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"]

def GetSignalHist(f,signals):
    h = f.Get(hist_name).Clone()
    h.Reset()
    for i in signals:
        t = f.Get("{}_{}".format(hist_name,i))
        if not t:
            print (hist_name,i)
        h.Add(t)
    return h

if __name__ == "__main__":
    Fout = ROOT.TFile.Open("emu_scale.root","RECREATE")
    fe = ROOT.TFile.Open(nuefile)
    fmu = ROOT.TFile.Open(numufile)
    he = GetSignalHist(fe,E_SIGNALS)
    hmu = GetSignalHist(fmu,MU_SIGNALS)
    he.AddMissingErrorBandsAndFillWithCV(hmu)
    hmu.AddMissingErrorBandsAndFillWithCV(he)
    Fout.cd()
    he.Divide(he,hmu)
    he.Write("signal_ratio")
    fe.Close()
    fmu.Close()
    Fout.Close()
