import ROOT
import os, time, sys, math
import datetime as dt
import PlotUtils
import UnfoldUtils
from multiprocessing import Process
from tools import PlotTools,Utilities

ROOT.TH1.AddDirectory(False)
numufile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme_t_muon_incmuon.root"
nuefile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme_t_electron_incmuon.root"
numuweightedfile =  {
    "mc":"/minerva/data/users/hsu/nu_e/kin_dist_mcme1A_nx_muon_weighted_incmuon.root",
    "data":"/minerva/data/users/hsu/nu_e/kin_dist_datame1A_nx_muon_weighted_incmuon.root"
}

MU_SIGNALS = ["CCQE","CCDelta","CC2p2h","CCDIS","CCOther"]
E_SIGNALS = ["CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"]

def GetSignalHist(f,signals,hist_name):
    h = f.Get(hist_name).Clone()
    h.Reset()
    for i in signals:
        t = f.Get("{}_{}".format(hist_name,i))
        if not t:
            print (hist_name,i)
        h.Add(t)
    return h

def MakeScaleFile():
    hist_name = "Enu"
    Fout = ROOT.TFile.Open("emu_scale.root","RECREATE")
    fe = ROOT.TFile.Open(nuefile)
    ePOT = Utilities.getPOTFromFile(nuefile)
    fmu = ROOT.TFile.Open(numufile)
    muPOT = Utilities.getPOTFromFile(nuefile)
    print(ePOT,muPOT)
    he = GetSignalHist(fe,E_SIGNALS,hist_name)
    hmu = GetSignalHist(fmu,MU_SIGNALS,hist_name)
    he.AddMissingErrorBandsAndFillWithCV(hmu)
    hmu.AddMissingErrorBandsAndFillWithCV(he)
    hmu.Scale(ePOT/muPOT)
    Fout.cd()
    he.Divide(he,hmu)
    he.Write("signal_ratio")
    fe.Close()
    fmu.Close()
    Fout.Close()

def MakeCompPlot():
    def getFileAndPOT(filename):
        POT = Utilities.getPOTFromFile(filename)
        f = ROOT.TFile.Open(filename)
        return f,POT
    def Draw(mnvplotter,data_hist,mc_hist1,mc_hist2):
        data_hist.Draw("E1 X0")
        mc_hist1.SetLineColor(ROOT.kRed)
        mc_hist1.Draw("HIST E1 SAME")
        mc_hist2.SetLineColor(ROOT.kGreen)
        mc_hist2.Draw("HIST SAME")

    nueMCfile,nueMCPOT = getFileAndPOT(nuefile)
    numuMCfile,numuMCPOT = getFileAndPOT(numuweightedfile["mc"])
    numuDatafile,numuDataPOT = getFileAndPOT(numuweightedfile["data"])
    print (nueMCPOT,numuDataPOT,numuMCPOT)
    for hist_name in [ "Eavail_q3", "Enu" ,"Eavail_Lepton_Pt"]:
        he = GetSignalHist(nueMCfile,E_SIGNALS,hist_name)
        hmu = GetSignalHist(numuMCfile,MU_SIGNALS,hist_name)
        hmu.Scale(numuDataPOT/numuMCPOT)
        he.Scale(numuDataPOT/nueMCPOT)
        hmudata = numuDatafile.Get(hist_name)
        Slicer = PlotTools.Make2DSlice if he.GetDimension()==2 else (lambda hist : [hist])
        PlotTools.MakeGridPlot(Slicer,Draw,[hmudata,he,hmu])
        PlotTools.Print("{}".format(hist_name))

if __name__ == "__main__":
    MakeCompPlot()
    #MakeScaleFile()
