import ROOT
import os, time, sys, math
import datetime as dt
import PlotUtils
import UnfoldUtils
from multiprocessing import Process
from tools import PlotTools,Utilities

ROOT.TH1.AddDirectory(False)
numufile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme_pre_muon_MAD.root"
nuefile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme_electron_MAD.root"
numuweightedfile =  {
    "mc":"/minerva/data/users/hsu/nu_e/kin_dist_mcme_muon_MAD.root",
    "data":"/minerva/data/users/hsu/nu_e/kin_dist_datame_muon_MAD.root"
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

def MakeScaleFile(hist_names):
    Fout = ROOT.TFile.Open("emu_scale.root","RECREATE")
    fe = ROOT.TFile.Open(nuefile)
    ePOT = Utilities.getPOTFromFile(nuefile)
    fmu = ROOT.TFile.Open(numufile)
    muPOT = Utilities.getPOTFromFile(nuefile)
    print(ePOT,muPOT)
    for hist_name in hist_names:
        he = GetSignalHist(fe,E_SIGNALS,hist_name)
        hmu = GetSignalHist(fmu,MU_SIGNALS,hist_name)
        he.AddMissingErrorBandsAndFillWithCV(hmu)
        hmu.AddMissingErrorBandsAndFillWithCV(he)
        hmu.Scale(ePOT/muPOT)
        Fout.cd()
        he.Divide(he,hmu)
        he.Write(hist_name)
    fe.Close()
    fmu.Close()
    Fout.Close()

def MakeCompPlot():
    def getFileAndPOT(filename):
        POT = Utilities.getPOTFromFile(filename)
        f = ROOT.TFile.Open(filename)
        return f,POT
    def Draw(mnvplotter,*args):
        mc_hist1 = args[0]
        mc_hist2 = args[1]
        data_hist = args[2] if len(args)>2 else None
        if data_hist:
            data_hist.Draw("E1 X0")
        if mc_hist2:
            mc_hist2.SetLineColor(ROOT.kGreen)
            mc_hist2.Draw("HIST SAME")
        if mc_hist1:
            mc_hist1.SetLineColor(ROOT.kRed)
            mc_hist1.Draw("HIST E1 SAME")

    nueMCfile,nueMCPOT = getFileAndPOT(nuefile)
    numuMCfile,numuMCPOT = getFileAndPOT(numuweightedfile["mc"])
    numuDatafile,numuDataPOT = getFileAndPOT(numuweightedfile["data"])
    print (nueMCPOT,numuDataPOT,numuMCPOT)
    for hist_name in [ "Eavail_q3", "Enu" ,"Eavail_Lepton_Pt","tEnu", "Enu_qe"]:
        he = GetSignalHist(nueMCfile,E_SIGNALS,hist_name)
        hmu = GetSignalHist(numuMCfile,MU_SIGNALS,hist_name)
        hmu.Scale(numuDataPOT/numuMCPOT)
        he.Scale(numuDataPOT/nueMCPOT)
        hmudata = numuDatafile.Get(hist_name)
        Slicer = PlotTools.Make2DSlice if he.GetDimension()==2 else (lambda hist : [hist])
        if hmudata:
            PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu,hmudata])
        else:
            PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu])
        PlotTools.Print("{}".format(hist_name))

if __name__ == "__main__":
    MakeCompPlot()
    #MakeScaleFile(["Enu","tEnu"])
