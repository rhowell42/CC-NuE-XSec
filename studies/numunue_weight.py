import ROOT
import os, time, sys, math
import datetime as dt
import PlotUtils
import UnfoldUtils
from multiprocessing import Process
from tools import PlotTools,Utilities


ROOT.TH1.AddDirectory(False)
numufile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme1A-E_muon_MAD.root"
nuefile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme1A-E_electron_MAD.root"
numuweightedfile =  {
    "mc":"/minerva/data/users/hsu/nu_e/kin_dist_mcme1A-E_weighted_muon_MAD.root",
    "data":"/minerva/data/users/hsu/nu_e/kin_dist_datame1A-E_weighted_muon_MAD.root"
}

MU_SIGNALS = ["CCQE","CCDelta","CC2p2h","CCDIS","CCOther"]
E_SIGNALS = ["CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"]
#E_SIGNALS = MU_SIGNALS

def repeat(f,n):
    def rf(x):
        for i in range(n):
            x = f(x)
        return x
    return rf

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
    muPOT = Utilities.getPOTFromFile(numufile)
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

def MakeScaleFile2():
    Fout = ROOT.TFile.Open("emu_scale2.root","RECREATE")
    fe = ROOT.TFile.Open(nuefile)
    ePOT = Utilities.getPOTFromFile(nuefile)
    fmu = ROOT.TFile.Open(numufile)
    muPOT = Utilities.getPOTFromFile(numufile)
    print(ePOT,muPOT)
    he = GetSignalHist(fe,E_SIGNALS,"tEnu")
    migration = fmu.Get("Enu_migration").Clone()
    he.Scale(muPOT/ePOT)
    tmp = FindScale2(he,migration)
    scale_hist = fmu.Get("Enu").Clone()
    for i in range(scale_hist.GetNbinsX()+2):
        scale_hist.SetBinContent(i,tmp[i][0])
    Fout.cd()
    scale_hist.Write("Enu")
    fe.Close()
    fmu.Close()
    Fout.Close()


def FindScale2(target,migration):
    m = ROOT.TMatrixD(migration.GetNbinsY()+2,migration.GetNbinsX()+2)
    t = ROOT.TMatrixD(target.GetNbinsX()+2,1)

    for x in range(migration.GetNbinsX()+2):
        for y in range(migration.GetNbinsY()+2):
            m[y][x]=migration.GetBinContent(x,y)

    for y in range(target.GetNbinsX()):
        t[y][0] = target.GetBinContent(y)

    SVD = ROOT.TDecompSVD(m)
    tmp = SVD.Invert()
    r = ROOT.TMatrixD(migration.GetNbinsX()+2,1)
    r.Mult(tmp,t)
    return repeat(smooth,4)(r)

def smooth(iput):
    tmp = iput[12][0]
    for i in range(13,iput.GetNrows()-1):
        new = (tmp+iput[i][0]+iput[i+1][0])/3
        tmp = iput[i][0]
        iput[i][0]=new
    iput[iput.GetNrows()-1][0]=(tmp+iput[iput.GetNrows()-1][0])/2
    print (list(iput[i][0] for i in range(iput.GetNrows())))
    return iput

def MakeCompPlot():
    def getFileAndPOT(filename):
        POT = Utilities.getPOTFromFile(filename)
        f = ROOT.TFile.Open(filename)
        return f,POT
    def Draw(mnvplotter,*args,canvas=PlotTools.CANVAS):
        mc_hist1 = args[0]
        mc_hist2 = args[1]
        data_hist = args[2] if len(args)>2 else None
        # if data_hist:
        #     data_hist.Draw("E1 X0")
        if mc_hist1:
            mc_hist1.SetTitle("nue MC")
            mc_hist1.SetLineColor(ROOT.kRed)
            mc_hist1.Draw("HIST E1 SAME")
        if mc_hist2:
            mc_hist2.SetTitle("weighted numu MC")
            mc_hist2.SetLineColor(ROOT.kGreen)
            mc_hist2.Draw("HIST SAME")

        
        leg = ROOT.TLegend(0.5,0.7,0.9,0.9).Clone()
        leg.AddEntry(mc_hist1,"nue MC")
        leg.AddEntry(mc_hist2,"weighted numu MC")
        leg.Draw()

       
    nueMCfile,nueMCPOT = getFileAndPOT(nuefile)
    numuMCfile,numuMCPOT = getFileAndPOT(numuweightedfile["mc"])
    numuDatafile,numuDataPOT = getFileAndPOT(numuweightedfile["data"])
    print (nueMCPOT,numuDataPOT,numuMCPOT)
    for hist_name in [ "Eavail_q3", "Enu" ,"Eavail_Lepton_Pt","Q3","Q0","Enu_qe","tQ0","tQ3","tLepton_Pt","tQ0_tQ3","tQ0_tLepton_Pt","Eavail","tEavail","tEavail_tQ3","tEavail_tLepton_Pt","Lepton_Pt","tEnu"]:
        try:
            he = GetSignalHist(nueMCfile,E_SIGNALS,hist_name)
            hmu = GetSignalHist(numuMCfile,MU_SIGNALS,hist_name)
        except :
            print ("hist: {} not found".format(hist_name))
            continue
        hmu.Scale(numuDataPOT/numuMCPOT)
        he.Scale(numuDataPOT/nueMCPOT)
        hmudata = None# numuDatafile.Get(hist_name)
        Slicer = PlotTools.Make2DSlice if he.GetDimension()==2 else (lambda hist : [hist])
        if hmudata:
            PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu])
        else:
            PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
        print(he.Integral(),hmu.Integral())
        PlotTools.CANVAS.Print("{}.png".format(hist_name))
        

if __name__ == "__main__":
    MakeCompPlot()
    #MakeScaleFile(["Enu","tEnu"])
    #MakeScaleFile2()
