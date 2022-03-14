import ROOT
import os, time, sys, math
import datetime as dt
import PlotUtils
import UnfoldUtils
from multiprocessing import Process
from tools import PlotTools,Utilities


ROOT.TH1.AddDirectory(False)
numufile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme1LMOP_muon_MAD.root"
nuefile = "/minerva/data/users/hsu/nu_e/kin_dist_mcme1LMOP_electron_MAD.root"
nuebkgtunedfile =  "/minerva/data/users/hsu/nu_e/kin_dist_mcme1D_nx_electron_MAD.root"
nuedatafile = "/minerva/data/users/hsu/nu_e/kin_dist_datame1D_nx_electron_MAD.root"
numuweightedfile =  {
    "mc":"/minerva/data/users/hsu/nu_e/kin_dist_mcme1D_nx_weighted_muon_MAD.root",
    "data":"/minerva/data/users/hsu/nu_e/kin_dist_datame1D_nx_weighted_muon_MAD.root"
}

MU_SIGNALS = ["CCQE","CCDelta","CC2p2h","CCDIS","CCOther"]
E_SIGNALS = ["CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"]
E_BACKGROUNDS = ["NuEElastic","NonFiducial","NonPhaseSpace","CCNuEAntiNu",
                 "ExcessModel","CCDIS","CCOther","NCCOH","NCDIS","NCRES","NCOther"]
MU_BACKGROUNDS = ["NC","CCNuE","CCAntiNuMu","NonFiducial","NonPhaseSpace","Other"]
PLOTPATH = "/minerva/data/users/hsu/numunueRatioPlots/"

SIGNAL_CHANNEL_PAIR = [
    ("CCQE","CCNuEQE"),
    ("CCDelta","CCNuEDelta"),
    ("CCDIS","CCNuEDIS"),
    ("CCOther","CCNuE")
]

#E_SIGNALS = MU_SIGNALS

def SubtractPoissonHistograms(h,h1):
    h.AddMissingErrorBandsAndFillWithCV(h1)
    errors = []
    for i in range(h.GetSize()):
        errors.append(math.sqrt(h.GetBinError(i)**2 + h1.GetBinError(i)**2))
    h.Add(h1,-1)
    for i in range(h.GetSize()):
        if errors[i]==0:
            continue
        h.SetBinError(i,errors[i])
    return h

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
        scale_hist.SetBinError(i,0)
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
    return repeat(smooth,0)(r)

def smooth(iput):
    tmp = iput[11][0]
    for i in range(12,iput.GetNrows()-1):
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
        data_hist1 = args[2] if len(args)>2 else None
        data_hist2 = args[3] if len(args)>3 else None
        leg = ROOT.TLegend(0.5,0.7,0.9,0.9).Clone()
        # if data_hist:
        #     data_hist.Draw("E1 X0")
        if data_hist1:
            #data_hist1.SetTitle("nue Data")
            data_hist1.SetLineColor(ROOT.kRed+2)
            #data_hist1.SetMarkerColor(ROOT.kRed)
            data_hist1.Draw("SAME")
            leg.AddEntry(data_hist1,"nue data")
        if data_hist2:
            #data_hist2.SetTitle("numu Data")
            data_hist2.SetLineColor(ROOT.kGreen+2)
            #data_hist2.SetMarkerColor(ROOT.kGreen)
            data_hist2.Draw("SAME")
            leg.AddEntry(data_hist2,"weighted numu data")
        if mc_hist1:
            #mc_hist1.SetTitle("nue MC")
            mc_hist1.SetLineColor(ROOT.kRed)
            mc_hist1.Draw("HIST SAME")
            leg.AddEntry(mc_hist1,"nue MC")
        if mc_hist2:
            #mc_hist2.SetTitle("weighted numu MC")
            mc_hist2.SetLineColor(ROOT.kGreen)
            mc_hist2.Draw("HIST SAME")
            leg.AddEntry(mc_hist2,"weighted numu MC")
        leg.Draw()

    def DrawRatio(mnvplotter,h1,h2,include_systematics = False):
        cast = (lambda x:x.GetCVHistoWithError()) if include_systematics else (lambda x:x.GetCVHistoWithStatError()) 
        mnvplotter.DrawDataMCRatio(cast(h1),cast(h2), 1.0 ,True,0,2)

    def DrawDoubleRatio(mnvplotter,h1,h2,h3,h4):
        h_r1 = h1.Clone("{}_ratio1".format(h1.GetName))
        h_r2 = h3.Clone("{}_ratio2".format(h3.GetName))
        h_r1.Divide(h_r1,h2)
        h_r2.Divide(h_r2,h4)
        mnvplotter.DrawDataMCRatio(h_r1.GetCVHistoWithError(), h_r2.GetCVHistoWithError(), 1.0 ,True,0,2)

    nueMCfile,nueMCPOT = getFileAndPOT(nuefile)
    nueBkgfile,nueBkgPOT = getFileAndPOT(nuebkgtunedfile)
    numuMCfile,numuMCPOT = getFileAndPOT(numuweightedfile["mc"])
    numuDatafile,numuDataPOT = getFileAndPOT(numuweightedfile["data"])
    nueDatafile,nueDataPOT = getFileAndPOT(nuedatafile)
    print (nueMCPOT,numuDataPOT,numuMCPOT,nueDataPOT)

    for hist_name in [ "Eavail_q3", "Enu" ,"Eavail_Lepton_Pt","Q3","Q0","tQ0","tQ3","tLepton_Pt","tQ0_tQ3","tQ0_tLepton_Pt","Eavail","tEavail","tEavail_tQ3","tEavail_tLepton_Pt","Lepton_Pt","tEnu","Eavail_q3_true_signal", "Eavail_Lepton_Pt_true_signal", "tEnu_true_signal"]:
        try:
            he = GetSignalHist(nueMCfile,E_SIGNALS,hist_name)
            hmu = GetSignalHist(numuMCfile,MU_SIGNALS,hist_name)
            hebkg = GetSignalHist(nueBkgfile,E_BACKGROUNDS,hist_name)
            hmubkg = GetSignalHist(numuMCfile,MU_BACKGROUNDS,hist_name)
        except :
            print ("hist: {} not found".format(hist_name))
            continue
        hmu.Scale(numuDataPOT/numuMCPOT)
        he.Scale(numuDataPOT/nueMCPOT)
        hebkg.Scale(numuDataPOT/nueBkgPOT)
        hmubkg.Scale(numuDataPOT/numuMCPOT)
        hmudata = numuDatafile.Get(hist_name)
        if hmudata:
            SubtractPoissonHistograms(hmudata,hmubkg)
        hedata = nueDatafile.Get(hist_name)
        if hedata:
            SubtractPoissonHistograms(hedata,hebkg)

        Slicer = PlotTools.Make2DSlice if he.GetDimension()==2 else (lambda hist : [hist])
        if hedata:
            PlotTools.MakeGridPlot(Slicer,DrawRatio,[hedata,he],draw_seperate_legend = he.GetDimension()==2)
            PlotTools.CANVAS.Print("{}{}_heratio.png".format(PLOTPATH,hist_name))
        if hmudata:
            PlotTools.MakeGridPlot(Slicer,DrawRatio,[hmudata,hmu],draw_seperate_legend = hmu.GetDimension()==2)
            PlotTools.CANVAS.Print("{}{}_hmuratio.png".format(PLOTPATH,hist_name))

        if hedata and hmudata:
            PlotTools.MakeGridPlot(Slicer,DrawDoubleRatio,[hedata,he,hmudata,hmu],draw_seperate_legend = he.GetDimension()==2)
            PlotTools.CANVAS.Print("{}{}_doubleratio.png".format(PLOTPATH,hist_name))

            PlotTools.MakeGridPlot(Slicer,DrawRatio,[hedata,hmudata],draw_seperate_legend = he.GetDimension()==2)
            PlotTools.CANVAS.Print("{}{}_dataratio.png".format(PLOTPATH,hist_name))
            PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu,hedata,hmudata],draw_seperate_legend = he.GetDimension()==2)
        else:
            PlotTools.MakeGridPlot(Slicer,Draw,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
        PlotTools.CANVAS.Print("{}{}.png".format(PLOTPATH,hist_name))
        PlotTools.MakeGridPlot(Slicer,DrawRatio,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
        PlotTools.CANVAS.Print("{}{}_MCratio.png".format(PLOTPATH,hist_name))

        for numu_sig,nue_sig in SIGNAL_CHANNEL_PAIR:
            he = GetSignalHist(nueMCfile,[nue_sig],hist_name)
            hmu = GetSignalHist(numuMCfile,[numu_sig],hist_name)
            hmu.Scale(numuDataPOT/numuMCPOT)
            he.Scale(numuDataPOT/nueMCPOT)
            PlotTools.MakeGridPlot(Slicer,DrawRatio,[he,hmu],draw_seperate_legend = he.GetDimension()==2)
            PlotTools.CANVAS.Print("{}{}_{}MCratio.png".format(PLOTPATH,hist_name,numu_sig))


if __name__ == "__main__":
    #MakeCompPlot()
    # MakeScaleFile(["Enu","tEnu","tEnu_true_signal"])
    MakeScaleFile2()
