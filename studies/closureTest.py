import ROOT
import PlotUtils
from config.AnalysisConfig import AnalysisConfig
from tools import Utilities, PlotTools
from tools.PlotLibrary import PLOT_SETTINGS
import subprocess
import os

ClosureFileDIR= "/minerva/data/users/hsu/Closure/bigNuE"
MYEVT = "/minerva/data/users/hsu/nu_e/kin_dist_mcme1D_nx_col13_Closure_MAD.root"
MYXSEC = "/minerva/data/users/hsu/nu_e/xsec_me1D_nx_col13_Closure_MAD.root"

CHANNEL_MAPPING = {
    "qe": "_CCNuEQE",
    "res": "_CCNuEDelta",
    "2p2h":"_CCNuE2p2h"
}
tolerance = 0.1

def ratio(h1,h2,suffix=None):
    h1.ClearAllErrorBands()
    htmp = h1.Clone("{}_Clone".format(h1.GetName()))
    htmp.Divide(h1,h2)
    htmp.SetMaximum(1+tolerance)
    htmp.SetMinimum(1-tolerance)
    can = ROOT.TCanvas("c2","c2")
    htmp.Draw("COLZ")
    name = h1.GetName() if suffix is None else "{}{}".format(h1.GetName(),suffix)
    can.Print("{}.png".format(name))


def EfficiencyDenominator(channel=None):
    MyEVTFile = ROOT.TFile.Open(MYEVT)
    ptname = "cc_visEPtlep_xsec" if channel is None else "cc_visEPtlep_{}_xsec".format(channel)
    q3name = "cc_visEq3_xsec" if channel is None else "cc_visEq3_{}_xsec".format(channel)
    GENIEXSecFilePt = ROOT.TFile.Open("{}/{}_2Devtrate.root".format(ClosureFileDIR,ptname))
    GENIEXSecFileQ3 = ROOT.TFile.Open("{}/{}_2Devtrate.root".format(ClosureFileDIR,q3name))
    suffix = CHANNEL_MAPPING[channel] if channel is not None else ""
    ratio(MyEVTFile.Get("Eavail_Lepton_Pt_true_signal"+suffix),GENIEXSecFilePt.Get(ptname))
    ratio(MyEVTFile.Get("Eavail_q3_true_signal"+suffix),GENIEXSecFileQ3.Get(q3name))
    MyEVTFile.Close()
    GENIEXSecFilePt.Close()
    GENIEXSecFileQ3.Close()

def Xsec():
    MyXsecFile= ROOT.TFile.Open(MYXSEC)
    GENIEXsecFile = ROOT.TFile.Open("{}/nue_MEC_xsec.root".format(ClosureFileDIR))
    Pt = (MyXsecFile.Get("Eavail_Lepton_Pt_dataxsec"),GENIEXsecFile.Get("cc_visEPtlep_xsec"))
    ratio(Pt[0],Pt[1],"xsec")
    q3 = (MyXsecFile.Get("Eavail_q3_dataxsec"),GENIEXsecFile.Get("cc_visEq3_xsec"))
    ratio(q3[0],q3[1],"xsec")
    MyXsecFile.Close()
    GENIEXsecFile.Close()

def FBFTest():
    with open("/minerva/data/users/hsu/Closure/me1D.txt") as f:
        s = 0
        for i,line in enumerate(f):
            print("file:{}, {}".format(i,line))
            f = ROOT.TFile.Open("/pnfs/minerva/scratch/users/hsu/CCNUE_selection_2022-03-19-150821_hists/kin_dist_mcme1D_nx_col13_Closure_MAD_{}.root".format(i))
            h=f.Get("Eavail_q3_true_signal")
            if h:
                my_evt=h.GetEntries()
                s+=my_evt
                print(my_evt)
                # os.system("runCCIncForNuEMEC {} | grep signal".format(line.rstrip("\n")))
        print(s)

if __name__ =="__main__":
    #FBFTest()
    Xsec()
    EfficiencyDenominator()
    EfficiencyDenominator("qe")
    EfficiencyDenominator("res")
