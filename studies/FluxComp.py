import ROOT
import PlotUtils
import os
from tools import Utilities


prefix = "/minerva/data/users/{}/nu_e/".format(os.environ["USER"])
fRHC = "{}{}".format(prefix,"kin_dist_mcmeRHC_NCDIF_RHCt_Had.root")
fFHC = "{}{}".format(prefix,"kin_dist_mcme1D_NCDIF_RHCt_MAD.root")

hlist = ["tEnu","tEnu_Excess_High_Inline"]

POT_FHC = Utilities.getPOTFromFile(fFHC)
POT_RHC = Utilities.getPOTFromFile(fRHC)
print("POT of FHC: {}, RHC: {}.".format(POT_FHC,POT_RHC))
TFFHC = ROOT.TFile.Open(fFHC)
TFRHC = ROOT.TFile.Open(fRHC)

for h in hlist:
    h_FHC = TFFHC.Get(h).Clone("tEnu_FHC")
    h_RHC = TFRHC.Get(h).Clone("tEnu_RHC")

    c1 = ROOT.TCanvas("C1","C1")
    #c1.SetLogy(1)
    ROOT.gStyle.SetOptStat(0);
    POT = 4e21

    h_FHC.Scale(POT/POT_FHC)
    h_RHC.Scale(POT/POT_RHC)


    h_FHC.SetLineColor(ROOT.kRed)
    h_FHC.Draw("HIST SAME")
    h_RHC.Draw("HIST SAME")

    leg = ROOT.TLegend(0.5,0.7,0.9,0.9)
    leg.AddEntry(h_FHC,"FHC flux")
    leg.AddEntry(h_RHC,"RHC flux")
    leg.Draw()
    c1.Print("test.png")
