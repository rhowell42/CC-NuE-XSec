import os
import sys
import ROOT
import PlotUtils

PREFIX = "/minerva/data/users/rhowell"
ROOT.TH1.AddDirectory(False)
c1 = ROOT.TCanvas("t","t",1920,1080)

f1 = ROOT.TFile.Open("{}/MINERvA101/cc_visEPtlep_xsec_2Devtrate.root".format(PREFIX))
f2 = ROOT.TFile.Open("{}/nu_e/unfolded_me1A_nx_CLOSURE_Had.root".format(PREFIX))
f3 = ROOT.TFile.Open("{}/nu_e/kin_dist_mcme1A_nx_CLOSURE_Had.root".format(PREFIX))

h1 = f1.Get("cc_visEPtlep_xsec")
h2 = f2.Get("Eavail_Lepton_Pt_bkg_unfolding")
h3 = f3.Get("Eavail_Lepton_Pt_true_signal")
h4 = f3.Get("Eavail_Lepton_Pt_migration_truth")

print((h1.Integral(),h3.Integral()))
print((h1.GetBinContent(10,10),h3.GetBinContent(10,10)))

print((h2.Integral(),h4.Integral()))
print((h2.GetBinContent(10,10),h4.GetBinContent(10,10)))


t1 = h1.GetCVHistoWithStatError()
t2 = h3.GetCVHistoWithStatError()
t1.Divide(t2)
resolution = 0.01
t1.GetZaxis().SetRangeUser(1-resolution,1+resolution)

t1.Draw("COLZ")
c1.Print("closure_effdept.png")

t1 = h2.GetCVHistoWithStatError()
t2 = h4.GetCVHistoWithStatError()
t1.Divide(t2)
resolution = 0.01
t1.GetZaxis().SetRangeUser(1-resolution,1+resolution)

t1.Draw("COLZ")
c1.Print("closure_effnupt.png")


f1 = ROOT.TFile.Open("{}/MINERvA101/cc_visEq3_xsec_2Devtrate.root".format(PREFIX))

h1 = f1.Get("cc_visEq3_xsec")
h2 = f2.Get("Eavail_q3_bkg_unfolding")
h3 = f3.Get("Eavail_q3_true_signal")
h4 = f3.Get("Eavail_q3_migration_truth")

print((h1.Integral(),h3.Integral()))
print((h1.GetBinContent(10,10),h3.GetBinContent(10,10)))

print((h2.Integral(),h4.Integral()))
print((h2.GetBinContent(10,10),h4.GetBinContent(10,10)))


t1 = h1.GetCVHistoWithStatError()
t2 = h3.GetCVHistoWithStatError()
t1.Divide(t2)
resolution = 0.01
t1.GetZaxis().SetRangeUser(1-resolution,1+resolution)

t1.Draw("COLZ")
c1.Print("closure_effdeq3.png")

t1 = h2.GetCVHistoWithStatError()
t2 = h4.GetCVHistoWithStatError()
t1.Divide(t2)
resolution = 0.01
t1.GetZaxis().SetRangeUser(1-resolution,1+resolution)

t1.Draw("COLZ")
c1.Print("closure_effnuq3.png")


#xsec level
f1 = ROOT.TFile.Open("{}/MINERvA101/nue_MEC_xsec.root".format(PREFIX))
f2 = ROOT.TFile.Open("{}/nu_e/xsec_me1A_nx_CLOSURE_Had.root".format(PREFIX))


h1 = f1.Get("cc_visEPtlep_xsec")
h2 = f2.Get("Eavail_Lepton_Pt_dataxsec")
h3 = f2.Get("Eavail_Lepton_Pt_mcxsec")

print((h1.Integral(),h2.Integral(),h3.Integral()))
print((h1.GetBinContent(10,10),h2.GetBinContent(10,10),h3.GetBinContent(10,10)))

t1 = h1.GetCVHistoWithStatError()
t2 = h2.GetCVHistoWithStatError()
t1.Divide(t2)
resolution = 0.01
t1.GetZaxis().SetRangeUser(1-resolution,1+resolution)

t1.Draw("COLZ")
c1.Print("closure_pt.png")

h1 = f1.Get("cc_visEq3_xsec")
h2 = f2.Get("Eavail_q3_dataxsec")

print((h1.Integral(),h2.Integral()))
print((h1.GetBinContent(10,10),h2.GetBinContent(10,10)))

t1 = h1.GetCVHistoWithStatError()
t2 = h2.GetCVHistoWithStatError()
t1.Divide(t2)
resolution = 0.01
t1.GetZaxis().SetRangeUser(1-resolution,1+resolution)

t1.Draw("COLZ")
c1.Print("closure_q3.png")

