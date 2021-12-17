import ROOT
import os, time, sys, math
import PlotUtils
from array import array

#f1 = ROOT.TFile.Open("/minerva/data/users/hsu/nu_e/kin_dist_mcme1A-E_col12.1_MAD.root")
f2 = ROOT.TFile.Open("/minerva/data/users/hsu/nu_e/kin_dist_mcme1A-F-BigNuE_col12.2_MAD.root")

hists = ["Eavail_Lepton_Pt_migration",
#          "Eavail_Lepton_Pt_migration_truth",
#          "Eavail_Lepton_Pt_migration_reco",
    "Eavail_q3_migration"]
    #"Lepton_Pt_migration"]
threshold = 10
p_value = 0.1


def exclude(h1,h2):
    if not threshold:
        return False
    for i in range(h1.GetSize()):
        if h1.GetBinContent(i) + h2.GetBinContent(i)<threshold:
            h1.SetBinContent(i,0)
            h2.SetBinContent(i,0)
    return True

for hname in hists:
    h = f2.Get(hname)
    n = h.GetNbinsX()+2
    for i in range(n):
        for j in [i-9,i+9]:
            if not 0<j<n:
                continue
            print(i,j)
            h1 = h.ProjectionX("_Px1",i,i)
            h2 = h.ProjectionX("_Px2",j,j)
            if h1.Integral(0,-1)==0 or h2.Integral(0,-1)==0:
                continue
            exclude(h1,h2)
            
            p = h1.Chi2Test(h2,"WW P")
            if p>p_value:
                print ("i and j are same distribution: {},{}. ".format(i,j))


#likely rejects

# for h in hists:

#     h1= f1.Get(h).Clone("h1")
#     h2= f2.Get(h).Clone("h2")

#     c1= ROOT.TCanvas("c1","c1")

#     h2.Scale(h1.Integral()/h2.Integral())
#     h1.SetBarOffset(0.2)
#     h1.Draw("COLZ TEXT E")
#     h2.SetBarOffset(-0.2)
#     h2.Add(h1,-1)
#     h2.Draw("TEXT SAME")
#     ROOT.gStyle.SetOptStat(0)
#     c1.Print("{}_linear.png".format(h))
