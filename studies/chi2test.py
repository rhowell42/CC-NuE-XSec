import ROOT
import os, time, sys, math
import PlotUtils
from array import array

f1 = ROOT.TFile.Open("/minerva/data/users/hsu/nu_e/kin_dist_mcme1A-E_col12.1_MAD.root")
f2 = ROOT.TFile.Open("/minerva/data/users/hsu/nu_e/kin_dist_mcme1A-E-BigNuE_col12.1_MAD.root")

hists = [#"Eavail_Lepton_Pt_migration",
#          "Eavail_Lepton_Pt_migration_truth",
#          "Eavail_Lepton_Pt_migration_reco",
    "Eavail_q3_migration",
    "Eavail_q3_migration_truth",
    "Eavail_q3_migration_reco",
    #"Eavail_migration",
    #"q3_migration",
    "Lepton_Pt_migration"]

def linearize(h):
    return h 
    # hr = ROOT.TH1D(h.GetName(),h.GetTitle(),h.GetSize(),0,h.GetSize())
    # for i in range(h.GetSize()):
    #     hr.SetBinContent(i,h.GetBinContent(i))
    # return hr

def exclude(h1,h2):
    for i in range(h1.GetSize()):
        if h1.GetBinContent(i)<0:
            h1.SetBinContent(i,0)
            h2.SetBinContent(i,0)
    return h1,h2



for h in hists:
    print (h,f1.Get(h),f2.Get(h))
    h1= linearize(f1.Get(h).Clone("h1"))
    h2= linearize(f2.Get(h).Clone("h2"))
    h1,h2 = exclude(h1,h2)
    #print(h1,h2)
    #h1.Scale(h2.Integral()/h1.Integral())
    
    #arr = array("d",[0]*200)
    chi2 = h2.Chi2Test(h1,"WW P")
    #l = [(abs(arr[i]),h2.GetBinContent(i)) for i in range(len(arr))]
    #l.sort()
    #print (len(arr))
    #print(l)


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
