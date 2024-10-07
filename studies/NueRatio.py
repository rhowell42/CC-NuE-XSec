import ROOT
import math
from ROOT import *
import glob
from array import *
import PlotUtils
from studies import NueWeight
from NueWeight import GetNueWeightHist, GetNueBarWeightHist

ROOT.gROOT.SetBatch(True)
CANVAS = ROOT.TCanvas("c2","c2",1600,1000)

MNVPLOTTER = PlotUtils.MnvPlotter()
f = ROOT.TFile.Open('/exp/minerva/data/users/rhowell/antinu_e/tEnuRatio_RhcFhcNue.root', 'recreate')

tEnuRatio_RhcFhcNue = GetNueWeightHist()
tEnuRatio_FhcRhcNueBar = GetNueBarWeightHist()
CANVAS.Update()

tEnuRatio_RhcFhcNue.Draw()

f.cd()
tEnuRatio_RhcFhcNue.Write()
f.Close()
