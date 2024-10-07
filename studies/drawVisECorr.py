import ROOT
import os
from ROOT import *
import glob
from PlotUtils.HistWrapper import HistWrapper
from tools.PlotLibrary import HistHolder
from config.PlotConfig import LOW_RECOIL_BIN_Q0
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
fitPath = "{0}/tools/".format(os.environ["CCNUEROOT"])
file_correction = 'outVisE_AuditSys.txt'.format(fitPath)

xbins = LOW_RECOIL_BIN_Q0

c1 = TCanvas("c1","c1",900,700)

datlist = []
#f1 = ROOT.TFile.Open('drawVisE.root', 'recreate')
#data = open('outVisE_AuditSys.txt','r')
data = open("{}/{}".format(fitPath,file_correction))

for value in data:
    value = value.rstrip('\n')
    datlist.append(double(value))

avg = []
for i in range(len(xbins)-1):
    v1 = xbins[i]
    v2 = xbins[i+1]
    v3 = (v1+v2)/2
    avg.append(v3)

#h2 = ROOT.TH2D('h2','',15, 0.0,1.2,15,-0.008,0.003)
h2 = ROOT.TProfile('h2','',len(avg),0.0,1.2)
for i in range(0,len(avg)):

    #h2.Fill(xbins[i],datlist[i])
    h2.Fill(xbins[i],datlist[i])

data.close()
fit_func = TF1("f1", "pol3",len(avg))


fit = h2.Fit('f1',"WS")

fitter = (fit_func.Eval(1.1))

#f1.cd()
#h2.Draw()
#h2.Write()
#h2.Write()
#h2.Draw()

#spline.Draw()
#spline.Write()
