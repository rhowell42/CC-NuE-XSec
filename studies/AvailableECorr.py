#Studying PC of tshower width

import ROOT
import math
from ROOT import *
import glob
from PlotUtils.HistWrapper import HistWrapper
from tools.PlotLibrary import HistHolder
ROOT.gROOT.SetBatch(True)

from config.PlotConfig import LOW_RECOIL_BIN_Q0


fitPath  = '/exp/minerva/data/users/rhowell/antinu_e'
#file1 = 'kin_dist_mcAvailDiffNoCorr_nx_new_fspline.root'.format(fitPath)
file1 = 'kin_dist_mcAuditSystematics_nx_collab1_fspline.root'.format(fitPath)

fRHC = ROOT.TFile.Open("{}/{}".format(fitPath,file1))
#histnotCCQElike = HistHolder("True Visible Energy vs Visible Difference", fRHC, "notCCQElike", True)
#histCCQElike = HistHolder("True Visible Energy vs Visible Difference", fRHC, "CCQElike", True)
histCCNuE = HistHolder("True Visible Energy vs Visible Difference", fRHC, "CCNuE", True)
histQE = HistHolder("True Visible Energy vs Visible Difference", fRHC, "CCNuEQE", True)
histDelta = HistHolder("True Visible Energy vs Visible Difference", fRHC, "CCNuEDelta", True)
histDIS = HistHolder("True Visible Energy vs Visible Difference", fRHC, "CCNuEDIS", True)
hist2p2h = HistHolder("True Visible Energy vs Visible Difference", fRHC, "CCNuE2p2h", True)

CCNuE = histCCNuE.GetHist()
QE = histQE.GetHist()
Delta = histDelta.GetHist()
DIS = histDIS.GetHist()
h2p2h = hist2p2h.GetHist()

#notCCQElike = histnotCCQElike.GetHist() #grab from histholder
#CCQElike = histCCQElike.GetHist()
#fRHC.Close()

#CCQElike.Add(notCCQElike) #I only want signal events
h2p2h.Add(CCNuE)
h2p2h.Add(QE)
h2p2h.Add(Delta)
h2p2h.Add(DIS)

with open('outVisE_Draw.txt', 'w') as f:
  for i in range(0,len(LOW_RECOIL_BIN_Q0)-1):

    h2p2h.GetYaxis().SetRange(i,i+1)
    mean = h2p2h.GetMean()
    rms = h2p2h.GetRMS()
 
    #f.write(str(mean)+'\n')
    f.write(mean) 

