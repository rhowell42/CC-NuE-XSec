import ROOT
import math
from functools import partial
import PlotUtils
import UnfoldUtils
ROOT.TH1.AddDirectory(False)

hist_names = {
    "Eavail_q3_bkgSubed",
    "Eavail_Lepton_Pt_bkg_unfolding",
}

dir_path = "/minerva/data/users/hsu/nu_e/"

file_names = {
    "bkgfit_me_q3_tune_RHC_Scale_Sarah.root",
    "bkgfit_me_pt_tune_RHC_Scale_Sarah.root",
    "bkgfit_me_q3_tune_FHC_Scale_Sarah.root",
    "bkgfit_me_pt_tune_FHC_Scale_Sarah.root",
    "bkgfit_me_3factors_tune_no_Scale_Sarah.root",
    "bkgfit_me_visE_tune_no_Scale_Sarah.root"
}



errorband_vector = []
for hist in hist_names:
    for fn in file_names:
        f = ROOT.TFile.Open(dir_path+fn)
        if fn is cv_file:
            base = f.Get(hist)
        else:
            errorband_vector.append(f.Get(hist))
ErrorBand = ROOT.PlotUtils.MnvVertErrorBand2D("bkgfit",base,len(errorband_vector))

for i in range(len(errorband_vector)):
    ErrorBand.GetHists()[i]=errorband_vector[i].Clone()

c1 = ROOT.TCanvas("c1","c1")
ErrorBand.GetErrroband(True).Draw("COLZ")
c1.Print("test.png")


