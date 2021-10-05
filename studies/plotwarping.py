import ROOT
import PlotUtils

hists = ["CCNUE_Warping_2021-10-04-225958"]
target = ["Eavail_Lepton_Pt","Eavail_q3"]
model = ["CV","MK_Model", "FSI_Weight0", "FSI_Weight1", "FSI_Weight2","SuSA2p2h", "LowQ2Pi0","LowQ2Pi2","LowQ2Pi1","LowQ2Pi3","2p2h0","2p2h1","2p2h2","RPA_highq20","RPA_lowq20","RPA_highq21","RPA_lowq21","GenieMaCCQE_UP","GenieMaCCQE_DOWN"]
Iterations = [i for i in range(0,15)]
c1 = ROOT.TCanvas("c1","c1",1920,1080)
c1.SetLogy(1)
for i in hists:
    for j in target:
        for k in model:
            path = "/pnfs/minerva/scratch/users/hsu/{}_hists/transWrap_{}{}_0.root".format(i,k,j)
            f = ROOT.TFile.Open(path)
            if not f:
                continue
            d = f.Get("Chi2_Iteration_Dists")
            if not d:
                continue
            h = d.Get("h_median_chi2_modelData_trueData_iter_chi2")
            h.Draw()
            c1.Print("{}{}{}.png".format(i,j,k))


#path = "/pnfs/minerva/persistent/users/hsu/{}_hists/transWrap_{}{}_0.root".format("CCNUE_Waring5","FSI_Weight0","Eavail_Lepton_Pt")
#path = "/minerva/app/users/hsu/cmtuser/Minerva_CCNuEAna/Ana/CCNuE/macros/Low_Recoil/studies/transWrapT_FSI_Weight0Eavail_Lepton_Pt.root"
# c1.SetLogy(0)
#for i in Iterations:
#     f = ROOT.TFile.Open(path)
#     if not f:
#         continue
#     d = f.Get("Unfolded_Data")
#     h1 = d.Get("Stat_0_Iter_{}_Linear".format(i))
#     if not h1:
#         continue
#     d=f.Get("Input_Hists")
#     h2 = d.Get("h_data_truth_Linear")
#    h2.Draw("HIST")
#    h1.Draw("SAME")
#    c1.Print("comp_{}.png".format(i))
