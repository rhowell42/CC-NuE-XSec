import ROOT
import PlotUtils
import sys

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
inputDir = "Input_Hists"
fakeDataTruthName = "h_data_truth"
unfoldedDataDir = "Unfolded_Data"
smearedDataDir = "Stat_Varied_Smeared_Data"
chi2SummaryDir = "Chi2_Iteration_Dists"
chi2SummaryName = "h_chi2_modelData_trueData_iter_chi2"
medianHistName = "h_median_chi2_modelData_trueData_iter_chi2"
meanChi2ProfileName = "m_avg_chi2_modelData_trueData_iter_chi2_truncated"

lineWidth = 3
iterChosen = None

def draw(fileName,univName,yNDF=None):
    myFile = ROOT.TFile.Open(fileName)
    if not myFile:
        return False
    #Try to infer a useful universe name from the file name
    #univName = fileName[fileName.find("merged") + len("merged") + 1:fileName.find(".root")]
    #if fileName.find("SuSA") != -1: #The SuSA warp is a stand-alone CV, so it needs special treatment
    #  univName = "SuSA"

    #calculate ndf on the fly
    if yNDF is None:
        yNDF = 0
        truth = myFile.Get(inputDir).Get(fakeDataTruthName)
        for i in range(truth.GetSize()):
            if truth.GetBinContent(i)>0:
                yNDF+=1

    spread = myFile.Get(chi2SummaryDir).Get(chi2SummaryName)
    if spread:
        spread.SetTitle("Universe: " + univName)
        spread.SetTitleOffset(0.75, "X")
        spread.SetTitleOffset(0.65, "Y")
        spread.Draw("colz")
  
    profile = myFile.Get(chi2SummaryDir).Get(meanChi2ProfileName)
    if profile:
        profile.SetTitle("Mean Chi2")
        profile.SetLineWidth(lineWidth)
        profile.SetLineColor(ROOT.kBlue)
        profile.SetMarkerStyle(0)
        profile.Draw("SAME")
  
    median = myFile.Get(chi2SummaryDir).Get(medianHistName)
    if median:
        median.SetTitle("Median Chi2")
        median.SetLineWidth(lineWidth)
        median.SetLineColor(ROOT.kBlack)
        median.Draw("HIST SAME")
    #Draw lines at number of degrees of freedom and 2x NDF
    ndfLine = ROOT.TLine(1, yNDF, spread.GetXaxis().GetXmax(), yNDF)
    ndfLine.SetLineWidth(lineWidth)
    ndfLine.SetLineStyle(ROOT.kDashed)
    ndfLine.Draw()
  
    doubleNDFLine = ROOT.TLine(1, 2*yNDF, spread.GetXaxis().GetXmax(), 2*yNDF)
    doubleNDFLine.SetLineColor(ROOT.kRed)
    doubleNDFLine.SetLineWidth(lineWidth)
    doubleNDFLine.SetLineStyle(ROOT.kDashed)
    doubleNDFLine.Draw()
  
    #Draw a line at the chosen number of iterations.
    if iterChosen:
        iterLine = ROOT.TLine(iterChosen + 0.5, 0, iterChosen + 0.5, spread.GetYaxis().GetXmax())
        iterLine.SetLineWidth(lineWidth)
        iterLine.SetLineStyle(ROOT.kDotted)
        iterLine.Draw()
    #Make a custom legend because I don't want to include the 2D histogram
    #while I must Draw() it first to set the right axis limits.
    leg = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
    leg.AddEntry(profile)
    leg.AddEntry(median)
    leg.AddEntry(ndfLine, "Number of Bins", "l")
    leg.AddEntry(doubleNDFLine, "2x Number of Bins", "l")
    if iterChosen:
        leg.AddEntry(iterLine, str(iterChosen) + " iterations", "l")
    leg.Draw()
    can.Print("{}.png".format("".join(fileName.split("/")[-2:])))
    myFile.Close()
    return True

def draw2(fileName,univName):
    f = ROOT.TFile.Open(path)
    if not f:
        return False
    h = f.Get(inputDir).Get("h_mc_reco_Linear")
    h.SetTitle("nue sample scaled")
    h.SetLineColor(ROOT.kRed)
    hdata = f.Get(inputDir).Get("h_data_Linear")
    hdata.SetTitle("standard MC")
    h.Scale(hdata.Integral()/h.Integral())
   
    hheat = ROOT.TH2D("recohists","reco histogram",h.GetNbinsX(),0,h.GetNbinsX(),100,0,h.GetMaximum()*1.1)
   
    for i in range(1,100):
        hthis = f.Get(smearedDataDir).Get("Input_Stat_{}_Iter_1_Linear".format(i))
        for j in range(0,hthis.GetSize()):
            hheat.Fill(j-0.1,hthis.GetBinContent(j))
    hheat.Draw("COLZ")
    h.Draw("HIST SAME")
    hdata.Draw("SAME")
    leg = ROOT.TLegend(0.1, 0.8, 0.3, 0.9)
    leg.AddEntry(h)
    leg.AddEntry(hdata)
    leg.Draw()
    can.Print("{}-Input.png".format("".join(fileName.split("/")[-2:])))
    f.Close()
    return True

def draw3(fileName,univName):
    f = ROOT.TFile.Open(path)
    if not f:
        return False
    #h = f.Get(inputDir).Get("h_mc_truth_Linear")
    #h.SetTitle("nue sample scaled")
    h = f.Get(unfoldedDataDir).Get("Stat_0_Iter_40_Linear")
    h.SetLineColor(ROOT.kRed)
    hdata = f.Get(inputDir).Get("h_data_truth_Linear")
    hdata.SetTitle("standard MC truth")
    h.SetTitle("standard MC Unfolded")
    #h.Scale(hdata.Integral()/h.Integral())
    hheat = ROOT.TH2D("truthhists","truth histogram",h.GetNbinsX(),0,h.GetNbinsX(),100,0,h.GetMaximum()*1.1)
    for i in range(1,100):
        hthis = f.Get(unfoldedDataDir).Get("Stat_{}_Iter_40_Linear".format(i))
        for j in range(0,hthis.GetSize()):
            hheat.Fill(j-0.1,hthis.GetBinContent(j))
    hheat.Draw("COLZ")
    h.Draw("HIST SAME")
    hdata.Draw("SAME")
    leg = ROOT.TLegend(0.1, 0.8, 0.3, 0.9)
    leg.AddEntry(h)
    leg.AddEntry(hdata)
    leg.Draw()
    can.Print("{}-Unfolded.png".format("".join(fileName.split("/")[-2:])))
    f.Close()
    return True

def draw4(fileName,univName):
    f = ROOT.TFile.Open(path)
    if not f:
        return False
    h = f.Get(unfoldedDataDir).Get("Stat_0_Iter_1")
    hrawchi2 = ROOT.TH1D("chi2contrib_0","chi2_contribution",h.GetSize(),0,h.GetSize())
    hheat = ROOT.TH2D("chi2contrib","chi2_contribution",h.GetSize(),0,h.GetSize(),100,-100,100)
    for i in range(0,100):
        TMatrix = f.Get("Chi2Maps").Get("Chi2Map_Stat_{}_Iter_1".format(i))
        projector = ROOT.TMatrixD(TMatrix.GetNrows(),1)
        for _ in range(TMatrix.GetNrows()):
            projector[_][0]=1
        chi2 = ROOT.TMatrixD(TMatrix,ROOT.TMatrixD.kMult,projector)

        for j in range(0,chi2.GetNrows()):
            nx = j%(h.GetNbinsX())
            ny = j//(h.GetNbinsX())
            binn = (h.GetNbinsX()+2)*(ny+1)+nx+2
            hheat.Fill(binn-0.1,chi2[j][0])
            if i==0:
                hrawchi2.Fill(binn-0.1,chi2[j][0])
    hheat.Draw("COLZ")
    hrawchi2.SetLineColor(ROOT.kRed)
    hrawchi2.Draw("HIST SAME")
    can.Print("{}-chi2.png".format("".join(fileName.split("/")[-2:])))
    f.Close()
    return True


if __name__ =="__main__":
    if len(sys.argv) > 1:
        print(sys.argv[1])
        hists = sys.argv[1:]
    else:
        hists = ["CCNUE_Warping_2021-12-15-123747_hists"]
    
    plotpathroot = "/minerva/data/users/hsu/nu_e/plot/warping/"
    target = ["Eavail_Lepton_Pt","Eavail_q3"]
    #model = ["CV","SuSA2p2h"]
    model = ["CV","MK_Model", "FSI_Weight0", "FSI_Weight1", "FSI_Weight2","SuSA2p2h", "LowQ2Pi0","LowQ2Pi2","LowQ2Pi1","LowQ2Pi3","2p2h0","2p2h1","2p2h2","RPA_highq20","RPA_lowq20","RPA_highq21","RPA_lowq21","GenieMaCCQE_UP","GenieMaCCQE_DOWN"]
    Iterations = [i for i in range(0,15)]
    can = ROOT.TCanvas("c1","c1",960,720)
    ROOT.gStyle.SetOptStat(0)
    #c1.SetLogy(1)
    for i in hists:
        for j in target:
            for k in model:
                path = "/pnfs/minerva/scratch/users/hsu/{}/transWrap_{}{}_-1.root".format(i,k,j)
                draw(path,k)
                # draw2(path,k)
                # draw3(path,k)
                # draw4(path,k)
