import ROOT
import math
from tools.unfoldingwrapper import DeOverflowWrapper,WithBackgroundWrapper,IdentityWrapper

do_fake = True
do_eff = False
NStat = 4
RandomGenerator = ROOT.TRandom(0)
DefaultOptions = {
    "do_fake":False,
    "do_eff":False,
    "name":"test",
    "iter":8,
    "do_sys":0
}
ROOT.TH1.AddDirectory(False)


def CalcChi2(h1,h2,cov):
    l1 = h1.GetSize()
    l2 = h2.GetSize()

    # if not (l1 == l2 == cov.GetNcols() == cov.GetNrows()):
    #     print (l1,l2,cov.GetNcols(),cov.GetNrows())
    #     raise ValueError("cov matrix/histograms dont have same dimention to calculate correlation.")

    diff = [(i,h1.GetBinContent(i)-h2.GetBinContent(i)) for i in range(l1) if (h1.GetBinContent(i)-h2.GetBinContent(i))!=0]
    v1 = ROOT.TMatrixD(1,len(diff))
    v2 = ROOT.TMatrixD(len(diff),1)
    this_cov = ROOT.TMatrixD(len(diff),len(diff))
    for i,entry in enumerate(diff):
        v1[0][i]=entry[1]
        v2[i][0]=entry[1]
        for j,entry2 in enumerate(diff):
            this_cov[i][j] = cov[entry[0]][entry2[0]]

    SVD =ROOT.TDecompSVD(this_cov)
    tmp = SVD.Invert()

    tmp = tmp*v2
    tmp = v1*tmp
    return tmp[0][0]
    #return sum([v1[0][i]**2/this_cov[i][i] for i in range(len(diff)) if v1[0][i]>0])

def DrawUnfoldingResult(h_unfolded,h_truth,h_cov,name):
    #return None
    c = ROOT.TCanvas("c2","c2")
    tmp = h_unfolded.Clone()
    tmp.Divide(h_truth)
    tmp.GetZaxis().SetRangeUser(0.5,1.5)
    tmp.Draw("COLZ")
    c.Print("/minerva/data/users/rhowell/tests/"+name+"ratio.png")
    h_cov.Draw("COLZ")
    c.Print("/minerva/data/users/rhowell/tests/"+name+"cov.png")
    

def MyUnfoldWrapper(warpper,data,options):
    response,reco,truth = warpper.GetMigration()
    miggg = ROOT.RooUnfoldResponse(reco,truth,response)
    hlin_measured = warpper.GetToUnfold(data,1)
    unfolder = ROOT.RooUnfold.New(ROOT.RooUnfold.kBayes,miggg,hlin_measured,options["iter"],options["name"])
    unfolder.IncludeSystematics(options["do_sys"])
    h_unfolded = unfolder.Hreco()
    h_unfolded_cov = unfolder.Ereco()
    #h_unfolded_cov.ResizeTo(h_unfolded.fN,h_unfolded.fN)
    return h_unfolded,h_unfolded_cov

def SubtractPoissonDist(h1,h2):
    for i in range(h1.GetSize()):
        mu = h1.GetBinContent(i)-h2.GetBinContent(i)
        sig = h1.GetBinContent(i)+h2.GetBinContent(i)
        h1.SetBinContent(i,mu)
        h1.SetBinError(i,math.sqrt(sig))

def StatShiftedReco(hist):
    hist2 = hist.Clone("new_{}".format(hist.GetName()))
    for i in range(hist.GetSize()):
        mu = RandomGenerator.Poisson(hist.GetBinContent(i))
        #print mu,hist.GetBinContent(i)
        hist2.SetBinContent(i,mu)
        hist2.SetBinError(i,math.sqrt(mu))
    return hist2


#statshifts background subtracted,Wrong
def Approach1(mig,mc_reco,mc_truth,fakedata_reco,fakedata_bkg,fakedata_truth):
    wrapper = DeOverflowWrapper()
    wrapper.setup(mig,mc_reco,fakedata_bkg)
    options = DefaultOptions.copy()
    chi2s = []
    for i in range(NStat):
        unfolded,cov = MyUnfoldWrapper(wrapper,StatShiftedReco(fakedata_reco),options)
        unfolded = wrapper.Delinearlizer(fakedata_truth,unfolded)
        DrawUnfoldingResult(unfolded,fakedata_truth,cov,"approach1{}".format(i))
        chi2 = CalcChi2(unfolded,fakedata_truth,cov)
        chi2s.append(chi2)
    return chi2s

#stat shifts measured, unfolding takes care of fakes. Should be right
def Approach2(mig,mc_reco,mc_truth, fakedata_reco, fakedata_bkg, fakedata_truth):
    wrapper = WithBackgroundWrapper()
    wrapper.setup(mig,mc_reco,fakedata_bkg)
    options = DefaultOptions.copy()
    chi2s = []
    for i in range(NStat):
        unfolded,cov = MyUnfoldWrapper(wrapper,StatShiftedReco(fakedata_reco),options)
        unfolded = wrapper.Delinearlizer(fakedata_truth,unfolded)
        DrawUnfoldingResult(unfolded,fakedata_truth,cov,"approach2{}".format(i))
        chi2 = CalcChi2(unfolded,fakedata_truth,cov)
        chi2s.append(chi2)
    return chi2s
    # chi2s = []
    # options = DefaultOptions.copy()
    # options["do_fake"]=True
    # #toUnfold = StatShiftedReco(fakedata_reco)
    # for i in range(NStat):
    #     toUnfold = StatShiftedReco(fakedata_reco)
    #     unfolded,cov = MyUnfoldWrapper(mig,mc_reco,mc_truth,toUnfold,options)
    #     unfolded = MyDelinerizer(fakedata_truth,unfolded)
    #     DrawUnfoldingResult(unfolded,fakedata_truth,cov,"approach2{}".format(i))
    #     chi2 = CalcChi2(unfolded,fakedata_truth,cov)
    #     chi2s.append(chi2)
    #     del toUnfold
    # return chi2s

#def stat shifts measured, and background, then pass to unfolding. should be right too.
def Approach3(mig,mc_reco,mc_truth, fakedata_reco, fakedata_bkg, fakedata_truth):
    wrapper = IdentityWrapper()
    wrapper.setup(mig,mc_reco,fakedata_bkg)
    options = DefaultOptions.copy()
    chi2s = []
    for i in range(NStat):
        unfolded,cov = MyUnfoldWrapper(wrapper,StatShiftedReco(fakedata_reco),options)
        #unfolded = wrapper.Delinearlizer(fakedata_truth,unfolded)
        #DrawUnfoldingResult(unfolded,fakedata_truth,cov,"approach2{}".format(i))
        chi2 = CalcChi2(unfolded,fakedata_truth,cov)
        chi2s.append(chi2)
    return chi2s

if __name__ == '__main__':
    f = ROOT.TFile.Open('/minerva/data/users/rhowell/nu_e/kin_dist_mcme_RHC_Scale_Sarah.root')
    # mig = ROOT.TH2D(f.Get("Eavail_migration"))
    # mc_reco = ROOT.TH1D(f.Get("Eavail"))
    # mc_truth = ROOT.TH1D(f.Get("Eavail_migration").ProjectionY())
    # fakedata_reco = ROOT.TH1D(f.Get("Eavail"))
    # fakedata_truth = ROOT.TH1D(f.Get("Eavail_migration").ProjectionY())
    # fakedata_signal = ROOT.TH1D(f.Get("Eavail_migration").ProjectionX())

    mig = ROOT.TH2D(f.Get("Eavail_q3_migration"))
    mc_reco = ROOT.TH2D(f.Get("Eavail_q3"))
    mc_truth = ROOT.TH2D(f.Get("Eavail_q3_migration_truth"))
    fakedata_reco = ROOT.TH2D(f.Get("Eavail_q3"))
    fakedata_truth = ROOT.TH2D(f.Get("Eavail_q3_migration_truth"))
    fakedata_signal = ROOT.TH2D(f.Get("Eavail_q3_migration_reco"))
    fakedata_bkg = fakedata_reco.Clone()
    fakedata_bkg.Add(fakedata_signal,-1)

    methods = [Approach1,Approach2]
    for i in methods:
        chi2s = i(mig,mc_reco,mc_truth,fakedata_reco,fakedata_bkg,fakedata_truth)
        print(chi2s)
        print(sum(chi2s)/NStat)
