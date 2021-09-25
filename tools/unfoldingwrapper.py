
import ROOT
import math
import PlotUtils
from array import array

class MyResponseMaker(object):
    def setup(self,response, mc_signal = None, mc_background=None,mc_truth = None):
        self.response=response
        self.mc_signal = mc_signal
        self.mc_truth = mc_truth
        self.mc_background = mc_background

    @staticmethod
    def ResponseCoverter(hist,nm,fakes=None):
        if hist.GetNbinsX()+2 == nm or hist.GetNbinsX() == nm:
            name = "new_{}".format(hist.GetName())
            new_hist =PlotUtils.MnvH2D(name,name,nm,0,nm,nm+(fakes is not None),0,nm+(fakes is not None))

            for i in range(nm):
                for j in range(nm):
                    new_hist.SetBinContent(i+1,j+1,hist.GetBinContent(i,j))
                    new_hist.SetBinError(i+1,j+1,hist.GetBinError(i,j))
            if fakes:
                truth_tmp = hist.ProjectionX()
                for i in range(nm):
                    new_hist.SetBinContent(i+1,nm+1,fakes.GetBinContent(i)-truth_tmp.GetBinContent(i))
                    new_hist.SetBinError(i+1,nm+1,math.sqrt(fakes.GetBinError(i)**2+truth_tmp.GetBinError(i)**2))
            return new_hist
        else:
            raise ValueError("index mismatch, res x bins :{}, mes x bins: {}".format(hist.GetNbinsX(),nm))

    @staticmethod
    def Linearlizer(hist):
        nd = hist.GetSize()
        name = "new_{}".format(hist.GetName())
        new_hist =PlotUtils.MnvH1D(name,name,nd,0,nd)

        for i in range(nd):
            new_hist.SetBinContent(i+1,hist.GetBinContent(i))
            new_hist.SetBinError(i+1,hist.GetBinError(i))
        return new_hist,nd

    @staticmethod
    def Delinearlizer(hist_sample,hist_input):
        return_hist=hist_sample.Clone("new_{}".format(hist_sample.GetName()))
        nd = hist_sample.GetSize()
        for i in range(nd):
            return_hist.SetBinContent(i,hist_input.GetBinContent(i+1))
            return_hist.SetBinError(i,hist_input.GetBinError(i+1))
        return return_hist

    @staticmethod
    def SubtractPoissonHistograms(h,h1):
        try:
            h.AddMissingErrorBandsAndFillWithCV (h1)
        except:
            print("Not MnvHXD, dont mind")
        errors = []
        for i in range(h.GetSize()):
            errors.append(math.sqrt(h.GetBinError(i)**2 + h1.GetBinError(i)**2))
        h.Add(h1,-1)
        for i in range(h.GetSize()):
            h.SetBinError(i,errors[i])

        return h

class IdentityWrapper(MyResponseMaker):
    def GetMigration(self):
        return self.response,self.mc_signal,self.mc_truth

    def GetToUnfold(self,data,pot_scale,background = None):
        if background is None:
            background = self.mc_background.Clone()
        background.Scale(pot_scale)
        MyResponseMaker.SubtractPoissonHistograms(data,background)
        return data

    def PostProcessing(self,data):
        return data

    def PrepareTruth(self,data_truth):
        return data_truth


class DeOverflowWrapper(MyResponseMaker):
    def GetMigration(self):
        nm = self.mc_background.GetSize()
        new_migration = MyResponseMaker.ResponseCoverter(self.response,nm)
        return new_migration,new_migration.ProjectionX("reco",0,-1),new_migration.ProjectionY("truth",0,-1)

    def GetToUnfold(self,data,pot_scale,background = None):
        if background is None:
            background = self.mc_background.Clone()
        background.Scale(pot_scale)
        MyResponseMaker.SubtractPoissonHistograms(data,background)
        return MyResponseMaker.Linearlizer(data)[0]

    def PostProcessing(self,data):
        return MyResponseMaker.DeLinearlizer(self.truth,data)

    def PrepareTruth(self,data_truth):
        return MyResponseMaker.Linearlizer(data_truth)[0]


class WithBackgroundWrapper(MyResponseMaker):
    def GetMigration(self):
        mc_reco = self.mc_signal.Clone()
        mc_reco.Add(self.mc_background)
        new_reco,nm = MyResponseMaker.Linearlizer(mc_reco)
        new_migration = MyResponseMaker.ResponseCoverter(self.response,nm)
        return new_migration,new_reco,new_migration.ProjectionY()

    def GetToUnfold(self,data,pot_scale,background=None):
        return MyResponseMaker.Linearlizer(data)[0]

    def PostProcessing(self,data):
        return MyResponseMaker.DeLinearlizer(self.truth,data)

    def PrepareTruth(self,data_truth):
        tmp = MyResponseMaker.Linearlizer(data_truth)[0]
        return tmp




class TruthHistogramMapper(object):
    def __init__(self,mig,truth):
        self.mig = mig
        self.truth = truth
        self.new_mig = None
        self.new_truth = None
        self.binnings = None

    def truncating1D(self,hist,target):
        if not isinstance(hist,TH1D):
            raise ValueError("put in 1d histogram please")
        segment = [hist.GetBinLowEdge(1)]
        accu = 0
        for i in range(1,hist.GetNbinsX()+1):
            if accu >= target:
                accu = 0
                segment.append(hist.GetBinLowEdge(i))
            accu += hist.GetBinContent(i)
        segment.append(hist.GetBinLowEdge(hist.GetNbinsX()+1))
        return array('d',segment)

    def truncating(self,hist,target):
        binnings = []
        if hist.GetNbinsY()==1:
            binnings.append(hist,target)
            return binnings

        for i in range(hist.GetNbinsY()+2):
            binnings.append(self.truncating1D(hist.ProjectionY("",i,i),target))
        return binnings

    def rebin(self,hist):
        if hist.GetNbinsY()==1:
            return hist.Rebin(len(self.binnings[0])-1,"{}_rebined".format(hist.GetName()),self.binnings[0])

        hists = []
        for i in range(hist.GetNbinsY()+2):
            hists.append(hist.ProjectionY("_px_{}".format(i),i,i).Rebin(len(self.binnings[i]-1),"",self.binnings[i]))
        return hists

    def match(self,hist,hists,Nbin):
        ybin = Nbin//(hist.GetNbinsX()+2)
        xbin = Nbin %(hist.GetNbinsX()+2)
        xcenter = hist.GetXaxis().GetBinCenter(xbin)
        offset = sum(hists[i].GetSize() for i in range(ybin))
        return offset + hists[ybin].GetBin(xcenter)

    def mergeMigbins(self):
        Nbins = sum(hist.GetSize() for hist in self.new_truth)
        hist = self.mig
        new_hist = ROOT.TH2D("{}_merged".format(hist.GetName()),hist.GetNbinsX(),0,hist.GetNbinsX(),Nbins,0,Nbins)
        for i in range(hist.GetNbinsY()+2):
            bin_new = self.match(self.truth,self.new_truth,i)
            for j in range(0,hist.GetNbinsX()+2):
                new_hist.AddBinContent(bin_new,hist.GetBinContent(j,i))
        return new_hist

    def flatten(self,hists):
        size = sum(hist.GetSize() for hist in hists)
        new_hist = ROOT.TH1D("flatten_{}".format(hists[0].GetName()),size-2,0,size-2)
        j = 0
        for h in hists:
            for i in range(h.GetSize()):
                new_hist.AddBinContent(j,h.GetBinContent(i))
                j+=1
        return new_hist

    def GetNewMigration(self,target):
        self.binnings = self.truncating(self.truth,target)
        self.new_truth = self.rebin(self.truth)
        self.new_mig = self.mergeMigbins()
        
