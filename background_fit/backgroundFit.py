import os
import sys
import ROOT
import PlotUtils
import math
import copy
from array import array
from collections import OrderedDict

from tools.PlotLibrary import HistHolder
from config.AnalysisConfig import AnalysisConfig
from config import BackgroundFitConfig
from tools import Utilities,PlotTools
from config.UnfoldingConfig import HISTOGRAMS_TO_UNFOLD
from config.DrawingConfig import Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer,PLOTS_TO_MAKE
from config.SignalDef import SIGNAL_DEFINATION

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)


#fit for single mnvHXD
# class Fitter( ROOT.TPyMultiGenFunction ):
#     def __init__( self, data_hist, mc_hists, previous_scale, NLL = False):
#         super(Fitter, self).__init__(self)
#         self.fitfn = data_hist # {E_e_High_Psi : mnvH2D, theta_e_High_Psi:mnvH2D}
#         self.mcfn = mc_hists #{E_e : { Signal:mnvH2d, Pi0:mnvH2D,Other:mnvH2D}, theta_e: ...}
#         self.prescl = previous_scale # { Signal: [scale_for_each_q3_bin,...],Pi0: [...],... }
#         self.Nybins = next(data_hist.itervalues()).GetNbinsY() #Number of q3 bins,excluding overflow and underflow bins
#         self.ybin = 1 #start from first bin
#         self.useNLL = NLL # use negative log likelihood
#         self.ndf = 0

#     def SetNthQ3bin(self,i):
#         if i < 1 or i>self.Nybins:
#             print ("unphysical q3 bin setting")
#             return False
#         else:
#             self.ybin = i
#             self.ndf= 0
#             return True

#     def NDim( self ):
#         return len(self.prescl)

#     def NDF(self):
#         return self.ndf

#     def DoEval( self, args ): #arguments are scale facotrs of each histogram 0: signal, 1:
#         first =  self.ndf == 0
#         chi2 = {hist : 0.0 for hist in self.fitfn}
#         NLL  = {hist : 0.0 for hist in self.fitfn} #negative log likelihood
#         for hist in self.fitfn:
#             for i in range(1,self.fitfn[hist].GetXaxis().GetNbins()+1):
#                 data_value = self.fitfn[hist].GetBinContent(i,self.ybin)
#                 mc_value = 0.0
#                 for j ,group in enumerate(BackgroundFitConfig.BACKGROUND_GROUP.keys()):
#                     mc_value += args[j] * self.mcfn[hist][group].GetBinContent(i,self.ybin)*self.prescl[group][self.ybin-1]

#                 chi2[hist] += math.pow((data_value-mc_value),2)/data_value if data_value >0.0 else 0.0
                #NLL[hist] += mc_value-data_value*math.log(mc_value)+math.lgamma(data_value+1) if mc_value>0.0 else 0.0
        #         NLL[hist] += mc_value-data_value + data_value * (math.log(data_value/mc_value) if ( data_value > 0 and mc_value > 0 ) else 0.0)
        #         if first:
        #             self.ndf += 1 if data_value >0.0 else 0

        # #print chi2, NLL
        # self.chi2 = chi2
        # self.NLL = NLL
        # return sum(NLL.itervalues()) if self.useNLL else sum(chi2.itervalues())

class Fitter(ROOT.TPyMultiGenFunction):
    def __init__(self):
        super(Fitter, self).__init__(self)
        self.data_values=[]
        self.mc_values={} # {"NCPI0":[2.0,3.9,4.5,...]}
        self.factorsmap = {} # {"NCPi0":[3,3,3,4,...]}
        self.constraint = BackgroundFitConfig.REGULATION_PARAMETER

    def DoEval(self, args):
        return self.CalChi2(args)+self.constraint*self.CalPenalty(args)

    def CalChi2(self,args):
        chi2 =0 
        for i in range(len(self.data_values)):
            data = self.data_values[i]
            mc = 0
            for cate in self.mc_values:
                try:
                    mc_factor = args[self.factorsmap[BackgroundFitConfig.CATEGORY_FACTORS[cate]][i%self.data_hist_nbins]]

                except KeyError:
                    mc_factor = 1.0
                mc+= self.mc_values[cate][i]*mc_factor
            if BackgroundFitConfig.USE_NLL:
                chi2 += mc - data*math.log(mc) + ROOT.Math.lgamma(data+1) if mc>0 else 0
            else:
                if data>0:
                    sig2 = data
                elif mc>0:
                    sig2 = mc
                else:
                    sig2 = None
                chi2 += math.pow((data-mc),2)/sig2 if sig2 else 0

        return chi2

    def CalPenalty(self,args):
        pen = 0
        for i in range(self.NDim()):
             pen += (args[i]-1.0)**2
        return pen

    def NDim( self ):
        #return 5
        return sum([len(li)+1 for li in iter(BackgroundFitConfig.SCALE_FACTORS.values())])

    def NDF(self):
        return len(data_values)-self.NDim()

    #input th1D, ie post choosing errorband
    def AddSideband(self,data_hist, mc_ints):
        for i in range(0,data_hist.GetSize()):
            self.data_values.append(data_hist.GetBinContent(i))
            for cate in mc_ints:
                if mc_ints[cate] is None:
                    continue
                self.mc_values.setdefault(cate,[]).append(mc_ints[cate].GetBinContent(i))

    def SetScaleFactorMap(self, hist, yaxis):
        i = 0
        for cate,factor_binning in BackgroundFitConfig.SCALE_FACTORS.items():
            for j in range(0,hist.GetSize()):
                ny = j//(hist.GetNbinsX()+2)
                nx = j%(hist.GetNbinsX()+2)
                bin_low_edge = hist.GetYaxis().GetBinLowEdge(ny) if yaxis else hist.GetXaxis().GetBinLowEdge(nx)
                set = False
                for k,v in enumerate(factor_binning):
                    if bin_low_edge<v:
                        self.factorsmap.setdefault(cate,[]).append(i+k)
                        set = True
                        break
                if not set:
                    self.factorsmap.setdefault(cate,[]).append(i+len(factor_binning))
            i+=(len(factor_binning)+1)
        print((self.factorsmap))
        self.data_hist_nbins = hist.GetSize()

    def Reset(self):
        del self.data_values[:]
        self.mc_values.clear()

    def Clone(self):
        return copy.copy(self)


def Get1DScaleFactor(variable_hists,scale_hists):
    scale_dict = {}
    comparable_scale = MakeComparableMnvHXD(variable_hists.GetHist(), scale_hists, False)
    for cate in variable_hists.hists:
        if variable_hists.hists[cate] is None:
            continue
        scaled =variable_hists.hists[cate].Clone()
        try:
            scale = comparable_scale[BackgroundFitConfig.CATEGORY_FACTORS[cate]]
            scaled.Multiply(scaled,scale)
            scale_dict[BackgroundFitConfig.CATEGORY_FACTORS[cate]] = scaled.Integral()/variable_hists.hists[cate].Integral()
        except KeyError:
            pass
        del scaled
    return scale_dict

def MakeComparableMnvHXD(hist, scale_hist, y_axis=False):
    new_scale = {}
    #print hist.GetName()
    for cate in scale_hist:
        new_scale[cate] = hist.Clone()
    print("y_axis?{}".format(y_axis))

    xbins = hist.GetNbinsX()+2 #including under/overflows.
    #ybins = hist.GetNbinsY()+2
    for i in range(0,hist.GetSize()):
        nx = i%xbins
        ny = i//xbins
        scale_bin_entry = hist.GetYaxis().GetBinCenter(ny) if y_axis else hist.GetXaxis().GetBinCenter(nx)
        #k = next(scale_hist.itervalues()).FindBin(scale_bin_entry)
        for cate in scale_hist:
            k = scale_hist[cate].FindBin(scale_bin_entry)
            new_scale[cate].SetBinContent(i,scale_hist[cate].GetBinContent(k))
            new_scale[cate].SetBinError(i,scale_hist[cate].GetBinError(k))
            for bandname in new_scale[cate].GetErrorBandNames():
                errorband = new_scale[cate].GetVertErrorBand(bandname)
                for ith in range(errorband.GetNHists()):
                    errorband.GetHist(ith).SetBinContent(i,scale_hist[cate].GetVertErrorBand(bandname).GetHist(ith).GetBinContent(k))
                    errorband.GetHist(ith).SetBinError(i,scale_hist[cate].GetVertErrorBand(bandname).GetHist(ith).GetBinError(k))

    # for cate in scale_hist:
    #     h = new_scale[cate]
    #     for i in range(h.GetNbinsX()+2):
    #         for j in range(h.GetNbinsY()+2):
    #             print (i,j,h.GetBinContent(i,j))

    return new_scale

def WriteScaleToMnvH1D(hist, scale, scale_err = None,  errorband=None,i=None):
    for group in scale:
        if errorband is None:
            universe_hist = hist[group]
        else:
            universe_hist= hist[group].GetVertErrorBand(errorband).GetHist(i)

        for j,value in enumerate(scale[group]):
            universe_hist.SetBinContent(j,value)
            if scale_err is not None:
                universe_hist.SetBinError(j,scale_err[group][j])
            else:
                universe_hist.SetBinError(j,0)

def RunUniverseMinimizer(minimizer, fitter, data_histholders, mc_histholders, error_band = None, i = None):
    fitter.Reset()
    for index in range(len(data_histholders)):
        data_input = data_histholders[index].GetHist()
        mc_input={}
        for cate in mc_histholders[index].hists:
            if cate == "Total" or mc_histholders[index].hists[cate] is None:
                continue
            mc_input[cate] = mc_histholders[index].hists[cate]
            if error_band is not None:
                mc_input[cate] = mc_input[cate].GetVertErrorBand(error_band).GetHist(i)
        fitter.AddSideband(data_input,mc_input)
    minimizer.Minimize()
    result = minimizer.X()
    errors = minimizer.Errors()
    scale ={}
    scale_error = {}
    j = 0
    # fill what the scale factor will be for each bin of scale factor histogram.
    # adjust for overflow bins
    for cate in BackgroundFitConfig.SCALE_FACTORS:
        scale[cate]=[]
        scale_error[cate]=[]
        offset = len(BackgroundFitConfig.SCALE_FACTORS[cate])+1
        for k in range(j,j+offset):
            scale[cate].append(result[k])
            scale_error[cate].append(errors[k])
        j += offset
    return scale,scale_error

def RunMinimizer(data_histholders, mc_histholders,scale_hists):
    minimizer = ROOT.Math.Factory.CreateMinimizer("Minuit2")
    minimizer.SetTolerance(0.0001)
    minimizer.SetPrintLevel(1)
    step = 0.1
    start = 1.0
    fitter = Fitter()
    minimizer.SetFunction(fitter)
    
    fitter.SetScaleFactorMap(data_histholders[0].GetHist(),BackgroundFitConfig.Yaxis)
    for i in range(fitter.NDim()):
        minimizer.SetVariable(i,str(i),start,step)
        minimizer.SetVariableLowerLimit(i,0.0)
    scale,scale_err = RunUniverseMinimizer(minimizer,fitter,data_histholders,mc_histholders)
    WriteScaleToMnvH1D(scale_hists,scale,scale_err)

    #then errorbands:
    for error_band in (mc_histholders[0].GetHist().GetErrorBandNames()):
        for i in range(mc_histholders[0].GetHist().GetVertErrorBand(error_band).GetNHists()):
            scale,_ = RunUniverseMinimizer(minimizer,fitter,data_histholders,mc_histholders,error_band,i)
            WriteScaleToMnvH1D(scale_hists,scale,None,error_band,i)

def TuneMC(hist_holder, scale_hists, x_axis=False, y_axis=False):
    if (x_axis and y_axis):
        return None # shouldnt happend
    elif not (x_axis or y_axis):
        try:
            ScaleCategories1D(hist_holder,scale_hists) #scale_dict is a global variable
        except AttributeError:
            return False
    else:
        comparable_scale = MakeComparableMnvHXD(hist_holder.GetHist(),scale_hists,y_axis)
        try:
            ScaleCategories(hist_holder,comparable_scale)
        except AttributeError:
            return False
    hist_holder.ResumTotal()
    return True

def ScaleCategories(hist_holder,scale_hists):
    for cate in hist_holder.hists:
        try:
            scale = scale_hists[BackgroundFitConfig.CATEGORY_FACTORS[cate]]
            hist_holder.hists[cate].Multiply(hist_holder.hists[cate],scale)
        except KeyError:
            continue
    

def ScaleCategories1D(hist_holder,scale_dict):
    for cate in hist_holder.hists:
        try:
            scale = scale_dict[BackgroundFitConfig.CATEGORY_FACTORS[cate]]
            hist_holder.hists[cate].Scale(scale)
        except KeyError:
            continue

def BackgroundSubtraction( data_hists, mc_hists):
    data_hists.POTScale(False)
    mc_hists.POTScale(False)
    out_data = data_hists.GetHist().Clone()
    out_mc = mc_hists.hists["Total"].Clone()
    out_data.AddMissingErrorBandsAndFillWithCV(out_mc)

    for group in mc_hists.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINATION:
            SubtractPoissonHistograms(out_data,mc_hists.hists[group])
            SubtractPoissonHistograms(out_mc,mc_hists.hists[group])
    return out_data,out_mc

def GetBackground(mc_hists):
    out_bkg = mc_hists.hists["Total"].Clone("bkgTotal")
    out_bkg.Reset()

    for group in mc_hists.hists:
        if group == "Total":
                continue
        elif group not in SIGNAL_DEFINATION:
            out_bkg.Add(mc_hists.hists[group])
    return out_bkg

def SubtractPoissonHistograms(h,h1):
    errors = []
    for i in range(h.GetSize()):
        errors.append(math.sqrt(h.GetBinError(i)**2 + h1.GetBinError(i)**2))
    h.Add(h1,-1)
    for i in range(h.GetSize()):
        h.SetBinError(i,errors[i])

    return h

def GetScaledDataMC(hist,datafile,mcfile,region):
    data_hist = HistHolder(hist,datafile,region,False,pot_scale)
    mc_hist = HistHolder(hist,mcfile,region,True,pot_scale)
    fit_on_axis = scaled_hist_name.upper() in data_hist.plot_name.upper()
    #print scaled_hist_name.upper(),data_hist.plot_name.upper()
    if fit_on_axis:
        fit_on_yaxis = ("_"+scaled_hist_name).upper() in data_hist.plot_name.upper()
        print(("fit {} on {} axis".format(data_hist.plot_name, "y" if fit_on_yaxis else "x")))
        TuneMC(mc_hist, scale_hists, not fit_on_yaxis , fit_on_yaxis )
    else:
        print(("fit {} on weighted ave".format(data_hist.plot_name)))
        variable_hist = HistHolder(BackgroundFitConfig.HIST_TO_FIT,mcfile,region,True,pot_scale)
        scale_dict = Get1DScaleFactor(variable_hist,scale_hists)
        TuneMC(mc_hist, scale_dict, False , False )
    return data_hist,mc_hist

def MakePlot(data_hists,mc_hists,config):
    if not (data_hists.valid and mc_hists.valid):
        return False
    #scale
    #mc_hists.ResumTotal()
    if "scale" in config:
        config["scale"](data_hists)
        config["scale"](mc_hists)
    else:
        Default_Scale(data_hists)
        Default_Scale(mc_hists)
    CanvasConfig = config.setdefault("canvasconfig",lambda x:True)
    PlotType = config.setdefault("plot_type",Default_Plot_Type)
    slicer = config.setdefault("slicer", DefaultSlicer(data_hists))
    draw_seperate_legend = config.setdefault("draw_seperate_legend",data_hists.dimension!=1 and PlotType != "migration")
    try:
        custom_tag = config["tag"]+PlotType if "tag" in config else PlotType+AnalysisConfig.bkgTune_tag
        if PlotType == "custom":
            plotfunction,hists=config["getplotters"](data_hists,mc_hists)
        else:
            if "args" in config:
                args = config["args"]
            elif "args" in DefaultPlotters[PlotType]:
                args = DefaultPlotters[PlotType]["args"]
            else:
                args = None
            if args is None:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists)
            else:
                plotfunction,hists = DefaultPlotters[PlotType]["func"](data_hists,mc_hists,*args)
            PlotTools.MakeGridPlot(slicer,plotfunction,hists,CanvasConfig,draw_seperate_legend)
            PlotTools.Print(AnalysisConfig.PlotPath(data_hists.plot_name,sideband,custom_tag))
            print("plot {} made.".format(data_hists.plot_name))
    except KeyError as e:
        print("plot {} not made.".format(data_hists.plot_name))
        print(e)
        return False
    return True


if __name__ == "__main__":
    #input knobs
    playlist=AnalysisConfig.playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    datafile,mcfile,pot_scale = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag)
    #output knobs:
    background_fit_tag = AnalysisConfig.bkgTune_tag
    scalefile=ROOT.TFile.Open(AnalysisConfig.BackgroundFitPath(playlist,background_fit_tag),"RECREATE")
    BackgroundFitConfig.SetGlobalParameter(background_fit_tag)

    data_histholders = []
    mc_histholders = []
    scaled_hist_name = None

    #fit scale histograms
    scale_hists = {}
    for factor,binning in BackgroundFitConfig.SCALE_FACTORS.items():
        mc_histholder = HistHolder(BackgroundFitConfig.HIST_TO_FIT,mcfile,"Signal",True,pot_scale)
        scale_hists[factor] = mc_histholder.GetHist().Rebin(len(binning)-1,factor,array('d',binning))
        scale_hists[factor].GetYaxis().SetTitle("Scale Factor")
        if scaled_hist_name is None:
            scaled_hist_name = mc_histholder.plot_name

    for region in AnalysisConfig.sidebands:
        data_histholders.append(HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,datafile,region,False))
        mc_histholders.append(HistHolder(BackgroundFitConfig.HIST_OBSERVABLE,mcfile,region,True,pot_scale))
        mc_histholders[-1].POTScale(False)

    RunMinimizer(data_histholders,mc_histholders,scale_hists)

    for hist in scale_hists.values():
        hist.Write()


    region = "Signal"

    for hist in HISTOGRAMS_TO_UNFOLD:
        data_hist,mc_hist = GetScaledDataMC(hist,datafile,mcfile,region)
        #data_hist.GetHist().Write(data_hist.plot_name)
        mc_hist.GetHist().Write(data_hist.plot_name)
        bkg = GetBackground(mc_hist)
        bkg.Write(data_hist.plot_name+"_mcbkg")
        #data.Write(data_hist.plot_name+"_bkgSubed")
        #mc.Write(data_hist.plot_name+"_mc")

    #change playlist
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    datafile,mcfile,pot_scale = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag)

    for config in PLOTS_TO_MAKE:
        sideband_group =  config.setdefault("sideband_group",["Signal"]+AnalysisConfig.sidebands)
        if isinstance(sideband_group,list):
            for sideband in sideband_group:
                data_hist,mc_hist = GetScaledDataMC(config["name"] if "name" in config else config,datafile,mcfile,sideband)
                MakePlot(data_hist,mc_hist,config)
        else:
            #assuing sideband_group is a tuple of name, and list of sidebands
            sideband = sideband_group[0]
            sidebands = sideband_group[1]
            data_hist,mc_hist = GetScaledDataMC(config["name"] if "name" in config else config,datafile,mcfile,sidebands[0])
            for _ in range(1,len(sidebands)):
                data_hist_tmp,mc_hist_tmp = GetScaledDataMC(config["name"] if "name" in config else config,datafile,mcfile,sidebands[_])
                data_hist.Add(data_hist_tmp)
                mc_hist.Add(mc_hist_tmp)
            MakePlot(data_hist,mc_hist,config)

    #make bkg subtracted data histogram

    datafile.Close()
    mcfile.Close()
    scalefile.Close()
