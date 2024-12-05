import os
import copy
from collections import OrderedDict
import argparse
import logging, sys
import ROOT
import PlotUtils
from Tools.FitTools import *
from Tools.PlotTools import *
from Tools.StitchTools import *
import numpy as np

import math
from array import array

from config.SignalDef import SWAP_SIGNAL_DEFINITION, SIGNAL_DEFINITION
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS 
from tools import Utilities
from tools.PlotLibrary import HistHolder
ccnueroot = os.environ.get('CCNUEROOT')

minBinCont = 1
errorbandDict = { #keys are the errorbands that need to be renamed, values are what to rename them to
        "GENIE_D2_MaRES":"GENIE_MaRES",
        "GENIE_D2_NormCCRES":"GENIE_NormCCRES",
        "EnergyScale":"ElectronScale",
        "HCALScale":"elE_HCAL",
        "Other_Response":"response_other",
        "Pion_Response":"response_meson",
        "Proton_Response":"response_p",
        "Crosstalk":"response_xtalk",
        "Numu Coh Scale 1":"numucoh1",
        "Numu Coh Scale 2":"numucoh2",
        "Numu Coh Scale 3":"numucoh3",
        "Numu Coh Scale 4":"numucoh4",
        "Numu Coh Scale 5":"numucoh5",
        "Numu Coh Scale 6":"numucoh6",
        "Numu Scale":"Numu CC Scale",
        "RPAHigh":"RPA_HighQ2",
        "RPALow":"RPA_LowQ2",
        "Target_Mass":"Target_Mass_CH",
        }

class StitchedHistogram:
    def __init__(self,name,pseudodata=False):
        self.name = name

        self.data_hists = OrderedDict()
        self.mc_hists = OrderedDict()
        self.templates = OrderedDict()

        self.keys = []

        self.mc_hist = None
        self.data_hist = None
        self.template = None

        self.is_pseudodata = pseudodata
        self.is_processed = False

    def AddHistogram(self,name,mc_hist,data_hist,template=None):
        if len(self.data_hists.keys()) > 0 and type(data_hist) != type(list(self.data_hists.values())[0]):
            raise ValueError("Cannot add {} to data histogram dictionary of {}".format(type(hist),type(list(self.data_hists.values())[0])))
        if len(self.mc_hists.keys()) > 0 and type(mc_hist) != type(list(self.data_hists.values())[0]):
            raise ValueError("Cannot add {} to data histogram dictionary of {}".format(type(hist),type(list(self.mc_hists.values())[0])))
        if template is not None and len(templates.keys()) > 0 and type(template) != type(templates.values()[0]):
            raise ValueError("Cannot add {} to template dictionary of {}".format(type(template),type(list(templates.values())[0])))

        self.data_hists[name] = data_hist
        self.mc_hists[name] = mc_hist

        if template is not None:
            self.templates[name] = template

        self.UpdateKeys()
        self.is_processed = False

    def RemoveHistogram(self,name):
        if name not in list(self.data_hists.keys()):
            raise ValueError("{} not in data histogram dictionary".format(name))
        if name not in list(self.mc_hists.keys()):
            raise ValueError("{} not in mc histogram dictionary".format(name))

        del self.data_hists[name]
        del self.mc_hists[name]
        
        if len(list(self.templates.keys())) > 0:
            if name not in list(self.templates.keys()):
                raise ValueError("{} not in templates dictionary".format(name))
            del self.templates[name]

        self.UpdateKeys()

    def UpdateKeys(self):
        self.keys = list(self.mc_hists.keys())
        if self.keys != list(self.data_hists.keys()):
            raise ValueError("MC dictionary incompatable with data dictionary")

    def MakeRatio(self,beam):
        beam = beam.lower()
        muonname = beam+"_muselection"
        elecname = beam+"_selection"

        if muonname in self.keys and elecname in self.keys:
            mc_clone = self.mc_hists[muonname].Clone()
            mc_clone.Divide(mc_clone,self.mc_hists[elecname])
            data_clone = self.data_hists[muonname].Clone()
            data_clone.Divide(data_clone,self.data_hists[elecname])

            self.RemoveHistogram(elecname)
            self.AddHistogram('{}_ratio'.format(beam),mc_clone,data_clone)

    def ApplyExclusion(self,exclude):
        if "fhc" in exclude:
            self.RemoveHistogram('fhc_muselection')
            self.RemoveHistogram('fhc_selection')
        if "rhc" in exclude:
            self.RemoveHistogram('rhc_muselection')
            self.RemoveHistogram('rhc_selection')
        if "numu" in exclude:
            self.RemoveHistogram('fhc_muselection')
            self.RemoveHistogram('rhc_muselection')
        if "nue" in exclude:
            self.RemoveHistogram('fhc_selection')
            self.RemoveHistogram('rhc_selection')
        if "elastic" in exclude:
            self.RemoveHistogram('fhc_nueel')
            self.RemoveHistogram('rhc_nueel')
        if "imd" in exclude:
            self.RemoveHistogram('fhc_imd')
            self.RemoveHistogram('rhc_imd')

    def CleanErrorBands(self,names=[]):
        self.LateralToVertical()
        self.SeparateBeamAngle()
        self.SyncErrorBands()

        for errname in names:
            for h in self.keys:
                if errname in self.data_hists[h].GetVertErrorBandNames():
                    self.data_hists[h].PopVertErrorBand(errname)
            for h in self.keys:
                if errname in self.mc_hists[h].GetVertErrorBandNames():
                    self.mc_hists[h].PopVertErrorBand(errname)

        self.is_processed = True

    def SeparateBeamAngle(self):
        for h in self.keys:
            if "beam_angle" not in self.data_hists[h].GetVertErrorBandNames():
                continue
            data_hists = self.data_hists[h].GetVertErrorBand("beam_angle").GetHists()
            self.data_hists[h].PopVertErrorBand("beam_angle")
            if "BeamAngleX" not in self.data_hists[h].GetVertErrorBandNames() and "BeamAngleY" not in self.data_hists[h].GetVertErrorBandNames():
                xdata_hists = data_hists[:2]
                ydata_hists = data_hists[2:]
                self.data_hists[h].AddVertErrorBand("BeamAngleX",xdata_hists)
                self.data_hists[h].AddVertErrorBand("BeamAngleY",ydata_hists)

            if "beam_angle" not in self.mc_hists[h].GetVertErrorBandNames():
                continue
            mc_hists = self.mc_hists[h].GetVertErrorBand("beam_angle").GetHists()
            self.mc_hists[h].PopVertErrorBand("beam_angle")
            if "BeamAngleX" not in self.mc_hists[h].GetVertErrorBandNames() and "BeamAngleY" not in self.mc_hists[h].GetVertErrorBandNames():
                xmc_hists = mc_hists[:2]
                ymc_hists = mc_hists[2:]
                self.mc_hists[h].AddVertErrorBand("BeamAngleX",xmc_hists)
                self.mc_hists[h].AddVertErrorBand("BeamAngleY",ymc_hists)

    def LateralToVertical(self):
        for h in self.keys:
            for name in self.data_hists[h].GetLatErrorBandNames():
                data_hists = self.data_hists[h].GetLatErrorBand(name).GetHists()
                self.data_hists[h].PopLatErrorBand(name)
                if h not in self.data_hists[h].GetVertErrorBandNames(): 
                    self.data_hists[h].AddVertErrorBand(name,data_hists)
            for name in self.mc_hists[h].GetLatErrorBandNames():
                mc_hists = self.mc_hists[h].GetLatErrorBand(name).GetHists()
                self.mc_hists[h].PopLatErrorBand(name)
                if h not in self.mc_hists[h].GetVertErrorBandNames(): 
                    self.mc_hists[h].AddVertErrorBand(name,mc_hists)

    def Stitch(self):
        # ----- Create empty ROOT histograms to fill with stitched content ----- #
        n_bins_new = 0
        for h in self.keys:
            for i in range(0,self.mc_hists[h].GetNbinsX()+1):
                if self.mc_hists[h].GetBinContent(i) > minBinCont:
                    n_bins_new+=1

        self.mc_hist   = ROOT.TH1D(self.name+"_mc",self.name+"_mc",n_bins_new,0,n_bins_new)
        self.data_hist = ROOT.TH1D(self.name+"_data",self.name+"_data",n_bins_new,0,n_bins_new)
        self.template  = ROOT.TH2D(self.name+"_template",self.name+"_template",n_bins_new,0,n_bins_new,34,0,0.495)

        # ----- Do some errorband cleaning ----- #
        if not self.is_processed:
            self.SeparateBeamAngle() # make beam angle systematics consistent across samples
            self.LateralToVertical() # convert lateral errorbands (deprecated) from old samples to vertical
            self.SyncErrorBands()  # make sure all samples have the same errorbands, reduce flux universes to 100
            self.is_processed = True

        # ----- Combine samples to one histogram ----- #
        self.StitchThis()   # combine 1D histograms
        if len(list(self.templates.keys())) > 0:
            self.StitchThis2D() # combine 2D templates
        
        # ----- Convert ROOT.TH1 to PlotUtils.Mnv1D and fill errorbands ----- #
        self.mc_hist = PlotUtils.MnvH1D(self.mc_hist)
        self.data_hist = PlotUtils.MnvH1D(self.data_hist)
        self.FillErrorBandsFromDict()

    def Write(self,filename):
        f = ROOT.TFile(filename,"UPDATE")
        self.mc_hist.Write()
        self.data_hist.Write()
        self.template.Write()
        f.Close()

    def StitchThis2D(self):
        i_new = 0
        for h in self.templates:
            h_old = self.templates[h]
            for x in range(1, self.mc_hists[h].GetNbinsX()+1):
                if self.mc_hists[h].GetBinContent(x) <= minBinCont:
                    continue

                i_new += 1
                colInt = 0
                for c in range(1,h_old.GetNbinsX()+1):
                    if isinstance(h_old,ROOT.TH1D):
                        colInt+=h_old.GetBinContent(c)
                    else:
                        colInt+=h_old.GetBinContent(c,x)

                for c in range(1,h_old.GetNbinsX()+1):
                    if isinstance(h_old,ROOT.TH1D):
                        bin_c = h_old.GetBinContent(c)
                    else:
                        bin_c = h_old.GetBinContent(c,x)
                    self.template.SetBinContent(i_new, c, bin_c/colInt if colInt > 0 else 0)

    def StitchThis(self):
        i_new = 0
        for h in self.keys:
            h_mc = self.mc_hists[h]
            h_data = self.data_hists[h]

            for i in range(1,h_mc.GetNbinsX()+1):
                if h_mc.GetBinContent(i) <= minBinCont:
                    continue # skip empty MC bins
                i_new += 1

                # ----- do MC stitching ----- #
                bin_c = h_mc.GetBinContent(i)
                self.mc_hist.SetBinContent(i_new,bin_c)
                # no statistical error on elastic scattering special production
                if 'nueel' in h:
                    self.mc_hist.SetBinError(i_new,0)

                # ----- do data stitching ----- #
                bin_c = h_data.GetBinContent(i)
                if self.is_pseudodata:
                    bin_c = h_mc.GetBinContent(i)                
                self.data_hist.SetBinContent(i_new,bin_c)

    def FillErrorBandsFromDict(self):
        offset = 1
        for h in self.keys:
            h_mc = self.mc_hists[h]
            h_data = self.data_hists[h]
            Nbins = h_mc.GetNbinsX()+1

            if not self.is_processed:
                raise ValueError("Histograms have not been synced")

            errorband_names_vert = h_mc.GetVertErrorBandNames()

            n_univ = 0
            sys_mc = 0.0
            sys_data = 0.0

            for error_band in errorband_names_vert:
                n_univ = h_mc.GetVertErrorBand(error_band).GetNHists()
                if not self.mc_hist.HasVertErrorBand(error_band) and h_mc.HasVertErrorBand(error_band):
                    self.mc_hist.AddVertErrorBandAndFillWithCV(error_band, n_univ)
                    self.mc_hist.GetVertErrorBand(error_band).SetUseSpreadError(h_mc.GetVertErrorBand(error_band).GetUseSpreadError())
                if not self.data_hist.HasVertErrorBand(error_band) and h_data.HasVertErrorBand(error_band):
                    self.data_hist.AddVertErrorBandAndFillWithCV(error_band, n_univ)
                    self.data_hist.GetVertErrorBand(error_band).SetUseSpreadError(h_data.GetVertErrorBand(error_band).GetUseSpreadError())

                for universe in range(0, n_univ):
                    bin_offset = offset
                    for b in range(1, Nbins):
                        bin_mc = h_mc.GetBinContent(b)
                        bin_data = h_data.GetBinContent(b)

                        if bin_mc <= minBinCont:
                            bin_offset += -1
                            continue

                        sys_mc = h_mc.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)
                        sys_data = h_data.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)

                        if self.is_pseudodata:
                            ratio = sys_data/bin_data if bin_data != 0 else 0
                            sys_data = bin_mc*ratio

                        bin_new = b + bin_offset - 1
                        self.mc_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(bin_new,sys_mc)
                        self.data_hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent(bin_new,sys_data)

            for i in range(1,Nbins):
                if self.mc_hists[h].GetBinContent(i) <= minBinCont:
                    continue
                offset+=1
        if False:
            for error in self.mc_hist.GetVertErrorBandNames():
                for n in range(self.mc_hist.GetVertErrorBand(error).GetNHists()):
                    for i in range(self.mc_hist.GetNbinsX()+1):
                        if self.mc_hist.GetBinContent(i) == 0:
                            self.mc_hist.GetVertErrorBand(error).SetBinContent(i,0)
                            self.mc_hist.GetVertErrorBand(error).GetHist(n).SetBinContent(i,0)

    def SyncHistograms(self,h_sync):
        if type(h_sync) == StitchedHistogram:
            mc_sync = h_sync.mc_hist
            data_sync = h_sync.data_hist
            if mc_sync.GetVertErrorBandNames() != self.mc_hist.GetVertErrorBandNames():
                self.Sync(mc_sync,self.mc_hist)
            if data_sync.GetVertErrorBandNames() != self.data_hist.GetVertErrorBandNames():
                self.Sync(data_sync,self.data_hist)
        elif type(h_sync) == PlotUtils.MnvH1D:
            if h_sync.GetVertErrorBandNames() != self.mc_hist.GetVertErrorBandNames():
                self.Sync(h_sync,self.mc_hist)
            if h_sync.GetVertErrorBandNames() != self.data_hist.GetVertErrorBandNames():
                self.Sync(h_sync,self.mc_hist)
        else:
            raise ValueError("Cannot sync {} to MnvH1D or StitchedHistogram".format(type(h_sync)))

    def SyncErrorBands(self,rename_bands=True):
        if rename_bands:
            for h in self.keys:
                self.RenameBands(self.mc_hists[h])
                self.RenameBands(self.data_hists[h])

        for _h1 in self.keys:
            for _h2 in self.keys:
                h1 = self.data_hists[_h1]
                h2 = self.data_hists[_h2]
                if h1 != h2:
                    self.Sync(h1,h2)

                h1 = self.mc_hists[_h1]
                h2 = self.mc_hists[_h2]
                if h1 != h2:
                    self.Sync(h1,h2)

            h1 = self.mc_hists[_h1]
            h2 = self.data_hists[_h1]
            self.Sync(h1,h2)

    def Sync(self,h1,h2):
        for name in h1.GetVertErrorBandNames():
            if name == "Flux":
                if h1.GetVertErrorBand(name).GetNHists() > 100:
                    h1_hists = h1.GetVertErrorBand(name).GetHists()
                    h1_hists = [h1_hists[i] for i in range(100)]
                    useSpread = h1.GetVertErrorBand(name).GetUseSpreadError()
                    errband = h1.GetVertErrorBand(name)
                    h1.PopVertErrorBand(name)
                    h1.AddVertErrorBand(name,h1_hists)
                    h1.GetVertErrorBand(name).SetUseSpreadError(useSpread)
                    for i in range(h1.GetNbinsX()+1):
                        h1.GetVertErrorBand(name).SetBinContent(i,errband.GetBinContent(i))
            if(name not in h2.GetVertErrorBandNames()):
                n_universes = h1.GetVertErrorBand(name).GetNHists()
                h2.AddVertErrorBandAndFillWithCV(name, n_universes)
                
        for name in h2.GetVertErrorBandNames():
            if name == "Flux":
                if h2.GetVertErrorBand(name).GetNHists() > 100:
                    h2_hists = h2.GetVertErrorBand(name).GetHists()
                    h2_hists = [h2_hists[i] for i in range(100)]
                    useSpread = h2.GetVertErrorBand(name).GetUseSpreadError()
                    errband = h2.GetVertErrorBand(name)
                    h2.PopVertErrorBand(name)
                    h2.AddVertErrorBand(name,h2_hists)
                    h2.GetVertErrorBand(name).SetUseSpreadError(useSpread)
                    for i in range(h2.GetNbinsX()+1):
                        h2.GetVertErrorBand(name).SetBinContent(i,errband.GetBinContent(i))
            if(name not in h1.GetVertErrorBandNames()):
                n_universes = h2.GetVertErrorBand(name).GetNHists()
                h1.AddVertErrorBandAndFillWithCV(name, n_universes)

    def RenameBands(self,h1):
        for name in h1.GetVertErrorBandNames():
            if str(name) in errorbandDict.keys():
                universes = h1.GetVertErrorBand(name).GetHists()
                useSpread = h1.GetVertErrorBand(name).GetUseSpreadError()
                h1.AddVertErrorBand(errorbandDict[str(name)],universes)
                h1.GetVertErrorBand(errorbandDict[str(name)]).SetUseSpreadError(useSpread)
                for i in range(h1.GetNbinsX()+1):
                    h1.GetVertErrorBand(errorbandDict[str(name)]).SetBinContent(i,h1.GetVertErrorBand(name).GetBinContent(i))
                h1.PopVertErrorBand(name)
