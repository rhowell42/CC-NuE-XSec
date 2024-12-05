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
    def __init__(self,name,is_mc,pseudodata=False):
        self.name = name
        self.is_mc = is_mc
        self.hists = OrderedDict()
        self.ref_hists = OrderedDict()
        self.templates = OrderedDict()
        self.is_pseudodata = pseudodata

        self.hist = None
        self.template = None

    def AddHistogram(self,name,hist,ref_hist,template = None):
        if len(self.hists.keys()) > 0 and type(hist) != type(list(self.hists.values())[0]):
            raise ValueError("Cannot add {} to histogram dictionary of {}".format(type(hist),type(list(self.hists.values())[0])))
        if template is not None and len(templates.keys()) > 0 and type(template) != type(templates.values()[0]):
            raise ValueError("Cannot add {} to template dictionary of {}".format(type(template),type(list(templates.values())[0])))

        self.hists[name] = hist
        self.ref_hists[name] = ref_hist

        if template is not None:
            self.templates[name] = template

    def RemoveHistogram(self,name):
        if name in list(self.hists.keys()):
            del self.hists[name]
            del self.ref_hists[name]
        else:
            raise ValueError("{} not in histogram dictionary".format(name))

        if name in list(self.templates.keys()):
            del self.templates[name]
        elif len(list(self.templates.keys())) > 0:
            raise ValueError("{} not in templates dictionary".format(name))

    def MakeRatio(self,beam):
        beam = beam.lower()
        muonname = beam+"_muselection"
        elecname = beam+"_selection"
        if muonname in self.hists.keys() and elecname in self.hists.keys():
            toClone = self.hists[muonname].Clone()
            toClone.Divide(toClone,self.hists[elecname])
            refClone = self.ref_hists[muonname].Clone()
            refClone.Divide(refClone,self.ref_hists[elecname])
            self.RemoveHistogram(elecname)
            self.AddHistogram('{}_ratio'.format(beam),toClone,refClone)

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
            for h in self.hists:
                if errname in self.hists[h].GetVertErrorBandNames():
                    self.hists[h].PopVertErrorBand(errname)

    def SeparateBeamAngle(self):
        for h in self.hists:
            if "beam_angle" not in self.hists[h].GetVertErrorBandNames():
                continue
            hists = self.hists[h].GetVertErrorBand("beam_angle").GetHists()
            self.hists[h].PopVertErrorBand("beam_angle")
            xhists = hists[:2]
            yhists = hists[2:]
            self.hists[h].AddVertErrorBand("BeamAngleX",xhists)
            self.hists[h].AddVertErrorBand("BeamAngleY",yhists)

    def LateralToVertical(self):
        for h in self.hists:
            for name in self.hists[h].GetLatErrorBandNames():
                hists = self.hists[h].GetLatErrorBand(name).GetHists()
                self.hists[h].PopLatErrorBand(name)
                self.hists[h].AddVertErrorBand(name,hists)

    def Stitch(self):
        # ----- Create empty ROOT histograms to fill with stitched content ----- #
        n_bins_new = 0
        for h in self.ref_hists:
            for i in range(0,self.ref_hists[h].GetNbinsX()+1):
                if self.ref_hists[h].GetBinContent(i) > minBinCont:
                    n_bins_new+=1

        self.hist     = ROOT.TH1D(self.name,self.name,n_bins_new,0,n_bins_new)
        self.template = ROOT.TH2D(self.name+"_template",self.name+"_template", n_bins_new, 0,  n_bins_new,34,0,0.495)

        # ----- Do some errorband cleaning ----- #
        self.SeparateBeamAngle() # make beam angle systematics consistent across samples
        self.LateralToVertical() # convert lateral errorbands (deprecated) from old samples to vertical
        self.SyncErrorBands()  # make sure all samples have the same errorbands, reduce flux universes to 100

        # ----- Combine samples to one histogram ----- #
        self.StitchThis()   # combine 1D histograms

        if len(list(self.templates.keys())) > 0:
            self.StitchThis2D() # combine 2D templates
        
        # ----- Convert ROOT.TH1 to PlotUtils.Mnv1D and fill errorbands ----- #
        self.hist = PlotUtils.MnvH1D(self.hist)
        self.FillErrorBandfromHist2()

    def Write(self,filename):
        f = ROOT.TFile(filename,"UPDATE")
        self.hist.Write()
        self.template.Write()
        f.Close()

    def StitchThis2D(self):
        i_new = 0
        for h in self.hists:
            h_old = self.templates[h]
            for x in range(1, self.ref_hist[h].GetNbinsX()+1):
                if self.ref_hists[h].GetBinContent(x) <= minBinCont:
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
        for h in self.hists:
            h_old = self.hists[h]
            for i in range(1,h_old.GetNbinsX()+1):
                bin_c = h_old.GetBinContent(i)
                if self.is_pseudodata:
                    bin_c = self.ref_hists[h].GetBinContent(i)
                
                if self.ref_hists[h].GetBinContent(i) <= minBinCont:
                    continue

                i_new += 1
                self.hist.SetBinContent(i_new,bin_c)

                # no statistical error on elastic scattering special production
                if 'nueel' in h and self.is_mc:
                    self.hist.SetBinError(i_new,0)

    def FillErrorBandfromHist2(self):
        offset = 1
        for h in self.hists:
            h_old = self.hists[h]
            Nbins = h_old.GetNbinsX()+1 if not self.ref_hists else self.ref_hists[h].GetNbinsX()+1

            errorband_names_vert = h_old.GetVertErrorBandNames()
            n_univ = 0
            sys_bc = 0.0

            for error_band in errorband_names_vert:
                if error_band == 'Flux':
                    n_univ = 100
                else:
                    n_univ = h_old.GetVertErrorBand(error_band).GetNHists()
                    if error_band == "GENIE_MaZExpCCQE":
                        n_univ = 100
                    elif error_band == "GENIE_MaCCQE":
                        n_univ = 2

                if not self.hist.HasVertErrorBand( error_band ) and h_old.HasVertErrorBand( error_band ):
                    self.hist.AddVertErrorBandAndFillWithCV( error_band, n_univ )
                    self.hist.GetVertErrorBand(error_band).SetUseSpreadError( h_old.GetVertErrorBand(error_band).GetUseSpreadError())

                for universe in range(0, n_univ):
                    bin_offset = offset
                    for b in range(1, Nbins):
                        bin_c = h_old.GetBinContent(b)
                        mcbin_c = self.ref_hists[h].GetBinContent(b)
                        sys_bc = h_old.GetVertErrorBand(error_band).GetHist(universe).GetBinContent(b)
                        ratio = sys_bc/bin_c if bin_c != 0 else 0
                        if self.is_pseudodata:
                            sys_bc = mcbin_c*ratio

                        if mcbin_c <= minBinCont:
                            bin_offset += -1
                            continue
                        bin_new = b + bin_offset - 1
                        self.hist.GetVertErrorBand(error_band).GetHist(universe).SetBinContent( bin_new, sys_bc )
            for i in range(1,Nbins):
                if self.ref_hists[h].GetBinContent(i) <= minBinCont:
                    continue
                offset+=1

        for error in self.hist.GetVertErrorBandNames():
            for n in range(self.hist.GetVertErrorBand(error).GetNHists()):
                for i in range(self.hist.GetNbinsX()+1):
                    if self.hist.GetBinContent(i) == 0:
                        self.hist.GetVertErrorBand(error).SetBinContent(i,0)
                        self.hist.GetVertErrorBand(error).GetHist(n).SetBinContent(i,0)

    def SyncHistograms(self,h_sync):
        if type(h_sync) == StitchedHistogram:
            h_sync = h_sync.hist
        elif type(h_sync) != PlotUtils.MnvH1D:
            raise ValueError("Cannot sync {} to MnvH1D".format(type(h_sync)))

        for name in self.hist.GetVertErrorBandNames():
            if name == "Flux":
                if self.hist.GetVertErrorBand(name).GetNHists() > 100:
                    h1_hists = self.hist.GetVertErrorBand(name).GetHists()
                    h1_hists = [h1_hists[i] for i in range(100)]
                    useSpread = self.hist.GetVertErrorBand(name).GetUseSpreadError()
                    errband = self.hist.GetVertErrorBand(name)
                    self.hist.PopVertErrorBand(name)
                    self.hist.AddVertErrorBand(name,h1_hists)
                    self.hist.GetVertErrorBand(name).SetUseSpreadError(useSpread)
                    for i in range(self.hist.GetNbinsX()+1):
                        self.hist.GetVertErrorBand(name).SetBinContent(i,errband.GetBinContent(i))
            elif(name not in h_sync.GetVertErrorBandNames()):
                n_universes = self.hist.GetVertErrorBand(name).GetNHists()
                h_sync.AddVertErrorBandAndFillWithCV(name, n_universes)
                
        for name in h_sync.GetVertErrorBandNames():
            if name == "Flux":
                if h_sync.GetVertErrorBand(name).GetNHists() > 100:
                    h2_hists = h_sync.GetVertErrorBand(name).GetHists()
                    h2_hists = [h2_hists[i] for i in range(100)]
                    useSpread = h_sync.GetVertErrorBand(name).GetUseSpreadError()
                    errband = h_sync.GetVertErrorBand(name)
                    h_sync.PopVertErrorBand(name)
                    h_sync.AddVertErrorBand(name,h2_hists)
                    h_sync.GetVertErrorBand(name).SetUseSpreadError(useSpread)
                    for i in range(h_sync.GetNbinsX()+1):
                        h_sync.GetVertErrorBand(name).SetBinContent(i,errband.GetBinContent(i))
            elif(name not in self.hist.GetVertErrorBandNames()):
                n_universes = h_sync.GetVertErrorBand(name).GetNHists()
                self.hist.AddVertErrorBandAndFillWithCV(name, n_universes)

    def SyncErrorBands(self,rename_bands=True):
        if rename_bands:
            for h in self.hists:
                for name in self.hists[h].GetVertErrorBandNames():
                    if str(name) in errorbandDict.keys():
                        universes = self.hists[h].GetVertErrorBand(name).GetHists()
                        useSpread = self.hists[h].GetVertErrorBand(name).GetUseSpreadError()
                        self.hists[h].AddVertErrorBand(errorbandDict[str(name)],universes)
                        self.hists[h].GetVertErrorBand(errorbandDict[str(name)]).SetUseSpreadError(useSpread)
                        for i in range(self.hists[h].GetNbinsX()+1):
                            self.hists[h].GetVertErrorBand(errorbandDict[str(name)]).SetBinContent(i,self.hists[h].GetVertErrorBand(name).GetBinContent(i))
                        self.hists[h].PopVertErrorBand(name)

        for _h1 in self.hists:
            for _h2 in self.hists:
                h1 = self.hists[_h1]
                h2 = self.hists[_h2]
                if h1 == h2:
                    continue

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
                    elif(name not in h2.GetVertErrorBandNames()):
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
                    elif(name not in h1.GetVertErrorBandNames()):
                        n_universes = h2.GetVertErrorBand(name).GetNHists()
                        h1.AddVertErrorBandAndFillWithCV(name, n_universes)
