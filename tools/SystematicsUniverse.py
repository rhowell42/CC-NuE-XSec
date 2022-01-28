"""
    SystematicsUniverse.py:
    Define systematics shift universe. 
"""

import math
import os,sys
import random
import itertools
import operator
import random
import re

import ROOT
import PlotUtils
import ctypes

from config import SystematicsConfig
from tools import Utilities
from collections import OrderedDict
from tools import pcweight
from .pcweight import GetModelWeight,MyWeighter

OneSigmaShift = [-1.0,1.0]
M_e = 0.511 # MeV
M_mu =  105.658 #MeV
M_e_sqr = M_e**2
M_mu_sqr = M_mu**2

#Alias:
#ROOT.PythonMinervaUniverse = PlotUtils.DefaultCVUniverse

#MethodName: UpperCaseCamel
#universeShortName: lower_case_underscore

# The base universe, define functions to be used for all universes.
class CVUniverse(ROOT.PythonMinervaUniverse, object):
    is_pc = False
    def __init__(self, chain, nsigma = None): #nsigma is None for data because we don't shift data
        super(CVUniverse,self).__init__(chain,0 if nsigma is None else nsigma)
        self.weight = None
        self.tuning_weight = None
        self.InitWithoutSuper(chain,nsigma)

    def InitWithoutSuper(self,chain,nsigma):
        self.nEntries = chain.GetEntries()
        self.mc = nsigma is not None
        self.chain=chain

    # magic function that allow direct access to TBranch by Universe.TBranch, but prohibit any python magic variable being generated.
    def __getattr__(self,attrName):
        def makelist(x):
            return list(makelist(_) for _ in x) if hasattr(x, "__getitem__") else x

        if attrName == "kin_cal" or attrName == "classifier":
            #the two attrs can't be in chainwrapper 
            raise AttributeError()

        # determine the dimension of leaf.
        try:
            return makelist(self.LeafGetters[attrName](attrName))
        except KeyError:
            #havn't see this function. Got to figure out how to get it.
            i = ctypes.c_int(1)
            leaf = self.chain.GetLeaf(attrName,-1)
            type = leaf.GetTypeName()
            if leaf.GetLeafCounter(i) or i.value > 1:
                dimension=1
            else:
                dimension=type.count("vector")
            # Assuming type names are Double_t,Int_t,Bool_t,vector<vector<Double_t>>, use regex to match the type
            regex="(<|^)([A-Za-z]+)(_t|>|$)"
            type = re.search(regex,type).group(2).capitalize()

            self.LeafGetters[attrName] = getattr(self, "Get{}{}".format("" if dimension==0 else ("Vec" if dimension==1 else "VecOfVec"),type))

        return makelist(self.LeafGetters[attrName](attrName))

    # return text when you print this class object. must reimplement or __getattr__ will hijack the method resolution
    def __repr__(self):
        return self.LatexName()+" Universe at sigma: "+ str(self.GetSigma())

    @property
    def n_prongs(self):
        return self.GetInt("n_prongs")

    # keep a reference to kinematicsCalulator and eventClassifier such that other function can pull info from universe
    def LoadTools(self,kinematicCalculator,eventClassifier):
        self.kin_cal = kinematicCalculator
        self.classifier = eventClassifier

    def SetEntry(self,n_entry):
        #go to another event, discard weight calculated.
        self.weight = None
        self.tuning_weight = None
        CVUniverse.LLR = None
        super(CVUniverse,self).SetEntry(n_entry)

    def SetLeptonType(self):
        if self.GetAnaToolName() != "MasterAnaDev" or self.HasNoBackExitingTracks:
            self.LeptonTheta = self.ElectronTheta
            self.LeptonEnergy = self.ElectronEnergy
            self.LeptonP3D = self.ElectronP3D
            self.M_lep_sqr =  M_e_sqr
            self.GetCorrection = self.GetLeakageCorrection
            return True
        else:
            self.LeptonTheta = self.GetThetamu
            self.LeptonEnergy  = self.GetEmu
            self.LeptonP3D = self.GetPmu
            self.M_lep_sqr =  M_mu_sqr
            self.GetCorrection = self.GetNuEFuzz
            return False

    @property
    def nsigma(self):
        return self.GetSigma() if self.mc else None

    def GetWeight(self,bkgtuning=True):
        if self.weight is None:
            if self.nsigma is None:
                self.weight = 1
            else:
                self.weight = self.GetModelWeight() if self.is_pc else self.GetStandardWeight()
        if bkgtuning or self.nsigma is None:
            if self.tuning_weight is None:
                self.tuning_weight = MyWeighter.GetWeight(self)
            return self.weight *self.tuning_weight
        else:
            return self.weight

    def GetStandardWeight(self):
        weight = 1.0
        weight *=self.GetGenieWeight()
        weight *=self.GetFluxAndCVWeight()
        weight *=self.GetLowRecoil2p2hWeight()
        weight *=self.GetRPAWeight()
        weight *=self.GetMyLowQ2PiWeight()
        weight *= self.GetGeantHadronWeight()
      
        # for _ in self.weighters:
        #     weight *= _()
        return weight

    #def GetPCWeight(self):
    #    return 1.0

    def GetMyLowQ2PiWeight(self):
        channel = SystematicsConfig.LowQ2PiWeightChannel
        if channel is None:
            return 1
        else:
            return super(CVUniverse,self).GetLowQ2PiWeight(channel.upper())

    def GetModelWeight(self):
        w = GetModelWeight(self, 'Eel',"Pi0") # (Ee, Theta, ETh), (PCElectron, PCPhoton, PCPi0)
        return w



    #################### functions to be overide in systematics shift universes #################

    def ElectronEnergyRaw(self):
        return self.GetVecElem("prong_part_E",0,3)

    @Utilities.decorator_ReLU
    def ElectronEnergy(self):
        return self.ElectronEnergyRaw() + self.GetEMEnergyShift() +self.GetLeakageCorrection()

    def ElectronP3D(self):
        allp = self.GetVecOfVecDouble("prong_part_E")
        #print([allp[0][i] for i in range(4)])
        scale = self.GetEMEnergyShift()/allp[0][3] if (allp[0][3]>0) else 0
        p = ROOT.Math.XYZVector(*tuple(list(allp[0])[:3]))*(1+scale)
        r = ROOT.Math.RotationX(SystematicsConfig.BEAM_ANGLE)
        s=r(p)
        if s is None:
            print (p,r,tuple(list(allp[0])[:3]))
        return s

    def ElectronTheta(self):
        return self.ElectronP3D().Theta()

    def ConeInsideE(self): 
        e = self.GetVecElem("prong_TrimmedCaloEnergy",0)/1e3
        psi = self.GetVecElem("Psi",0)
        #print e, 'inside', psi, 'psi'
        return e*psi

    def ConeOutsideE(self):
        localVec = list(self.GetVecDouble("ExtraEnergy"))
        e = sum(localVec)
	#print e, 'out'
        return e

    @staticmethod
    def Scale4VectorMomentum(vet,scale):
        return vet

    def PrimaryProtonScore(self):
        score = self.MasterAnaDev_proton_score 
        return score if score>=0 else -1

    def PrimaryProtonTheta(self):
        theta = self.MasterAnaDev_proton_theta
        return math.degrees(theta) if theta>=0 else -1

    @property
    def UpstreamInlineEnergy(self):
        uie = self.GetDouble("UpstreamInlineEnergy")
        vtx = 0
        #vtx = sum(self.GetVecElem("Bined_UIE",i) for i in range(2))
        vtx = self.GetDouble("recoil_energy_nonmuon_vtx50mm")
        return max(0,uie-vtx)

    # =============== collcetion of all recoil energy stuff.========
    def GetEAvail(self):
        if self.GetAnaToolName() == "MasterAnaDev":
            return self.blob_recoil_E_tracker+self.blob_recoil_E_ecal
        else:
            return self.recoile_passive_tracker+self.recoile_passive_ecal

    def GetNuEFuzz(self):
        try:
            fuzz = self.blob_nuefuzz_E_tracker+self.blob_nuefuzz_E_ecal
        except:
            fuzz = 0
        #return 0
        return max(0,fuzz) * SystematicsConfig.AVAILABLE_E_CORRECTION

    def GetLeakageCorrection(self):
        return SystematicsConfig.LEAKAGE_CORRECTION(self.ElectronEnergyRaw())

    def Pi0_Additional_Leakage(self):
        return random.uniform(0,20) if self.mc and (self.mc_intType == 4 or (self.mc_current==2 and self.mc_intType==10)) else 0

    # A function wrapper that maps negative values to zero
    # returning MeV
    @Utilities.decorator_ReLU
    def AvailableEnergy(self):
        #return self.GetEAvail()* SystematicsConfig.AVAILABLE_E_CORRECTION
        return (self.GetEAvail() -self.GetCorrection()) * SystematicsConfig.AVAILABLE_E_CORRECTION

    @Utilities.decorator_ReLU
    def RecoilEnergy(self):
        result = None
        if self.GetAnaToolName() == "MasterAnaDev":
            result = sum(self.__getattr__("blob_recoil_E_{}".format(i)) for i in ["tracker","ecal","hcal","od","nucl"])
        else:
            result = sum(self.__getattr__("recoile_passive_{}".format(i)) for i in ["tracker","ecal","hcal","odclus","nucl"])
        return result - self.GetCorrection()

    def GetFuzzCorrection(self):
        return 0

    def GetEMEnergyShift(self):
        return self.prong_ECALCalibE[0]*SystematicsConfig.EM_ENERGY_SCALE_SHIFT_ECAL if self.nsigma is None else 0

    def dedxDelta(self):
        dedx = self.GetVecDouble("dEdXs")
        dedx_dz  = self.GetVecDouble("dEdXs_dz")
        dedx_proj  = self.GetVecDouble("dEdXs_proj")
        total = 0
        outer = 0
        if dedx.size() == 0:
            return 0.9999
        for i in range(0,dedx.size()):
            if dedx_proj[i]>400:
                outer+=dedx[i]
            total+=dedx[i]
        return outer/total 

    # def GetLLR(self):
    #     if CVUniverse.LLR is not None:
    #         return CVUniverse.LLR
    #     else:
    #         CVUniverse.LLR = self.CalculateLLR()
    #         return CVUniverse.LLR

    # def CalculateLLR(self):
    #     dedx = self.GetVecOfVecDouble("prong_dEdXs")
    #     dedx_dz  = self.GetVecOfVecDouble("prong_dEdXs_dz")
    #     dEdX_sum = {}
    #     LLR = 0
    #     for i in range(min(dedx[0].size(),dedx_dz[0].size())):
    #         try: 
    #             dEdX_sum[dedx_dz[0][i]-dedx_dz[0][0]] += dedx[0][i]
    #         except KeyError:
    #             dEdX_sum[dedx_dz[0][i]-dedx_dz[0][0]] = dedx[0][i]

    #     #print dEdX_sum
    #     for i in range(0,20):
    #         dEdXstep = dEdX_sum.setdefault(i,0)*math.cos(self.prong_axis_vector[0][0])
    #         LLR += self.GetLLRStep(i,dEdXstep,self.ElectronEnergy()/1e3)

    #     return LLR

    # @staticmethod
    # def GetLLRStep(dplane,dEdX,energy):
    #     energy_index = max(2,min(11,int(math.ceil(energy))))
    #     try:
    #         prob_e = CVUniverse.GetLLRStep.emap[dplane*12+energy_index]
    #         prob_g = CVUniverse.GetLLRStep.gmap[dplane*12+energy_index]
    #     except AttributeError:
    #         print "get LLR histograms."
    #         CVUniverse.GetLLRStep.f = ROOT.TFile.Open("{}/macros/Low_Recoil/studies/hist_IA.root".format(os.environ['CCNUEROOT']))
    #         CVUniverse.GetLLRStep.f2 = ROOT.TFile.Open("{}/macros/Low_Recoil/studies/hist.root".format(os.environ['CCNUEROOT']))
    #         CVUniverse.GetLLRStep.emap ={}
    #         CVUniverse.GetLLRStep.gmap ={}

    #         for i in range(0,30):
    #             for j in range(1,12):
    #                 CVUniverse.GetLLRStep.emap[i*12+j]= CVUniverse.GetLLRStep.f2.Get("elikelihood{}_{}".format(i,j))
    #                 CVUniverse.GetLLRStep.gmap[i*12+j]= CVUniverse.GetLLRStep.f2.Get("glikelihood{}_{}".format(i,j))

    #         return CVUniverse.GetLLRStep(dplane,dEdX,energy)

    #     except KeyError:
    #         #print dplane, energy_index, energy
    #         #print self.event.ShortName()
    #         return 0

    #     pe = max(prob_e.GetBinContent(prob_e.FindBin(dEdX)),1e-5)
    #     pg = max(prob_g.GetBinContent(prob_g.FindBin(dEdX)),1e-5)
    #     #print pe,pg

    #     return -math.log(pe)+math.log(pg)

# The order of inheritance decide the method resolution order, thus PlotUtils.XXXXUniverse must be the first.
class FluxUniverse( ROOT.PlotUtils.FluxUniverse(ROOT.PythonMinervaUniverse),CVUniverse, object):

    def __init__(self,chain,universe_number):
        super(FluxUniverse,self).__init__(chain,1,universe_number)
        super(ROOT.PlotUtils.FluxUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,1)

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [FluxUniverse(chain,  i ) for i in range(0,SystematicsConfig.NUM_FLUX_UNIVERSE)]

class GenieUniverse(ROOT.PlotUtils.GenieUniverse(ROOT.PythonMinervaUniverse), CVUniverse, object):
    def __init__(self,chain,nsigma,universe_name):
        super(GenieUniverse,self).__init__(chain,nsigma,universe_name)
        super(ROOT.PlotUtils.GenieUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieUniverse(chain,i,j) for j in SystematicsConfig.GENIE_UNIVERSES for i in OneSigmaShift]

class GenieRvx1piUniverse(ROOT.PlotUtils.GenieRvx1piUniverse(ROOT.PythonMinervaUniverse),CVUniverse,object):
    def __init__(self,chain,nsigma,universe_name):
        super(GenieRvx1piUniverse,self).__init__(chain,nsigma,universe_name)
        super(ROOT.PlotUtils.GenieRvx1piUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieRvx1piUniverse(chain,i,j) for j in ["Rvn1pi","Rvp1pi"] for i in OneSigmaShift]

class GenieFaCCQEUniverse(ROOT.PlotUtils.GenieFaCCQEUniverse(ROOT.PythonMinervaUniverse),CVUniverse,object):
    def __init__(self,chain,nsigma,universe_number):
        super(GenieFaCCQEUniverse,self).__init__(chain,nsigma,universe_number)
        super(ROOT.PlotUtils.GenieFaCCQEUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        if ROOT.PythonMinervaUniverse.UseZExpansionFaReweight():
            return [GenieFaCCQEUniverse(chain,1.0,i) for i in range(0,SystematicsConfig.NumZExpansionUniverses)]
        else:
            return [GenieFaCCQEUniverse(chain,i,-1) for i in OneSigmaShift]


# I don't think nu_e analysis need minos shift
# class MinosUniverse(ROOT.PlotUtils.GenieUniverse(ROOT.PythonMinervaUniverse), CVUniverse, object):

class Universe2p2h(ROOT.PlotUtils.Universe2p2h(ROOT.PythonMinervaUniverse), CVUniverse, object):
    def __init__(self,chain,universe_number):
        super(Universe2p2h,self).__init__(chain,1,universe_number)
        super(ROOT.PlotUtils.Universe2p2h(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,1)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [Universe2p2h(chain, i) for i in SystematicsConfig.UNIVERSES_2P2H]

class RPAUniverse(ROOT.PlotUtils.RPAUniverse(ROOT.PythonMinervaUniverse), CVUniverse, object):

    def __init__(self,chain, universe_number,q2_region):
        super(RPAUniverse,self).__init__(chain,1,universe_number,q2_region)
        super(ROOT.PlotUtils.RPAUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain, 1)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [RPAUniverse(chain, i,j) for j in SystematicsConfig.RPA_UNIVERSES for i in SystematicsConfig.RPA_UNIVERSES[j] ]

class ResponseUniverse(ROOT.PlotUtils.ResponseUniverse(ROOT.PythonMinervaUniverse), CVUniverse, object):
    def __init__(self,chain,nsigma,name):
        super(ResponseUniverse,self).__init__(chain,nsigma,name)
        super(ROOT.PlotUtils.ResponseUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,1)
        self.re = "^blob_.*_E_(tracker|ecal|od|nucl|hcal)$"


    def __getattr__(self,attrName):
        if re.match(self.re,attrName) is None:
            return super(ResponseUniverse,self).__getattr__(attrName)
        else:
            r = super(ResponseUniverse,self).__getattr__(attrName)
            shift = super(ResponseUniverse,self).__getattr__("{}_{}".format(attrName,self.m_name))* self.m_frac_unc
            return r+self.nsigma*shift


    def RecoilClustersCalorimetry(self,splined):
        return super(ResponseUniverse,self).RecoilClustersCalorimetry(splined)+self.GetRecoilShift()

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [ResponseUniverse(chain,i,j) for j in SystematicsConfig.RESPONSE_BRANCHES for i in OneSigmaShift]

class LowQ2PionUniverse(ROOT.PlotUtils.LowQ2PionUniverse(ROOT.PythonMinervaUniverse),CVUniverse, object):
    def __init__(self,chain,nsigma,channel):
        super(LowQ2PionUniverse,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.LowQ2PionUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)
        self.channel = channel.upper() if channel is not None else None

    def GetMyLowQ2PiWeight(self):
        if self.nsigma == 0 or self.channel is None:
            return 1
        else :
            return super(LowQ2PionUniverse,self).GetLowQ2PiWeight(self.channel)

    @staticmethod
    def GetSystematicsUniverses(chain):
        cvshifts = [LowQ2PionUniverse(chain,i,SystematicsConfig.LowQ2PiWeightChannel) for i in OneSigmaShift] if SystematicsConfig.LowQ2PiWeightChannel is not None else []
        cvshifts.extend([LowQ2PionUniverse(chain,i,j) for j in SystematicsConfig.LowQ2PiWeightSysChannel for i in ([-1,1] if j is not None else [0])])
        return cvshifts

class MuonUniverseMinerva(ROOT.PlotUtils.MuonUniverseMinerva(ROOT.PythonMinervaUniverse),CVUniverse, object):
    def __init__(self,chain,nsigma):
        super(MuonUniverseMinerva,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.MuonUniverseMinerva(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)


    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonUniverseMinerva(chain,i) for i in OneSigmaShift]


class MuonUniverseMinos(ROOT.PlotUtils.MuonUniverseMinos(ROOT.PythonMinervaUniverse),CVUniverse, object):
    def __init__(self,chain,nsigma):
        super(MuonUniverseMinos,self).__init__(chain,nsigma)
        super(ROOT.PlotUtils.MuonUniverseMinos(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)


    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonUniverseMinos(chain,i) for i in OneSigmaShift]

###########################################################################
class ElectronEnergyShiftUniverse(CVUniverse):
    def __init__(self,chain, nsigma,region):
        super(ElectronEnergyShiftUniverse,self).__init__(chain, nsigma)
        self.region = region

    def GetEMEnergyShift(self):
        return self.nsigma*SystematicsConfig.EM_ENERGY_SCALE_UNCERTAINTY[self.region]*self.GetVecElem("prong_"+self.region+"CalibE",0)

    def ShortName(self):
        return "elE_"+self.region

    def LatexName(self):
        return "EM Energy Scale in "+ self.region

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [ElectronEnergyShiftUniverse(chain, i,region) for region in SystematicsConfig.EM_ENERGY_SCALE_UNCERTAINTY for i in OneSigmaShift]


###########################################################################
class ElectronAngleShiftUniverse(CVUniverse):
    def __init__(self,chain, nsigma):
        super(ElectronAngleShiftUniverse,self).__init__(chain, nsigma)
        self.axis_angle = random.random()*math.pi
        self.shift_angle = random.gauss(0,self.nsigma*SystematicsConfig.ELECTRON_ANGLE_UNCERTAINTY)

    def ShortName(self):
        return "eltheta"

    def LatexName(self):
        return "Election Candidate Angle"

    def ElectronP3D(self):
        p = super(ElectronAngleShiftUniverse,self).ElectronP3D()
        unit_3_vector = p/p.R()
        normal_vector1 = p.Cross(ROOT.Math.XYZVector(1.0,0,0))
        if normal_vector1.R() != 0:
            normal_vector1 *= 1.0/normal_vector1.R()
        r1 = ROOT.Math.AxisAngle(unit_3_vector,self.axis_angle)
        rotation_axis = r1(normal_vector1)
        r2 = ROOT.Math.AxisAngle(rotation_axis,self.shift_angle)
        return r2(p)

    def ElectronTheta2(self):
        theta = super(ElectronAngleShiftUniverse,self).ElectronTheta()
        return theta+self.shift_angle

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [ElectronAngleShiftUniverse(chain,  1) for i in range(SystematicsConfig.NUM_UNIVERSE)]

    

###########################################################################
class BirksShiftUniverse(CVUniverse):
    def __init__(self,chain, nsigma):
        super(BirksShiftUniverse,self).__init__(chain, nsigma)

    @property
    def prong_part_score(self):
        return self.prong_part_score_BirkUP if self.nsigma>0 else self.prong_part_score_BirkDOWN

    @property
    def prong_dEdXMeanFrontTracker(self):
        return self.prong_dEdXMeanFrontTracker_BirkUP if self.nsigma>0 else self.prong_dEdXMeanFrontTracker_BirkDOWN

    @property
    def prong_NonMIPClusFrac(self):
        return self.prong_NonMIPClusFrac_BirkUP if self.nsigma>0 else self.prong_NonMIPClusFrac_BirkDOWN

    def ShortName(self):
        return "birks"

    def LatexName(self):
        return "Birk's Constant"

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [BirksShiftUniverse(chain, i) for i in OneSigmaShift]

###########################################################################

class BeamAngleShiftUniverse(CVUniverse):
    def __init__(self,chain, nsigma, x):
        super(BeamAngleShiftUniverse,self).__init__(chain, nsigma)
        self.rotation =  ROOT.Math.RotationX(SystematicsConfig.BEAM_XANGLE_UNCERTAINTY*nsigma) if x else ROOT.Math.RotationY(SystematicsConfig.BEAM_YANGLE_UNCERTAINTY*nsigma)

    def ShortName(self):
        return "beam_angle"

    def LatexName(self):
        return "Beam Angle"

    def ElectronP3D(self):
        p = super(BeamAngleShiftUniverse,self).ElectronP3D()
        return self.rotation(p)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [BeamAngleShiftUniverse(chain, i, j) for j in [True,False] for i in OneSigmaShift]


class GeantHadronUniverse(ROOT.PlotUtils.GeantHadronUniverse(ROOT.PythonMinervaUniverse), CVUniverse, object):
   def __init__(self,chain,nsigma,pdg):
       super(GeantHadronUniverse,self).__init__(chain,nsigma,pdg)
       super(ROOT.PlotUtils.GeantHadronUniverse(ROOT.PythonMinervaUniverse),self).InitWithoutSuper(chain,nsigma)

   @staticmethod
   def GetSystematicsUniverses(chain):
       return [GeantHadronUniverse(chain,i,j) for j in SystematicsConfig.GEANT_PARTICLES for i in OneSigmaShift]

class MKModelUniverse(CVUniverse,object):
    def __init__(self,chain,nsigma):
        super(MKModelUniverse,self).__init__(chain,nsigma)
        self.reweighter = ROOT.PlotUtils.MKReweighter(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)()

    def IsVerticalOnly(self):
        return True

    def GetStandardWeight(self):
        weight = super(MKModelUniverse,self).GetStandardWeight()
        weight*= self.reweighter.GetWeight(self,ROOT.PlotUtils.detail.empty())
        return weight

    def ShortName(self):
        return "MK_model"

    def LatexName(self):
        return "MK Model"

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MKModelUniverse(chain,1)]


class FSIWeightUniverse(CVUniverse,object):
    def __init__(self,chain,nsigma,iweight):
        super(FSIWeightUniverse,self).__init__(chain,nsigma)
        self.reweighter = ROOT.PlotUtils.FSIReweighter(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)((iweight+1)//2,(iweight+1)%2)

    def IsVerticalOnly(self):
        return True

    def GetStandardWeight(self):
        weight = super(FSIWeightUniverse,self).GetStandardWeight()
        weight*= self.reweighter.GetWeight(self,ROOT.PlotUtils.detail.empty())
        return weight

    def ShortName(self):
        return "fsi_weight"

    def LatexName(self):
        return "FSI weight"

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [FSIWeightUniverse(chain,1,i) for i in range(3)]

class SusaValenciaUniverse(CVUniverse,object):
    def __init__(self,chain,nsigma):
        super(SusaValenciaUniverse,self).__init__(chain,nsigma)
        self.reweighter = ROOT.PlotUtils.SuSAFromValencia2p2hReweighter(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)()

    def IsVerticalOnly(self):
        return True

    def GetLowRecoil2p2hWeight(self):
        return 1.0

    def GetRPAWeight(self):
        return 1.0

    def ShortName(self):
        return "SuSA_Valencia_Weight"

    def LatexName(self):
        return "SuSA Valencia Weight"

    def GetStandardWeight(self):
        weight = super(SusaValenciaUniverse,self).GetStandardWeight()
        weight*= self.reweighter.GetWeight(self,ROOT.PlotUtils.detail.empty())
        return weight

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [SusaValenciaUniverse(chain,1)]

class LeakageUniverse(CVUniverse,object):
    def __init__(self,chain,nsigma):
        super(LeakageUniverse,self).__init__(chain,nsigma)

    def GetLeakageCorrection(self):
        #return super(LeakageUniverse,self).GetLeakageCorrection() * (1+self.nsigma*SystematicsConfig.LEAKAGE_SYSTEMATICS)

        return super(LeakageUniverse,self).GetLeakageCorrection() + ( self.nsigma*10 if abs(self.mc_primaryLepton) ==11 else 0)

    def ShortName(self):
        return "Leakage_Uncertainty"

    def LatexName(self):
        return "Leakage Unvertainty"

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [LeakageUniverse(chain,i) for i in OneSigmaShift]


def GetAllSystematicsUniverses(chain,is_data,is_pc =False,exclude=None,playlist=None):
    #use list because I want to control the order, process cv first, then vertical shifts, then lateral shifts
    universes = []
    CVUniverse.is_pc = is_pc
    CVUniverse.LeafGetters= {}
    if is_data:
        #data has only cv universe
        universes.append(CVUniverse(chain, None))
    else:
        #append cv universe
        universes.append(CVUniverse(chain,0))
        
        #Set Playlist, only MC cares about playlist
        if playlist is None:
            #attempt guess playlist from the chain
            playlist = Utilities.PlaylistLookup(chain.GetValue("mc_run",0))

        if not is_pc:
            CVUniverse.SetPlaylist(playlist)
            #Set NuE constraint
            CVUniverse.SetNuEConstraint(SystematicsConfig.USE_NUE_CONSTRAINT)
        else:
            CVUniverse.SetPlaylist("dummy")
            CVUniverse.SetNuEConstraint(False)

        #Set nu_type
        CVUniverse.SetAnalysisNuPDG(SystematicsConfig.AnaNuPDG)

        #Set NonResPi weight
        CVUniverse.SetNonResPiReweight(SystematicsConfig.NonResPi)
        CVUniverse.SetDeuteriumGeniePiTune(False)
        CVUniverse.SetZExpansionFaReweight(SystematicsConfig.NumZExpansionUniverses)
        CVUniverse.SetNFluxUniverses(SystematicsConfig.NUM_FLUX_UNIVERSE)

        if chain.GetTree().GetName() == "Truth":
            CVUniverse.SetTruth(True)

        if exclude is None or "all" not in exclude:
            # Vertical shift first to skip some cut calculation

            #Electron momentum universe
            if abs(SystematicsConfig.AnaNuPDG)==12:
                universes.extend(ElectronEnergyShiftUniverse.GetSystematicsUniverses(chain ))
            elif abs(SystematicsConfig.AnaNuPDG)==14:
                universes.extend(MuonUniverseMinerva.GetSystematicsUniverses(chain ))
                universes.extend(MuonUniverseMinos.GetSystematicsUniverses(chain ))
            else:
                raise ValueError ("AnaNuPDG should be \pm 12 or 14, but you set {}".format(SystematicsConfig.AnaNuPDG))

            # #Electron angle universe
            universes.extend(ElectronAngleShiftUniverse.GetSystematicsUniverses(chain ))
            # #beam angle shift universe
            universes.extend(BeamAngleShiftUniverse.GetSystematicsUniverses(chain ))
            #particle response shift universe
            universes.extend(ResponseUniverse.GetSystematicsUniverses(chain ))

            #Flux universe
            universes.extend(FluxUniverse.GetSystematicsUniverses(chain ))

            #Genie universe
            universes.extend(GenieUniverse.GetSystematicsUniverses(chain ))
            universes.extend(GenieRvx1piUniverse.GetSystematicsUniverses(chain ))
            universes.extend(GenieFaCCQEUniverse.GetSystematicsUniverses(chain ))

            # #2p2h universes
            universes.extend(Universe2p2h.GetSystematicsUniverses(chain ))

            # #RPA universe:
            universes.extend(RPAUniverse.GetSystematicsUniverses(chain ))


            # # #LowQ2PionUniverse
            universes.extend(LowQ2PionUniverse.GetSystematicsUniverses(chain ))

            # #birk shift universe
            # #universes.extend(BirksShiftUniverse.GetSystematicsUniverses(chain ))

            # MKModelUniverse
            universes.extend(MKModelUniverse.GetSystematicsUniverses(chain ))

            #FSIWeighUniverse
            universes.extend(FSIWeightUniverse.GetSystematicsUniverses(chain ))

            #SuSAValenciaUniverse
            universes.extend(SusaValenciaUniverse.GetSystematicsUniverses(chain ))

            #hadron reweight shifting universe
            universes.extend(GeantHadronUniverse.GetSystematicsUniverses(chain ))

            #leakage universe
            universes.extend(LeakageUniverse.GetSystematicsUniverses(chain ))


    # Group universes in dict.
    univ_dict = OrderedDict()
    for univ in universes:
        univ_dict.setdefault(univ.ShortName(), []).append(univ)

    if exclude is not None and len(univ_dict)>1:
        for i in exclude:
            if i in univ_dict:
                del univ_dict[i]

    return univ_dict

#code for testing new cv universe.
if __name__ == "__main__":

    #path = "/minerva/data/users/hsu/merged_files/merged-1558038217.root"
    path = "/pnfs/minerva/persistent/users/jyhan/NuECCQE/v21r1p1_recoilE_final/mcme1A/grid/central_value/minerva/ana/v21r1p1/00/11/00/17/SIM_minerva_00110017_Subruns_0325_NuECCQE_Ana_Tuple_v21r1p1.root"
    
    chain = ROOT.TChain("NuECCQE")
    chainWrapper = PlotUtils.ChainWrapper("NuECCQE")
    chain.Add(path)
    chainWrapper.Add(path)
    chainWrapper.GetValue("prong_part_score",0)
    univs = GetAllSystematicsUniverses(chainWrapper,False)
    for i in range(0,1):
        univ = univs["Genie"][0]
        univ.SetEntry(i)
        chain.GetEntry(i)
        print(univ.GetWeight())
        #print univ.GetVecElem("prong_part_E",0,3)
        e =  univ.GetVecOfVecDouble("prong_part_E")
        print(e[0].size())
        print(e[0][0],e[0][1],e[0][2])
        print(chainWrapper.GetValue("prong_part_E",i,0))
        print(chain.prong_part_E[0][0],chain.prong_part_E[0][1],chain.prong_part_E[0][2])
        univ2 = univs["Flux"][0]
        univ2.SetEntry(i)
        print(univ.UseNuEConstraint())
        print(univ2.UseNuEConstraint())
        print(univ2.GetWeight())
        print(univ.mc_FSPartPDG)
        #for _ in univs["eltheta"]:
        #    print _.shift_angle


    #hist = PlotUtils.HistWrapper(ROOT.PythonMinervaUniverse)("name","title",5,0,10,map)
