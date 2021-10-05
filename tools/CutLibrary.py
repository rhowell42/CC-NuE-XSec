"""
    Analysis Cut Library that provide functions/class to be used in EventClassification.
    Cut value is configured in config/CutConfig.py

    Adapted from Jeremy's CutConfig.

    author: H. Su
    date: Jul 2019
"""


import array
import functools
import math
import ROOT

from config import CutConfig


class SelectionCut(object):
    def __init__(self, **kwargs):
        self.name = kwargs["name"]
        try:
            self._value_getter = kwargs["value_getter"]
        except KeyError:
            self._value_getter = functools.partial(lambda name, event, nprong: getattr(event, name), self.name)

        self._cut_fn = kwargs["cut_fn"]

        try:
            self._variable_range = kwargs["variable_range"]
        except KeyError:
            self._variable_range = None

    def Values(self, event, nprong=0):
        return self._value_getter(event, nprong)

    def DoesEventPass(self, event, nprong=0):
        return self.DoValuesPass( self.Values(event, nprong) )

    def DoValuesPass(self, vals):
        return self._cut_fn(vals)

# these are used often.  just consolidate them here.
REQUIRE_POSITIVE_INT = lambda val: val >= 1
REQUIRE_UNITY_INT = lambda val: int(val) == 1  #val >= 1 and val < 2

VARIABLE_RANGE_01 = [-0.1,0.9,1.9] # for true/false varialbe, avoiding ambiguity


"""
  How to define a new cut:
1) add an entry to CUT_CONFIGS, format is :
  "name of cut" : {
   "value_getter" : lambada function for getting the cut variable,
   "cut_fn" : lambda function for determing if the cut varialbe passes.
  }
2) add the cut to the list of sample cuts in CutConfig.py 
"""

def CalcApothem(x,y):
    x=abs(x)
    y=abs(y)
    if ( x == 0 or y/x > 1/math.sqrt(3)):
        return (y+x/math.sqrt(3))/2*math.sqrt(3)
    else:
        return x

CUT_CONFIGS = {
    "NoCut": {
        "value_getter": lambda event, nprong: None,
        "cut_fn": lambda val: True,
        "variable_range": VARIABLE_RANGE_01,
    },
    
    "HasFiducialVertex": {
        "cut_fn": REQUIRE_POSITIVE_INT,
        "variable_range": VARIABLE_RANGE_01,
    },
    "Low UIE": {
        "value_getter": lambda event,nprong: event.UpstreamInlineEnergy,
        "cut_fn":lambda val : val<10,
    },
    "High UIE": {
        "value_getter": lambda event,nprong: event.UpstreamInlineEnergy,
        "cut_fn":lambda val : val>10,
    },
    "Vertex_Z": {
        "value_getter": lambda event, nprong: event.vtx[2],
        "cut_fn": lambda val: CutConfig.FIDUCIAL_Z_RANGE[0] <= val <= CutConfig.FIDUCIAL_Z_RANGE[1],
        "variable_range": [100 * i for i in range(100)]
    },
    "Vertex_Z_Centralized" : {
        "value_getter": lambda event,nprong: event.classifier.ZRange.passesCut(event,ROOT.PlotUtils.detail.empty()),
        "cut_fn" : lambda val:val
    },
    "Vertex_Apothem_Centralized" : {
        "value_getter": lambda event,nprong: event.classifier.Apothem.passesCut(event,ROOT.PlotUtils.detail.empty()),
        "cut_fn" : lambda val:val
    },
    "Vertex_Apothem": {
        "value_getter":lambda event,nprong: CalcApothem(event.vtx[0],event.vtx[1]),
        "cut_fn": lambda val: val <= CutConfig.FIDUCIAL_APOTHEM,
        "variable_range": [10*i for i in range(120)]
    },

    "Truth_Vertex_Z": {
        "value_getter": lambda event, nprong: event.mc_vtx[2],
        "cut_fn": lambda val: CutConfig.FIDUCIAL_Z_RANGE[0] <= val <= CutConfig.FIDUCIAL_Z_RANGE[1],
        "variable_range": [100 * i for i in range(100)]
    },
    "Truth_Vertex_Apothem": {
        "value_getter":lambda event,nprong: CalcApothem(event.mc_vtx[0],event.mc_vtx[1]),
        "cut_fn": lambda val: val <= CutConfig.FIDUCIAL_APOTHEM,
        "variable_range": [10*i for i in range(120)]
    },

    "HasTracks": {
        "value_getter": lambda event, nprong: event.n_prongs>0,
        "cut_fn": REQUIRE_POSITIVE_INT,
        "variable_range": VARIABLE_RANGE_01,
    },
    
    "HasNoBackExitingTracks": {
        "cut_fn": REQUIRE_POSITIVE_INT,
        "variable_range": VARIABLE_RANGE_01,
    },

    "EMLikeTrackScore": {
        "value_getter": lambda event,n_prongs: event.prong_part_score,
        "cut_fn": lambda val: [_>=CutConfig.PID_SCORE_CUT for _ in val].count(True)>0,
        "variable_range": [0.1*i for i in range(5,11)]
    },

    "UniqueEMLikeTrackScore": {
        "value_getter": lambda event,n_prongs: event.prong_part_score,
        "cut_fn": lambda val: [_>=CutConfig.PID_SCORE_CUT for _ in val].count(True)==1,
        "variable_range": [0.1*i for i in range(5,11)]
    },
    
    "DSCalVisE" : {
        "value_getter": lambda event, nprong:
            # return None if there's HCAL energy but not ECAL energy...
            event.prong_HCALVisE[nprong]/event.prong_ECALVisE[nprong] if event.prong_ECALVisE[nprong] > 0.1 else (1e10 if event.prong_HCALVisE[nprong] > 0 else 0),
        "cut_fn": lambda val: val <= CutConfig.DS_CAL_VISE_CUT,
        "variable_range": [0.1* i for i in range(0,11)]
    },
    
    "ODCalVisE": {
        "value_getter": lambda event, nprong:
            # return None if there's OD energy but not side ECAL energy...
            event.prong_ODVisE[nprong]/event.prong_SideECALVisE[nprong] if event.prong_SideECALVisE[nprong] > 0.1 else (1e10 if event.prong_ODVisE[nprong] > 0 else 0),
        "cut_fn": lambda val: val <= CutConfig.OD_CAL_VISE_CUT,
        "variable_range": [0.1* i for i in range(0,11)]
    },

    "StartPointVertexMultiplicity": {
        "cut_fn": REQUIRE_UNITY_INT,
        "variable_range": list(range(0,6))
    },


    "HasNoVertexMismatch": {
        "cut_fn": REQUIRE_UNITY_INT,
        "variable_range": VARIABLE_RANGE_01,
    },


    "VertexTrackMultiplicity": {
        "cut_fn": lambda val: CutConfig.MIN_VERTEX_TRACK_MULTIPLICITY <=  val <= CutConfig.MAX_VERTEX_TRACK_MULTIPLICITY,
        "variable_range": list(range(0,6))
    },

    "HasNoNonEMExitingTracks": {
        "cut_fn": REQUIRE_POSITIVE_INT,
        "variable_range": VARIABLE_RANGE_01,
    },

    "Psi": {
        "value_getter": lambda event, nprong: (event.prong_TotalVisE[nprong]/1e3, event.Psi) ,
        "cut_fn": lambda vals: vals[1] < CutConfig.PSI_FLAT_CUT,
        "variable_range": [0.1* i for i in range(0,11)]
    },

    "Wexp": {
        "value_getter": lambda event, nprong: event.kin_cal.reco_W,
        "cut_fn": lambda vals: vals>0 and vals <= CutConfig.WEXP_CUT,
        "variable_range": [0.5* i for i in range(0,11)]
    },
    "Eavail": {
        "value_getter": lambda event, nprong: event.kin_cal.reco_visE,
        "cut_fn": lambda vals: vals <= CutConfig.Reco_visEcut,
        "variable_range": [0.4* i for i in range(0,11)]
    },

    # "InversePsi": {
    #     "value_getter": lambda event, nprong: (event.prong_TotalVisE[nprong]/1e3, event.Psi) ,
    #     "cut_fn": lambda vals: vals[1] > CutConfig.PSI_FLAT_CUT,
    #     "variable_range": [0.1* i for i in range(0,11)]
    # },

    "MeanFrontdEdX": {
        "value_getter": lambda event, nprong: event.prong_dEdXMeanFrontTracker[nprong],
        #"cut_fn": lambda val: val > 0,
        "cut_fn": lambda val: 0 < val <= CutConfig.FRONT_DEDX_CUT,
#       "value_getter": lambda event: (event.prong_dEdXMeanFrontPositionTracker[0] if hasattr(event, "prong_dEdXMeanFrontPositionTracker") else 0, event.prong_dEdXMeanFrontTracker[0]),
#       "value_getter": _dedx_cut_getter,
#       "cut_fn": lambda vals: vals[1] > 0 and vals[1] <= FRONTDEDX_CUT_GRAPH.Eval(vals[0]),
        "variable_range": [0.1* i for i in range(0,51)]
    },
    
    #   "UpstreamODEnergy": {
    #       "value_getter": lambda event: event.UpstreamODEnergy,
    #       "cut_fn": lambda val: val < UPSTREAM_OD_ENERGY_CUT,
    #   },

    "MidMeanFrontdEdX" : {
        "value_getter": lambda event, nprong: event.prong_dEdXMeanFrontTracker[nprong],
        #"cut_fn": lambda val: val > 0,
        "cut_fn": lambda val: val > CutConfig.FRONT_DEDX_CUT and val<CutConfig.FRONT_DEDX_PI0_UPPERBOUND,
        "variable_range": [0.1* i for i in range(0,51)]
    },

    "HighMeanFrontdEdX" : {
        "value_getter": lambda event, nprong: event.prong_dEdXMeanFrontTracker[nprong],
        #"cut_fn": lambda val: val > 0,
        "cut_fn": lambda val: val > CutConfig.FRONT_DEDX_PI0_UPPERBOUND,
        "variable_range": [0.1* i for i in range(0,51)]
    },

    "NonMIPClusFrac": {
        "value_getter": lambda event, nprong: event.prong_NonMIPClusFrac[nprong],
        "cut_fn": lambda val: val > CutConfig.NONMIP_CLUS_FRAC_CUT,
        "variable_range": [0.1* i for i in range(0,11)]
    },

    "TransverseGapScore": {
        "value_getter": lambda event, nprong: event.prong_TransverseGapScore[nprong],
        "cut_fn": lambda val: val > CutConfig.TRANSVERSE_GAP_SCORE_CUT,
        "variable_range": [1.5* i for i in range(0,21)]
    },
    
    #"HasNoMichelElectrons": {
    #"value_getter": lambda event, nprong: (event.michel_digits, event.michel_energy, event.michel_slice_energy),
    #   "cut_fn": lambda vals: all([not(vals[0][i] < 35 and vals[1][i] < 55 and vals[2][i] < 100 and vals[1][i]/vals[0][i] > 0.8) for i in range(len(vals[0]))])
    #},#Better way for this?
    
    "DeadTime": {
        "value_getter": lambda event, nprong: event.phys_n_dead_discr_pair_upstream_prim_track_proj,
        "cut_fn": lambda val: val >= 0 and val <= 1,
        "variable_range": [0.1* i for i in range(0,11)]
    },
    
    "Afterpulsing": {
        "value_getter": lambda event, nprong: event.prong_FirstFireFraction[nprong],
        "cut_fn": lambda val: val >= CutConfig.FIRST_FIRE_FRACTION_CUT,
        "variable_range": [0.1* i for i in range(0,11)]
    },

    "Exuv" : {
        "value_getter": lambda event, nprong: event.kin_cal.Exuv,
        "cut_fn": lambda val:abs(val) <= CutConfig.EXUV_CUT,
        "variable_range": [0.1* i -0.5 for i in range(0,11)]
    },

    "Euv" : {
        "value_getter": lambda event, nprong: event.kin_cal.Euv,
        "cut_fn": lambda val: abs(val) <= CutConfig.EUV_CUT,
        "variable_range": [0.1* i -0.5 for i in range(0,11)]
    },

    "LLR" : {
        "value_getter": lambda event, nprong: max(-20,min(19.9, event.kin_cal.LLR)) if event.kin_cal.LLR is not None else None,
        "cut_fn": lambda val: val is not None and val <= 0,
        "variable_range": [i for i in range(-20,20)]
    },

    "Neutrino Helicity" : {
        "value_getter": lambda event, nprong: -1 if event.MasterAnaDev_nuHelicity ==1 else 1,
        "cut_fn": lambda val: CutConfig.HELICITY is None or val * CutConfig.HELICITY>0,
        "variable_range": [-1.1,0,1.1]
    },
    "HasMINOSMatch" : {
        "value_getter": lambda event,nprong: event.IsMinosMatchMuon(),
        "cut_fn": lambda val: val,
        "variable_range": [-1.1,0,1.1]
    }

}  # CUT_CONFIGS


KINEMATICS_CUT_CONFIGS = {
    "RecoLeptonEnergy": {
        "value_getter": lambda event,nprong: event.kin_cal.reco_E_lep,
        "cut_fn": lambda val: CutConfig.ELECTRON_ENERGY_RANGE[0] <= val < CutConfig.ELECTRON_ENERGY_RANGE[1],
        "variable_range": [0.5*i for i in range(0,21)]
    },
    
    "RecoLeptonAngle": {
        "value_getter": lambda event,nprong: event.kin_cal.reco_theta_lep,
        "cut_fn": lambda val: CutConfig.ELECTRON_ANGLE_RANGE[0] <= val < CutConfig.ELECTRON_ANGLE_RANGE[1],
    },

    "RecoNeutrinoEnergy": {
        "value_getter": lambda event,nprong: event.kin_cal.reco_E_nu_cal,
        "cut_fn": lambda val: CutConfig.NEUTRINO_ENERGY_RANGE[0] <= val < CutConfig.NEUTRINO_ENERGY_RANGE[1],
    },
    "RecoQ3": {
        "value_getter": lambda event,nprong: event.kin_cal.reco_q3,
        "cut_fn": lambda val: CutConfig.RECO_Q3_RANGE[0] <= val < CutConfig.RECO_Q3_RANGE[1],
        "variable_range": [0.1*i for i in range(0,21)]
    },

    "TrueLeptonEnergy": {
        "value_getter": lambda event,nprong: event.kin_cal.true_E_lep,
        "cut_fn": lambda val: CutConfig.ELECTRON_ENERGY_RANGE[0] <= val < CutConfig.ELECTRON_ENERGY_RANGE[1],
    },
    
    "TrueLeptonAngle": {
        "value_getter": lambda event,nprong: event.kin_cal.true_theta_lep,
        "cut_fn": lambda val: CutConfig.ELECTRON_ANGLE_RANGE[0] <= val < CutConfig.ELECTRON_ANGLE_RANGE[1],
    },

    "TrueNeutrinoEnergy": {
        "value_getter": lambda event,nprong: event.kin_cal.true_enu_genie,
        "cut_fn": lambda val: CutConfig.NEUTRINO_ENERGY_RANGE[0] <= val < CutConfig.NEUTRINO_ENERGY_RANGE[1],
    },
    "TrueQ3": {
        "value_getter": lambda event,nprong: event.kin_cal.true_q3,
        "cut_fn": lambda val: CutConfig.TRUE_Q3_RANGE[0] <= val < CutConfig.TRUE_Q3_RANGE[1],
        "variable_range": [0.1*i for i in range(0,21)]
    },
}

if CutConfig.HACK_R2:
    print("\n\n=====================================================================")
    print("WARNING: using R < 850mm and z cuts instead of ntuple fiducial volume cut...")
    print("=====================================================================\n\n")
    
    CUT_CONFIGS["HasFiducialVertex"] = {
        "value_getter": lambda event, nprong: (math.sqrt(event.vtx[0]**2 + event.vtx[1]**2), event.vtx[2]),
        "cut_fn": lambda vals: vals[0] <= 850 and vals[1] > 5990 and vals[1] < 8340,
    }

# make SelectionCut objects now

class Cuts(object):
    def __init__(self):
        self.cutmap = {}
    def __getitem__(self,key):
        if key in self.cutmap:
            return self.cutmap[key]
        elif key in CUT_CONFIGS:
            config = CUT_CONFIGS[key]
            config["name"]=key
            self.cutmap[key]=SelectionCut(**config)
            return self.cutmap[key]
        elif key in KINEMATICS_CUT_CONFIGS:
            config = KINEMATICS_CUT_CONFIGS[key]
            config["name"]=key
            self.cutmap[key]=SelectionCut(**config)
            return self.cutmap[key]
        else:
            raise KeyError("Cut {} not defined")


CUTS = Cuts()
# for name, config in CUT_CONFIGS.items():
#     config["name"] = name
#     CUTS[name] = SelectionCut(**config)

# for name, config in KINEMATICS_CUT_CONFIGS.items():
#     config["name"] = name
#     CUTS[name] = SelectionCut(**config)

