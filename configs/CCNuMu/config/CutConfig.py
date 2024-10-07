"""
  Analysis cut values and the functions that represent them.
  Tune them using studies/CutTuning.py
  
   author: J. Wolcott <jwolcott@fnal.gov>
   date: December 2013

"""

############################################################################

# hack the fiducial vertex cut,
HACK_R2 = False

# tuned cut values
DS_CAL_VISE_CUT = 0.2
OD_CAL_VISE_CUT = 0.05
PSI_FLAT_CUT = 0.1
FRONT_DEDX_CUT = 2.4  # in MeV/cm
PID_SCORE_CUT = 0.7
MIN_VERTEX_TRACK_MULTIPLICITY = 2
MAX_VERTEX_TRACK_MULTIPLICITY = 6 

NONMIP_CLUS_FRAC_CUT = 0.4
TRANSVERSE_GAP_SCORE_CUT = 15
FIRST_FIRE_FRACTION_CUT = 0.25
UPSTREAM_OD_ENERGY_CUT = 5000
EXUV_CUT = 0.2
EUV_CUT = 0.3

DEDX_PLANES = [5,9]
HELICITY= 1

FIDUCIAL_APOTHEM = 850
FIDUCIAL_Z_RANGE = [5980,8422]

# Kinematics cutoffs
ELECTRON_ENERGY_RANGE = [2.5, float('inf')] # in GeV
NEUTRINO_ENERGY_RANGE = [0, 100] # in GeV.
ELECTRON_ANGLE_RANGE = [0, 20] # in deg
RECO_Q3_RANGE = [0,4]
TRUE_Q3_RANGE = [0,4]

RECO_PT_RANGE= [.2,1.6]
RECO_VISE_RANGE= [0,1.2]
RECO_ZDIFF_RANGE= [-30,50]
RECO_VTX_ZDIFF_RANGE= [-20,20]
RECO_PROTON_END = 8500
WEXP_CUT = 1.4
visE_RANGE = [0.0,0.7]
visE_FHCRANGE = [0.0,0.7]
RECO_ANGLE = [20,90]
FRONT_DEDX_PI0_UPPERBOUND = 5
PsiEe_CUT = 0.5
Ethetasquared_CUT = .003

############################################################################
# choose the cuts you want from cut library

SAMPLE_CUTS = {
    "Signal" : [
        "InverseHasNoBackExitingTracks",
        "Vertex_Z",
        "Vertex_Apothem",
        "Eavail",
        "Pt",
        "Etheta",
        #"ZDifference",
        #"ZDifference_vtx",
        #"ProtonEnd",
        #"ContainedProton",
     ]
}

KINEMATICS_CUTS = [
    "LeptonAngle",
]

#######################################
