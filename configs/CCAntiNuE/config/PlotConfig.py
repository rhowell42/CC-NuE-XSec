
"""
  PlotConfig.py:
   Centralization of common plotting configurations.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    May 2014
"""

import math
import ROOT
#from collections import OrderedDict

#from AnalysisConfig import AnalysisConfig

# constants needed in calculations.
M_e = 0.511 # MeV
M_p = 938.272 # MeV
M_n = 939.565 # MeV
QE_binding_E = 34 # MeV  (same as nu_mu PRL for FHC)

# these will be used repeatedly in calculations.
# better just do the arithmetic once here.
M_n_star = M_n - QE_binding_E
M_e_sqr = M_e**2
M_n_star_sqr = M_n_star**2
M_p_sqr = M_p**2

NQ3 = 8
NQ0 = 19

LOW_RECOIL_BIN_Q3 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
LOW_RECOIL_BIN_Q3_Truth = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]

#LOW_RECOIL_BIN_Q0_Truth = [0.0, 0.2, 0.4, 0.8, 1.0, 1.2]
#LOW_RECOIL_BIN_Q0 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2] 

## Hang Bins ##
#LOW_RECOIL_BIN_Q0_Truth = [0.,0.05]+[ 0.1*i for i in range(1,7)]+[0.2*i for i in range(4,7)]
#LOW_RECOIL_BIN_Q0 = [0,0.05,0.1,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.20,0.24,0.28,0.32,0.4,0.5]#,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
LOW_RECOIL_BIN_Q0 = [0,0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.3,0.4,0.5,0.7,1.0,1.5,2.0]
LOW_RECOIL_BIN_Q0_Truth = [0,0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.3,0.4,0.5,0.7,1.0,1.5,2.0]
#LOW_RECOIL_BIN_Q0 = [0,0.01,.02,.03,.04,.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.3,0.4,0.5]
#LOW_RECOIL_BIN_Q0_Truth = [0,0.01,.02,.03,.04,.05,0.1,0.12,0.14,0.16,0.18,0.20,0.25,0.3,0.4,0.5]

#LOW_RECOIL_BIN_Q0_Truth = [ 0.1*i for i in range(7)]+[0.2*i for i in range(4,7)]
#LOW_RECOIL_BIN_Q0 = [0.0,0.04,0.08,0.12,0.16,0.24,0.32,0.4,0.6,0.8,1.0,1.2]

#LOW_RECOIL_BIN_Q0 = [0.0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80, 1.00, 1.2]
#LOW_RECOIL_BIN_Q0_Truth = [0.0, 0.08, 0.16, 0.25, 0.35, 0.50, 0.8, 1.2]

#LOW_RECOIL_BIN_Q3 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2]
#LOW_RECOIL_BIN_Q3_Truth = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2] #, 1.6, 2]
#LOW_RECOIL_BIN_Q0 = [0.0, 0.04, 0.08,
#                             0.12, 0.16, 0.2,
#                             0.25, 0.30, 0.35, 0.40,
#                             0.50, 0.60, 0.80, 1.00, 1.2, 2.0]
#LOW_RECOIL_BIN_Q0_Truth = [0.0, 0.04, 0.08,
#                             0.12, 0.16, 0.2,
#                             0.25, 0.30, 0.35, 0.40,
#                             0.50, 0.60, 0.80, 1.00, 1.2, 2.0]

#LOW_RECOIL_BIN_Q3 = [0.0, 0.2, 0.4, 0.6, 0.9, 1.2]
#LOW_RECOIL_BIN_Q0 = [0.0, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.4, 0.6, 0.8, 1.0, 1.2]
#LOW_RECOIL_BIN_Q3_Truth = LOW_RECOIL_BIN_Q3
#LOW_RECOIL_BIN_Q0_Truth = LOW_RECOIL_BIN_Q0

#PT_BINNING = [0.2,0.4,0.6,0.8,1]
#PT_BINNING = [0.2 * i for i in range(9)]
PT_BINNING = [0,0.07,0.15,0.25,0.33,0.4,0.47,0.55,0.7,0.85,1.0,1.25,1.5,2.5]
PT_BINNING_AARON = [0,0.15,0.3,0.45,0.6,0.75,0.9,1.25]
PT_BINNING_Truth = [0.2 * i for i in range(9)]
PSI_EE_BINNING = [0.15*i for i in range(17)]
EE_PSI_BINNING = [2.5,  3.5, 4.5, 6, 8, 12, 16, 22]
BACKGROUND_FIT_Q3_BIN = [0.0, 0.6, 0.8, 1.0, 1.2, 1.6, 2]
#LOW_RECOIL_BIN_Q0 = [0.0, 0.04, 0.08,
#		             0.12, 0.16, 0.2,
#		             0.25, 0.30, 0.35, 0.40,
#		             0.50, 0.60, 0.80, 1.00, 1.2, 2.0]
#LOW_RECOIL_BIN_Q0_Truth = [0.0, 0.04, 0.08,
#                             0.12, 0.16, 0.2,
#                             0.25, 0.30, 0.35, 0.40,
#                             0.50, 0.60, 0.80, 1.00, 1.2, 2.0]

RESOLUTION_BINNING = [-1.0+0.1*i for i in range(41)]
MC_EM_ENERGY_SCALE = 1.05  # this from Trung: Doc 9370, slides 2-5, 15; Doc 9911; Doc 10102, slides 18-20
MC_MICHEL_ELECTRON_ENERGY_SCALE = 1.03  # from the NIM
MC_DEDX_ENERGY_SCALE = 1.015  # Based on eyeballing the dE/dx distribution for this analysis in the extra energy sideband.


W_RECOIL_BIN_Q3_Truth = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2] 
# WARNING -- change this at your peril.  it affects the accepted region for true events
# via the 'truth' hacks in SelectedSamplePlots,
# which consequently affects the range that the efficiency correction corrects into!
"""
if AnalysisConfig.signal_defn == AnalysisConfig.SIGNAL_DEFNS.CCQE_LIKE:
	NEUTRINO_ENERGY_RANGE = [0, 10] # in GeV.  
elif AnalysisConfig.signal_defn == AnalysisConfig.SIGNAL_DEFNS.GENIE_CCQE:
	NEUTRINO_ENERGY_RANGE = [1.5, 10]
elif AnalysisConfig.signal_defn == AnalysisConfig.SIGNAL_DEFNS.EXCESS:
	NEUTRINO_ENERGY_RANGE = [0, float("inf")]  # this is presumptively an NC process, so anything's fair game
"""

# used in making plots of various particle types
PID_LABELS = {
	0: "other",
	11: "e^{#pm}",
	13: "#mu^{#pm}",
	22: "#gamma",
	111: "#pi^{0}",
	211: "#pi^{#pm}",
	2212: "p^{+}",
}

# for plotting atomic number.
# want bins for H, C, O, Si, Cl, Ti, Fe, Pb 
A_LABELS = {
	0: "other",
	1: "H",
	12: "C",
#	14: "N",
	16: "O",
	27: "Al",
	28: "Si",
	35: "Cl",
	48: "Ti",
	56: "Fe",
	207: "Pb", 
}

#ELECTRON_ANGLE_BINNING =  range(10) + [10, 12, 15, 20, 27, 35]
ELECTRON_ANGLE_BINNING =  [1 * i for i in range(41)]
ELECTRON_ANGLE_RESIDUAL_BINNING =  [-1.+0.1* i for i in range(0,21)]
PROTON_ANGLE_RESIDUAL_BINNING = [-20 +i for i in range(41)]
#ELECTRON_ANGLE_2D_BINNING = [-x for x in reversed(ELECTRON_ANGLE_BINNING)] + ELECTRON_ANGLE_BINNING[1:]
#EXCESS_ANGLE_BINNING = range(0, 15, 3) + [15, 20, 27, 35]
EXCESS_ANGLE_BINNING = list(range(0, 15, 3)) + [20]
LOW_RECOIL_BIN_HIGH_Q0 = [0.2, 0.25, 0.30, 0.35, 0.40,
                     0.50, 0.60, 0.80, 1.00, 1.5]
LOW_RECOIL_BIN_LOW_Q0 = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1,
                     0.12, 0.14, 0.16, 0.2]
PROTON_ANGLE_BINNING = [2*i for i in range(51)]

NEUTRINO4_EE_BINNING = [0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9,10,12.5,15,17.5,20]
NEUTRINO4_P_BINNING = [i*.025 for i in range(40)]
NEUTRINO4_LE_BINNING = [i*.015 for i in range(34)]
NEUTRINO4_EE_THETA_BINNING = [.003 * i for i in range(len(NEUTRINO4_EE_BINNING))]
NEUTRINO4_LENGTH_BINNING = [.309, .339, .369, .399, .429, .459, .489, .519, .549, .579, .609, .639, .670, .700, .730, .760, .790, .820, .850, .880, .910, .940, .970, 1.000]

NEUTRINO4_EE_BINNING_INV = [0.0,1/20,1/15,1/12.5,1/10,1/8,1/7.5,1/7,1/6.5,1/6,1/5.5,1/5,1/4.5,1/4,1/3.75,1/3.5,1/3.25,1/3,1/2.75,1/2.5,1/2,1/1.75,1/1.5]
NEUTRINO4_L_OVER_E_BINNING = [0.0, 2., 5., 10., 17., 23., 29., 37., 46., 56., 69., 85., 104., 120., 139., 160., 185., 215., 250., 293., 347., 417., 510.,600]
NEUTRINO4_L_OVER_E_BINNING = [i/1000 for i in NEUTRINO4_L_OVER_E_BINNING]
ELECTRON_ENERGY_BINNING = [0.0,2.5,5,7.5,10,12.5,15,20]
#ELECTRON_ENERGY_BINNING = [0.5 * i for i in range(31)]+[20]
SUM_VISIBLE_ENERGY_BINNING = [1 * i for i in range(1,11)]
VISIBLE_ENERGY_RESIDUAL_BINNING = [-1+0.04* i for i in range(0,51)]
ELECTRON_ENERGY_RESIDUAL_BINNING = [-1+0.05* i for i in range(0,41)]
#EXCESS_ENERGY_BINNING = [0, 3, 6, 9, 12, 15, 20]  # ELECTRON_ENERGY_BINNING + [12, 15, 20]
#ELECTRON_ENERGY_BINNING = [0.75, 2, 3, 5, 7, 9, 20]  # Jaewon's bins
#NEUTRINO_ENERGY_BINNING = [i for i in range(6)] + [7, 10] # , 13, 18, 25]
NEUTRINO_ENERGY_BINNING = [10,11,12,13,14,15,17,19,21,23,25,30,35,40,45,50]
#NEUTRINO_ENERGY_BINNING_BIGGER = [i for i in range(6)] + [7, 10, 13, 18, 25]
#OD_ENERGY_BINNING = [0, 0.5, 1, 2, 3, 4, 5] 
#VISIBLE_ENERGY_BINNING = [ 0.1*i for i in range(20) ] + [0.2*i for i in range(10, 20)] + [0.5*i for i in range(8, 14)] + range(7, 10)
VISIBLE_ENERGY_BINNING = [ 0.02*i for i in range(9) ] + [0.05*i for i in range(4,9)]+ [0.1*i for i in range(5, 7)] + [0.2*i for i in range(4, 8)] 
APOTHEM_BINNING = [50*i for i in range(17)]
VERTEX_Y_BINNING = [100*i -1000 for i in range(21)]
#ETH2_BINNING = [0, 0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05]
ETH2_BINNING = [0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.5]

QSQUARED_BINNING_CCQE_LIKE = [0, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 1.0, 1.25, 2.0]
QSQUARED_BINNING_GENIE_CCQE = [0, 0.1, 0.2, 0.4, 0.8, 1.2, 2.0]
QSQUARED_BINNING_INC = [0.1 *i for i in range(101) ]
EXCESS_ENERGY_BINNING = [0, 3, 6, 9, 12, 15, 20]
NUECONE_BINNING = [10* i for i in range(0,51)]
NUE_SCATTERING_BINNING = [0.8,2,3,5,7,9,20]


#QSQUARED_BINNING = QSQUARED_BINNING_CCQE_LIKE if AnalysisConfig.signal_defn != AnalysisConfig.SIGNAL_DEFNS.GENIE_CCQE else QSQUARED_BINNING_GENIE_CCQE
QSQUARED_BINNING = QSQUARED_BINNING_INC

QSQUARED_SLICE_EDGES = [0.2, 0.75, 2,]

W_BINNING = [0.1*i for i in range(21)]

PSI_BINNING = [0.05*i for i in range(41)]
PSI_TAIL_BINNING = [0.05*i for i in range(4,15)]
PSI_FRONT_BINNING = [0.01*i for i in range(0,11)]
BIN_AVAIL_OPTIMIZATION = [0.0, 0.080,0.260,0.540,1.2,2.0]
TRANS_BINNING = [0.1*i for i in range(0,21)]
LONG_BINNING = [0.1*i for i in range(-20,51)]
LONG_2D_BINNING = [0,0.25,0.5,0.75,1.0]
LONG_NEG_BINNING = [0.1*i for i in range(-20,1)]
LONG_POS_BINNING = [0.1*i for i in range(0,51)]

UPSTREAM_INLINE_ENERGY_BINS = [0.4, 1, 2, 4, 6.5, 10, 20, 40, 65, 100,150,250]  # designed to look ok on log10 axis

#VERTEX_ENERGY_BINS = [0, 5, 10, 15, 20, 30, 40,] + list(range(50, 350,30))
VERTEX_ENERGY_BINS = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,  0.09, 0.1]

WHOLE_DET_Z_BINS = list(range(4000, 10000, 100))

VERTEX_Z_RESOLUTION_BINNING = [0.1*i for i in range(41)]


BIRK_MIGRATION_PID_CV_BINNING=[0.49 + 0.05*i for i in range(0,11)]
BIRK_MIGRATION_PID_SHIFTED_BINNING=[0.49 + 0.05*i for i in range(0,11)]

BIRK_MIGRATION_DEDX_BINNING=[2+0.05* i for i in range(0,21)]

FRACTION_BINNING = [0.02 * i for i in range(51) ]

LLR_BINNING = [i for i in range(-20,21)]

ENERGY_BALANCE_BINNING = [0.02* i for i in range(26)]

LOG_ENERGY_BINNING = [0.1*i for i in range(41)]
FRONTDEDX_BINNING= [0.4 * i for i in range(26)]
FRONTDEDX_TAIL_BINNING = [0.4 * i for i in range(6,26)]
FRONTDEDX_FRONT_BINNING = [0.4 * i for i in range(7)]
#FRONTDEDX_BINNING = [0.2 * i for i in range(51)]
FRONTDEDX_POSITION_BINNING = [5 * i for i in range(-3, 83)]
DEDX_BINNING = [0.2 * i for i in range(52)]
FRONTDEDX_MEDIAN_BINNING = [0.4 * i for i in range(26)]
FRONTDEDX_2D_BINNING = [0,2.401,10]

PION_ENERGY_BINNING = [0.1 * i for i in range(51)]
PROTON_ENERGY_BINNING = [0.1 * i for i in range(21)]
PROTON_KE_BINNING = [0.05 * i for i in range(21)]
#VERTEX_DIFF_BINNING = range(-20,40)
VERTEX_DIFF_BINNING = [20*i for i in range(-3,31)]
VERTEX_Z_BINNING = [5000+ 200*i for i in range(21)]
PI0_ENERGY_BINNING = [0.5 * i for i in range(21)]
#EXTRA_ENERGY_BINNING = [0.05*i for i in range(41)]
EXTRA_ENERGY_BINNING = [0.05*i for i in range(15)]
OUTSIDE_ENERGY_BINNING = [0.05*i for i in range(21)]

PARTICLE_PDG_BINNING = list(range(-11,10))
GENIE_EVENT_BINNING = [0.5 * i for i in range(0,17)]
FRACTION_BINNING = [0.1 * i for i in range(-0,12)]
CLUSTER_PDG_FRAC = [0.05 * i for i in range(0,23)]
CLUSTER_PDG = [5* i for i in range(-45, 45)]
CLUSTER_PDG_HIGH = [5* i for i in range(400, 460)]
#CLUSTER_PDG_FULL = [5*i for i in range(-6,460)]
CLUSTER_PDG_FULL = [5*i for i in range(-50,460)]
PDG_TEMP_BINNING = [0.2*i for i in range(50,81)]

POSITION_BINNING_LOW = [0, 25.01, 50.01, 75.01, 100.01, 125.01, 150.01, 175.01, 200.01, 225.01]
POSITION_BINNING_HIGH = [225.01, 275.01, 300.01, 325.01, 350.01, 375.01, 400.01]
POSITION_BINNING_AVG = [0, 100.01, 200.01, 300.01, 400.01]
DIFFERENCE_BINNING = [0.1 * i for i in range(-20, 50)]
SHOWER_WIDTH_BINNING = [0.1 * i for i in range(0, 21)]

#Jeremy comparison plots
ZCOORDINATE_BINNING = [10*i for i in range(50, 90)]
INLINE_BINNING = [1*i for i in range(0,26)]
INLINE_LOW_BINNING = [0.01*i for i in range(0,11)]
INLINE_HIGH_BINNING = [0.1*i for i in range(1,51)]
INLINE_W_BINNING = [0.1*i for i in range(0,51)]
VERTEXR_BINNING = [0.5*i for i in range(0,21)]
RADIAL_BINNING = [0.5*i for i in range(-10,11)]
#TRANSVERSE_SHOWER_BINNING = [0.05*i for i in range(0,41)]
TRANSVERSE_SHOWER_BINNING = [0.1*i for i in range(0,21)]
#EXCESS_ENERGY_FIT = [1.5, 3, 4.5, 6, 9, 12, 15, 20]
EETHETA_BINNING = [0.2*i for i in range(0,7)]
EXCESS_ENERGY_FIT = [1.5, 2, 3.5, 5, 7, 10]
LEAKAGE_BINNING  = [0.01 * i for i in range(-20,20)]

RECO_W_BINNING = [0.1 * i for i in range(51)]

def FindQ2Bin(q2_val):
	q2_bin = None
	if q2_val >= QSQUARED_BINNING[-1]:
		return len(QSQUARED_BINNING)-1
	
	if q2_val < QSQUARED_BINNING[0]:
		print("underflow")
		return 1
	
	for bin_num in range(1, len(QSQUARED_BINNING)):
		if q2_val >= QSQUARED_BINNING[bin_num-1] and q2_val < QSQUARED_BINNING[bin_num]:
			q2_bin = bin_num
			break
	assert q2_bin is not None
	return q2_bin

def Printvar(event):
    print(event.mc_FSPartPDG)
    print(event.mc_FSPartE)
    print((event.kin_cal.true_visE, event.kin_cal.reco_visE))

#make sure all categories have a color and a name
#assert set(INT_COLORS.keys())==set(INT_NAMES.keys())

HISTS_TO_MAKE = [
### Hists Needed for XSec ######
   # "Visible Energy Migration",
   # "Q3 Migration",
   # "Visible Energy vs q3",
   # "True Signal Visible Energy vs q3",
   # "Visible Energy vs q3 Migration",
    
#### Hists Needed for pT XSec ####
   # "Lepton Pt Migration",  
   # {"variables":["Visible Energy","Lepton Pt"],
   #  "tags": {"sideband","truth_class"},
   #  },
   # "True Signal Visible Energy vs Lepton Pt",
   # "Visible Energy vs Lepton Pt Migration",
 
### Else ############    
   # "Visible Energy",
   # "Q0",
   # "Q3",
   # "Lepton Pt", 
   # "True Visible Energy",
   # {"variables":["True Visible Energy","True Q3"],
   #  "tags": {"sideband","truth_class"},
   # },
   # {"variables":["PsiEe","Lepton Energy"],
   #  "tags":{"sideband","truth_class"},
   # },
   # "True Signal Neutrino Energy",
   # {"variables":["Neutrino Energy"],
   #  "tags":{"truth_class"}},
   # "True Neutrino Energy",
   # "Neutrino Energy Migration",
   # "Lepton Energy High Inline",
   # "Lepton Energy Low Inline",
   # "Lepton Energy",


#### Thesis Plots #######
    #"Psi",
    #"Psi*Ee vs Ee",
    #"Vertex Z coordinate",
    #"Vertex R",
    #"Vertex Energy",
    #"Inline Upstream Energy",
    #"Extra Energy",
    #{"variables":["Inline Upstream Energy"],
    # "tags":{"sideband","truth_class"},
    #},
    #"Extra Energy Long",
    #"Extra Energy Trans",
    #"Transverse Shower Asymmetry",
    #"EeTheta2",

    #"Truth Visible Energy vs q3",
    #"Truth Visible Energy vs Lepton Pt",  
    #{"variables":["Visible Energy","Q3"],
    # "tags": {"sideband","truth_class"},
    #},
    #"True Visible Energy vs Visible Difference",
    # {"variables":["PsiEe","Lepton Theta"],
    #  "tags":{"sideband","truth_class","mc_only"},
    #  },
    #{"variables":["Neutrino Energy"],
    # "tags":{"truth_class"}},
    #{"variables":["Neutrino Energy QE"],
    # "tags":{"truth_class"}},
    #"Lepton Theta",
    
    ##### Neutrino4 Plots ########
    #{"variables":["Lepton Energy","Q3"],
    # "tags": {"truth_class","sideband"},
    #},
    #{"variables":["Biased Neutrino Energy","Q3"],
    #"tags": {"truth_class","sideband"},
    #},
    #{"variables":["Visible Energy","Lepton Energy"],
    #"tags": {"truth_class","sideband"},
    #},
    #{"variables":["Lepton Energy","Lepton Pt"],
    # "tags": {"truth_class","sideband"},
    #},
    #{"variables":["Biased Neutrino Energy","Lepton Pt"],
    #"tags": {"truth_class","sideband"},
    #},
    #{"variables":["Neutrino Length Travelled"],
    #"tags": {"truth_class","sideband"},
    #},
    #{"variables":["Lepton Energy"],
    # "tags": {"truth_class","sideband"},
    #},
    #"True Energy Inverse vs Biased Neutrino Energy",
    #"True Energy vs Neutrino Length Travelled",
    #"Biased Neutrino 4 Energy vs q3",
    #"Biased Neutrino 4 Energy vs Lepton Pt",
    #"Lepton Energy vs Available Energy",
    #"True Lepton Energy vs Available Energy",

    ##### Kinematic Variable Plots #####
    #{"variables":["Neutrino Length Travelled"],
    #"tags":{"truth_class","sideband"},
    #},
    #{"variables":["Lepton Pt"],
    #"tags":{"truth_class","sideband"},
    #},
    #{"variables":["Visible Energy"],
    #"tags":{"truth_class","sideband"},
    #},
    #"Available Energy vs Lepton Energy",
    #"Pion Energy vs Length Travelled",
    #"Biased Neutrino Energy",
    #{"variables":["Neutrino Length Travelled"],
    #"tags":{"truth_class","sideband"},
    #},

    #{"variables":["Lepton Pt"],
    #"tags":{"truth_class","sideband"},
    #},
    #{"variables":["Visible Energy"],
    #"tags":{"truth_class","sideband"},
    #},
    #"Available Energy vs Lepton Energy",
    #"Available Energy vs Lepton Pt",
    #"Lepton Pt vs Lepton Energy",
    #"Lepton Pt vs Available Energy",

    #"Reco Energy vs L/E",
    #"Biased Neutrino Energy",
    #"Signal Biased Neutrino Energy",
    #"True Signal Biased Neutrino Energy",
    ##### Neutrino4 Truth Plots ########
    #"Signal Reco Energy vs L/E",
    #"Signal Reco Energy vs sin",
    #"template 0",
    #"template 1",
    #"template 2",
    #"template 3",
    #"template 4",
    #"template 5",
    #"template 6",
    #"template 7",
    #"template 8",
    #"template 9",
    #"template 10",
    #"template 11",
    #"template 12",
    #"template 13",
    #"template 14",
    #"template 15",
    #"template 20",
    #"template 30",
    #"template 50",


    #{"variables": ["Inline Upstream Energy","Lepton Energy"],
   # # "tags":{"truth_class","sideband"}},
    #  {"variables":["Epi(1-cos(pi))"],
    #  "tags":{"sideband","truth_class"},
    #  },
    # {"variables":["Delta"],
    #  "tags":{"sideband","truth_class"},
    #  },
     # {"variables":["Shower Width"],
     # "tags":{"sideband","truth_class"},
     # },

     # {"variables":["Visible Zoomin","Lepton Energy"],
     # "tags": {"truth_class","sideband"},
     # "cuts":[lambda event: event.kin_cal.reco_q2_cal<0.02, lambda event: event.kin_cal.reco_Etheta2<0.0032]
     # },

     # {"variables":["Front dEdX", "Visible Energy"],
     # "tags": {"sideband","truth_class"}
     # },
    # {"variables":["Front dEdX", "Q3"],
    #  "tags": {"sideband","truth_class"}
    #  },
    # {"variables":["Front dEdX", "Lepton Pt"],
    #  "tags": {"sideband","truth_class"}
    #  },
    # {"variables":["Front dEdX", "Sum Visible Energy"],
    #  "tags": {"sideband","truth_class"}
    #  },
    #  {"variables":["Vertex Apothem","Visible Energy"],
    #  "tags": {"sideband","truth_class"}
    #  },
    #  {"variables":["Vertex Y","Visible Energy"],
    #  "tags": {"sideband","truth_class"}
    #  },
    #  {"variables":["Vertex Z","Visible Energy"],
    #  "tags": {"sideband","truth_class"}
     # },
   
    # "Q0 Migration",
    # {"variables":["Q0","True Visible Energy"],
    #  "tags": {"mc_only","signal_only"}
    #  },
    # {"variables":["Visible Energy","True Q0"],
    #  "tags": {"mc_only","signal_only"}
    #  },
    # "Lepton Pt Migration",
    # {"variables": ["Euv"],
    #  "tags":{"truth_class","sideband"}},
    # {"variables": ["Exuv"],
    #  "tags":{"truth_class","sideband"}},
    # {"variables":["W"],
    #  "tags": {"sideband","truth_class"}
    #  },
    # "Sum Visible Energy low inline",
    # "Sum Visible Energy high inline",

    # {"variables":["Print Arachne"],
    #  "tags": {"mc_only"},
    #  "cuts": [(lambda event: event.classifier.is_true_signal),(lambda universe: universe.ShortName() == "cv"),(lambda event: event.classifier.side_band == "Signal"),(lambda event :event.kin_cal.true_visE - event.kin_cal.reco_visE>0.4)]
    #  },
    #  {"variables":["Print Var"],
    #   "tags": {"mc_only"},
    #   "cuts": [(lambda event: event.classifier.is_true_signal),(lambda universe: universe.ShortName() == "cv"),(lambda event: event.classifier.side_band == "Signal"),(lambda event :event.kin_cal.true_visE - event.kin_cal.reco_visE>0.4)]
    #   },
    # {"variables":["Start Multiplicity"],
    #  "tags": {"sideband","truth_class"}
    #  },
    # {"variables":["Vertex Multiplicity"],
    #  "tags": {"sideband","truth_class"}
    #  },
    #"nue_EL",
    #"electron_energy",

    "Nu Parent Energy vs Length Travelled",
    "Nu Parent Energy vs Nu Energy",
    "Reco Energy vs Longitudinal Distance",
    "Reco Energy vs Transverse Distance",

    #"Reco Energy vs L/E",
    #"True Energy vs Biased Neutrino Energy",
    #"Estimator vs Front dEdX",
    #"Biased Neutrino Energy",


    #"Q2",
    #"Available Energy vs True W",
    #"Available Energy vs Lepton Pt",
    #"Front dEdX",
    #"True Energy vs Biased Neutrino Energy",



    #"True Signal Lepton Pt CCQE",
    #"True Signal Lepton Pt Aaron",
]

#for i in LOW_RECOIL_BIN_Q0:
#    HISTS_TO_MAKE.append("Eavail bin "+str(i))
#    HISTS_TO_MAKE.append("E Estimator Eavail bin "+str(i))
