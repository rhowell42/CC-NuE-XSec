"""
   Define your signal and background here.

   author: H. Su
"""
from collections import OrderedDict
from config.CutConfig import KINEMATICS_CUTS
from tools.CutLibrary import CUTS
from tools import TruthTools

#helper functions for signal defination

def IsFiducial(event):
    return all(map((lambda cutname:CUTS[cutname].DoesEventPass(event)),["Truth_Vertex_Z","Truth_Vertex_Apothem"] ))

IsCC = lambda event: event.mc_current == 1
IsNC = lambda event: event.mc_current == 2

IsQE =  lambda event : event.mc_intType==1 and event.mc_charm!=1
IsDelta = lambda event : event.mc_intType == 2 and event.mc_resID==0
IsDIS = lambda event : event.mc_intType == 3
IsCoherent = lambda event: event.mc_intType == 4
IsElastic = lambda event: event.mc_intType == 7
Is2p2h = lambda event: event.mc_intType == 8


IsPC = lambda event: event.mc_processType ==5
IsUnknown  = lambda event : event.mc_intType == 10

def countPDG(fsparticles,pdgs):
    return sum([1 if x in pdgs else 0 for x in fsparticles])

IsNuE = lambda event: abs(event.mc_incoming) == 12
IsAntiNu = lambda event: event.mc_incoming < 0 # assuming only neutrino incomes
IsPi0InFinalState = lambda event: 111 in event.mc_FSPartPDG
IsProtonInFinalState = lambda event: 2212 in event.mc_FSPartPDG
IsMultiMeson = lambda event: countPDG(event.mc_FSPartPDG, [211,-211, 321,-321,323,-323,111,130,310,311])>1


IsHeavyBaryon = lambda event: 3112 in event.mc_FSPartPDG or 3122 in event.mc_FSPartPDG or 3212 in event.mc_FSPartPDG or 3222 in event.mc_FSPartPDG or 4112 in event.mc_FSPartPDG or 4122 in event.mc_FSPartPDG or 4212 in event.mc_FSPartPDG or 4222 in event.mc_FSPartPDG
IsMeson = lambda event: 211 in event.mc_FSPartPDG or -211 in event.mc_FSPartPDG or 321 in event.mc_FSPartPDG or -321 in event.mc_FSPartPDG or 323 in event.mc_FSPartPDG or -323 in event.mc_FSPartPDG  or 111 in event.mc_FSPartPDG or 130 in event.mc_FSPartPDG or 310 in event.mc_FSPartPDG or 311 in event.mc_FSPartPDG
IsDeexcitationPhoton =  lambda event: event.mc_FSPartPDG[0] == 22 and event.mc_FSPartE[0] < 10
IsPhoton = lambda event: 22 in event.mc_FSPartPDG and event.mc_FSPartPDG[0] != 22

def IsInKinematicPhaseSpace(event):
    return all(CUTS["True{}".format(cut)].DoesEventPass(event) for cut in KINEMATICS_CUTS)


# In case a event satisfy multiple definations, the first takes priority.
TRUTH_CATEGORIES = OrderedDict()
TRUTH_CATEGORIES["ExcessModel"] = lambda event: IsPC(event) or IsUnknown(event)
TRUTH_CATEGORIES["NuEElastic"] = lambda event: IsElastic(event)
TRUTH_CATEGORIES["NonFiducial"] = lambda event: IsCC(event) and IsNuE(event) and not IsFiducial(event)
TRUTH_CATEGORIES["NonPhaseSpace"] = lambda event: IsCC(event) and IsNuE(event) and not IsInKinematicPhaseSpace(event)
TRUTH_CATEGORIES["CCNuEAntiNu"] = lambda event: IsCC(event) and IsNuE(event) and IsAntiNu(event)

TRUTH_CATEGORIES["CCNuEQE"] = lambda event: IsCC(event) and IsNuE(event) and IsQE(event)
TRUTH_CATEGORIES["CCNuEDelta"] = lambda event: IsCC(event) and IsNuE(event) and IsDelta(event)
TRUTH_CATEGORIES["CCNuEDIS"] = lambda event: IsCC(event) and IsNuE(event) and IsDIS(event)
TRUTH_CATEGORIES["CCNuE2p2h"] = lambda event: IsCC(event) and IsNuE(event) and Is2p2h(event)
TRUTH_CATEGORIES["CCNuE"] = lambda event: IsCC(event) and IsNuE(event)

TRUTH_CATEGORIES["CCDIS"] = lambda event: IsCC(event) and IsDIS(event)
TRUTH_CATEGORIES["CCOther"] = lambda event: IsCC(event)
TRUTH_CATEGORIES["NCCOH"] = lambda event: IsNC(event) and IsCoherent(event)
TRUTH_CATEGORIES["NCDIS"] = lambda event: IsNC(event) and IsDIS(event)
TRUTH_CATEGORIES["NCRES"] = lambda event: IsNC(event) and IsDelta(event)
TRUTH_CATEGORIES["NCOther"] = lambda event: IsNC(event) 

# My signal is one or more of the listed categories.
SIGNAL_DEFINATION = [
    "CCNuEQE",
    "CCNuEDelta",
    "CCNuEDIS",
    "CCNuE",
    "CCNuE2p2h",
]

EXTRA_OTHER = [
]
