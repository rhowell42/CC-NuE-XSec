"""
  PlotLibrary.py:
   Define plots to be handled by MakePlotProcessors.
   Configure by PlotConfig

   Original author: J. Wolcott (jwolcott@fnal.gov)
                    May 2014
"""
import ROOT
import math
from config import PlotConfig,SignalDef
from tools import TruthTools,Utilities


#tags are treated in MyHistograms.MakePlotProcessors
# sideband: make the plot in sideband sample as well. Default is only signal sample
# truth_class : make the plot for different truth_categroy seperately, only works for mc of course
# mc_only: only make this plot for mc sample
# signal_only: only make this plot for true signal, only works for mc of course.

#tags for histograms that will be used to produce cross section.
reco_tags={"sideband","truth_class"}
migration_tags={"mc_only","signal_only"} 
truth_signal_tags={"mc_only","signal_only","ignore_selection","truth_class"}

#tags for histograms for various study
resolution_tags={"sideband","mc_only"}
truth_tags={"sideband","truth_class","mc_only"}

skip_sys = lambda universe: universe.ShortName() == "cv"

"""
How to define a new plot:
1) add a entry to PLOT_SETTINGS. format:
"name of plot":
{
"name" : the name of mnvHXD,
"title": the title of mnvHXD, "title;x-axis-lable;y-axis-lable",
"binning": the binning(s) of mnvHXD. [ binning for x-axis, binning for y-axis if 2D histogram ]
           binning is a list of bin edges for example: [0,1,2,] means two bins, 0-1 and 1-2.
"value_getter" : the function to get variable(s) to be filled to mnvHXD. [lambda function for x-variable, lambda function for y-variable if 2D,]
"tags" : giving tags to plot such that it will be treated specially in MakePlotProcessors.
         available tags at the time of writing are:
         "sideband" : make this plot for sideband region ( default is for signal region only)
         "mc_only" : only make this plot for mc samples.
         "truth_class" : make variants of this plot for each truth_class defined in SignalDef.py
"cuts" : list of lambda functions decide if a entry will be fill per universe/event. can be used to make additional cut for individual plot rather all plots.
}

2) add the name of this plot to HIST_TO_MAKE in PlotConfig.py

"""



DecMap = {
    "low" : lambda func, binning: Utilities.decorator_ReLU(func,binning[0]),
    "high" : lambda func, binning: Utilities.decorator_Cap(func,binning[-2]),
}


def VariantPlotsNamingScheme(*args):
    return "_".join(args)

def TranslateSettings(key): 
    if isinstance(key,str):
        settings = PLOT_SETTINGS[key].copy()  
    elif isinstance(key,dict):
        #make a new settings out of the dict 
        settings = key.copy() 
        for _ in settings["variables"]:
            var = VARIABLE_DICT[_] 
            settings.setdefault("binning",[]).append(var["binning"])
            tmp = var["value_getter"]
            if "dec" in var:
                for key in var["dec"]:
                    tmp = DecMap[key](tmp,var["binning"])

            settings.setdefault("value_getter",[]).append(tmp)
        settings["title"] = ";"+";".join(VARIABLE_DICT[_]["title"] for _ in settings["variables"])
        settings["name"] = "_".join(VARIABLE_DICT[_]["name"] for _ in settings["variables"])
        del settings["variables"]

    return settings


VARIABLE_DICT = {
    "Leakage" : {
        "name":"LeakgageE",
        "title": "Outside Cone Energy (GeV) - TrueE*Corr",
        "binning": [0.01 * i for i in range(-20,20)],
        "value_getter" : Utilities.decorator_ReLU(lambda event: sum(event.ExtraEnergy)/1e3-0.008*event.mc_FSPartE[0]/1e3)
    },

    "Leading Pi0 E" : {
        "name":"Leading_Pi0_E",
        "title":"Leading E_{pi0} (GeV)",
        "binning":PlotConfig.PI0_ENERGY_BINNING,
        "value_getter": lambda event: (TruthTools.LeadingParticleEnergy(event,111,False) or -1.0)/1e3
    },

    
    "NuEFuzz Muon" : {
        "name":"NuECone_E_muon",
        "title":"Muon Contribution in Cone(MeV)",
        "binning": PlotConfig.NUECONE_BINNING,
        "value_getter": lambda event: event.blob_nuefuzz_E_tracker_mu+event.blob_nuefuzz_E_ecal_mu
    },

     "NuEFuzz Proton" : {
        "name":"NuECone_E_proton",
        "title":"Proton Contribution in Cone(MeV)",
        "binning": PlotConfig.NUECONE_BINNING,
        "value_getter": lambda event: event.blob_nuefuzz_E_tracker_p+event.blob_nuefuzz_E_ecal_p
    },

    "NuEFuzz E" : {
        "name":"NuECone_E",
        "title":"Energy in Cone(MeV)",
        "binning": PlotConfig.NUECONE_BINNING,
        "value_getter": lambda event: event.blob_nuefuzz_E_tracker+event.blob_nuefuzz_E_ecal
    },

    "NuEFuzz Meson" : {
        "name":"NuECone_E_meson",
        "title":"Energy in Cone(MeV)",
        "binning": PlotConfig.NUECONE_BINNING,
        "value_getter": lambda event: event.blob_nuefuzz_E_tracker_meson+event.blob_nuefuzz_E_ecal_meson
    },

     "NuEFuzz Neutron" : {
        "name":"NuECone_E_neutron",
        "title":"Energy in Cone(MeV)",
        "binning": PlotConfig.NUECONE_BINNING,
        "value_getter": lambda event: event.blob_nuefuzz_E_tracker_lown+event.blob_nuefuzz_E_ecal_lown+event.blob_nuefuzz_E_tracker_midn+event.blob_nuefuzz_E_ecal_midn+event.blob_nuefuzz_E_tracker_highn+event.blob_nuefuzz_E_ecal_highn
    },

     "NuEFuzz Xtalk" : {
        "name":"NuECone_E_xtalk",
        "title":"Energy in Cone(MeV)",
        "binning": PlotConfig.NUECONE_BINNING,
        "value_getter": lambda event: event.blob_nuefuzz_E_tracker_xtalk+event.blob_nuefuzz_E_ecal_xtalk
    },

     "NuEFuzz Other" : {
        "name":"NuECone_E_other",
        "title":"Energy in Cone(MeV)",
        "binning": PlotConfig.NUECONE_BINNING,
        "value_getter": lambda event: event.blob_nuefuzz_E_tracker_other+event.blob_nuefuzz_E_ecal_other
    },
     "NuEFuzz EM" : {
        "name":"NuECone_E_em",
        "title":"Energy in Cone(MeV)",
        "binning": PlotConfig.NUECONE_BINNING,
        "value_getter": lambda event: event.blob_nuefuzz_E_tracker_em+event.blob_nuefuzz_E_ecal_em
    },

    "Epi(1-cos(pi))" : {
        "name":"Epi_cospi",
        "title":"E_#pi(1-cos(#theta_#pi)) (GeV)",
        "binning": [0.005* i for i in range(0,41)],
        "value_getter": lambda event: event.kin_cal.reco_E_lep*(1-event.kin_cal.reco_cos_theta_lep)
    },

    "PsiEe" : {
        "name":"PsiEe",
        "title":"Psi * E_e (GeV)",
        "binning": [0.2*i for i in range(21)],
        "value_getter": lambda event: event.Psi*event.kin_cal.reco_E_lep
    },

    "Delta" : {
        "name":"Delta",
        "title":"energy ratio outer/total",
        "binning": [0.02*i for i in range(51)],
        "value_getter": lambda event: event.dedxDelta()
    },

    "Shower Width" : {
        "name":"Shower_width",
        "title":"median shower width (cm)",
        "binning": [0.1*i for i in range(31)],
        "value_getter": lambda event: event.prong_MedianPlaneShowerWidth[0]
    },

    "Visible Zoomin" : {
        "name":"Visible_zoomin",
        "title": "Visible Energy (MeV)",
        "binning": [10 * i for i in range(30)],
        "value_getter" : lambda event:event.kin_cal.reco_visE*1000,
        "dec":["low"]
    },

    "CorrCalc" : {
        "name":"CalcCorr",
        "title": "Outside Cone Energy (GeV) / TrueE",
        "binning": [0.01 * i for i in range(0,20)],
        "value_getter" : lambda event: (sum(event.ExtraEnergy)/1e3)/(event.mc_FSPartE[0]/1e3)
    },
    "True q3" :
    {
        "name" : "tQ3",
        "title": "q3 (GeV)",
        "binning" : PlotConfig.LOW_RECOIL_BIN_Q3,
        "value_getter" : lambda event: event.kin_cal.true_q3,
    },
    "True Lepton Theta" :
    {
        "name" : "teltheta",
        "title" : "#theta_{lep} (degree)", 
        "binning": [1 * i for i in range(41)],
        "value_getter": lambda event: event.kin_cal.true_theta_lep,
    },
    "Xtalk Energy in eProng":
    {
        "name" : "xtalk",
        "title" : "Xtalk dE (MeV)",
        "binning" : [5*i for i in range (51)],
        "value_getter" : lambda event: min(event.xtalkEnergy_ID_eProng,249),
    },
    "EM Energy out of eProng":
    {
        "name" : "EM_leak",
        "title" : "EM particle dE (MeV)",
        "binning" : [5*i for i in range (51)],
        "value_getter" : lambda event: min(event.emEnergy_ID_passive,249),
    },
    "EM Energy out of eProng around vtx":
    {
        "name" : "EM_leak_vtx",
        "title" : "EM particle dE (MeV)",
        "binning" : [5*i for i in range (51)],
        "value_getter" : lambda event: min(event.emEnergy_vtx_passive,249),
    },
    "Other Energy in eProng around vtx":
    {
        "name" : "other_vtx",
        "title" : "Other particle dE (MeV)",
        "binning" : [5*i for i in range (51)],
        "value_getter" : lambda event: min(event.otherEnergy_vtx_eProng,249),
    },
     "EM Energy out of eProng downstream":
    {
        "name" : "EM_leak_ds",
        "title" : "EM particle dE (MeV)",
        "binning" : [5*i for i in range (51)],
        "value_getter" : lambda event: min(event.emEnergy_ID_passive-event.emEnergy_vtx_passive,249),
    },
    "Other Energy in eProng downstream":
    {
        "name" : "other_ds",
        "title" : "Other particle dE (MeV)",
        "binning" : [5*i for i in range (51)],
        "value_getter" : lambda event: min(event.otherEnergy_ID_eProng-event.otherEnergy_vtx_eProng,249),
    },
    "Exuv":
    {
        "name" : "Exuv",
        "title" : "X-UV Energy balance",
        "binning" : PlotConfig.ENERGY_BALANCE_BINNING,
        "value_getter" : lambda event: abs(event.kin_cal.Exuv),
        "dec":["high"]
    },
     "Euv":
    {
        "name" : "Euv",
        "title" : "U-V Energy balance",
        "binning" : PlotConfig.ENERGY_BALANCE_BINNING,
        "value_getter" : lambda event: abs(event.kin_cal.Euv),
        "dec":["high"]
    },
    "Inline Upstream Energy":
    {
        "name": "inline_upstream_E",
        "title": "Log (Inline Upstream Energy) ( Log (MeV))",
        "binning": [0.25 * i for i in range(-1,11)],
        "value_getter" : lambda event: math.log10(max(0.9,event.UpstreamInlineEnergy))
    },
    "OutsideConeE" : {
        "name":"OutsideConeE",
        "title": "Outside Cone Energy (GeV)",
        "binning": [0.01 * i for i in range(-15,15)],
        "value_getter" : lambda event: sum(event.ExtraEnergy)/1e3#-0.016*event.mc_FSPartE[0]/1e3
    },
     "Psi":
    {
        "name" : "Psi",
        "title" : "Psi",
        "binning" : PlotConfig.PSI_BINNING,
        "value_getter" : lambda event: event.Psi,
    },
    "Lepton Pt":
    {
        "name" : "Lepton_Pt",
        "title" : "Pt_lepton (GeV)",
        "binning" : PlotConfig.PT_BINNING,
        "value_getter" : lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad),
        #"dec":["high"]
    },
    "Vertex Apothem":
    {
        "name":"vtx_apothem",
        "title":"Vertex Apothem(mm)",
        "binning": PlotConfig.APOTHEM_BINNING,
        "value_getter": lambda event: TruthTools.CalcApothem(event.vtx[0],event.vtx[1])
    },
    "True Lepton Pt":
    {
        "name" : "tLepton_Pt",
        "title" : "tPt_lepton",
        "binning" : PlotConfig.PT_BINNING,
        "value_getter" : lambda event: event.kin_cal.true_P_lep*math.sin(event.kin_cal.true_theta_lep_rad)
    },
    "Visible Energy":
    {
        "name" : "Eavail",
        "title" : "E_{avail} (GeV)",
        #"binning" : [PlotConfig.VISIBLE_ENERGY_BINNING],
        "binning" : PlotConfig.LOW_RECOIL_BIN_Q0,
        "value_getter" : lambda event: event.kin_cal.reco_visE
    },
    "Sum Visible Energy":
    {
        "name" : "sum_visE",
        "title" : "Evis_{event} (GeV)",
        #"binning" : [PlotConfig.VISIBLE_ENERGY_BINNING],
        "binning" : PlotConfig.SUM_VISIBLE_ENERGY_BINNING,
        "value_getter" : lambda event: event.prong_passive_idclus[0]/1e3+event.kin_cal.reco_visE,
        "dec":["high"]
    },
    "Lepton Visible Energy":
    {
        "name" : "lepton_visE",
        "title" : "Evis_{lepton} (GeV)",
        #"binning" : [PlotConfig.VISIBLE_ENERGY_BINNING],
        "binning" : PlotConfig.ELECTRON_ENERGY_BINNING,
        "value_getter" : lambda event: event.prong_passive_idclus[0]/1e3
    },
    "Lepton Energy Course Binning":
    {
        "name" : "Eel",
        "title" : "Reconstructed E_e (GeV)", 
        "binning" : PlotConfig.EXCESS_ENERGY_FIT,
        "value_getter" : lambda event: event.kin_cal.reco_E_lep
    },
     "Lepton Theta":
    {
        "name" : "eltheta",
        "title" : "#theta_{lep} (degree)",
        "binning" : PlotConfig.ELECTRON_ANGLE_BINNING,
        "value_getter" : lambda event: event.kin_cal.reco_theta_lep,
    },
    "Lepton Energy":
    {
        "name" : "Eel",
        "title" : "Reconstructed E_lep (GeV)", 
        "binning" : PlotConfig.ELECTRON_ENERGY_BINNING,  
        "value_getter" : lambda event: event.kin_cal.reco_E_lep,
    },
    "Electron Candidate True PID":
    {
        "name": "TruePID",
        "title": "Log10(TruePID of Electron Prong)",
        "binning" : [x * 0.1 for x in range(-10, 41)],
        "value_getter" : lambda event : math.log10(max(0.1,abs(event.prong_TruePID[0])))
    },
    "N_EMLike_prongs":
    {
        "name": "n_prongs",
        "title": "n_prongs",
        "binning" : [i for i in range(1,10)],
        "value_getter" : lambda event : event.n_prongs
    },
     "Q3" :
    {
        "name" : "Q3",
        "title": "q3 (GeV)",
        "binning" : PlotConfig.LOW_RECOIL_BIN_Q3,
        "value_getter" : lambda event: event.kin_cal.reco_q3,
    },
    "Print Arachne":
    {
        "name":"dummy",
        "title":"dummy",
        "binning":[0,1,2],
        "value_getter": lambda event: Utilities.ArachneLink(event) or 1
    },
    "Print Var":
    {
        "name":"dummy",
        "title":"dummy",
        "binning":[0,1,2],
        "value_getter": lambda event: PlotConfig.Printvar(event)
    },
    "Proton Angle":
    {
        "name":"p_theta",
        "title": "Proton Theta (degree)",
        "binning": PlotConfig.PROTON_ANGLE_BINNING,
        "value_getter": lambda event:event.PrimaryProtonTheta()
    },
    "True Proton Angle":
    {
        "name":"tp_theta",
        "title": "True Proton Theta (degree)",
        "binning": PlotConfig.PROTON_ANGLE_BINNING,
        "value_getter": lambda event:TruthTools.FindThetaFromIndex(event,TruthTools.MostEnergeticParticle(event,2212,False))
    },
    "Proton Angle Residual":
    {
        "name":"dtp_theta",
        "title": "True Proton Theta Residual (reco-truth)",
        "binning": PlotConfig.PROTON_ANGLE_RESIDUAL_BINNING,
        "value_getter": lambda event: event.PrimaryProtonTheta() - (TruthTools.FindThetaFromIndex(event,TruthTools.MostEnergeticParticle(event,2212,False)) or  -9999),
        "dec" : ["low","high"]
    },
    "W" :
    {
        "name":"Wexp",
        "title": "Reconstructed W (GeV)",
        "binning": PlotConfig.RECO_W_BINNING,
        "value_getter": lambda event: event.kin_cal.reco_W,
        "dec" : ["high"]
    },
     "Qsquare":
    {
        "name" : "Q2",
        "title" : "Q^{2}_{cal} (GeV^2)",
        "binning" : PlotConfig.QSQUARED_BINNING,
        "value_getter" : lambda event: event.kin_cal.reco_q2_cal,
        "dec" : ["high"]
    },
    "Vertex Y" :
    {
        "name": "Vertex_y",
        "title": "Vertex y position",
        "binning":PlotConfig.VERTEX_Y_BINNING,
        "value_getter": lambda event: event.vtx[1],
    },
    "Neutrino Energy":
    {
        "name" : "Enu",
        "title" : "E_{#nu} (GeV)",
        "binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "value_getter" : lambda event: event.kin_cal.reco_E_nu_cal,
    },
     "True Neutrino Energy":
    {
        "name" : "tEnu",
        "title" : "E_{#nu} (GeV)",
        "binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "value_getter" : lambda event: event.kin_cal.true_enu_genie,
    },
     "Neutrino Energy QE":
    {
        "name" : "Enu_qe",
        "title" : "E_{#nu}^{QE} (GeV)",

        "binning" : PlotConfig.NEUTRINO_ENERGY_BINNING,
        "value_getter" : lambda event: event.kin_cal.reco_E_nu_QE,
    },
    "Start Multiplicity":
    {
        "name" : "Mul_start",
        "title" : "Startpoint Multiplicity",
        "binning" : [i for i in range(11)],
        "value_getter" : lambda event: event.StartPointVertexMultiplicity,
    },
    "Vertex Multiplicity":
    {
        "name" : "Mul_vtx",
        "title" : "Vtx Multiplicity",

        "binning" : [i for i in range(11)],
        "value_getter" : lambda event: event.VertexTrackMultiplicity,
    },
    "VertexDiff" : {
        "name" : "diff_vtx_z",
        "title": "Diff evt vtx and prong vtx ",
        "binning" : [-200+ 20*i for i in range(21)],
        "value_getter" : lambda event: event.prong_axis_vertex[0][2]-event.vtx[2],
        "dec":["high","low"]
    },
    "Vertex Difference" : {
        "name" : "RecoTruth_vtx_z",
        "title": "Vtx Diff ; Reco - True Vtx; NEvents",
        "binning" : PlotConfig.VERTEX_DIFF_BINNING,
        "value_getter" : lambda event: TruthTools.vtxdiff(event),
        "dec":["high","low"]
    },
    "Vertex Z" : {
        "name" : "vtx_z",
        "title": "Vertex Z",
        "binning" : PlotConfig.VERTEX_Z_BINNING,
        "value_getter" : lambda event: event.vtx[2],
        "dec":["high","low"]
    },
     "Front dEdX":
    {
        "name": "frontdedx",
        "title": "Mean front dE/dx",
        "binning": PlotConfig.DEDX_BINNING,
        "value_getter" : lambda event: event.prong_dEdXMeanFrontTracker[0],
        "dec":["high"],
    },
    "True Visible Energy":
    {
        "name" : "tEavail",
        "title" : "True Available Energy",
        "binning" : PlotConfig.LOW_RECOIL_BIN_Q0,
        "value_getter" : lambda event: event.kin_cal.true_visE,
    },
    "Q0" :
    {
        "name" : "Q0",
        "title": "q0 (GeV)",
        "binning" : PlotConfig.LOW_RECOIL_BIN_Q0,
        "value_getter" : lambda event: event.kin_cal.reco_q0,
    },
    "True Q0" :
    {
        "name" : "tQ0",
        "title": "True q0 (GeV)",
        "binning" : PlotConfig.LOW_RECOIL_BIN_Q0,
        "value_getter" : lambda event: event.kin_cal.true_q0,
    },
     "UIE Mean Pos" :
    {
        "name" : "UIE_Mean_Pos",
        "title": "Mean UIE Position (rad. length)",
        "binning" : [0.2*i for i in range(21)],
        "value_getter" : lambda event: event.UpstreamInlineEnergyWgtdPosMean,
        "dec":["high"],
    },
    "UIE Pos RMS" :
    {
        "name" : "UIE_Pos_RMS",
        "title": "UIE Position RMS (rad. length)",
        "binning" : [0.1*i for i in range(21)],
        "value_getter" : lambda event: event.UpstreamInlineEnergyWgtdPosRMS,
        "dec":["high"],
    },
     "Rock Muon Removed" :
    {
        "name" : "rockmu_rmed",
        "title": "Rock Muon Removed",
        "binning" : [0,1,2],
        "value_getter" : lambda event: int(event.rock_muons_removed),
    },
    "Vertex Energy":
    {
        "name":"vertex_blob_E_50mm",
        "title":"Vertex Energy in 50mm ",
        "binning": [20* i for i in range(26)],
        "value_getter": lambda event: event.recoil_energy_nonmuon_vtx50mm,
    }
    
}

PLOT_SETTINGS= {
    "Pseudo q3":
    {
        "name": "pseudo_q3",
        "title": "pseudo q3 (Gev); Nevents/GeV",
        "binning": [PlotConfig.RECO_W_BINNING],
        "value_getter":[lambda event: event.kin_cal.pseudo_q3],
        "tags":reco_tags
    },

    "FrontdEdX":
    {
        "name": "frontdedx",
        "title": "Mean front dE/dx; FrontdEdX; Nevents/GeV",
        "binning": [PlotConfig.FRONTDEDX_BINNING],
        "value_getter" : [lambda event: event.prong_dEdXMeanFrontTracker[0]],
        "tags" : reco_tags,
        #"cuts": [skip_sys, lambda event: event.classifier.truth_class == "notCCQElike"]
    },
    "Leakage":
    {
        "name": "LeakageE",
        "title": "Outside Cone Energy (GeV) - TrueE*Corr",
        "binning": [PlotConfig.LEAKAGE_BINNING],
        "value_getter" : [lambda event: TruthTools.Leakage(event)],
        "tags" : reco_tags,
        #"cuts": [skip_sys, lambda event: event.classifier.truth_class == "notCCQElike"]
    },
    "True Lepton Theta" :
    {
        "name" : "teltheta",
        "title" : "#theta_{e} (degree)",
        "binning": [PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter": [lambda event: event.kin_cal.true_theta_lep],
        "tags" : truth_tags,
    },
    "FrontdEdX_position":
    {
        "name": "frontdedx_position",
        "title": "Position of Min front dE/dx; Position (mm); Nevents/GeV",
        "binning": [PlotConfig.FRONTDEDX_POSITION_BINNING],
        "value_getter" : [lambda event: event.prong_dEdXMeanFrontPositionTracker[0]],
        "tags" : reco_tags
    },
    "FrontdEdX_front":
    {
        "name": "frontdedx_front",
        "title": "Mean front dE/dx front; FrontdEdX; Nevents/GeV",
        "binning": [PlotConfig.ELECTRON_ANGLE_BINNING],
        #"binning": [PlotConfig.FRONTDEDX_FRONT_BINNING],
        "value_getter" : [lambda event: event.prong_dEdXMeanFrontTracker[0]],
        "tags" : reco_tags
    },

    "FrontdEdX_tail":
    {
        "name": "frontdedx_tail",
        "title": "Mean front dE/dx tail; FrontdEdX; Nevents/GeV",
        "binning": [PlotConfig.FRONTDEDX_TAIL_BINNING],
        "value_getter" : [lambda event: event.prong_dEdXMeanFrontTracker[0]],
        "tags" : reco_tags
    },

    "Lepton Theta":
    {
        "name" : "eltheta",
        "title" : "Lepton Theta; #theta_{lep} (degree); dNEvents/d#theta",
        "binning" : [PlotConfig.ELECTRON_ANGLE_BINNING],
       # "binning" : [PlotConfig.EXCESS_ANGLE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_theta_lep],
        "tags": reco_tags,
        #"cuts": [skip_sys, lambda event: event.classifier.truth_class == "CCQElike"]
    },
    "Lepton Energy":
    {
        "name" : "Eel",
        "title" : "Lepton Energy; Reconstructed E_lep (GeV); dNEvents/dE_lep", 
        "binning" : [PlotConfig.ELECTRON_ENERGY_BINNING],  
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep],    
        #"value_getter" : [lambda event: TruthTools.Print(event)], 
        "tags": reco_tags, 
    },
    "Lepton Energy Low Inline":
    {
        "name" : "Eel_low_inline",
        "title" : "Lepton Energy; Reconstructed E_lep (GeV); dNEvents/dE_lep", 
        "binning" : [PlotConfig.ELECTRON_ENERGY_BINNING],  
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep],    
        "cuts": [lambda event: TruthTools.Rebin(event)<10,lambda event: event.kin_cal.reco_q3<1.2],
        "tags": reco_tags, 
    },
    "Lepton Energy High Inline":
    {
        "name" : "Eel_high_inline",
        "title" : "Lepton Energy; Reconstructed E_lep (GeV); dNEvents/dE_lep", 
        "binning" : [PlotConfig.ELECTRON_ENERGY_BINNING],  
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep],
        "cuts": [lambda event: TruthTools.Rebin(event)>10,lambda event: event.kin_cal.reco_q3<1.2],
        "tags": reco_tags, 
    },

     "Low nu Lepton Energy":
    {
        "name" : "Eel_low_nu",
        "title" : "Lepton Energy; E_lep (GeV); dNEvents/dE_lep",
        "binning" : [PlotConfig.ELECTRON_ENERGY_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep],
        "tags": reco_tags,
        "cuts": [ lambda event: event.kin_cal.reco_q0<0.8]
    },
    "Neutrino Energy":
    {
        "name" : "Enu",
        "title" : "Neutrino Energy; E_{#nu} (GeV); dNEvents/dE_{#nu}",

        "binning" : [PlotConfig.NEUTRINO_ENERGY_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_nu_cal],
        "tags": reco_tags
    },
    "True Lepton Energy":
    {
        "name" : "tEel",
        "title" : "True Lepton Energy; E_e (GeV); dNEvents/dE_e",
        "binning" : [PlotConfig.EXCESS_ENERGY_FIT],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep],
        "tags": reco_tags

    },
    "True Neutrino Energy":
    {
        "name" : "tEnu",
        "title" : "True Neutrino Energy; E_{#nu} (GeV); NEvents/GeV",
        "binning" : [PlotConfig.NEUTRINO_ENERGY_BINNING],
        "value_getter" : [lambda event:event.kin_cal.true_enu_genie],
        "tags": truth_tags
    },

    "True Lepton Energy":
    {
        "name" : "tEl",
        "title" : "True Lepton Energy; True Lepton E (GeV); NEvents",
        "binning" : [PlotConfig.ELECTRON_ENERGY_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep], 
        "tags": truth_tags
    },

    "GENIE Event Information":
    {
        "name" : "tEvent",
        "title" : "GENIE Event Type; mcIntType; NEvents",
        "binning" : [PlotConfig.GENIE_EVENT_BINNING],
        "value_getter" : [lambda event: TruthTools.HighEnergyPi0(event, TruthTools.MostEnergeticParticle(event, 111), event.mc_intType)],
        "tags": reco_tags
    },

    "Fractional Pi0 E/Eevent":
    {
        "name" : "tFractionPi0E",
        "title" : "Fractional Pi0 Energy; Pi0E / TotalE; NEvents",
        "binning" : [PlotConfig.FRACTION_BINNING],
        "value_getter" : [lambda event: TruthTools.LeadingPi0Fraction(event, TruthTools.MostEnergeticParticle(event, 111))],
        "tags": reco_tags, 
    },

    "Sum Visible Energy high inline":
    {
        "name" : "Sum_visE_high_inline",
        "title" : "Evis_{lepton} (GeV)",
        "binning" : [PlotConfig.ELECTRON_ENERGY_BINNING[:15]],
        "value_getter" : [lambda event: event.prong_TotalVisE[0]/1e3+event.kin_cal.reco_visE],
        "tags": reco_tags,
        "cuts": [lambda event: TruthTools.Rebin(event)>10],
    },


    "Sum Visible Energy low inline":
    {
        "name" : "Sum_visE_low_inline",
        "title" : "Evis_{lepton} (GeV)",
        "binning" : [PlotConfig.ELECTRON_ENERGY_BINNING[:15]],
        "value_getter" : [lambda event: event.prong_TotalVisE[0]/1e3+event.kin_cal.reco_visE],
        "tags": reco_tags,
        "cuts": [lambda event: TruthTools.Rebin(event)<10],
    },

    "True W":
    {
        "name" : "tW",
        "title" : "True W; W (GeV); NEvents",
        "binning" : [PlotConfig.W_BINNING],
        "value_getter" : [lambda event: TruthTools.HighEnergyPi0(event, TruthTools.MostEnergeticParticle(event, 111), event.mc_intType)], 
        "tags": reco_tags,
    },

    "Visible Energy":
    {
        "name" : "Eavail",
        "title" : "Available Energy; E_{avail} (GeV); dNEvents/dE_{avail}",
        #"binning" : [PlotConfig.VISIBLE_ENERGY_BINNING],
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0],
        #"value_getter" : [lambda event: TruthTools.EvailStudy(event)],
        "value_getter" : [lambda event: event.kin_cal.reco_visE],
        "tags": reco_tags
    },

    "True Visible Energy":
    {
        "name" : "tEavail",
        "title" : "Available Energy; E_{avail} (GeV); dNEvents/dE_{avail}",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0],
        "value_getter" : [lambda event: event.kin_cal.true_visE],
        "tags": truth_tags
    },

    "Qsquare":
    {
        "name" : "Q2",
        "title" : "QSquared; Q^{2}_{cal} (GeV^2); dNEvents/dQ^{2}_{cal}",

        "binning" : [PlotConfig.QSQUARED_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_q2_cal],
        "tags": reco_tags
    },
    "True Qsquare":
    {
        "name" : "tQ2",
        "title" : "QSquared; Q^{2}_{cal} (GeV^2); dNEvents/dQ^{2}_{cal}",

        "binning" : [PlotConfig.QSQUARED_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_q2],
        "tags": truth_tags
    },
    "Psi":
    {
        "name" : "Psi",
        "title" : "Psi; Psi; dNEvents/dPsi", 
        "binning" : [PlotConfig.PSI_BINNING],
        "value_getter" : [lambda event: event.Psi],
        "tags": reco_tags, 
    },

    "Psi_front":
    {
        "name" : "Psi_front",
        "title" : "Psi; Psi; NEvents/GeV",
        "binning" : [PlotConfig.PSI_FRONT_BINNING],
        "value_getter" : [lambda event: event.Psi],
        "tags": reco_tags
    },


    "Psi_tail":
    {
        "name" : "Psi_tail",
        "title" : "Psi; Psi; NEvents/GeV",
        "binning" : [PlotConfig.PSI_TAIL_BINNING],
        "value_getter" : [lambda event: event.Psi],
        "tags": reco_tags
    },

    "Q0" :
    {
        "name" : "Q0",
        "title": "q0 ; q0 (GeV); NEvents/GeV",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0],
        "value_getter" : [lambda event: event.kin_cal.reco_q0],
        "tags": reco_tags,
    },

    "Q3" :
    {
        "name" : "Q3",
        "title": "q3 ; q3 (GeV); dNEvents/dq3",

        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q3],
        "value_getter" : [lambda event: event.kin_cal.reco_q3],
        "tags": reco_tags
    },

    "BkgFit_q3" :
    {
        "name" : "Q3_bkgfit",
        "title": "q3 ; q3 (GeV); dNEvents/dq3",
        "binning" : [PlotConfig.BACKGROUND_FIT_Q3_BIN],
        "value_getter" : [lambda event: event.kin_cal.reco_q3],
        "tags": reco_tags
    },

    "True Q0" :
    {
        "name" : "tQ0",
        "title": "q0 ; q0 (GeV); dNEvents/dq0",

        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0],
        "value_getter" : [lambda event: event.kin_cal.true_q0],
        "tags": truth_tags,
    },
    "True q3" :
    {
        "name" : "tQ3",
        "title": "q3 ; q3 (GeV); dNEvents/dq3",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q3],
        "value_getter" : [lambda event: event.kin_cal.true_q3],
        "tags": truth_tags
    },
    "True Pi0 Energy" :
    {
        "name" : "tEpi0",
        "title": "Epi0 ; Epi0 (GeV); NEvents",
        "binning" : [PlotConfig.PI0_ENERGY_BINNING],
        "value_getter" : [lambda event: TruthTools.EnergyofParticle(event, 111)],
        "tags": reco_tags
    },
    # "True Lepton Pt" :
    # {   
    #     "name" : "tEePt",
    #     "title": "ETruePt ; True_E_Pt (GeV); NEvents/GeV",
    #     "binning" : [PlotConfig.PT_BINNING],
    #     "value_getter" : [lambda event: TruthTools.TransverseMomentum(event, 11)],
    #     "tags": reco_tags
    # },

    "True PDG" :
    {
        "name" : "PDG",
        "title": "PDG ; PDG Code; NEvents",
        "binning" : [PlotConfig.PARTICLE_PDG_BINNING],
        "value_getter" : [lambda event: TruthTools.ParticlePDG(event, TruthTools.MostEnergeticParticle(event, 111))],
        "tags": reco_tags
    },

    "Lepton Pt":
    {
        "name" : "Lepton_Pt",
        "title" : "Pt_lepton (GeV); NEvents",
        "binning" : [PlotConfig.PT_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad)],
        "tags":reco_tags
    },
    "Energy Inside Cone" :
    {
        "name" : "InsideConeE",
        "title": "InsideConeE ; InsideCone_E (GeV); NEvents/GeV",
        "binning" : [PlotConfig.EXTRA_ENERGY_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_InsideCone_E],
        "tags": reco_tags
    },

    "Energy Outside Cone" :
    {
        "name" : "OutsideConeE",
        "title": "OutsideConeE ; Psi*Ee (GeV); NEvents/GeV",
        "binning" : [PlotConfig.OUTSIDE_ENERGY_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_OutsideCone_E],
        "tags": reco_tags
    },

    #"Cluster dE PDG" :
    #{   
    #    "name" : "ClusterdEPDG high",
    #    "title": "ClusterdEPDG ; de Particle PDG; NEvents",
    #    "binning" : [PlotConfig.CLUSTER_PDG],
    #    "value_getter" : [lambda event: TruthTools.dE(event)],
    #    "tags": reco_tags,
    #    #"cuts": [skip_sys, lambda event: event.classifier.truth_class == "notCCQElike"]
    #},

    #"Fraction of Bin Electron" :
    #{
    #    "name" : "BinFraction_electron",
    #    "title": "BinFraction ; Electron Fraction per Event; NEvents",
    #    "binning" : [PlotConfig.FRACTION_BINNING],
    #    "value_getter" : [lambda event: TruthTools.dE(event)],
    #    "tags": reco_tags,
    #    #"cuts": [skip_sys, lambda event: event.classifier.truth_class == "notCCQElike"]
    #},
    "Vertex Difference":
    {
        "name" : "RecoTruth_vtx_z",
        "title": "Vtx Diff ; Reco - True Vtx; NEvents",
        "binning" : [PlotConfig.VERTEX_DIFF_BINNING],
        "value_getter" : [lambda event: TruthTools.vtxdiff(event)],
        "tags": reco_tags,
    },
    "Median Plane Shower Width" :
    {
        "name": "Shower_width",
        "title" : "Shower Width; Median Shower Width; NEvents",
        "binning": [PlotConfig.SHOWER_WIDTH_BINNING],
        "value_getter" :[ lambda event: event.prong_MedianPlaneShowerWidth[0]],
        "tags" : reco_tags,
    },
    "Width Ee" :
    {
        "name": "Shower_width_Ee",
        "title" : "Shower Width Ee; Median Shower Width; NEvents",
        "binning": [PlotConfig.SHOWER_WIDTH_BINNING],
        "value_getter" :[ lambda event: TruthTools.binEe(event)] ,
        "tags" : reco_tags,
    },

################Jeremy comparison plots ###########################

    "Vertex Z coordinate" :
    {
        "name": "Vertex_zcoord",
        "title" : "Vertex_zcoord; Vertex Z Coordinate (cm); NEvents/cm",
        "binning": [PlotConfig.ZCOORDINATE_BINNING],
        "value_getter" :[ lambda event: TruthTools.Zcoordinate(event)],
        "tags" : reco_tags,
    },
    "Inline Upstream Energy" :
    {
        "name": "InlineUpstream_Energy",
        "title" : "InlineUpstream; In-line energy upstream of vertex (MeV); NEvents ",
        "binning": [PlotConfig.UPSTREAM_INLINE_ENERGY_BINS],
        #"value_getter" :[ lambda event: event.UpstreamInlineEnergy],
        "value_getter" :[ lambda event: TruthTools.Rebin(event)],
        "tags" : reco_tags,
    },
    "Inline Upstream Energy Weighted Position" :
    {
        "name": "InlineUpstream_WeightPos",
        "title" : "InlineUpstream; In-line energy weight position (rad. lengths); NEvents ",
        "binning": [PlotConfig.INLINE_W_BINNING],
        "value_getter" :[ lambda event: event.UpstreamInlineEnergyWgtdPosMean],
        "tags" : reco_tags,
    },
    "Inline Upstream Energy Low" :
    {
        "name": "InlineUpstream_LowEnergy",
        "title" : "InlineUpstreamLow; In-line energy upstream of vertex (MeV); NEvents ",
        "binning": [PlotConfig.INLINE_LOW_BINNING],
        "value_getter" :[ lambda event: event.UpstreamInlineEnergy],
        "tags" : reco_tags,
    },
    "Inline Upstream Energy High" :
    {
        "name": "InlineUpstream_HighEnergy",
        "title" : "InlineUpstreamHigh; In-line energy upstream of vertex (MeV); NEvents ",
        "binning": [PlotConfig.INLINE_HIGH_BINNING],
        "value_getter" :[ lambda event: event.UpstreamInlineEnergy],
        "tags" : reco_tags,
    },
    "Vertex R" :
    {
        "name": "Vertex_R",
        "title" : "VertexR2; Vertex R^2 (cm^2); NEvents ",
        "binning": [PlotConfig.VERTEXR_BINNING],
        "value_getter" :[ lambda event: TruthTools.VertexR(event)],
        "tags" : reco_tags,
    },
    "Extra Energy":
    {
        "name" : "Extrael",
        "title" : "Electron Extra Energy; Extra Energy (GeV); dNEvents/dE_e",
        "binning" : [PlotConfig.EXTRA_ENERGY_BINNING],
        "value_getter" : [lambda event: event.ExtraEnergy[0]/1e3], 
        "tags": reco_tags,
    },
    "e Radial":
    {
        "name" : "Etraj_radial",
        "title" : "Shower Axis; e traj. radial porjection; Events",
        "binning" : [PlotConfig.RADIAL_BINNING],
        "value_getter" : [lambda event: TruthTools.Print(event)],
        "tags": reco_tags,
    },
    "EeTheta2":
    {
        "name" : "EeTheta2",
        "title" : "EeTheta2; Reconstructed EeTheta2 (GeV x radian2); Events/(GeV rad2)",
        "binning" : [PlotConfig.ETH2_BINNING],
        "value_getter" : [lambda event: TruthTools.Eth2(event)],
        "tags": reco_tags,
    },
    "EeTheta":
    {
        "name" : "EeTheta",
        "title" : "EeTheta; Reconstructed EeTheta (GeV x radian); dEvents/(GeV rad)",
        "binning" : [PlotConfig.EETHETA_BINNING], 
        "value_getter" : [lambda event: TruthTools.Eth(event)],
        "tags": reco_tags,
    },
    "True EeTheta":
    {
        "name" : "tEeTheta",
        "title" : "EeTheta; Reconstructed EeTheta (GeV x radian); dEvents/(GeV rad)",
        "binning" : [PlotConfig.EETHETA_BINNING],
        "value_getter" : [lambda event: TruthTools.tEth(event)],
        "tags": reco_tags,
    },
    "Extra energy ratio":
    {
        "name" : "EeRatio",
        "title" : "Psi; Extra energy ratio; dNEvents/dPsi",
        "binning" : [PlotConfig.EXTRA_ENERGY_BINNING],
        "value_getter" : [lambda event: event.Psi],
        "tags": reco_tags,
    },
    "Transverse Shower Asymmetry":
    {
        "name" : "TransAsymm",
        "title" : "Transverse Shower Asymmetry; Transverse Shower Asymmetry; dNEvents/dA",
        "binning" : [PlotConfig.TRANSVERSE_SHOWER_BINNING],
        "value_getter" : [lambda event: TruthTools.TransverseShower(event)],
        "tags": reco_tags,
    },
     "Vertex Energy":
    {
        "name" : "VertexEnergy",
        "title" : "Vertex Energy; Vertex Attached Energy (MeV); dNEvents/dA",
        "binning" : [PlotConfig.VERTEX_ENERGY_BINS],
        "value_getter" : [lambda event: event.VertexEnergy[2]*1e3],
        "tags": reco_tags,
    },

################## 2D Reco Plots: ##################################


     "True Z Vertex":
    {
        "name" : "tzVertex",
        "title": "zVertex; zVertex(mm); NEvents",
        "binning" : [PlotConfig.WHOLE_DET_Z_BINS],
        "value_getter" : [lambda event: event.mc_vtx[2]],
        "tags": truth_signal_tags,
        'cuts': [ lambda event: event.classifier.reco_cuts_passed["NoCut"], lambda event: event.classifier.reco_cuts_passed["EMLikeTrackScore"]]
    },

    "Front dEdX":
    {
        "name" : "frontdEdX",
        "title" : "FrontdEdX; dEdX (MeV/cm);  NEvents",
        "binning" : [PlotConfig.DEDX_BINNING],
        "value_getter" : [Utilities.decorator_Cap(lambda event: event.prong_dEdXMeanFrontTracker[0],PlotConfig.DEDX_BINNING[-2])],
        "tags":reco_tags
    },

    ################## 2D Reco Plots: ##################################
    "q3 vs Visible Energy":
    {
        "name" : "q3_Eavail",
        "title": "q3 vs Available Energy; q3 (GeV); Available Energy (GeV) ;d^{2}NEvents/dE_{avail}dq3",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q3,
                     PlotConfig.VISIBLE_ENERGY_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_q3,
                          lambda event: event.kin_cal.reco_visE],
        "tags": reco_tags
    },

    "Visible Energy vs q3":
    {
        "name" : "Eavail_q3",
        "title": "Available Energy v.s. q3; Available Energy (GeV); q3 (GeV); d^{2}NEvents/dE_{avail}dq3",

        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0,
                     PlotConfig.LOW_RECOIL_BIN_Q3],
        "value_getter" : [lambda event: event.kin_cal.reco_visE,
                          lambda event: event.kin_cal.reco_q3],
        "tags": reco_tags
    },

    "q0 vs q3":
    {
        "name" : "q0_q3",
        "title": "q0 v.s. q3; q0 (GeV); q3 (GeV); d^{2}NEvents/dq0dq3",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0,
                     PlotConfig.LOW_RECOIL_BIN_Q3],
        "value_getter" : [lambda event: event.kin_cal.reco_q0,
                          lambda event: event.kin_cal.reco_q3],
        "tags": reco_tags
    },
    "Long vs Trans":
    {
        "name" : "Long_Trans",
        "title": "Long v.s. Trans; Trans (m); Long (m); d^{2}ClusterE/dLdT",
        "binning" : [PlotConfig.TRANS_BINNING,
                     PlotConfig.LONG_2D_BINNING],
        "value_getter" : [lambda event: TruthTools.Trans(event),
                          lambda event: TruthTools.Long(event)],
        #"weight_function": lambda event: TruthTools.GetWeight(event),    
        "tags": reco_tags
    },
    "Lepton Theta vs q3":
    {
        "name" : "eltheta_q3",
        "title" : "Lepton Theta v.s. q3; #theta_{lep} (degree); q3 (GeV); d^{2}NEvents/dq3d#theta_{e}",
        "binning" : [PlotConfig.EXCESS_ANGLE_BINNING,
                     PlotConfig.BACKGROUND_FIT_Q3_BIN],
        "value_getter" : [lambda event: event.kin_cal.reco_theta_lep,
                          lambda event: event.kin_cal.reco_q3],
        "tags":reco_tags
    },

    "Lepton Energy vs q3":
    {
        "name" : "Eel_q3",
        "title" : "Lepton Energy v.s. q3; E_{lep} (GeV); q3 (GeV); d^{2}NEvents/dq3dE_{e}",
        "binning" : [PlotConfig.EXCESS_ENERGY_FIT,
                     PlotConfig.BACKGROUND_FIT_Q3_BIN],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep,
                          lambda event: event.kin_cal.reco_q3],
        "tags":reco_tags
    },

    "Lepton Energy vs Theta":
    {   
        "name" : "Eel_eltheta",
        "title" : "Lepton Energy v.s. Theta;  #theta_{lep} (degree);Lepton Energy (GeV); d^{2}NEvents/dE_{e}d#theta_{e}",
        "binning" : [PlotConfig.EXCESS_ANGLE_BINNING,
                     PlotConfig.EXCESS_ENERGY_FIT],
        "value_getter" : [lambda event:event.kin_cal.reco_theta_lep,
                          lambda event: event.kin_cal.reco_E_lep],
        "tags":reco_tags
    },
    "Lepton Energy vs EeTheta":
    {
        "name" : "Eel_Eeeltheta",
        "title" : "Lepton Energy v.s. EeTheta;  Ee#theta_{e} (GeV Radians);Lepton Energy (GeV); d^{2}NEvents/dE_{e}dEe#theta_{e}",
        "binning" : [PlotConfig.EETHETA_BINNING,
                     PlotConfig.EXCESS_ENERGY_FIT],
        "value_getter" : [lambda event: TruthTools.Eth(event),
                          lambda event: event.kin_cal.reco_E_lep],
        "tags":reco_tags
    },
    "Lepton Energy vs q0":
    {
        "name": "Eel_highq0",
        "title": "Lepton Energy vs q0; E_{e} (GeV); Eavailable (GeV); NEvents",
        "binning" : [PlotConfig.EXCESS_ENERGY_FIT,
                     PlotConfig.LOW_RECOIL_BIN_HIGH_Q0],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep, 
                          lambda event: event.kin_cal.reco_q0],
        "tags":reco_tags,
    #    "cuts" : [skip_sys, lambda event: event.classifier.truth_class == "NuEElastic"] 
    },

    "Visible Energy vs Lepton Pt":
    {
        "name" : "Eavail_Lepton_Pt",

        "title": "Available Energy v.s. q3; E_{avail} (GeV); Lepton Pt (GeV); NEvents",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0,
                     PlotConfig.PT_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_visE,
                          lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad)],
        "tags": reco_tags   
    },


    "Lepton Energy vs low q0":
    {
        "name": "Eel_lowq0",
        "title": "Lepton Energy vs low q0; E_{e} (GeV); Eavailable (GeV); NEvents",
        "binning" :  [PlotConfig.EXCESS_ENERGY_FIT,
                     PlotConfig.LOW_RECOIL_BIN_LOW_Q0],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep,
                          lambda event: event.kin_cal.reco_q0],
        "tags":reco_tags,
    #    "cuts" : [skip_sys, lambda event: event.classifier.truth_class == "NuEElastic"]
    },
    "Front dEdX vs q3":
    {
        "name" : "frontdEdX_q3",
        "title" : "Front dEdX v.s. q3; dEdX (MeV/cm); q3 (GeV); d^{2}NEvents/d(dEdX)dq3",
        "binning" : [PlotConfig.DEDX_BINNING,
                     PlotConfig.BACKGROUND_FIT_Q3_BIN],
        "value_getter" : [lambda event: event.prong_dEdXMeanFrontTracker[0],
                          lambda event: event.kin_cal.reco_q3],
        "tags":reco_tags
    },
    "Psi vs q3":
    {
        "name" : "Psi_q3",
        "title" : "Psi v.s. q3; Psi; q3 (GeV); d^{2}NEvents/dPsidq3",
        "binning" : [PlotConfig.PSI_BINNING,
                     PlotConfig.BACKGROUND_FIT_Q3_BIN],
        "value_getter" : [lambda event: event.Psi,
                          lambda event: event.kin_cal.reco_q3],
        "tags":reco_tags
    },

    "Pi0 theta" :
    {
        "name": "Pi0_theta",
        "title" : "Pi0 Kinematics: theta, Pi0_theta(degree), reco_theta(degree)",
        "binning": [PlotConfig.ELECTRON_ANGLE_BINNING,
                    PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter" :[ lambda event: TruthTools.FindThetaFromIndex(event,TruthTools.MostEnergeticParticle(event,111)),
                          lambda event: event.kin_cal.reco_theta_e],
        "tags" : resolution_tags,
        "cuts" : [skip_sys, lambda event: event.classifier.truth_class == "CCPi0" or event.classifier.truth_class == "NCPi0"]
    },

    "Pi0 Energy" :
    {
        "name": "Pi0_energy",
        "title" : "Pi0 Kinematics: E, Pi0_E(degree), reco_E(degree)",
        "binning": [PlotConfig.ELECTRON_ENERGY_BINNING,
                    PlotConfig.ELECTRON_ENERGY_BINNING],
        "value_getter" :[ lambda event: event.mc_FSPartE[TruthTools.MostEnergeticParticle(event,111)]/1e3,
                          lambda event: event.kin_cal.reco_E_e],
        "tags" : resolution_tags,
        "cuts" : [skip_sys, lambda event: event.classifier.truth_class == "CCPi0" or event.classifier.truth_class == "NCPi0"] 
    },
    "Shower Width strip" :
    {
        "name": "ShowerWidth_energy",
        "title" : "Median Shower Width (strip) vs RecoEe; Width; reco_E (GeV); d^{2}NEvents/dStrip",
        "binning": [PlotConfig.SHOWER_WIDTH_BINNING,
                    PlotConfig.ELECTRON_ENERGY_BINNING],
        "value_getter" :[ lambda event: event.prong_MedianPlaneShowerWidth[0],
                          lambda event: event.kin_cal.reco_E_e],
        "tags" : resolution_tags,
    },
    "Inline Vertex Study" :
    {
        "name": "InlineVertexStudy",
        "title" : "Reco-True VertexZ vs InlineUpstreamEnergy; Reco VertexZ-True (mm); InlineUpstream E (GeV)",
        "binning": [PlotConfig.VERTEX_DIFF_BINNING,
                    PlotConfig.UPSTREAM_INLINE_ENERGY_BINS],
        "value_getter" :[ lambda event: TruthTools.InlineVertexStudyVtx(event),
                          lambda event: TruthTools.InlineVertexStudyE(event)],
        "tags" : resolution_tags,
    },

    "Lepton Energy vs q0":
    {
        "name": "Eel_q0",
        "title": "Lepton Energy vs q0; E_{e} (GeV); Eavailable (GeV); NEvents",
        "binning" : [PlotConfig.ELECTRON_ENERGY_BINNING,
                     PlotConfig.LOW_RECOIL_BIN_Q0],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep,
                          lambda event: event.kin_cal.reco_q0],
        "tags":reco_tags,
        "cuts" : [skip_sys, lambda event: event.classifier.truth_class == "NuEElastic"]
    },

    ################## 1D Migration Plots: ##################################

    "Q3 Migration":
    {
        "name" : "q3_migration",
        "title" : "q3 Migration; Reco q3 (GeV); True q3 (GeV)",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q3,
                     PlotConfig.LOW_RECOIL_BIN_Q3_Truth],
        #"binning" : [PlotConfig.BIN_Q3_OPTIMIZATION,
        #             PlotConfig.BIN_Q3_OPTIMIZATION],
        "value_getter" : [lambda event: event.kin_cal.reco_q3,lambda event: event.kin_cal.true_q3],
        "tags":migration_tags,
    },

    "Visible Energy Migration":
    {
        "name" : "Eavail_migration",

        "title" : "Available Energy Migration; Reco E_{avail} (GeV); True E_{avail} (GeV)",
        #"binning" : [PlotConfig.LOW_RECOIL_BIN_AVAIL,
        #            PlotConfig.LOW_RECOIL_BIN_AVAIL],
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0,
                     PlotConfig.LOW_RECOIL_BIN_Q0_Truth],
        "value_getter" : [lambda event: event.kin_cal.reco_visE,lambda event: event.kin_cal.true_visE],
        "tags":migration_tags,
    },

    "Q0 Migration":
    {
        "name" : "q0_migration",
        "title" : "q0 Migration; Reco q0 (GeV); True q0 (GeV)",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0,
                     PlotConfig.LOW_RECOIL_BIN_Q0_Truth],
        "value_getter" : [lambda event: event.kin_cal.reco_q0,lambda event: event.kin_cal.true_q0],
        "tags":migration_tags,
    },

     "Lepton Energy Migration":
    {
        "name" : "Eel_migration",
        "title" : "Lepton Energy Migration; Reco E (GeV); True E (GeV)",
        "binning" : [PlotConfig.ELECTRON_ENERGY_BINNING,
                     PlotConfig.ELECTRON_ENERGY_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep,lambda event: event.kin_cal.true_E_lep],
        "tags":migration_tags,
    },

      "Neutrino Energy Migration":
    {
        "name" : "Enu_migration",
        "title" : "Neutrino Energy Migration; Reco E (GeV); True E (GeV)",
        "binning" : [PlotConfig.NEUTRINO_ENERGY_BINNING,
                     PlotConfig.NEUTRINO_ENERGY_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_nu_cal,lambda event: event.kin_cal.true_enu_genie],
        "tags":migration_tags,
    },

     "Lepton Pt Migration":
    {
        "name" : "Lepton_Pt_migration",
        "title" : "Lepton Pt Migration; Reco bins; Truth bins",
        "binning" : [PlotConfig.PT_BINNING,
                     PlotConfig.PT_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad),lambda event: event.kin_cal.true_P_lep*math.sin(event.kin_cal.true_theta_lep_rad)],
        "tags":migration_tags,
    },

     "Lepton Energy Migration High Angle":
    {
        "name" : "Eel_migration_high_theta",
        "title" : "Lepton Energy Migration; Reco E (GeV); True E (GeV)",
        "binning" : [PlotConfig.ELECTRON_ENERGY_BINNING,
                     PlotConfig.ELECTRON_ENERGY_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_lep,lambda event: event.kin_cal.true_E_lep],
        "tags":migration_tags,
        "cuts": [lambda event:event.kin_cal_true_theta_e>20]
    },

     "Lepton Theta Migration":
    {
        "name" : "eltheta_migration",
        "title" : "Lepton Theta Migration; Reco #theta (degree); True #theta (degree)",
        "binning" : [PlotConfig.EXCESS_ANGLE_BINNING,
                     PlotConfig.EXCESS_ANGLE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_theta_e,lambda event: event.kin_cal.true_theta_e],
        "tags":reco_tags,
    },
     "Lepton Energy Theta":
    {
        "name" : "Ee_eltheta_migration",
        "title" : "Lepton Energy Theta; Reco Energy (GeV); Reco #theta (degree)",
        "binning" : [PlotConfig.EXCESS_ENERGY_FIT,
                     PlotConfig.EXCESS_ANGLE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.reco_E_e,lambda event: event.kin_cal.reco_theta_e],
        "tags":reco_tags,
    },
    "Extra Energy Difference":
    {
        "name" : "ExtraEnergyDiff",
        "title": "Pos Diff ; Longitudinal Difference; Transverse Difference",
        "binning" : [PlotConfig.LONG_BINNING,
                     PlotConfig.TRANS_BINNING],
        "value_getter" : [lambda event: TruthTools.Long(event), lambda event: TruthTools.Trans(event)],
        "tags": reco_tags,
    },
    "Extra Energy Long":
    {
        "name" : "ExtraEnergyLong",
        "title": "Long Extra Energy Clusters, Energy Weighted ; Longitudinal Difference (m); ClusterE",
        "binning" : [PlotConfig.LONG_BINNING],
        "value_getter" : [lambda event: TruthTools.Long(event)],
        #"weight_function": lambda event: TruthTools.GetWeight(event),
        "tags": reco_tags,
    },
    "Extra Energy Trans":
    {
        "name" : "ExtraEnergyTrans",
        "title": "Trans Distance Extra Energy Clusters ; Transverse Difference (m); ClusterE",
        "binning" : [PlotConfig.TRANS_BINNING],
        "value_getter" : [lambda event: TruthTools.Trans(event)],
        #"weight_function": lambda event: TruthTools.GetWeight(event),
        "tags": reco_tags,
    },
    "Extra Energy Long Negative":
    {
        "name" : "NegExtraEnergyLong",
        "title": "Long Extra Energy Clusters; Longitudinal Difference (m); ClusterE",
        "binning" : [PlotConfig.LONG_NEG_BINNING],
        "value_getter" : [lambda event: TruthTools.Long(event)],
        #"weight_function": lambda event: TruthTools.GetWeight(event),
        "tags": reco_tags,
    },
    "Extra Energy Long Positive":
    {
        "name" : "PosExtraEnergyLong",
        "title": "Long Extra Energy Clusters; Longitudinal Difference (m); ClusterE",
        "binning" : [PlotConfig.LONG_POS_BINNING],
        "value_getter" : [lambda event: TruthTools.Long(event)],
        #"weight_function": lambda event: TruthTools.GetWeight(event),
        "tags": reco_tags,
    },
    ################## Truth Signal Plots: ##################################
    "True Signal Lepton Theta":
    {
        "name" : "eltheta_true_signal",
        "title": "#theta; #theta; NEvents",
        "binning" : [PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_theta_lep],
        "tags": truth_signal_tags
    },

    "True Signal q0 vs q3":
    {
        "name" : "q0_q3_true_signal",
        "title": "q0 v.s. q3; q0 (GeV); q3 (GeV); NEvents",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0_Truth,
                     PlotConfig.LOW_RECOIL_BIN_Q3_Truth],
        "value_getter" : [lambda event: event.kin_cal.true_visE,
                          lambda event: event.kin_cal.true_q3],
        "tags": truth_signal_tags
    },

    "True Signal Visible Energy vs q3":
    {
        "name" : "Eavail_q3_true_signal",

        "title": "Available Energy v.s. q3; E_{avail} (GeV); q3 (GeV); NEvents",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0_Truth,
                     PlotConfig.LOW_RECOIL_BIN_Q3_Truth],
        "value_getter" : [lambda event: event.kin_cal.true_visE,
                          lambda event: event.kin_cal.true_q3],
        #"tags": truth_signal_tags
        "tags": truth_signal_tags
    },

     "True Lepton Pt":
    {
        "name" : "tLepton_Pt",
        "title" : "tPt_lepton",
        "binning" : [PlotConfig.PT_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_P_lep*math.sin(event.kin_cal.true_theta_lep_rad)],
        "tags":truth_tags
    },
     "True Signal Visible Energy vs Lepton Pt":
    {
        "name" : "Eavail_Lepton_Pt_true_signal",

        "title": "Available Energy v.s. q3; E_{avail} (GeV); Lepton Pt (GeV); NEvents",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0_Truth,
                     PlotConfig.PT_BINNING_Truth],
        "value_getter" : [lambda event: event.kin_cal.true_visE,
                          lambda event: event.kin_cal.true_P_lep*math.sin(event.kin_cal.true_theta_lep_rad)],
        "tags": truth_signal_tags   
    },
    "True Signal Lepton Energy vs q3":
    {
        "name" : "Ee_q3_true_signal",
        "title": "E_e v.s. q3; E_e (GeV); q3 (GeV); NEvents",
        "binning" : [PlotConfig.EXCESS_ENERGY_FIT,
                     PlotConfig.LOW_RECOIL_BIN_Q3],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep,
                          lambda event: event.kin_cal.true_q3],
        "tags": truth_signal_tags
    },
    "True Signal Lepton Theta vs q3":
    {
        "name" : "Theta_e_q3_true_signal",
        "title": "Theta_e v.s. q3; Theta_e (degree); q3 (GeV); NEvents",
        "binning" : [PlotConfig.EXCESS_ANGLE_BINNING,
                     PlotConfig.LOW_RECOIL_BIN_Q3],
        "value_getter" : [lambda event: event.kin_cal.true_theta_lep,
                          lambda event: event.kin_cal.true_q3],
        "tags": truth_signal_tags
    },

    "True Signal Z Vertex":
    {
        "name" : "zVertex_true_signal",
        "title": "zVertex; zVertex(mm); NEvents",
        "binning" : [PlotConfig.WHOLE_DET_Z_BINS],
        "value_getter" : [lambda event: event.mc_vtx[2]],
        "tags": truth_signal_tags
    },


    ################## 2D Migration Plots: ##################################

    "q0 vs q3 Migration":
    {
        "name" : "q0_q3_migration",
        "title" : "q0, q3 Migration; Reco bins; Truth bins",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0,
                     PlotConfig.LOW_RECOIL_BIN_Q3,
                     PlotConfig.LOW_RECOIL_BIN_Q0_Truth,
                     PlotConfig.LOW_RECOIL_BIN_Q3_Truth],
        "value_getter" : [lambda event: event.kin_cal.reco_q0,lambda event: event.kin_cal.reco_q3,
                          lambda event: event.kin_cal.true_visE,lambda event: event.kin_cal.true_q3],
        "tags":migration_tags,
    },

    "Visible Energy vs q3 Migration":
    {
        "name" : "Eavail_q3_migration",
        "title" : "Available Energy, q3 Migration; Reco bins; Truth bins",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0,
                     PlotConfig.LOW_RECOIL_BIN_Q3,
                     PlotConfig.LOW_RECOIL_BIN_Q0_Truth,
                     PlotConfig.LOW_RECOIL_BIN_Q3_Truth],
        "value_getter" : [lambda event: event.kin_cal.reco_visE,lambda event: event.kin_cal.reco_q3,
                          lambda event: event.kin_cal.true_visE,lambda event: event.kin_cal.true_q3],
        "tags":migration_tags,
    },

    "Visible Energy vs Lepton Pt Migration":
    {
        "name" : "Eavail_Lepton_Pt_migration",
        "title" : "Available Energy, Lepton Pt Migration; Reco bins; Truth bins",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q0,
                     PlotConfig.PT_BINNING,
                     PlotConfig.LOW_RECOIL_BIN_Q0_Truth,
                     PlotConfig.PT_BINNING_Truth],
        "value_getter" : [lambda event: event.kin_cal.reco_visE,lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad),
                          lambda event: event.kin_cal.true_visE,lambda event: event.kin_cal.true_P_lep*math.sin(event.kin_cal.true_theta_lep_rad)],
        "tags":migration_tags,
    },

    "Lepton Energy vs q3 Migration":
    {
        "name" : "Ee_q3_migration",
        "title" : "Lepton Energy, q3 Migration; Reco bins; Truth bins",
        "binning" : [PlotConfig.EXCESS_ENERGY_BINNING,
                     PlotConfig.LOW_RECOIL_BIN_Q3,
                     PlotConfig.EXCESS_ENERGY_BINNING,
                     PlotConfig.LOW_RECOIL_BIN_Q3],
        "value_getter" : [lambda event: event.kin_cal.reco_E_e,lambda event: event.kin_cal.reco_P_lep*math.sin(event.kin_cal.reco_theta_lep_rad),
                          lambda event: event.kin_cal.true_E_e,lambda event: event.kin_cal.true_P_lep*math.sin(event.kin_cal.true_theta_lep_rad)],
        "tags":migration_tags,
    },

    "Lepton Theta vs q3 Migration":
    {
        "name" : "Etheta_q3_migration",
        "title" : "Visible Energy, q3 Migration; Reco bins; Truth bins",
        "binning" : [PlotConfig.EXCESS_ANGLE_BINNING,
                     PlotConfig.LOW_RECOIL_BIN_Q3,
                     PlotConfig.EXCESS_ANGLE_BINNING,
                     PlotConfig.LOW_RECOIL_BIN_Q3],
        "value_getter" : [lambda event: event.kin_cal.reco_theta_e,lambda event: event.kin_cal.reco_q3,
                          lambda event: event.kin_cal.true_theta_e,lambda event: event.kin_cal.true_q3],
        "tags":migration_tags,
    },

    ################## 1D Resolution Plots: ############################

    "Vertex Resolution":
    {
        "name" : "zvertex_resolution",
        "title": "Z vertex Resolution; Log (#DeltaZ (mm)/10mm); NEvents",
        "binning": [PlotConfig.VERTEX_Z_RESOLUTION_BINNING],
        "value_getter": [lambda event: math.log10(max(abs(event.vtx[2]-event.mc_vtx[2]),1))],
        "tags": resolution_tags,
        "cuts": [lambda event: event.kin_cal.reco_q3<=1.2]
    },

    "Lepton Theta Resolution":
    {
        "name" : "eltheta_resolution",
        "title" : "Lepton Theta Resolution; Lepton Theta Resolution; NEvents",
        "binning" : [PlotConfig.RESOLUTION_BINNING],
        "value_getter" : [lambda event: (event.kin_cal.reco_theta_e-event.kin_cal.true_theta_e)/event.kin_cal.true_theta_e],
        "tags":resolution_tags,
        "cuts":[skip_sys]
    },

    "Lepton Theta Residual":
    {
        "name" : "del_eltheta",
        "title" : "True Lepton Theta; Lepton Theta Residual (degree); NEvents",
        "binning" : [PlotConfig.ELECTRON_ANGLE_BINNING, PlotConfig.ELECTRON_ANGLE_RESIDUAL_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_theta_e, lambda event: event.kin_cal.reco_theta_e - event.kin_cal.true_theta_e],
        "tags":resolution_tags,
        "cuts":[skip_sys]
    },

    "Lepton Energy Resolution":
    {
        "name" : "Eel_resolution",
        "title" : ";True Lepton Energy (GeV);Lepton Energy Resolution; Lepton Energy NEvents",
        "binning" : [PlotConfig.ELECTRON_ENERGY_BINNING,PlotConfig.RESOLUTION_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_E_lep, lambda event: event.kin_cal.reco_E_lep/event.kin_cal.true_E_lep -1],
        "tags":migration_tags,
        "cuts":[skip_sys]
    },

    "Neutrino Energy Resolution":
    {
        "name" : "Enu_resolution",
        "title" : "Neutrino Energy Resolution; Neutrino Energy Resolution; NEvents",
        "binning" : [PlotConfig.RESOLUTION_BINNING],
        "value_getter" : [lambda event: (event.kin_cal.reco_E_nu_cal-event.kin_cal.true_enu_genie)/event.kin_cal.true_enu_genie],
        "tags":resolution_tags,
        "cuts":[skip_sys]
    },

    "Q3 Resolution":
    {
        "name" : "q3_resolution",
        "title" : "q3 Resolution; q3(GeV), q3 Resolution; NEvents",
        "binning" : [PlotConfig.LOW_RECOIL_BIN_Q3,PlotConfig.RESOLUTION_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_q3,
                          lambda event: min((event.kin_cal.reco_q3-event.kin_cal.true_q3)/event.kin_cal.true_q3,PlotConfig.RESOLUTION_BINNING[-2])],
        "tags":migration_tags,
        "cuts":[skip_sys]
    },

    "Visible Energy Resolution":
    {
        "name" : "Eavail_resolution",
        "title" : "; True Available Energy; Available Energy Resolution",
        "binning" : [PlotConfig.VISIBLE_ENERGY_BINNING,PlotConfig.RESOLUTION_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_q0,
                          lambda event: PlotConfig.RESOLUTION_BINNING[-2] if event.kin_cal.true_visE<=0 else min((event.kin_cal.reco_visE-event.kin_cal.true_q0)/event.kin_cal.true_q0,PlotConfig.RESOLUTION_BINNING[-2])],
        "tags":migration_tags,
        "cuts":[skip_sys]
    },

    "Q0 Resolution":
    {
        "name" : "Q0_resolution",
        "title" : "Q0 Resolution; Q0; Q0 Resolution",
        "binning" : [PlotConfig.VISIBLE_ENERGY_BINNING,PlotConfig.RESOLUTION_BINNING],
        "value_getter" : [lambda event: event.kin_cal.true_q0,
                          lambda event: PlotConfig.RESOLUTION_BINNING[-2] if event.kin_cal.true_q0<=0 else min((event.kin_cal.reco_q0-event.kin_cal.true_q0)/event.kin_cal.true_q0,PlotConfig.RESOLUTION_BINNING[-2])],
        "tags":migration_tags,
        "cuts":[skip_sys]
    },

    "dEdX vs Theta":
    {
        "name" : "dEdX_vs_eltheta",
        "title": "dEdX v.s. Lepton Theta; dEdX;Theta",
        "binning" : [PlotConfig.DEDX_BINNING,
                     PlotConfig.ELECTRON_ANGLE_BINNING],
        "value_getter" : [lambda event: event.prong_dEdXMeanFrontTracker[0],
                          lambda event: event.kin_cal.reco_theta_e],
        "tags": reco_tags
    },

    "ProngPurity vs Available Energy Residual":
    {
        "name" : "prong_purity_vs_Eavail_residual",
        "title": "eProngPurity v.s. Visible Energy shift ;Prong Purity; Residual Visible Energy",
        "binning" : [PlotConfig.FRACTION_BINNING,
                     PlotConfig.RESOLUTION_BINNING],
        "value_getter" : [lambda event: event.prong_FracProngFromBestTrueTraj[0],
                          lambda event: (event.kin_cal.reco_visE-event.kin_cal.true_visE)/event.kin_cal.true_visE],
        "tags": reco_tags
    },
}

#like histfolio, connect related histograms together
class HistHolder:
    def __init__(self, name, f, sideband, is_mc, pot=1.0,data_pot=None):
        self.plot_name = TranslateSettings(name)["name"]
        self.dimension = len(TranslateSettings(name)["value_getter"])
        self.is_mc = is_mc
        self.hists= {}
        self.sideband=sideband
        self.pot_scale = data_pot/pot if (data_pot is not None and pot is not None) else pot
        self.valid = False
        self.bin_width_scaled = False
        self.POT_scaled =False
        if f is not None:
            self.valid = self.ReadHistograms(f)

    def ReadHistograms(self,f):
        variant_arg = [self.plot_name]
        if self.sideband != "Signal":
            variant_arg.append(self.sideband)
        namestring = VariantPlotsNamingScheme(*variant_arg)
        self.hists["Total"] = Utilities.GetHistogram(f, namestring)
        if self.is_mc:
            for cate in list(SignalDef.TRUTH_CATEGORIES.keys())+["Other"]:
                cate_namestring = VariantPlotsNamingScheme(*(variant_arg+[cate]))
                self.hists[cate]=Utilities.GetHistogram(f, cate_namestring)

        return self.hists["Total"] is not None


    def Scale(self, scale, bin_width_normalize):
        if not self.valid:
            return None
        if self.bin_width_scaled:
            raise ValueError("hist holder alreay scaled by bin_width")

        for k, v in self.hists.items():
            if v is not None:
                v.Scale(scale,"width" if bin_width_normalize else "")
        self.bin_width_scaled = bin_width_normalize

    def AreaScale(self,bin_width_normalize = False):
        if not self.valid:
            return None
        scale = 1.0/self.hists["Total"].Integral()
        self.Scale(scale,bin_width_normalize)


    def POTScale(self,bin_width_normalize = True):
        if not self.valid:
            return None

        if self.POT_scaled and self.pot_scale != 1.0:
            raise ValueError("hist holder alreay scaled by POT")

        for k, v in self.hists.items():
            if v is not None:
                v.Scale(self.pot_scale if self.is_mc else 1.0,"width" if bin_width_normalize else "")

        self.POT_scaled = True

    def GetCateList(self,grouping = None):
        if not self.valid:
            return None
        _mc_ints = []
        _colors= []
        _titles= []
        local_grouping = grouping
        for cate in list(local_grouping.keys())[::-1]:
            config = local_grouping[cate]
            hist = None
            if "cate" in config:
                for fine_cate in config["cate"]:
                    if not self.hists[fine_cate]:
                        continue
                    if hist is None:
                        hist = self.hists[fine_cate].Clone()
                    else:
                        hist.Add(self.hists[fine_cate])
            else:
                hist = self.hists[cate]

            if hist:
                hist.SetTitle(config["title"])
                _mc_ints.append(hist)
                _titles.append(config["title"])
                _colors.append(config["color"])
        return _mc_ints,_colors,_titles

    def GetHist(self):
        if not self.valid:
            return None
        return self.hists["Total"]

    def ResumTotal(self):
        self.hists["Total"].Reset("ICESM")
        for i in self.hists:
            if i=="Total":
                continue
            else:
                self.hists["Total"].Add(self.hists[i])

    def GetSumCate(self,cate_list):
        if len(cate_list) == 0:
            return None
        hnew = self.hists[cate_list[0]].Clone()
        for i in range(1,len(cate_list)):
            htmp = self.hists[cate_list[i]]
            if htmp:
                hnew.Add(htmp)
            else:
                continue
        return hnew

    def Add(self,hist_holder):
        if not isinstance(hist_holder,HistHolder):
            raise ValueError("Can only add histholder to histholder")
        for i in self.hists:
            if self.hists[i] and hist_holder.hists[i]:
                self.hists[i].Add(hist_holder.hists[i])
            else:
                print(i,self.hists[i],hist_holder.hists[i])
