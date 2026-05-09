import ROOT
from collections import OrderedDict
from config import PlotConfig
from tools.PlotLibrary import HistHolder
from tools import PlotTools
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS



Default_Plot_Type="stacked"
Default_Scale = lambda histHolder:histHolder.POTScale(True)
DefaultSlicer = PlotTools.PrepareSlicer


# the order and color of each band in stacked plots defined here.
COLORS=ROOT.MnvColors.GetColors(ROOT.MnvColors.k36Palette)
Categories = OrderedDict()

Categories["CCNuE"]= {
    "title": "CC Anti#nu_{e}",
    "cate":{"CCNu","CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"},
    "color" : COLORS[0]
}

#Categories["CCQElike"]={
#    "title" : "CCQElike #bar#nu_{e}",
#   #"cate":{"CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE","NonPhaseSpace"},
#    "color": COLORS[0]
#}
#Categories["notCCQElike"]= {
#    "title": "notCCQElike #bar{#nu}_{e}",
#    #"cate" : {"CCPi0","CCOther"},
#    "color" : COLORS[1]
#}
Categories["NCDiff"] = {
     "title":"NC Diffractive #pi^{0}",
     "color" : COLORS[2]
}
Categories["CCPi0"]= {
    "title": "Not #nu_{e} CC #pi^{0}",
    "color" : COLORS[3]
}
Categories["NCCohPi0"]= {
    "title": "NC Coh #pi^{0}",
    "color" : COLORS[4]
}
Categories["NCPi0"]= {
    "title": "NC #pi^{0}",
    "color" : COLORS[5]
}
Categories["NCPi"]= {
    "title": "NC #pi+-",
    "color" : COLORS[8]
}
Categories["CCPi"]= {
    "title": "CC #pi+-",
    "color" : COLORS[9]
}
Categories["NuEElastic"] = {
    "title":"#nu + e elastic",
    "color" : COLORS[6]
}
Categories["NCOther"] = {
    "title":"NCOthers",
    "color" : COLORS[7]
}
Categories["CCOther"] = {
    "title":"CCOthers",
    "color" : COLORS[10]
}

BackgroundDecomposition = {
    "NCDiff" :
    {
        "title" : "NC Diffractive",
        "color": COLORS[0]
    },
    "NuEElastic" : {
        "title" : "#nu + e",
        "color": COLORS[1]
    },
    "CCPi0" : {
        "title" : "CC #pi^{0}",
        "color": COLORS[3]
    },
    "NCPi0" : {
        "title" : "NC #pi^{0}",
        "color": COLORS[5]
    },
    "CCPi" : {
        "title" : "CC #pi",
        "color": COLORS[6]
    },
    "NCPi" : {
        "title" : "NC #pi",
        "color": COLORS[7]
    },
    "NCCohPi0": {
        "title" : "NC Coh #pi^{0}", 
        "color": COLORS[8]
    },
    "NCOther" : {
        "title" : "NCOther",
        "color": COLORS[9]
    },
    "CCOther" : {
        "title" : "CCOther",
        "color": COLORS[10]
    }
    }

BackgroundDecomposition = Categories.copy()
del BackgroundDecomposition["CCNuE"]

#Categories = OrderedDict()
#Categories["NCnumu"] = {
#     "title":"NC #nu_{#mu}",
#     "color" : COLORS[2]
#}
#Categories["NCnue"] = {
#     "title":"NC #nu_{e}",
#     "color" : COLORS[3]
#}
#Categories["CCnue"] = {
#     "title":"CC #nu_{e}",
#     "color" : COLORS[4]
#}
#Categories["CCantinue"] = {
#     "title":"CC anti #nu_{e}",
#     "color" : COLORS[5]
#}
#Categories["NueOther"] = {
#    "title":"NueOthers",
#    "color" : COLORS[0]
#}

#Categories["ExcessModel"] = {
#    "title":"Ad hoc excess model",
#    "color" : COLORS[2]
#}

SignalDecomposition = {
    "CCQE" :
    {
        "title" : "CC #nu_{#mu}-QE",
        "color": COLORS[5]
    },
    "CC2p2h" : {
        "title" : "CC #nu_{#mu}-2p2h",
        "color": COLORS[6]
    },
    "CCDelta" : {
        "title" : "CC #nu_{#mu}-Delta Res",
        "color": COLORS[1]
    },
    "CCOther": {
        "title" : "CC #nu_{#mu}-Res", 
        "color": COLORS[4]
    },
    "CCDIS" : {
        "title" : "CC #nu_{#mu}-DIS",
        "color": COLORS[2]
    },
    "CCWrongSign": {
        "title" : "CC #nu_{#mu} wrong sign", 
        "color": COLORS[7]
    }
    }

#SignalDecomposition = {
#    "NueBar" :
#    {
#        "title" : "Anti#nu_{e}", 
#        "color": COLORS[0]
#    },
#    "Nue" :
#    {
#        "title" : "#nu_{e}", 
#        "color": COLORS[1]
#    },
#    "Background" : {
#        "title" : "Backgrounds",
#        "cate": {"Other"},
#        "color": COLORS[2]
#    }
#}
#}


#SignalDecomposition = {
#    "CCNuE" :
#    {
#        "title" : "Signal",
#        "color": COLORS[0]
#    },
#    "Background" : {
#        "title" : "Backgrounds",
#        "cate": {"NCDiff", "NuEElastic", "NonPhaseSpace", "NonFiducial", "CCPi0", "NCCohPi0", "NCPi0", "Other"},
#        "color": COLORS[2]
#    }
#}

SignalBackground = {
    "Signal" :
    {
        "title" : "CC #nu_{e}",
        "cate":{"CCNu","CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"},
        "color": COLORS[0]
    },
    "Background" : {
        "title" : "Backgrounds",
        #"cate": {"CCPi0","Other","NCDiff","NuEElastic","NonPhaseSpace","NonFiducial","NCCohPi0","NCPi0","CCNu"},
        #"cate": {"CCPi0","Other","NCDiff","NuEElastic","NCCohPi0","NCPi0","CCNu"},
        "cate": {"CCPi0","CCPi","NCPi","NCOther","CCOther","NCDiff","NuEElastic","NCCohPi0","NCPi0"},
        "color": COLORS[4]
    }
}

SignalChargedBackground = {
    "Signal" :
    {
        "title" : "CC #nu_{e}",
        "cate":{"CCNu","CCNuEQE","CCNuEDelta","CCNuE2p2h","CCNuEDIS","CCNuE"},
        "color": COLORS[0]
    },
    "CCBackground" : {
        "title" : "Charged Current Backgrounds",
        "cate": {"CCPi0","CCPi0Proton","CCOther","CCPi"},
        "color": COLORS[1]
    },
    "NCBackground" : {
        "title" : "Neutral Current Backgrounds",
        "cate": {"Other","NCOther","NCDiff","NCCohPi0","NCPi0","NCPi0Proton","NCPi"},
        "color": COLORS[2]
    },
    "NuEEl" :{
        "title" : "#nu + e elastic",
        "cate": {"NuEElastic"},
        "color": COLORS[3]
    }
}
SignalChargedBackground = {
    "Signal" :
    {
        "title" : "CC #nu_{e}",
        "cate":{"CCQE","CCDelta","CC2p2h","CCDIS","CCOther","CCWrongSign"},
        "color": COLORS[0]
    }
}

ChargedBackground = {
    "CCBackground" : {
        "title" : "Charged Current Backgrounds",
        "cate": {"CCPi0","CCPi0Proton","CCOther","CCPi"},
        #"cate": {"CCPi0","CCPi"},
        "color": COLORS[1]
    },
    "NCBackground" : {
        "title" : "Neutral Current Backgrounds",
        "cate": {"Other","NCOther","NCDiff","NCCohPi0","NCPi0","NCPi0Proton","NCPi"},
        "color": COLORS[2]
    },
    "NuEEl" :{
        "title" : "#nu + e elastic",
        "cate": {"NuEElastic"},
        "color": COLORS[3]
    }
}

NuEElasticCategory = {
    "Nu+e" :{
        "title": "nu+e",
        "cate" :{"NuEElastic"},
        "color":COLORS[1]
    },
    "CCNuE":
    {
        "title" : "CC #nu_{e}",
        "cate": {
            "CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"
        },
        "color": COLORS[0]
    },
    "NCCOH":{
        "title" : "NCCoh",
        "color": COLORS[3]
    },
    "ExcessModel": {
        "title": "NC Diffractive",
        "color":COLORS[2]
    },
    "Other":{
        "title" : "Others",
        #"cate": {"ExcessModel","NonFiducial","CCDIS","CCOther","NCRES","NCDIS","NCOther","Other","CCNuEAntiNu"},
        "cate": {"ExcessModel","CCDIS","CCOther","NCRES","NCDIS","NCOther","Other","CCNuEAntiNu"},
        "color": COLORS[4]
    }
}


MODELS = {
    "MnvTune v1": {
        "errorband":(None,None),
        "color":COLORS[0]
    },
    "2p2h Tune (QE)": {
        "errorband":("Low_Recoil_2p2h_Tune",2),
        "color":COLORS[1]
    },
    "SuSA 2p2h" : {
        "errorband":("SuSA_Valencia_Weight",0),
        "color":COLORS[2]
    },
    "MK model": {
        "errorband": ("MK_model",0),
        "color":COLORS[7]
    },
    "Low Q2 Pion Joint": {
        "errorband" : ("LowQ2Pi_Joint",0),
        "color": COLORS[5]
    },
    "Low Q2 Pion NuPi0": {
        "errorband" : ("LowQ2Pi_NUPI0",0),
        "color": COLORS[6]
    }
}

DefaultPlotters={
    "comp":{"func": PlotTools.PrepareComp},
    "ratio":{"func": PlotTools.PrepareRatio},
    "bkgratio":{"func": PlotTools.PrepareBkgRatio},
    "err": {"func": PlotTools.PrepareErr,
            "args":(SignalDecomposition,True,False,CONSOLIDATED_ERROR_GROUPS)},
    "stacked":{"func":PlotTools.PrepareStack,
        "args": (SignalDecomposition,)},
    "diff":{"func":PlotTools.PrepareDiff},
    "migration":{"func":PlotTools.PrepareMigration,
        "args": (SignalDecomposition,)},
    "hist2d":{"func":PlotTools.PrepareHist2D,
        "args": (SignalDecomposition,)},
    "category_hist":{"func":PlotTools.CategoryHist,
                     "args": (SignalChargedBackground,)},
    "sigdep":{"func":PlotTools.PrepareSignalDecompose,
              "args": (ChargedBackground,True,False)},
    "sigdepratio":{"func":PlotTools.PrepareSignalDecomposeRatio,
              "args": (SignalDecomposition,True)},
    "errband":{"func":PlotTools.PrepareErrorBand},
    "model":{},
    "model_ratio" :{},
    #"sigdep":{},
}

PROB_PLOTS = [
    {"name":"True Energy vs Neutrino Length Travelled"},
    {"name":["Biased Neutrino Energy"]},
    {"name":"True Energy Inverse vs Biased Neutrino Energy",
      "plot_type" : "migration"},
]

PLOTS_TO_MAKE = [
    ###### ME STERILE PAPER PLOTS ########
    {"name":"True Energy vs Biased Neutrino Energy",
            "plot_type" : 'category_hist'},
    {"name":"True Energy vs L/E",
            "plot_type" : 'category_hist'},

    #{"name":"True Energy vs Biased Neutrino Energy",
    #        "plot_type" : "category_hist"},
    #{"name":"True Energy vs Biased Neutrino Energy",
    #        "plot_type" : 'migration'},
    #{"name":"Reco Energy vs L/E",
    #    "plot_type" : "migration"},
    #{"name":"True Energy vs Biased Neutrino Energy",
    #     "plot_type" : "category_hist"},
    #{"name":"True Energy vs Biased Neutrino Energy",
    #     "plot_type" : "migration"},
    #{"name":"Reco Energy vs L/E",
    #    "plot_type" : "category_hist"},
    
    #{"name":"Biased Neutrino Energy"},
    #{"name":"Biased Neutrino Energy",
    #     "plot_type" : "err"},
    ]
