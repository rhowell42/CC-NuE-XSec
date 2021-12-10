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
COLORS=ROOT.MnvColors.GetColors()
Categories = OrderedDict()
Categories["ExcessModel"] = {
    "title":"Ad hoc excess model",
    "color" : COLORS[2]
}
Categories["CCNuE"]={
    "title" : "CC #nu_{e}",
    "cate":{"CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE"},
    "color": COLORS[0]
}

Categories["2p2h"]={
    "title": "CC #nu_{e} 2p2h",
    "cate" :{"CCNuE2p2h"},
    "color": COLORS[1]
}
Categories["CCDIS"]= {
    "title": "CC #nu_{#mu} DIS",
    "color" : COLORS[3]
}
Categories["NCCOH"]= {
    "title": "NC coh",
    "color" : COLORS[4]
}
Categories["NCRES"]= {
    "title": "NC Res",
    "color" : COLORS[5]
}

Categories["NCDIS"]= {
    "title": "NC DIS",
    "color" : COLORS[6]
}

# Categories["NuEElastic"] = {
#     "title":"#nu + e elastic",
#     "color" : COLORS[6]
# }
Categories["Other"] = {
    "title":"Others",
    "cate" : {"Other","NonFiducial","NonPhaseSpace","NuEElastic","CCNuEAntiNu","NCOther","CCOther"},
    "color" : COLORS[7]
}


SignalDecomposition = {
    "CCNuEQE" :
    {
        "title" : "CC #nu_{e}-QE",
        "color": COLORS[0]
    },
    "CCNuEDelta" : {
        "title" : "CC #nu_{e}-Delta",
        "color": COLORS[1]
    },
    "CCNuEDIS" : {
        "title" : "CC #nu_{e}-DIS",
        "color": COLORS[2]
    },
    "CCNuE2p2h" : {
        "title" : "CC #nu_{e}-2p2h",
        "color": COLORS[3]
    },
    "Background" : {
        "title" : "Backgrounds",
        "cate": {"ExcessModel","NonFiducial","CCDIS","CCOther","NCCOH","NCRES","NCDIS","NCOther","NuEElastic","Other","CCNuEAntiNu"},
        "color": COLORS[4]
    }
    # "CCNuEAntiNu" : {
    #     "title" : "CC #bar{#nu_{e}}",
    #     "color": COLORS[4]
    # }
}

SignalBackground = {
    "Signal" :
    {
        "title" : "CC #nu_{e}",
        "cate": {
            "CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"
        },
        "color": COLORS[0]
    },
    "Background" : {
        "title" : "Backgrounds",
        "cate": {"ExcessModel","NonFiducial","CCDIS","CCOther","NCCOH","NCRES","NCDIS","NCOther","NuEElastic","Other","CCNuEAntiNu"},
        "color": COLORS[4]
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
        "cate": {"ExcessModel","NonFiducial","CCDIS","CCOther","NCRES","NCDIS","NCOther","Other","CCNuEAntiNu"},
        "color": COLORS[4]
    }
}


DefaultPlotters={
    "comp":{"func": PlotTools.PrepareComp},
    "ratio":{"func": PlotTools.PrepareRatio},
    "err": {"func": PlotTools.PrepareErr,
            "args":(True,False,CONSOLIDATED_ERROR_GROUPS)},
    "stacked":{"func":PlotTools.PrepareStack,
               "args": (Categories,)},
    "diff":{"func":PlotTools.PrepareDiff},
    "migration":{"func":PlotTools.PrepareMigration},
    "sigdep":{"func":PlotTools.PrepareSignalDecompose,
              "args": (SignalDecomposition,True,False)}

}

PLOTS_TO_MAKE = [
    # {"name":"Lepton Energy",
    #  "plot_type":"stacked"},
    # {"name":"Neutrino Energy",
    #  "plot_type":"stacked"},
    # {"name":"Neutrino Energy",
    #   "sideband_group": ("ex_and_pi0",("Excess","Pi0"))},
    # {"name":"Visible Energy",
    #   "sideband_group": ("ex_and_pi0",("Excess","Pi0"))},
    # {"name": "Lepton Theta"},
    # {"name": "Lepton Theta"},
    # {"name": "Psi"},
    # {"name":"Q0"},
    # {"name":"Q3"},
    # {"name":"Visible Energy",
    #   "plot_type" : "comp",},
    # {"name":"Q3",
    #  "plot_type" : "sigdep",},
    # {"name":"Visible Energy",
    #  "plot_type" : "err",},
    # {"variables":["Leading Pi0 E","PsiEe"],
    #  "slicer":lambda hist: PlotTools.Make2DSlice(hist,True,interval=2),
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
    #  },
    # {"variables":["PsiEe","Lepton Energy"],
    #  "slicer": lambda hist: PlotTools.Make2DSlice(hist,interval=5),
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
    #  },
    #  {"variables":["PsiEe","Lepton Energy"],
    #  "slicer": lambda hist: PlotTools.Make2DSlice(hist,interval=5),
    #  },
    #  {"name":"Lepton Theta",
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
    #  },
    
    # {"name":"Q0",
    #  "plot_type" : "err",},
    # {"name":"Q3",
    #  "plot_type" : "err",},
    # {"name":"Q3 Migration",
    #   "plot_type":"migration",
    #  "slicer": lambda x: [x.Clone()]},
    # {"name":"Q0 Migration",
    #  "plot_type":"migration",
    #  "slicer": lambda x: [x.Clone()]},
     # {"name":"Visible Energy Migration",
     #  "plot_type":"migration",
     #  "slicer": lambda x: [x.Clone()]},
    # {"name":"Visible Energy vs q3"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type":"comp"},
    #  {"name":"Visible Energy vs q3",
    #  "plot_type":"err"},
    # {"variables":["Q0","True Visible Energy"],
    #  "scale":lambda hist:hist.AreaScale(True),
    #  "slicer":lambda hist: PlotTools.Make2DSlice(hist,True),
    #  "plot_type":"comp",
    #  "args":(True,False,True),
    #  },
    #  {"name":"Visible Energy Migration",
    #   "slicer":lambda hist: PlotTools.Make2DSlice(hist,True),
    #   "scale":lambda hist:hist.AreaScale(True),
    #   "args":(True,False,True),
    #   "plot_type":"comp",
    #   },
    #  {"name":"Q0 Migration",
    #   "slicer":lambda hist: PlotTools.Make2DSlice(hist,True),
    #   "scale":lambda hist:hist.AreaScale(True),
    #   "args":(True,False,True),
    #   "plot_type":"comp",
    #   },
    # {"variables":["Visible Energy","True Q0"],
    #  "scale":lambda hist:hist.AreaScale(True),
    #  "slicer":lambda hist: PlotTools.Make2DSlice(hist,True),
    #  "args":(True,False,True),
    #  "plot_type":"comp",
    #  },
    # {"name":"Visible Energy Migration",
    #  "plot_type":"comp",
    #  "slicer":lambda hist: PlotTools.Make2DSlice(hist,True)},
     {"name":"Neutrino Energy Migration",
      "plot_type":"migration"},
    #  "slicer":lambda hist: PlotTools.Make2DSlice(hist,True)},
    # {"name":"True Signal Visible Energy vs q3"
    #  "plot_type":"err"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type" : "comp"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type" : "ratio"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type" : "err"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type" : "sigdep",},
    # {"name":"q0 vs q3",
    #  "plot_type" : "sigdep",},
    # {"name":"Visible Energy"},
    {"variables":["Visible Zoomin","Lepton Energy"],
     "slicer": lambda hist: PlotTools.Make2DSlice(hist,interval=5),
     "args":(NuEElasticCategory,),
     "scale":lambda histHolder:histHolder.POTScale(False)},
    # {"name":"FrontdEdX_position" },
    # {"name":"Front dEdX" },
    # {"name":"Inline Upstream Energy Weighted Position",
    #  "canvasconfig":PlotTools.Logy},
    # {"name":"Inline Upstream Energy",
    #  "canvasconfig":PlotTools.Logx,
    #  "scale" : lambda histHolder:histHolder.POTScale(False)},
    {"variables":["Visible Energy","Lepton Pt"]},
    {"variables":["Visible Energy","Lepton Pt"],
     "plot_type":"comp"},
     {"variables":["Visible Energy","Lepton Pt"],
     "plot_type":"err"},
    # {"variables":["Visible Energy","Lepton Pt"],
    #  "slicer": PlotTools.IdentitySlicer,
    #  "scale":lambda histHolder:histHolder.POTScale(False),
    #  "plot_type":"diff"},
    # {"variables":["Visible Energy","Lepton Pt"],
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
    #  "slicer": PlotTools.IdentitySlicer,
    #  "scale":lambda histHolder:histHolder.POTScale(False),
    #  "plot_type":"diff"},
    {"variables":["Visible Energy","Lepton Pt"],
     "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline"))},
    {"name":"Visible Energy vs q3",
     "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline"))},
    #  {"variables":["Lepton Energy"],
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline"))},
    # {"variables":["Sum Visible Energy","Lepton Pt"],
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
    #  "slicer":PlotTools.ProjectionX},
    # {"variables":["Sum Visible Energy","Lepton Pt"],
    #  "slicer":PlotTools.ProjectionX},
    {"variables":["Visible Energy","Lepton Pt"],
     "plot_type":"err"},
    {"variables":["True Neutrino Energy","Neutrino Energy"],
     "plot_type":"migration"}
]

