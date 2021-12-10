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
Categories["CCQE"]={
    "title" : "QE",
    "color": COLORS[0]
}

Categories["CCDelta"]={
    "title" : "Delta",
    "color": COLORS[1]
}

Categories["CC2p2h"]={
    "title": "2p2h",
    "color": COLORS[3]
}
Categories["CCDIS"]= {
    "title": "DIS",
    "color" : COLORS[2]
}
Categories["CCOther"]= {
    "title": "CC-Other",
    "color" : COLORS[4]
}

# Categories["NuEElastic"] = {
#     "title":"#nu + e elastic",
#     "color" : COLORS[6]
# }
Categories["Other"] = {
    "title":"Others",
    "cate" : {"Other","NonFiducial","NonPhaseSpace","CCNuE","NC"},
    "color" : COLORS[6]
}


SignalDecomposition = {
    "CCQE" :
    {
        "title" : "CC #nu_{mu}-QE",
        "color": COLORS[0]
    },
    "CCDelta" : {
        "title" : "CC #nu_{mu}-Delta",
        "color": COLORS[1]
    },
    "CCDIS" : {
        "title" : "CC #nu_{mu}-DIS",
        "color": COLORS[2]
    },
    "CC2p2h" : {
        "title" : "CC #nu_{mu}-2p2h",
        "color": COLORS[3]
    },
    "Background" : {
        "title" : "Backgrounds",
        "cate": {"Other","NonFiducial","NonPhaseSpace","CCNuE","NC"},
        "color": COLORS[4]
    }
    # "CCNuEAntiNu" : {
    #     "title" : "CC #bar{#nu_{e}}",
    #     "color": COLORS[4]
    # }
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
              "args": (SignalDecomposition,True,False)},
    "sigdepRatio":{"func":PlotTools.PrepareSignalDecomposeRatio,
                   "args": (SignalDecomposition,True,False)},

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
    #  {"name":"Visible Energy vs q3",
    #  "plot_type":"sigdepRatio"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type":"sigdep"},

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
    # },
    # {"name":"Visible Energy Migration",
    #  "plot_type":"comp",
    #  "slicer":lambda hist: PlotTools.Make2DSlice(hist,True)},
    #  {"name":"Q0 Migration",
    #  "plot_type":"comp",
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
    # {"variables":["Visible Zoomin","Lepton Energy"],
    #  "slicer": lambda hist: PlotTools.Make2DSlice(hist,interval=5),
    #  "args":(NuEElasticCategory,),
    #  "scale":lambda histHolder:histHolder.POTScale(False)},
    # {"name":"FrontdEdX_position" },
    # {"name":"Front dEdX" },
    # {"name":"Inline Upstream Energy Weighted Position",
    #  "canvasconfig":PlotTools.Logy},
    # {"name":"Inline Upstream Energy",
    #  "canvasconfig":PlotTools.Logx,
    #  "scale" : lambda histHolder:histHolder.POTScale(False)},
    # {"variables":["Visible Energy","Lepton Pt"]},
    # {"variables":["Visible Energy","Lepton Pt"],
    #  "plot_type":"comp"},
    # {"variables":["Visible Energy","Lepton Pt"],
    #  "plot_type":"sigdepRatio"},
    #  {"variables":["Visible Energy","Lepton Pt"],
    #  "plot_type":"err"},
    {"name":"Neutrino Energy Migration",
     "plot_type":"migration"},
    {"variables":["True Neutrino Energy","Neutrino Energy"],
     "plot_type":"migration"}
]

