"""
   Macros for studying efficiency, purity, etc of a cut variable
"""

import os
import ROOT
import string
#import cProfile
import signal
import sys
from itertools import chain

#start loading my modules
from tools import Utilities
from config.AnalysisConfig import AnalysisConfig
from tools.SystematicsUniverse import GetAllSystematicsUniverses
from tools.EventClassification import EventClassifier
from tools.KinematicsCalculator import KinematicsCalculator
from tools.MyHistograms import MakeSelectionPlotsByCut,MakeCutVariablePlot
#from selection.

ROOT.TH1.AddDirectory(False)


PLOTS_TO_MAKE = [
    #{"key":"True Electron Theta"},
    #{"key":"Electron Energy"},
    {"key":"True Visible Energy"},
    {"key":"Visible Energy"},
    # {"key":"Q3"},
    # {"key":"True q3"},
    # {"key": {
    #     "variables":["W"],
    # }},
]

NMinus1CutPlot = [
    #("Wexp")
    ("MeanFrontdEdX",),
    #("LLR",),("RecoQ3",),("RecoElectronEnergy",)
]

def plotRecoKin(mc, chainwrapper, outname):
    """ The main code of event selection """
    kin_cal = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False, calc_true = mc)
    eventClassifier = EventClassifier(classifiers=["Reco","Truth"] if mc else ["Reco"], use_kin_cuts=True, use_sideband=[])
    universes = GetAllSystematicsUniverses(chainwrapper, not mc, False, "all")
    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal,eventClassifier)

    output_file = ROOT.TFile.Open(outname,"RECREATE")
    Plots = preparePlots(universes,mc)
    nEvents = chainwrapper.GetEntries()
    if AnalysisConfig.testing and nEvents > 5000:
        nEvents = 5000

    for counter in range(nEvents):
        #half an hour for a event, should be much more than needed unless stuck in I/O
        if counter %100000 == 0:
            print(counter)

        for universe in chain.from_iterable(iter(universes.values())):
            universe.SetEntry(counter)

            #only update kin_cal & eventClassifier when universe in not vertical only.
            if not universe.IsVerticalOnly():
                kin_cal.CalculateKinematics(universe)
                eventClassifier.Classify(universe)

            for entry in Plots:
                entry.Process(universe)

    output_file.cd()
    for entry in Plots:
        entry.Finalize()

    print("counter for reco,truth: ", eventClassifier.counter)
    output_file.Close()


def preparePlots(universes,mc):
    # make a bunch of Plot Processor, grouped by signal/sideband
    plots=set([])
    #for region in ["Signal"]+AnalysisConfig.sidebands:

    for settings in PLOTS_TO_MAKE:
        plots.update(MakeSelectionPlotsByCut(**settings))

    for cut in NMinus1CutPlot:
        plots.update(MakeCutVariablePlot(*cut))

    # add the errorband map to plot processor.
    for plot in plots:
        plot.AddErrorBands(universes)

    return plots



if __name__ == "__main__":

    Reco = AnalysisConfig.run_reco
    POT_cal = AnalysisConfig.POT_cal
    print("playlist %s running ---------" % AnalysisConfig.playlist)
    for st in AnalysisConfig.data_types:
        outputSelectionHistogram = AnalysisConfig.CutStudyPath(AnalysisConfig.playlist,"data" in st,True)
        print(outputSelectionHistogram)
        #if AnalysisConfig.count[1] is None:
        #    used_POT,total_POT = Utilities.getPOT(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,POT_cal)
        #    print("POT for given playlist: " + AnalysisConfig.playlist + ", and data_type: " + st + ", is "+str(used_POT))

        if Reco :
            #cProfile.run('plotRecoKin(st=="mc", Utilities.fileChain(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,"NuECCQE",AnalysisConfig.count[0],AnalysisConfig.count[1]), outputSelectionHistogram)')
            plotRecoKin(st=="mc", Utilities.fileChain(AnalysisConfig.playlist,st,AnalysisConfig.ntuple_tag,None,AnalysisConfig.count[0],AnalysisConfig.count[1]), outputSelectionHistogram)

        print("selection is done for ", st, AnalysisConfig.playlist)

