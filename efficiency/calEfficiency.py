import os
import ROOT
import PlotUtils
from itertools import chain
#from config.PlotConfig import HISTS_TO_MAKE_
from config.AnalysisConfig import AnalysisConfig
from tools.EventClassification import EventClassifier
from tools.KinematicsCalculator import KinematicsCalculator
from tools import Utilities
from tools.SystematicsUniverse import GetAllSystematicsUniverses
from tools.MyHistograms import MakePlotProcessors

HISTS_TO_MAKE = [
    #"True Signal q0 vs q3",
    #"True Signal Visible Energy vs q3",
    "Biased Neutrino Energy",
    #"True Visible Energy",
    #"True Signal Electron Theta"
]

def MakeTruthHistogram(chainwrapper,outname):
    kin_cal = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False, calc_true = True, calc_reco = False)
    eventClassifier = EventClassifier(classifiers=["Truth"],use_kin_cuts=True, use_sideband=[])
    universes = GetAllSystematicsUniverses(chainwrapper, False)
    for univ in chain.from_iterable(iter(universes.values())):
        univ.LoadTools(kin_cal,eventClassifier)

    nEvents = chainwrapper.GetEntries()
    output_file = ROOT.TFile.Open(outname,"RECREATE")
    Plots = preparePlots(universes)
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

                if eventClassifier.is_true_signal:
                    for entry in Plots:
                        entry.Process(universe)

    output_file.cd()
    for entry in Plots:
        entry.Finalize()

    print("counter for reco,truth: ", eventClassifier.counter)
    output_file.Close()


def preparePlots(universes):
    # make a bunch of Plot Processor, grouped by signal/sideband
    plots=[]
    for entry in HISTS_TO_MAKE:
        settings = {"key":entry,"region":"Signal","mc":True}
        plots.extend(MakePlotProcessors(**settings))

    for entry in plots:
        entry.AddErrorBands(universes)

    return plots



if __name__ =="__main__":
    outputTruthHistogram = AnalysisConfig.TruthHistoPath(AnalysisConfig.playlist,True)
    Tchain = Utilities.fileChain(AnalysisConfig.playlist,"mc",AnalysisConfig.ntuple_tag,"Truth",AnalysisConfig.count[0],AnalysisConfig.count[1])
    MakeTruthHistogram(Tchain,outputTruthHistogram)
