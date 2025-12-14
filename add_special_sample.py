import os
import sys
import re
import multiprocessing as mpi
import subprocess
import ROOT
import PlotUtils
from itertools import chain
from config.AnalysisConfig import AnalysisConfig
from config.SignalDef import SIGNAL_DEFINATION
from tools import Utilities

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
SelectionFilesRegex = "(kin|truth)_dist_(data|mc)(.+)_"+ AnalysisConfig.selection_tag+"_"+AnalysisConfig.ntuple_tag+"(_[0-9]+)?\.root"

def AddOneFile(selection_string,special_string,bigGenie=False):
    print("Opening {} to get special sample hists".format(special_string))
    print("Opening {} to set signal histograms to special sample hists".format(selection_string))
    special_file = ROOT.TFile.Open(special_string)
    selection_file = ROOT.TFile.Open(selection_string,"UPDATE")
    keylist = special_file.GetListOfKeys()
    for key in keylist:
        special_hist = special_file.Get(key.GetName())
        if isinstance(special_hist,ROOT.TTree) or not special_hist:
            continue

        selection_hist = selection_file.Get(key.GetName())
        if selection_hist:
            if special_hist and bigGenie and any(signal in key.GetName() for signal in SIGNAL_DEFINATION):
                print("scaling {} integral from {} to {}".format(key.GetName(),special_hist.Integral(),selection_hist.Integral()))
                if special_hist.Integral() > 0:
                    special_hist.Scale(selection_hist.Integral()/special_hist.Integral())
                    selection_hist = special_hist.Clone()
                print("new integral is {}".format(selection_hist.Integral()))
            else:
                print("not doing anything with {}".format(key.GetName()))

            
            selection_hist.Write("", ROOT.TObject.kOverwrite)

    print("Done merging special hists... cleaning up")
    selection_file.Close()
    special_file.Close()
    print("done a file")

def MaddWrapper(output_playlist,input_files,is_data):
    args = ["madd", AnalysisConfig.SelectionHistoPath(output_playlist,is_data)]
    args.extend(input_files)
    #cmd = "madd {} {}".format(AnalysisConfig.SelectionHistoPath(output_playlist,is_data)," ".join(input_files))
    #print (cmd)
    #os.system(cmd)
    print(args)
    subprocess.run(args,stdout=subprocess.DEVNULL)
    print ("done")

def MergeHistograms():
    for sample_type in dict_of_files:
        MaddWrapper(AnalysisConfig.playlist+"temporary_special_hist",chain.from_iterable(iter(dict_of_files[sample_type]["kin"].values())),False)
        MergeTuples(AnalysisConfig.selection_file,AnalysisConfig.SelectionHistoPath(AnalysisConfig.playlist+"temporary_special_hist",False),None,None,True)

def MergeTuples(selTuple,specialTuple,pot1=None,pot2=None,bigGenie=False):
    """
    Wanna merge tuple two tuples, by direct merging or scale by reading Meta tree.
    """
    AddOneFile(selTuple,specialTuple,bigGenie)

def AddRegexMatchedFiles(dir_path,f = None):
    if f is None:
        files = os.listdir(dir_path)
    else:
        files = [f]
    filesmap = {}
    count = 0
    for string in files:
        match = re.match(SelectionFilesRegex,string)
        #print match.group(1,2,3)
        if match is not None:
            filesmap.setdefault(match.group(2),{}).setdefault(match.group(1),{}).setdefault(match.group(3),[]).append(dir_path+"/"+match.string)
            count+=1
    print(("added {} files".format(count)))
    return filesmap

if __name__ == '__main__':
    dict_of_files = AddRegexMatchedFiles(AnalysisConfig.input_dir)
    MergeHistograms()
