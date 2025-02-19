import os
import datetime as dt
import argparse
ccnueroot = os.environ.get('CCNUEROOT')
import logging, sys
import ROOT
import PlotUtils
import numpy as np
from scipy import optimize,integrate

import math
import psutil
import multiprocessing
import threading
nthreads = 1
from array import array

#insert path for modules of this package.
from config.SignalDef import SIGNAL_DEFINATION
from tools.PlotLibrary import HistHolder
from fit_tools.FitTools import *

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

MNVPLOTTER = PlotUtils.MnvPlotter()

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

if __name__ == "__main__":
    BLUEARC = "/exp/minerva/data/users/{}/chi2s".format(os.environ["USER"])
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--grid",
                        action="store_true",
                        default = False,
                        help = "Run macro on grid, Input/Output path must be updated to avoid direct access to BlueArc"
    )
    parser.add_argument("-o", "--output",
                        dest = "output_dir",
                        help="Use alternate location for output file.",
                        default=BLUEARC
    )
    parser.add_argument("-rs", "--rowstart",
                        dest = "row",
                        help="Row of the scale matrix to start at.",
                        default=0
    )
    parser.add_argument("-rtd", "--rowstodo",
                        dest = "rowstodo",
                        help="Rows of the scale matrix to do",
                        default=10
    )
    args = parser.parse_args()
    outdir_chi2s = args.output_dir
    rowstart = int(args.row)
    rowstodo = int(args.rowstodo)
    runongrid = args.grid
    try:
        logging.info("loading scale matrices...")
        scale_matrix = np.loadtxt("{}/FeldmanCousins/ScaleMatrix_0.txt".format(ccnueroot),dtype=np.float32)
        for i in range(1,10):
            matrix = np.loadtxt("{}/FeldmanCousins/ScaleMatrix_{}.txt".format(ccnueroot,i),dtype=np.float32)
            scale_matrix = np.concatenate((scale_matrix,matrix))
    except:
        logging.info("Couldn't find all of the scale matrices, you should make them.")

    num_rows, num_cols = scale_matrix.shape
    logging.info("Got scale matrix with shape {} x {}".format(num_rows,num_cols))

    stitched_data = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('data_stitched')
    stitched_mc = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched')
    
    stitched_nueTemp = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('LE_template_nue')
    stitched_anueTemp = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('LE_template_anue')
    stitched_numuTemp = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('LE_template_numu')
    stitched_anumuTemp = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('LE_template_anumu')

    stitched_nue_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_nue')
    stitched_anue_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_anue')
    stitched_numu_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_numu')
    stitched_anumu_energy = ROOT.TFile.Open("{}/FeldmanCousins/NuE_stitched_hists.root".format(ccnueroot)).Get('mc_stitched_anumu')


    for i in range(1,stitched_mc.GetNbinsX()):
        summed = stitched_nue_energy.GetBinContent(i) + stitched_anue_energy.GetBinContent(i) + stitched_numu_energy.GetBinContent(i) + stitched_anumu_energy.GetBinContent(i)
        ratio = stitched_nue_energy.GetBinContent(i) / summed
        new_nue = stitched_mc.GetBinContent(i) * ratio
        ratio = stitched_anue_energy.GetBinContent(i) / summed
        new_anue = stitched_mc.GetBinContent(i) * ratio
        ratio = stitched_numu_energy.GetBinContent(i) / summed
        new_numu = stitched_mc.GetBinContent(i) * ratio
        ratio = stitched_anumu_energy.GetBinContent(i) / summed
        new_anumu = stitched_mc.GetBinContent(i) * ratio
        if i <= 12:
            stitched_nue_energy.SetBinContent(i,new_nue)
            stitched_anue_energy.SetBinContent(i,new_anue)
            stitched_numu_energy.SetBinContent(i,new_numu)
            stitched_anumu_energy.SetBinContent(i,new_anumu)

    templates = {
            "nue":stitched_nueTemp,
            "anue":stitched_anueTemp,
            "numu":stitched_numuTemp,
            "anumu":stitched_anumuTemp,
            "nue_energy":stitched_nue_energy,
            "anue_energy":stitched_anue_energy,
            "numu_energy":stitched_numu_energy,
            "anumu_energy":stitched_anumu_energy
            }

    if not runongrid:
        logging.info("Making chi2 text files...")
        for p in range(rowstart,num_rows,rowstodo*nthreads):
            upto = p+rowstodo
            
            getFits(scale_matrix[p:upto,:],stitched_data,p,templates, stitched_mc)
            #print("going to loop through: {} -> {}".format(p,upto))
            #print("going to loop through: {} -> {}".format(upto,upto+rowstodo))
            #print("going to loop through: {} -> {}".format(upto+rowstodo,upto+2*rowstodo))
            #t1 = multiprocessing.Process(target=getFits, args=(scale_matrix[p:upto,:],stitched_mc,p), name='t1')
            #t2 = multiprocessing.Process(target=getFits, args=(scale_matrix[upto:upto+rowstodo,:],stitched_mc,upto), name='t2')
            #t3 = multiprocessing.Process(target=getFits, args=(scale_matrix[upto+rowstodo:upto+2*rowstodo,:],stitched_mc,upto+rowstodo), name='t3')

            #t1.start()
            #t2.start()
            #t3.start()

            #t1.join()
            #t2.join()
            #t3.join()
            rowstart+=1

        logging.info("Done making chi2 text files, exiting...")
    else:
        getFits(scale_matrix[rowstart:rowstart+rowstodo,:],stitched_mc,rowstart,templates,stitched_data)
