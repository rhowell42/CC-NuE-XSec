import os
import sys
import ROOT
import PlotUtils

import numpy as np
import math
from array import array

from tools.StitchedHistogram import *
from tools.Fitters import *
ccnueroot = os.environ.get('CCNUEROOT')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

def MergeAsimovs():
    while True:
        print("Type the path of the directory that contains the Asimov delta chi2 files")
        path = input()
        if len(path) == 0:
            print("Type a valid path")
            continue

        if os.path.isdir(path):
            results = []
            for f in os.listdir(path):
                f_path = path+'/'+f
                print(f_path)
                if os.path.isdir(f_path):
                    for f_ in os.listdir(f_path):
                        file = f_path+'/'+f_
                        res = np.loadtxt(file,delimiter=',')
                        if np.isnan(res).any() or np.isinf(res).any():
                            print("Bad delta chi2 computed in {}".format(file))
                            continue
                        results.append(res)
                else:
                    res = np.loadtxt(f_path,delimiter=',')
                    if np.isnan(res).any() or np.isinf(res).any():
                        print("Bad delta chi2 computed in {}".format(f_path))
                        continue
                    results.append(res)
            results = np.array(results).flatten()
            print('saving asimov_deltachi2s.npy with {} entries'.format(results.shape[0]))
            np.save("asimov_deltachi2s",results)
            return(results)
            break
        else:
            print("Not a valid path")
            break

def MergeChi2s():
    while True:
        print("Type the path of the directory that contains the chi2 contours you want to merge")
        path = input()
        if len(path) == 0:
            print("Enter a valid path")
            break

        if path[-1] != '/':
            path = path+'/'
            
        if os.path.isdir(path):
            data_files = [f for f in os.listdir(path) if "chi2_surface" in f]
            if len(data_files) == 0:
                print("No files found in {}".format(path))
                break
        else:
            print("directory does not exist")
            break

        # ----- sort file names to order PMNS parametesr ----- #
        start = '_m_'
        end = '_Ue4'
        m_names = [float(s[s.find(start)+len(start):s.rfind(end)]) for s in data_files]
        m_names = list(set(m_names))
        m_names.sort()

        start = 'Ue4_'
        end = '.dat.npy'
        e_names = [float(s[s.find(start)+len(start):s.rfind(end)]) for s in data_files]
        e_names = list(set(e_names))
        e_names.sort()

        datas = []
        asimovs = []
        datas_penalty = []
        asimovs_penalty = []

        for m in m_names:
            m_data   = []
            m_data_penalty = []

            m_asimov = []
            m_asimov_penalty = []

            for e in e_names:
                data = np.load(path+"chi2_surface_data_m_{}_Ue4_{}.dat.npy".format(m,e))
                data_penalty = np.load(path+"chi2_penalty_data_m_{}_Ue4_{}.dat.npy".format(m,e))

                asimov = np.load(path+"chi2_surface_pseudodata_m_{}_Ue4_{}.dat.npy".format(m,e))
                asimov_penalty = np.load(path+"chi2_penalty_pseudodata_m_{}_Ue4_{}.dat.npy".format(m,e))

                m_data.append(data)
                m_data_penalty.append(data_penalty)

                m_asimov.append(asimov)
                m_asimov_penalty.append(asimov_penalty)

            datas.append(m_data)
            datas_penalty.append(m_data_penalty)

            asimovs.append(m_asimov)
            asimovs_penalty.append(m_asimov_penalty)

        datas = np.array(datas)
        datas_penalty = np.array(datas_penalty)

        asimovs = np.array(asimovs)
        asimovs_penalty = np.array(asimovs_penalty)

        np.save("data_chi2s",datas)
        np.save("data_penalties",datas_penalty)
        
        np.save("asimov_chi2s",asimovs)
        np.save("asimov_penalties",asimovs_penalty)

        print("Done saving files")

        break

if __name__ in "__main__":
    print("Do you want to merge the chi2 contour files? (y/n)")
    ans = input().lower()
    if ans == 'y':
        MergeChi2s()

    print("Do you want to merge asimov delta chi2 files? (y/n)")
    ans = input().lower()
    if ans == 'y':
        asimovs = MergeAsimovs()

        if False:
            # Plot a histogram of the asimov delta chi2s with the critical levels marked
            # First fit to data and use that for the critical chi2s
            filename = "rootfiles/NuE_stitched_hists.root"
            
            file_path = "{}/oscillations/{}".format(ccnueroot,filename)

            sample_histogram = StitchedHistogram("sample")
            sample_histogram.Load(file_path)

            stat = Statistics(sample_histogram,lam=AnalysisConfig.lambdaValue,exclude=AnalysisConfig.exclude)
            chi2_null,penalty = stat.Chi2DataMC(marginalize=True)
            print("null chi2: {:.3f}".format(chi2_null))

            fitter = OscillationFitter(sample_histogram,lam=AnalysisConfig.lambdaValue,exclude=AnalysisConfig.exclude)
            chi2_fit,res = fitter.DoFit()

            delta_chi2 = chi2_null-chi2_fit

            print("Data fit: delta chi2 = {:.3f} = {:.3f} - {:.3f}".format(delta_chi2,chi2_null,chi2_fit))
            print("Best fit params:")
            print("   delta m^2 = {:.3f} eV^2 +- {:.4f}".format(res['m'],0))
            print("   U_e4^2    = {:.3f}      +- {:.4f}".format(res['ue4'],0))
            print("   U_mu4^2   = {:.5f}    +- {:.4f}".format(res['umu4'],0))
            print("   U_tau4^2  = {:.3f}      +- {:.4f}".format(res['utau4'],0))

        delta_chi2 = 1.097
        #plt.figure(figsize=(8,6))
        n, bins, patches = plt.hist(asimovs,bins=100,range=(0,100),alpha=0.5,label='Null Hypothesis\nPseudo Experiments',density=True)
        plt.plot((delta_chi2,delta_chi2),(0,np.max(n)),label='Data $\Delta\chi^2={:.2f}$'.format(delta_chi2))
        for i in [68,90,95,99]:
            for j in np.linspace(0,50,1000):
                c = asimovs[asimovs<j].shape[0]
                if 100*c/asimovs.shape[0] > i:
                    plt.plot((j,j),(0,np.max(n)),label="{}% at $\Delta\chi^2={:.2f}$".format(i,j),linestyle='--')
                    break
            
        plt.title("$\Delta\chi^2s$ for Null Hypothesis Pseudo Experiments\n"+"{} Experiments".format(asimovs.shape[0]),fontsize=20)
        plt.xlabel(r"$\Delta\chi^2=\chi^2(null) - \chi^2_{min}(best\ osc.)$",fontsize=16)
        plt.ylabel("Frequency",fontsize=16)
        plt.yscale('log')
        plt.xlim(0,100)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=14)
        plt.savefig('histogram_plot.png', dpi=300, bbox_inches='tight')
