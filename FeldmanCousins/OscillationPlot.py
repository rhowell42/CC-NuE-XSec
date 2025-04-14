import os
import sys
import ROOT
import PlotUtils

import numpy as np
import math
from array import array

from tools.PlotLibrary import HistHolder
from Tools.Histogram import *
from Tools.FitTools import *
from config.AnalysisConfig import AnalysisConfig
ccnueroot = os.environ.get('CCNUEROOT')

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FixedLocator
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT.TH1.AddDirectory(False)
ROOT.SetMemoryPolicy(ROOT.kMemoryStrict)

str_to_latex = {
        "dm2"   : r"$\Delta m^{2}$",
        "ue4"   : r"$|U_{e4}|^{2}$",
        "umu4"  : r"$|U_{\mu4}|^{2}$",
        "utau4" : r"$|U_{{\tau 4}}|^2"
        }

str_to_axis = {
        "dm2"   : np.logspace(-1,2,100),
        "ue4"   : 0.15*np.logspace(-4,0,100),
        "umu4"  : 0.41*np.logspace(-5,0,100)
        }

str_to_index = {
        "dm2"   : 0,
        "ue4"   : 2,
        "umu4"  : 1
        }

str_to_indices = {
        "dm2"   : [20,27,32,39,44,51,61,71,92],
        "ue4"   : [0,53,70,78,79,80,82,84,92],
        "umu4"  : [0,50,58,62,65,67,69,72,87]
        }

hatch_styles = ["//",r"\\","||","--"]
hatch_index = 0

class ExperimentContour:
    def __init__(self,title,fname,fill,color="black"):
        self.title = title

        self.color = color
        self.patch = Line2D([0], [0], color=self.color)
        self.style = "solid"
        self.fill = fill
        self.fname = fname

        self.xaxis = ""
        self.yaxis = ""
        self.data = {}
        self.LoadData()

    def interpolate(self,inArr,size=100):
        ret = []
        for i in range(len(inArr)-1):
            x = [inArr[i],inArr[i+1]]
            new = np.linspace(x[0],x[1],size)
            ret.append(new)
        ret = np.concatenate(ret,axis=None)
        return(ret)

    def LoadData(self,header=1):
        if isinstance(self.fname,str):
            head = open(self.fname).readline().strip('\n')
            self.xaxis = head.split(',')[0]
            self.yaxis = head.split(',')[1]
            result = np.loadtxt(self.fname,delimiter=',',skiprows=header)
            if "theta" in self.xaxis:
                xdata = np.array([np.roots([-4,4,-coeff])[1] for coeff in result[:,0]])
            else:
                xdata = result[:,0]

            xdata = self.interpolate(xdata)
            if "e4" in self.xaxis:
                self.xaxis = "ue4"
            elif "mu4" in self.xaxis:
                self.xaxis = "umu4"
            self.data[self.xaxis] = xdata

            if "theta" in self.yaxis:
                ydata = np.array([np.roots([-4,4,-coeff])[1] for coeff in result[:,1]])
            else:
                ydata = result[:,1]

            ydata = self.interpolate(ydata)
            if "e4" in self.yaxis:
                self.yaxis = "ue4"
            elif "mu4" in self.yaxis:
                self.yaxis = "umu4"
            self.data[self.yaxis] = ydata

        elif isinstance(self.fname,list):
            head = open(self.fname[0]).readline().strip('\n')
            self.xaxis = head.split(',')[0]
            self.yaxis = head.split(',')[1]
            self.data[self.xaxis] = []
            self.data[self.yaxis] = []
            for name in self.fname:
                result = np.loadtxt(name,delimiter=',',skiprows=header)

                if "theta" in self.xaxis:
                    xdata = np.array([np.roots([-4,4,-coeff])[1] for coeff in result[:,0]])
                else:
                    xdata = result[:,0]

                xdata = self.interpolate(xdata)
                self.data[self.xaxis].append(xdata)

                if "theta" in self.yaxis:
                    ydata = np.array([np.roots([-4,4,-coeff])[1] for coeff in result[:,1]])
                else:
                    ydata = result[:,1]

                ydata = self.interpolate(ydata)
                self.data[self.yaxis].append(ydata)
    
            xdata = self.data[self.xaxis]
            ydata = self.data[self.yaxis]

            if "e4" in self.xaxis:
                del self.data[self.xaxis]
                self.xaxis = "ue4"
                self.data[self.xaxis] = xdata
            elif "mu4" in self.xaxis:
                del self.data[self.xaxis]
                self.xaxis = "umu4"
                self.data[self.xaxis] = xdata

            if "e4" in self.yaxis:
                del self.data[self.yaxis]
                self.yaxis = "ue4"
                self.data[self.yaxis] = ydata
            elif "mu4" in self.yaxis:
                del self.data[self.yaxis]
                self.yaxis = "umu4"
                self.data[self.yaxis] = ydata

    def SetColor(self,color):
        self.color = color

    def SetLineStyle(self,style):
        self.style = style

    def SetPatch(self,graphic):
        global hatch_index
        graphic = graphic.lower()
        if graphic == "line":
            self.patch = Line2D([0], [0], color=self.color, linestyle=self.style)
        elif graphic == "patch":
            self.patch = Patch(color=self.color) 
        elif graphic == "hatch":
            self.patch = Patch(fill=False,hatch=hatch_styles[hatch_index], edgecolor=self.color)
        else:
            raise ValueError("Graphic type not supported for legend")

    def GetPatch(self):
        return(self.patch)

    def GetTitle(self):
        return(self.title)

    def GetIntersects(self,data_ref,data_comp,axis,line):
        intersections = []
        if line > data_comp[0]:
            intersections.append(data_ref[0])

        for i in range(1,len(data_ref)-1):
            if line > data_comp[i] and line < data_comp[i-1]:
                intersections.append(data_ref[i])
            elif line < data_comp[i] and line > data_comp[i-1]:
                intersections.append(data_ref[i])

        if line > data_comp[-1]:
            intersections.append(data_ref[-1])
        box1 = []
        box2 = []
        if len(intersections) > 0:
            for i in range(0,len(intersections),2):
                point1 = intersections[i]
                point2 = intersections[i+1]

                box1.append([axis[0],axis[0],axis[-1],axis[-1]])
                box2.append([point1,point2,point2,point1])
        return(box1,box2)

    def PlotBox(self,axis,boxx,boxy):
        if self.fill:
            self.SetPatch('Patch')
            for i in range(len(boxx)):
                axis.fill(boxx[i],boxy[i],alpha=0.4,color=self.color)
        else:
            self.SetPatch('Hatch')
            for i in range(len(boxx)):
                axis.fill(boxx[i],boxy[i],fill=False,hatch=hatch_styles[hatch_index],alpha=0.2,color=self.color)

    def Plot(self,axis,xaxis,yaxis,line,panel):
        axes = self.data.keys()
        t = self.data[xaxis] if xaxis in axes else self.data[yaxis]
        global hatch_index

        if isinstance(t,list):
                for i in range(len(t)):
                    if xaxis in axes and yaxis in axes:
                        x = self.data[xaxis][i]
                        y = self.data[yaxis][i]
                        if self.fill:
                            self.SetPatch("patch")
                            axis.fill(x,y,color=self.color,alpha=0.4)
                        else:
                            self.SetPatch("line")
                            axis.plot(x,y,color=self.color,linestyle=self.style)
                    else: # need to recreate an axis
                        if xaxis not in axes:
                            ydata = self.data[yaxis][i]
                            pdata = self.data[panel][i]
                            boxx,boxy = self.GetIntersects(ydata,pdata,str_to_axis[xaxis],line)
                            self.PlotBox(axis,boxx,boxy)
                        elif yaxis not in axes:
                            xdata = self.data[xaxis][i]
                            pdata = self.data[panel][i]
                            boxx,boxy = self.GetIntersects(xdata,pdata,str_to_axis[yaxis],line)
                            self.PlotBox(axis,boxy,boxx)
        else:
            if xaxis in axes and yaxis in axes:
                x = self.data[xaxis]
                y = self.data[yaxis]
                if self.fill:
                    self.SetPatch("patch")
                    axis.fill(x,y,color=self.color,alpha=0.4)
                else:
                    self.SetPatch("line")
                    axis.plot(x,y,color=self.color,linestyle=self.style)
            else: # need to recreate an axis
                if xaxis not in axes:
                    ydata = self.data[yaxis]
                    pdata = self.data[panel]
                    boxx,boxy = self.GetIntersects(ydata,pdata,str_to_axis[xaxis],line)
                    self.PlotBox(axis,boxx,boxy)
                elif yaxis not in axes:
                    xdata = self.data[xaxis]
                    pdata = self.data[panel]
                    boxx,boxy = self.GetIntersects(xdata,pdata,str_to_axis[yaxis],line)
                    self.PlotBox(axis,boxy,boxx)

        # increment hatch styles if previously used
        if xaxis not in axes or yaxis not in axes:
            hatch_index+=1

class PanelPlot:
    def __init__(self,title,xaxis,yaxis,panel):
        self.title = title
        self.xaxis = xaxis
        self.yaxis = yaxis
        self.panel = panel

        self.x = str_to_axis[self.xaxis]
        self.y = str_to_axis[self.yaxis]
        self.p = str_to_axis[self.panel]

        self.indices = str_to_indices[panel]

        self.exclusion_results = []
        self.fit_results = []
        self.artists = []

        self.sens = None
        self.excl = None

    def SetTitle(self,title):
        self.title = title

    def SetXaxis(self,axis):
        self.xaxis = axis
        self.x = str_to_axis[axis]

    def SetYaxis(self,axis):
        self.yaxis = axis
        self.y = str_to_axis[axis]

    def SetPanel(self,axis):
        self.panel = axis
        self.p = str_to_axis[axis]
        self.indices = str_to_indices[self.panel]

    def AddExclusions(self,exp):
        for e in exp:
            self.exclusion_results.append(e)

    def AddAlloweds(self,exp):
        for e in exp:
            self.fit_results.append(e)

    def CreateAxis(self):
        nrows = math.ceil(math.sqrt(len(self.indices)))
        ncols = nrows
        self.fig, self.axes = plt.subplots(nrows,ncols,sharex=True,sharey=True,figsize=(12,8),gridspec_kw={'wspace':0, 'hspace':0})
        art1 = self.fig.suptitle(self.title,y=0.97,size=20)
        art2 = self.fig.supxlabel(str_to_latex[self.xaxis],y=0.02,size=18)
        art3 = self.fig.supylabel(str_to_latex[self.yaxis],x=.07,size=18)

        self.artists.extend([art1,art2,art3])
        textprops = dict(facecolor='white',edgecolor='white', alpha=0.8)

        for i,ax in enumerate(self.axes.flatten()):
            ax.label_outer()
            ax.set_aspect("auto")
            ax.set_xscale("log")
            ax.set_yscale("log")
            if self.yaxis == "dm2":
                ax.yaxis.set_major_locator(FixedLocator([1e2,1e1,1e0]))
            elif self.yaxis == "umu4":
                ax.yaxis.set_major_locator(FixedLocator([1e-1,1e-2,1e-3,1e-4,1e-5]))
            elif self.yaxis == "ue4":
                ax.yaxis.set_major_locator(FixedLocator([1e-1,1e-2,1e-3,1e-4]))
                
            ax.set_xlim((self.x[0],self.x[-1]))
            ax.set_ylim((self.y[0],self.y[-1]))

            ax.text(self.x[7],self.y[7],s=str_to_latex[self.panel]+ "= {:.5f}".format(self.p[self.indices[i]]),bbox=textprops,c="black",fontweight="bold")

    def PlotFeldmanCousins(self, FC_excl, FC_sens, limits):
        contour_labels = [str(i)+'%' for i in limits]
        contour_colors = ['blue','red']

        X,Y = np.meshgrid(self.x,self.y)

        for i,ax in enumerate(self.axes.flatten()):
            self.sens = ax.contour(X,Y,FC_sens[i],levels=limits,colors=contour_colors,origin="lower",linestyles='dashed')
            self.excl = ax.contour(X,Y,FC_excl[i],levels=limits,colors=contour_colors,origin="lower")

    def PlotFeldmanCousinsLists(self, FC_excl, FC_sens, colors, limits):
        contour_labels = [str(i)+'%' for i in limits]

        X,Y = np.meshgrid(self.x,self.y)

        for i,ax in enumerate(self.axes.flatten()):
            for j in range(len(colors)):
                self.sens = ax.contour(X,Y,FC_sens[j][i],levels=limits,colors=colors[j],origin="lower",linestyles='dashed')
                self.excl = ax.contour(X,Y,FC_excl[j][i],levels=limits,colors=colors[j],origin="lower")

    def PlotExclusions(self):
        global hatch_index
        for i,ax in enumerate(self.axes.flatten()):
            panel = self.p[self.indices[i]]
            for exp in self.exclusion_results:
                exp.Plot(ax,self.xaxis,self.yaxis,panel,self.panel)
            hatch_index = 0 #reset hatch index for next plot

    def PlotAlloweds(self):
        global hatch_index
        for i,ax in enumerate(self.axes.flatten()):
            panel = self.p[self.indices[i]]
            for exp in self.fit_results:
                exp.Plot(ax,self.xaxis,self.yaxis,panel,self.panel)
            hatch_index = 0 #reset hatch index for next plot

    def PlotLegend(self,limits):
        handles = []
        contour_labels = [str(i)+"%" for i in limits]
        contour_colors = ['blue','red','green']

        for c in range(len(contour_labels)):
            line = plt.Line2D([0,0], [0,0], color=contour_colors[c], linestyle='dashed')
            handles.append(line)
        for c in range(len(contour_labels)):
            line = plt.Line2D([0,0], [0,0], color=contour_colors[c])
            handles.append(line)
            
        ph = [plt.plot([],marker="", ls="")[0]]*2
        contour_labels+=contour_labels

        half = len(limits)
        handles = ph[:1] + handles[:half] + ph[1:] + handles[half:]
        contour_labels = ["Sensitivity:"] + contour_labels[:half] + ["Exclusion:"] + contour_labels[half:]

        leg = self.fig.legend(handles,contour_labels,ncol=2,title='MINERvA Feldman Cousins',bbox_to_anchor=(1.2, 0.5),loc='upper right',fontsize=14)
        leg.get_title().set_fontsize(14)

        for vpack in leg._legend_handle_box.get_children():
            for hpack in vpack.get_children()[:1]:
                hpack.get_children()[0].set_width(0)
                
        excl_handles = [exp.GetPatch() for exp in self.exclusion_results]
        excl_handles.extend([exp.GetPatch() for exp in self.fit_results])

        excl_labels = [exp.GetTitle() for exp in self.exclusion_results]
        excl_labels.extend([exp.GetTitle() for exp in self.fit_results])
        leg2 = self.fig.legend(excl_handles,excl_labels,title="External Experiments",bbox_to_anchor=(1.2, 0.75),loc='upper right',fontsize=14)
        leg2.get_title().set_fontsize(14)

        self.artists.extend([leg,leg2])

    def PlotListLegend(self,limits,titles,colors):
        handles = []
        contour_labels = [title for title in titles]
        contour_colors = colors

        for c in range(len(contour_labels)):
            line = plt.Line2D([0,0], [0,0], color=contour_colors[c], linestyle='dashed')
            handles.append(line)
        for c in range(len(contour_labels)):
            line = plt.Line2D([0,0], [0,0], color=contour_colors[c])
            handles.append(line)
            
        ph = [plt.plot([],marker="", ls="")[0]]*2
        contour_labels+=contour_labels

        half = len(titles)
        handles = ph[:1] + handles[:half] + ph[1:] + handles[half:]
        contour_labels = ["95% Sensitivity:"] + contour_labels[:half] + ["95% Exclusion:"] + contour_labels[half:]

        leg = self.fig.legend(handles,contour_labels,ncol=2,title='MINERvA Feldman Cousins',bbox_to_anchor=(1.2, 0.5),loc='upper right',fontsize=12)
        leg.get_title().set_fontsize(14)

        for vpack in leg._legend_handle_box.get_children():
            for hpack in vpack.get_children()[:1]:
                hpack.get_children()[0].set_width(0)
                
        excl_handles = [exp.GetPatch() for exp in self.exclusion_results]
        excl_handles.extend([exp.GetPatch() for exp in self.fit_results])

        excl_labels = [exp.GetTitle() for exp in self.exclusion_results]
        excl_labels.extend([exp.GetTitle() for exp in self.fit_results])
        leg2 = self.fig.legend(excl_handles,excl_labels,title="External Experiments",bbox_to_anchor=(1.2, 0.75),loc='upper right',fontsize=14)
        leg2.get_title().set_fontsize(14)

        self.artists.extend([leg,leg2])

    def MakePlot(self,FC_excl,FC_sens,limits,name):
        self.CreateAxis()
        self.PlotFeldmanCousins(FC_excl,FC_sens,limits)
        self.PlotExclusions()
        self.PlotAlloweds()
        self.PlotLegend(limits)
        print("saving figure {}".format(name))
        plt.savefig(name,bbox_extra_artists=self.artists,bbox_inches='tight')

    def MakeCompPlot(self,excls,sens,titles,limits,name):
        colors = ["red","blue","green"]
        self.CreateAxis()
        self.PlotFeldmanCousinsLists(excls,sens,colors,limits)
        self.PlotExclusions()
        self.PlotAlloweds()
        self.PlotListLegend(limits,titles,colors)
        print("saving figure {}".format(name))
        plt.savefig(name,bbox_extra_artists=self.artists,bbox_inches='tight')

    def Animate(self,FC_excl,FC_sens,limits,name,index):
        self.fig, self.axes = plt.subplots(figsize=(12,8),gridspec_kw={'wspace':0, 'hspace':0})
        self.axes.set_title(self.title,size=20)
        self.axes.set_xlabel(str_to_latex[self.xaxis],size=18)
        self.axes.set_ylabel(str_to_latex[self.yaxis],size=18)

        textprops = dict(facecolor='white',edgecolor='white', alpha=0.8)

        self.axes.label_outer()
        self.axes.set_aspect("auto")
        self.axes.set_xscale("log")
        self.axes.set_yscale("log")
            
        self.axes.set_xlim((self.x[0],self.x[-1]))
        self.axes.set_ylim((self.y[0],self.y[-1]))

        self.axes.text(self.x[7],self.y[7],s=str_to_latex[self.panel]+ "= {:.5f}".format(self.p[index]),bbox=textprops,c="black",fontweight="bold")

        contour_labels = [str(i)+'%' for i in limits]
        contour_colors = ['blue','red']

        X,Y = np.meshgrid(self.x,self.y)

        self.axes.contour(X,Y,FC_sens,levels=limits,colors=contour_colors,origin="lower",linestyles='dashed')
        self.axes.contour(X,Y,FC_excl,levels=limits,colors=contour_colors,origin="lower")

        global hatch_index
        panel = self.p[index]
        for exp in self.exclusion_results:
            exp.Plot(self.axes,self.xaxis,self.yaxis,panel,self.panel)
        hatch_index = 0 #reset hatch index for next plot
        panel = self.p[index]
        for exp in self.fit_results:
            exp.Plot(self.axes,self.xaxis,self.yaxis,panel,self.panel)
        hatch_index = 0 #reset hatch index for next plot
        self.PlotLegend(limits)
        plt.savefig(name,bbox_extra_artists=self.artists,bbox_inches='tight')
        plt.close()


def GetFCSlices(dchi2s,achi2s,results,panel):
    indices = str_to_indices[panel]
    slc = str_to_index[panel]
    sens_list = []
    excl_list = []
    for i in indices:
        FC_sens = achi2s.take(indices=i,axis=slc)
        FC_excl = dchi2s.take(indices=i,axis=slc)
        for iy, ix in np.ndindex(FC_sens.shape):
            FC_sens[iy,ix] = 100*(results[FC_sens[iy,ix] > results].shape[0]/results.shape[0])
            FC_excl[iy,ix] = 100*(results[FC_excl[iy,ix] > results].shape[0]/results.shape[0])

        sens_list.append(FC_sens)
        excl_list.append(FC_excl)
    return(sens_list,excl_list)

def GetFCSlice(dchi2s,achi2s,results,index,slc):
    FC_sens = achi2s.take(indices=index,axis=slc)
    FC_excl = dchi2s.take(indices=index,axis=slc)
    for iy, ix in np.ndindex(FC_sens.shape):
        FC_sens[iy,ix] = 100*(results[FC_sens[iy,ix] > results].shape[0]/results.shape[0])
        FC_excl[iy,ix] = 100*(results[FC_excl[iy,ix] > results].shape[0]/results.shape[0])

    return(FC_sens,FC_excl)

if __name__ == "__main__":
    filename = "NuE_stitched_hists.root"
    file_path = "{}/FeldmanCousins/rootfiles/{}".format(ccnueroot,filename)

    lam = int(AnalysisConfig.lambdaValue)
    exclude = AnalysisConfig.exclude

    histogram = StitchedHistogram("sample")
    histogram.Load(file_path)

    #invCov = histogram.GetInverseCovarianceMatrix(sansFlux=True)
    #fitter = Fitter(sample_histogram,invCov=invCov,lam=lam,exclude=exclude)
    #best_fit,res = fitter.DoFit()

    #indexed like [dm2,umu4,ue4]

    pplot = PanelPlot(title,'ue4','dm2','umu4')

    stereo = ExperimentContour("STEREO 95% Excl.","exp_results/stereo_2Dexcl.csv",False)
    neutrino4 = ExperimentContour("Neutrino-4 $2\sigma$ Conf.",["exp_results/n4_c1.csv","exp_results/n4_c2.csv","exp_results/n4_c3.csv","exp_results/n4_c4.csv"],True,"pink")
    raa = ExperimentContour("RAA 90% Allowed","exp_results/RAA.csv",True,"gray")
    minos = ExperimentContour("MINOS 90% Excl.","exp_results/MINOS.csv",False,"black")

    stereo.SetPatch("Line")
    minos.SetPatch("Line")
    neutrino4.SetPatch("Patch")
    raa.SetPatch("Patch")

    limits = [95]
    pplot.AddExclusions([stereo,minos])
    pplot.AddAlloweds([neutrino4,raa])

    if AnalysisConfig.animate:
        for i in range(len(str_to_axis["dm2"])):
            pplot.SetXaxis("ue4")
            pplot.SetYaxis("umu4")
            pplot.SetPanel("dm2")
            sens,excl = GetFCSlice(data_chi2s,asimov_chi2s,results,i,str_to_index["dm2"])
            name = "plots/ue4_vs_umu4_%02d.png" % i
            pplot.Animate(sens,excl,limits,name,i)

            pplot.SetXaxis("ue4")
            pplot.SetYaxis("dm2")
            pplot.SetPanel("umu4")
            sens,excl = GetFCSlice(data_chi2s,asimov_chi2s,results,i,str_to_index["umu4"])
            name = "plots/ue4_vs_dm2_%02d.png" % i
            pplot.Animate(sens,excl,limits,name,i)
            
            pplot.SetXaxis("umu4")
            pplot.SetPanel("ue4")
            sens,excl = GetFCSlice(data_chi2s,asimov_chi2s,results,i,str_to_index["ue4"])
            name = "plots/umu4_vs_dm2_%02d.png" % i
            pplot.Animate(sens,excl,limits,name,i)
    elif AnalysisConfig.compare_profiles:
        lam = 1
        best_fit = 101.24
        exclude = "none"
        data_chi2s = np.load("chi2s/lambda{}_{}/data_chi2s.npy".format(lam,exclude)) - best_fit
        asimov_chi2s = np.load("chi2s/lambda{}_{}/asimov_chi2s.npy".format(lam,exclude))
        results = np.load("chi2s/lambda{}_{}/asimov_deltachi2s.npy".format(lam,exclude))
        sens_list1,excl_list1 = GetFCSlices(data_chi2s,asimov_chi2s,results,"umu4")

        lam = 12
        best_fit = 152.65
        exclude = "none"
        data_chi2s = np.load("chi2s/lambda{}_{}/data_chi2s.npy".format(lam,exclude)) - best_fit
        asimov_chi2s = np.load("chi2s/lambda{}_{}/asimov_chi2s.npy".format(lam,exclude))
        results = np.load("chi2s/lambda{}_{}/asimov_deltachi2s.npy".format(lam,exclude))
        sens_list2,excl_list2 = GetFCSlices(data_chi2s,asimov_chi2s,results,"umu4")

        lam = 1
        exclude = "ratio"
        best_fit = 139.336
        data_chi2s = np.load("chi2s/lambda{}_{}/data_chi2s.npy".format(lam,exclude)) - best_fit
        asimov_chi2s = np.load("chi2s/lambda{}_{}/asimov_chi2s.npy".format(lam,exclude))
        results = np.load("chi2s/lambda{}_{}/asimov_deltachi2s.npy".format(lam,exclude))
        sens_list3,excl_list3 = GetFCSlices(data_chi2s,asimov_chi2s,results,"umu4")

        pplot.SetTitle("MINERvA Sterile Neutrino Search\nFlux Profiling Comparison")
        pplot.SetXaxis("ue4")
        pplot.SetYaxis("dm2")
        pplot.SetPanel("umu4")
        pplot.MakeCompPlot([excl_list1,excl_list2,excl_list3],[sens_list1,sens_list2,sens_list3],["$\lambda=1$","$\lambda=12$","$\lambda=1$ Sans\nFlavor Ratio"],limits,"plots/FC_ue4_vs_dm2_profiling.png")

        lam = 1
        best_fit = 101.24
        exclude = "none"
        data_chi2s = np.load("chi2s/lambda{}_{}/data_chi2s.npy".format(lam,exclude)) - best_fit
        asimov_chi2s = np.load("chi2s/lambda{}_{}/asimov_chi2s.npy".format(lam,exclude))
        results = np.load("chi2s/lambda{}_{}/asimov_deltachi2s.npy".format(lam,exclude))
        sens_list1,excl_list1 = GetFCSlices(data_chi2s,asimov_chi2s,results,"ue4")

        lam = 12
        best_fit = 152.65
        exclude = "none"
        data_chi2s = np.load("chi2s/lambda{}_{}/data_chi2s.npy".format(lam,exclude)) - best_fit
        asimov_chi2s = np.load("chi2s/lambda{}_{}/asimov_chi2s.npy".format(lam,exclude))
        results = np.load("chi2s/lambda{}_{}/asimov_deltachi2s.npy".format(lam,exclude))
        sens_list2,excl_list2 = GetFCSlices(data_chi2s,asimov_chi2s,results,"ue4")

        lam = 1
        exclude = "ratio"
        best_fit = 139.336
        data_chi2s = np.load("chi2s/lambda{}_{}/data_chi2s.npy".format(lam,exclude)) - best_fit
        asimov_chi2s = np.load("chi2s/lambda{}_{}/asimov_chi2s.npy".format(lam,exclude))
        results = np.load("chi2s/lambda{}_{}/asimov_deltachi2s.npy".format(lam,exclude))
        sens_list3,excl_list3 = GetFCSlices(data_chi2s,asimov_chi2s,results,"ue4")

        pplot.SetXaxis("umu4")
        pplot.SetYaxis("dm2")
        pplot.SetPanel("ue4")
        pplot.MakeCompPlot([excl_list1,excl_list2,excl_list3],[sens_list1,sens_list2,sens_list3],["$\lambda=1$","$\lambda=12$","$\lambda=1$ Sans\nFlavor Ratio"],limits,"plots/FC_umu4_vs_dm2_profiling.png")
        
        lam = 1
        best_fit = 101.24
        exclude = "none"
        data_chi2s = np.load("chi2s/lambda{}_{}/data_chi2s.npy".format(lam,exclude)) - best_fit
        asimov_chi2s = np.load("chi2s/lambda{}_{}/asimov_chi2s.npy".format(lam,exclude))
        results = np.load("chi2s/lambda{}_{}/asimov_deltachi2s.npy".format(lam,exclude))
        sens_list1,excl_list1 = GetFCSlices(data_chi2s,asimov_chi2s,results,"dm2")

        lam = 12
        best_fit = 152.65
        exclude = "none"
        data_chi2s = np.load("chi2s/lambda{}_{}/data_chi2s.npy".format(lam,exclude)) - best_fit
        asimov_chi2s = np.load("chi2s/lambda{}_{}/asimov_chi2s.npy".format(lam,exclude))
        results = np.load("chi2s/lambda{}_{}/asimov_deltachi2s.npy".format(lam,exclude))
        sens_list2,excl_list2 = GetFCSlices(data_chi2s,asimov_chi2s,results,"dm2")

        lam = 1
        exclude = "ratio"
        best_fit = 139.336
        data_chi2s = np.load("chi2s/lambda{}_{}/data_chi2s.npy".format(lam,exclude)) - best_fit
        asimov_chi2s = np.load("chi2s/lambda{}_{}/asimov_chi2s.npy".format(lam,exclude))
        results = np.load("chi2s/lambda{}_{}/asimov_deltachi2s.npy".format(lam,exclude))
        sens_list3,excl_list3 = GetFCSlices(data_chi2s,asimov_chi2s,results,"dm2")

        pplot.SetXaxis("ue4")
        pplot.SetYaxis("umu4")
        pplot.SetPanel("dm2")
        pplot.MakeCompPlot([excl_list1,excl_list2,excl_list3],[sens_list1,sens_list2,sens_list3],["$\lambda=1$","$\lambda=12$","$\lambda=1$ Sans\nFlavor Ratio"],limits,"plots/FC_ue4_vs_umu4_profiling.png")

    else:
        title = 'MINERvA Sterile Neutrino Search\nFlux Profiling $\lambda={}$ '.format(lam)
        if exclude == 'ratio':
            title+='Sans Flavor Ratio'
            best_fit = 139.336
        elif lam == 1:
            best_fit = 101.24
        elif lam == 12:
            best_fit = 152.65
        data_chi2s = np.load("chi2s/lambda{}_{}/data_chi2s.npy".format(lam,exclude)) - best_fit
        asimov_chi2s = np.load("chi2s/lambda{}_{}/asimov_chi2s.npy".format(lam,exclude))
        results = np.load("chi2s/lambda{}_{}/asimov_deltachi2s.npy".format(lam,exclude))

        pplot.SetXaxis("ue4")
        pplot.SetYaxis("dm2")
        pplot.SetPanel("umu4")
        sens_list,excl_list = GetFCSlices(data_chi2s,asimov_chi2s,results,"umu4")
        pplot.MakePlot(excl_list,sens_list,limits,"plots/FC_ue4_vs_dm2_lambda{}_exclude_{}.png".format(lam,exclude))

        pplot.SetXaxis("umu4")
        pplot.SetYaxis("dm2")
        pplot.SetPanel("ue4")
        sens_list,excl_list = GetFCSlices(data_chi2s,asimov_chi2s,results,"ue4")
        pplot.MakePlot(excl_list,sens_list,limits,"plots/FC_umu4_vs_dm2_lambda{}_exclude_{}.png".format(lam,exclude))
        
        pplot.SetXaxis("ue4")
        pplot.SetYaxis("umu4")
        pplot.SetPanel("dm2")
        sens_list,excl_list = GetFCSlices(data_chi2s,asimov_chi2s,results,"dm2")
        pplot.MakePlot(excl_list,sens_list,limits,"plots/FC_ue4_vs_umu4_lambda{}_exclude_{}.png".format(lam,exclude))
