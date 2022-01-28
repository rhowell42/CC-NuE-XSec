import math
import ROOT
import PlotUtils
import heapq
from array import array
from collections import Iterable #for checking whether iterable.
from functools import partial
import ctypes
import code

from config.AnalysisConfig import AnalysisConfig

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.02
MNVPLOTTER.extra_top_margin = -.035# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46 
MNVPLOTTER.mc_bkgd_line_color = 46

MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .03

CANVAS = ROOT.TCanvas("c2","c2",1600,1000)

def Logx(canvas):
    canvas.SetLogx(1)
    return True

def Logy(canvas):
    canvas.SetLogy(1)
    return True

def PrepareSlicer(hist_holder):
    if hist_holder.dimension==1:
        return lambda x:[x.Clone()]
    elif hist_holder.dimension == 2:
        return Make2DSlice
    else:
        return None
        #raise KeyError("Only 1D,2D histograms are supported")

def IdentitySlicer(hist):
    return [hist.Clone()]

def ProjectionX(hist):
    return [hist.ProjectionX()]

def ProjectionY(hist):
    return [hist.ProjectionY()]

def PrepareSignalDecompose(data_hists,mc_hists,cates,sys_on_mc = True, sys_on_data = False,only_cated = False):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew

    if (data_hists.valid and mc_hists.valid):
        hists = [data_hists.GetHist().GetCVHistoWithError() if sys_on_data else data_hists.GetHist().GetCVHistoWithStatError()]
        cate_hists,colors,titles  = mc_hists.GetCateList(cates)
        if only_cated :
            totalHist = ReSumHists(cate_hists)
        else:
            totalHist = mc_hists.GetHist()
        hists.append(totalHist.GetCVHistoWithError() if sys_on_mc else totalHist.GetCVHistoWithStatError())
        hists.extend(cate_hists)
        plotfunction = lambda mnvplotter,data_hist, mc_hist, *mc_ints : partial(MakeSignalDecomposePlot,color=colors,title=titles)(data_hist,mc_hist,mc_ints)
    elif mc_hists.valid:
        cate_hists,colors,titles  = mc_hists.GetCateList(cates)
        if only_cated :
            totalHist = ReSumHists(cate_hists)
        else:
            totalHist = mc_hists.GetHist()
        hists = [totalHist.GetCVHistoWithError() if sys_on_mc else totalHist.GetCVHistoWithStatError()]
        hists.extend(cate_hists)
        plotfunction = lambda mnvplotter, mc_hist, *mc_ints : partial(MakeSignalDecomposePlot,color=colors,title=titles)(None,mc_hist,mc_ints)
    else:
        raise KeyError("Non sense making signal decomposition plots without mc")
    return plotfunction,hists

def PrepareSignalDecomposeRatio(data_hists,mc_hists,cates,sys_on_mc = True, sys_on_data = False,only_cated = False):
    def ReSumHists(h_list):
        if len(h_list) == 0:
            return None
        hnew = h_list[0].Clone()
        for i in range(1,len(h_list)):
            hnew.Add(h_list[i])
        return hnew

    if (data_hists.valid and mc_hists.valid):
        hists = [data_hists.GetHist().GetCVHistoWithError() if sys_on_data else data_hists.GetHist().GetCVHistoWithStatError()]
        cate_hists,colors,titles  = mc_hists.GetCateList(cates)
        if only_cated :
            totalHist = ReSumHists(cate_hists)
        else:
            totalHist = mc_hists.GetHist()
        hists.append(totalHist.GetCVHistoWithError() if sys_on_mc else totalHist.GetCVHistoWithStatError())
        hists.extend(cate_hists)
        plotfunction = lambda mnvplotter,data_hist, mc_hist, *mc_ints : partial(MakeSignalDecomposePlot,color=colors,title=titles)(data_hist,mc_hist,mc_ints)
    elif mc_hists.valid:
        cate_hists,colors,titles  = mc_hists.GetCateList(cates)
        if only_cated :
            totalHist = ReSumHists(cate_hists)
        else:
            totalHist = mc_hists.GetHist()
        hists = [totalHist.GetCVHistoWithError() if sys_on_mc else totalHist.GetCVHistoWithStatError()]
        hists.extend(cate_hists)
        plotfunction = lambda mnvplotter, mc_hist, *mc_ints : partial(MakeSignalDecomposePlot,color=colors,title=titles)(None,mc_hist,mc_ints)
    else:
        raise KeyError("Non sense making signal decomposition plots without mc")
    for i in hists:
        i.Divide(i,mc_hists.GetHist())
    return plotfunction,hists

def PrepareComp(data_hists,mc_hists,sys_on_mc = True, sys_on_data = False,as_frac=False):
    if not (data_hists.valid or mc_hists.valid):
        raise KeyError("both data and mc is None")
    if (data_hists.valid and mc_hists.valid):
        plotfunction = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCWithErrorBand(data_hist,mc_hist,1.0,"TR")
        hists = [data_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_data else data_hists.GetHist().GetCVHistoWithStatError(),
                 mc_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_mc else mc_hist.GetHist().GetCVHistoWithStatError()]
    elif data_hists.valid:
        plotfunction = lambda mnvplotter,data_hist: data_hist.Draw("E1X0")
        hists = [data_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_data else data_hists.GetHist().GetCVHistoWithStatError()]
    else:
        plotfunction = lambda mnvplotter,mc_hist: mnvplotter.DrawMCWithErrorBand(mc_hist)
        hists = [mc_hists.GetHist().GetCVHistoWithError(True,as_frac) if sys_on_mc else mc_hist.GetHist().GetCVHistoWithStatError()]

    return plotfunction,hists

def PrepareRatio(data_hists,mc_hists):
    if not (data_hists.valid and mc_hists.valid):
        raise KeyError("both data and mc is Required for ratio")
    plotfunction = lambda mnvplotter,data_hist, mc_hist: mnvplotter.DrawDataMCRatio(data_hist, mc_hist, 1.0 ,True,True,0,2)
    hists = [data_hists.GetHist(),
             mc_hists.GetHist()]
    return plotfunction,hists

def PrepareErr(data_hists,mc_hists,sys_on_mc=True,sys_on_data=False,grouping=None):
    plotfunction = lambda mnvplotter,data_hist: mnvplotter.DrawErrorSummary(data_hist)
    if data_hists.valid and sys_on_data:
        hists =[data_hists.GetHist()]
    elif mc_hists.valid and sys_on_mc:
        hists =[mc_hists.GetHist()]
    else:
        raise KeyError("Can't make error plot for systematics config mismatch")
    updatePlotterErrorGroup(grouping)
    #AdaptivePlotterErrorGroup(hists[0],7)
    return plotfunction,hists

def PrepareErrorBand(data_hists,mc_hists, name, sys_on_mc=True,sys_on_data=False):
    plotfunction = lambda mnvplotter,hist: MakeErrorBandPlot(hist,name,mnvplotter)
    if data_hists.valid and sys_on_data:
        hists =[data_hists.GetHist()]
    elif mc_hists.valid and sys_on_mc:
        hists =[mc_hists.GetHist()]
    else:
        raise KeyError("Can't make error plot for systematics config mismatch")
    return plotfunction,hists

def PrepareStack(data_hists,mc_hists,Grouping = None):
    if not mc_hists.valid:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")
    mc_list,color,title = mc_hists.GetCateList(Grouping)
    if data_hists.valid :
        plotfunction =  lambda mnvplotter, data_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend="TR")(data_hist,mc_ints)
        hists = [data_hists.GetHist()]
    else :
        plotfunction =  lambda mnvplotter, mc_hist, *mc_ints: partial(MakeDataMCStackedPlot, color=color,title=title, legend = "TR")(mc_hist,mc_ints)
        tmp = mc_hists.GetHist().Clone()
        tmp.Reset()
        hists =[tmp]
    hists.extend(mc_list)
    return plotfunction,hists

def PrepareDiff(data_hists,mc_hists):
    if not (data_hists.valid and mc_hists.valid):
        raise KeyError("both data and mc is Required for differece")
    hists = [data_hists.GetHist(),mc_hists.GetHist()]
    def plotDifference(mnvplotter,data_hist,mc_hist):
        #print data_hist
        ndf = ctypes.c_int()
        chi2mnv=mnvplotter.Chi2DataMC(data_hist,mc_hist,ndf)
        print(chi2mnv,ndf,data_hist.Integral("width"),mc_hist.Integral("width"))
        chi2,ndf=CalChi2(data_hist,mc_hist)
        tmp = data_hist.Clone()
        tmp.Add(mc_hist,-1)
        #print(data_hist.GetName(),tmp.Integral(0,-1,"width"))
        tmp.Draw("E1")
        size = 0.035
        align = ctypes.c_int()
        xLabel = ctypes.c_double()
        yLabel = ctypes.c_double()
        mnvplotter.DecodePosition("TR", size, align, xLabel, yLabel )
        mnvplotter.AddPlotLabel("chi2/ndf: {:.4f}/{:d}".format(chi2,ndf),xLabel.value, yLabel.value, size, 4, 112, align.value)
        print(("chi2/ndf: {:.4f}/{:d}".format(chi2,ndf)))

    return plotDifference,hists

def PrepareMigration(data_hists,mc_hists):
    if not(mc_hists.valid):
        raise KeyError("No MC histogram to plot migration")
    hists = [mc_hists.GetHist()]
    plotfunction = lambda mnvplotter,mc_hist: mnvplotter.DrawNormalizedMigrationHistogram(mc_hist,False,False,True,True)
    return plotfunction,hists

def updatePlotterErrorGroup(group,mnvplotter=MNVPLOTTER):
    mnvplotter.error_summary_group_map.clear();
    for k,v in group.items():
        vec = ROOT.vector("std::string")()
        for vs in v :
            vec.push_back(vs)
        mnvplotter.error_summary_group_map[k]= vec

def TopNErrorBand(hist,topN):
    #find N largest errorband given histogram hist
    heap = []
    for errorband_name in hist.GetVertErrorBandNames():
        sum_error = hist.GetVertErrorBand(errorband_name).GetErrorBand(False,False).Integral()
        #print(errorband_name,hist.GetVertErrorBand(errorband_name).GetErrorBand(False,False).Integral())
        if len(heap)<topN:
            heapq.heappush(heap, (sum_error,errorband_name))
        elif sum_error>heap[0][0]:
            heapq.heappushpop(heap,(sum_error,errorband_name))
    #print(heap)
    result = []
    while heap:
        result.append(heapq.heappop(heap)[1])
    return result


def AdaptivePlotterErrorGroup(hist,topN,mnvplotter=MNVPLOTTER):
    result = TopNErrorBand(hist,topN)
    rest = [i for i in hist.GetVertErrorBandNames() if i not in result]
    updatePlotterErrorGroup({"Rest":rest},mnvplotter)


def CalMXN(N_plots,max_horizontal_plot = 4):
    height = math.ceil(1.0*N_plots/max_horizontal_plot)  #hard code max 3 plots per line
    width = math.ceil(N_plots/height) #number of plots per line
    return int(width),int(height)

def Make2DSlice(hist2D_o, X_slice=False, bin_start = 1 , bin_end = 0,interval = 1):
    # x slice true produce 1d histograms of Y variable for each x variable bins.
    # bin_end = 0 : all except overflow, = -1: all including overflow, other: bins before(i.e. *not* including ) bin_end

    slicing_hists = []
    hist2D = hist2D_o.Clone()
    axis = hist2D.GetXaxis() if X_slice else hist2D.GetYaxis()
    Nbins = axis.GetNbins()
    start = max(0,bin_start)
    end = Nbins - bin_end + 1 if bin_end <=0 else bin_end
    for i in range(start,end,interval):
        slicing_hists.append(hist2D.ProjectionX(hist2D.GetName()+str(i),i,i+interval-1,"o") if not X_slice else hist2D.ProjectionY(hist2D.GetName()+str(i),i,i+interval-1,"o"))
        slicing_hists[-1].SetTitle("%.1f<%s<%.1f"%(axis.GetBinLowEdge(i),axis.GetTitle(),axis.GetBinUpEdge(i+interval-1)))
        slicing_hists[-1].GetYaxis().SetTitle(hist2D.GetZaxis().GetTitle())

    del hist2D
    return slicing_hists

def MakeDataMCPlot(data_hist, mc_hist, pot_scale=1.0, sys_on_mc = True, sys_on_data = False, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    local_mc = (mc_hist.GetCVHistoWithError() if sys_on_mc else mc_hist.GetCVHistoWithStatError()) if mc_hist else None
    local_data = (data_hist.GetCVHistoWithError() if sys_on_data else data_hist.GetCVHistoWithStatError()) if data_hist else None
    if not (local_mc or local_data):
        raise KeyError("both data and mc is None")
    if local_mc and local_data:
        mnvplotter.DrawDataMCWithErrorBand(local_data,local_mc,pot_scale,"TR")
    elif local_mc:
        mnvplotter.DrawMCWithErrorBand(local_mc,pot_scale)
    else:
        local_data.Draw("E1X0")

def MakeDataMCStackedPlot(data_hist, mc_hists, legend = "TR", pot_scale=1.0, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    if not mc_hists:
        raise KeyError("Doesn't make sense to plot stacked histogram without MC")
    TArray = ROOT.TObjArray()
    for i in range(len(mc_hists)):
        if color:
            mc_hists[i].SetFillColor(color[i])
        if title:
            mc_hists[i].SetTitle(title[i])
        TArray.Add(mc_hists[i])

    if data_hist is not None:
        mnvplotter.DrawDataStackedMC(data_hist,TArray,pot_scale,legend,"Data",0,0,1001)
    else:
        mnvplotter.DrawStackedMC(TArray,pot_scale,legend,0,0,1001)

def MakeSignalDecomposePlot(data_hist, mc_hist, mc_hists, title, color, pot_scale = 1.0, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    if data_hist is not None:
        mnvplotter.DrawDataMCWithErrorBand(data_hist,mc_hist,pot_scale,"TR")
    else:
        mnvplotter.DrawMCWithErrorBand(mc_hist,pot_scale)
    TLeg = GetTLegend(ROOT.gPad)
    for i in range(len(mc_hists)):
        mc_hists[i].SetLineColor(color[i])
        mc_hists[i].Scale(pot_scale)
        mc_hists[i].Draw("HIST SAME")
        if TLeg:
            TLeg.AddEntry(mc_hists[i],title[i])
    if TLeg:
        TLeg.Draw()

def MakeRatioPlot(data_hist,mc_hist,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    if not (data_hist and mc_hist):
        raise KeyError("both data and mc is Required for ratio")
    mnvplotter.DrawDataMCRatio(data_hist, mc_hist, 1.0 ,True,True,0,2)

def MakeErrPlot(hist,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    mnvplotter.axis_maximum = 0.3
    mnvplotter.DrawErrorSummary(hist)
    mnvplotter.axis_maximum = -1111

def MakeErrorBandPlot(hist,name,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    errorband = hist.GetVertErrorBand(name)
    errorband.DivideSingle(errorband,hist)
    errorband.DrawAll("",True)

def MakeGridPlot(MakingSlice,MakingEachPlot,input_hists,CanvasConfig=lambda canvas:True, draw_seperate_legend = False, mnvplotter=MNVPLOTTER,canvas=CANVAS,outname=None):
    if canvas is CANVAS:
        canvas.Clear()
    slices = list(map(lambda *args: args, *list(map(MakingSlice,input_hists))))
    N_plots = len(slices)
    canvas.Divide(*CalMXN(N_plots+int(draw_seperate_legend)))
    for i in range(N_plots):
        canvas.cd(i+1)
        tcanvas = canvas.GetPad(i+1)
        if not CanvasConfig(tcanvas):
            print("Warning: failed setting canvas.")
        if N_plots>1:
            SetMargin(tcanvas)
        MakingEachPlot(mnvplotter,*slices[i])
        mnvplotter.AddHistoTitle(slices[i][0].GetTitle())
        #code.interact(local=locals())
    if draw_seperate_legend:
        Tleg = None
        for i in range(N_plots):
            Tleg = GetTLegend(canvas.GetPad(i+1)) or Tleg
        if Tleg:
            canvas.cd(N_plots+1)
            Tleg.SetX1(0)
            Tleg.SetX2(1)
            Tleg.SetY1(0)
            Tleg.SetY2(1)
            Tleg.SetTextSize(2*MNVPLOTTER.legend_text_size);
            Tleg.Draw()
    if outname:
        
        mnvplotter.MultiPrint(canvas,outname,"png")

def Print(outname,mnvplotter=MNVPLOTTER,canvas=CANVAS):
    mnvplotter.MultiPrint(canvas,outname,"png")

def SetMargin(pad):
    pad.SetRightMargin(0)
    #pad.SetLeftMargin(0)
    pad.SetTopMargin(0.1)
    #pad.SetBottomMargin(0)

def GetTLegend(pad):
    tlist = pad.GetListOfPrimitives()
    for i in tlist:
        if (isinstance(i,ROOT.TLegend)) :
            pad.RecursiveRemove(i)
            pad.Update()
            return i
    return None

def MakeDataMCStackedPlot(data_hist, mc_hists, color=None, title=None, legend = "TR", pot_scale=1.0, mnvplotter=MNVPLOTTER,canvas=CANVAS):
    TArray = ROOT.TObjArray()
    for i in range(len(mc_hists)):
        if color is not None:
            mc_hists[i].SetFillColor(color[i])
        if title is not None:
            mc_hists[i].SetTitle(title[i])
        TArray.Add(mc_hists[i])
    if data_hist is not None:
        mnvplotter.DrawDataStackedMC(data_hist,TArray,pot_scale,legend,"Data",0,0,1001)
    else:
        mnvplotter.DrawStackedMC(TArray,pot_scale,legend,0,0,1001)


def MakeMigrationPlots(hist, output, no_text = False, fix_width = True, mnvplotter= MNVPLOTTER, canvas = CANVAS):
    if canvas is CANVAS:
        canvas.Clear()
    #make sure migration matrix has fix bin width
    if fix_width:
        hist.GetXaxis().Set(hist.GetNbinsX(),0,hist.GetNbinsX())
        hist.GetXaxis().SetTitle("reco bin number")
        hist.GetYaxis().Set(hist.GetNbinsY(),0,hist.GetNbinsY())
        hist.GetYaxis().SetTitle("truth bin number")
    mnvplotter.DrawNormalizedMigrationHistogram(hist,False,False,True,no_text)
    mnvplotter.MultiPrint(canvas,output)

def Make2DPlot(hist,output,mnvplotter= MNVPLOTTER, canvas = CANVAS):
    if canvas is CANVAS:
        canvas.Clear()

    hist.Draw("COLZ")
    mnvplotter.MultiPrint(canvas,output)

def CalChi2(data_hist,mc_hist,pot_scale=1.0):
    chi2 = 0
    ndf = 0
    for i in range (0,data_hist.GetSize()):
        data = data_hist.GetBinContent(i)
        mc = mc_hist.GetBinContent(i)*pot_scale
        sig = data_hist.GetBinError(i)
        #print (i,data,mc,chi2)
        chi2 += (data-mc)**2/sig/sig if data>0 else 0
        ndf += 1 if data>0 else 0
    print("data/mc integral: {}/{}".format(data_hist.Integral(),mc_hist.Integral()))
    print("chi2/ndf of histogram {} is {}/{}.".format(data_hist.GetName(),chi2,ndf))
    return chi2,ndf

