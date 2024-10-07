import os
import sys
import ROOT
import PlotUtils
from array import array

#insert path for modules of this package.
from config import PlotConfig
from config.AnalysisConfig import AnalysisConfig
from config.DrawingConfig import PLOTS_TO_MAKE,Default_Plot_Type,Default_Scale,DefaultPlotters,DefaultSlicer
from tools import Utilities,PlotTools
from tools.PlotLibrary import HistHolder

MNVPLOTTER = PlotUtils.MnvPlotter()
#config MNVPLOTTER:
MNVPLOTTER.draw_normalized_to_bin_width=False
MNVPLOTTER.legend_text_size = 0.02
#MNVPLOTTER.extra_top_margin = -.035# go slightly closer to top of pad
MNVPLOTTER.mc_bkgd_color = 46 
MNVPLOTTER.mc_bkgd_line_color = 46

MNVPLOTTER.data_bkgd_color = 12 #gray
MNVPLOTTER.data_bkgd_style = 24 #circle  

#legend entries are closer
MNVPLOTTER.height_nspaces_per_hist = 1.2
MNVPLOTTER.width_xspace_per_letter = .4
MNVPLOTTER.legend_text_size        = .03

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)
SELECTED_SIDEBANDS = AnalysisConfig.sidebands

def MakePlot(templates,mcsignal,datasignal,true_mc=False):
    #scale
    Default_Scale(datasignal)
    Default_Scale(mcsignal)
    Default_Scale(templates)
    
    CanvasConfig = template.setdefault("canvastemplate",lambda x:True)

    hist = templates.GetHist() 
    gr = hist.ProfileY()
    sinavg.append(gr)

    CanvasConfig = template.setdefault("canvasconfig", PlotTools.Logz)
    
    iterations = 100
    _h = ROOT.TH2D( 'chi2', 'sin2ee vs sin2ee', iterations, 0, 1, iterations, 0, 1)

    fakedata = mcsignal.GetHist().Clone()
    ### Only want to look at flux error
    #for errorband in fakedata.GetVertErrorBandNames():
        #if errorband != "Statistical":
        #        fakedata.PopVertErrorBand(errorband)
    
    errors = datasignal.GetHist().Clone()
    errors.PopVertErrorBand("SuSA_Valencia_Weight")
    errors.PopVertErrorBand("MK_model")
    errors.PopVertErrorBand("LowQ2Pi_None")
    fakedata.PopVertErrorBand("SuSA_Valencia_Weight")
    fakedata.PopVertErrorBand("MK_model")
    fakedata.PopVertErrorBand("LowQ2Pi_None")
    mcTotalWithSys = fakedata.Clone() # for MCRatio plot
    #errors.PopVertErrorBand("Statistical")
    #c = 0
    ### Only want to look at flux error
    #for errorband in errors.GetVertErrorBandNames():
    #    if errorband != "Statistical":# and errorband != "Leakage_Uncertainty":
    #        errors.PopVertErrorBand(errorband)
    #        MNVPLOTTER.DrawErrorSummary(errors)
    #        PlotTools.Print(AnalysisConfig.PlotPath(templates.plot_name,sideband,("data_errorsummary_without_"+str(c)+str(errorband))))            
    #        c+=1
    #MNVPLOTTER.DrawErrorSummary(errors)
    #PlotTools.Print(AnalysisConfig.PlotPath(templates.plot_name,sideband,("data_errorsummary_with_Leakage_Uncertainty")))            
    #print("GOT HERE*************")

    #c = 0
    ### Only want to look at flux error
    #for errorband in fakedata.GetVertErrorBandNames():
    #    if errorband != "Statistical":# and errorband != "Leakage_Uncertainty":
    #        fakedata.PopVertErrorBand(errorband)
    #        MNVPLOTTER.DrawErrorSummary(fakedata)
    #        PlotTools.Print(AnalysisConfig.PlotPath(templates.plot_name,sideband,("data_fakedataummary_without_"+str(c)+str(errorband))))            
    #        c+=1
    #MNVPLOTTER.DrawErrorSummary(fakedata)
    #PlotTools.Print(AnalysisConfig.PlotPath(templates.plot_name,sideband,("mc_fakedataummary_with_Leakage_Uncertainty")))            
    #print("GOT HERE*************")
    
    fakedata = fakedata.GetCVHistoWithError(False)
    
    for m in range(0,iterations):
        sinee = m/iterations
        _fakedata = fakedata.Clone()
        _test = fakedata.Clone()
        for i in range(gr.GetNbinsX()):
            if gr.GetBinContent(i) == 0.0125: # catch for the tprofile being fooled by binwidth in case of no oscillation
                gr.SetBinContent(i,0)
            if _fakedata.GetBinContent(i) < 0:
                _fakedata.SetBinContent(i,0)
            P_ee = 1 - sinee*gr.GetBinContent(i)
            before = _fakedata.GetBinContent(i)
            after = P_ee*before
            try:
                newerr = ((errors.GetCVHistoWithError(False).GetBinError(i)*after/before)**2 + (after**.5)**2)**.5
            except:
                newerr = ((errors.GetCVHistoWithError(False).GetBinError(i))**2 + (after**.5)**2)**.5
            _fakedata.SetBinContent(i,after)
            _fakedata.SetBinError(i,newerr)
        #_fakedata.SetBinError(i,after**.5)


        for n in range(0,iterations):
            testee = n/iterations
            _mcfit = _test.Clone()
            for j in range(gr.GetNbinsX()):
                if _mcfit.GetBinContent(j) < 0:
                    _mcfit.SetBinContent(j,0)
                P_ee = 1 - testee*gr.GetBinContent(j)
                before = _mcfit.GetBinContent(j)
                after = P_ee*before
                try:
                    newerr = ((_mcfit.GetBinError(j)*before/after)**2 + (after**.5)**2)**.5
                except:
                    newerr = ((_mcfit.GetBinError(j))**2 + (after**.5)**2)**.5
                _mcfit.SetBinContent(j,after)
                _mcfit.SetBinError(j,newerr)
        #_mcfit.SetBinError(j,after**.5)

            original = _fakedata.Clone()
            chi2 = MNVPLOTTER.Chi2DataMC(original,_mcfit)/28
            _h.Fill(m,n)
            _h.SetBinContent(m,n,chi2)
    _h.GetZaxis().SetRangeUser(0,100)
    _h.GetXaxis().SetTitle("pseudo data sin^{2}(2#theta_{ee})")
    _h.GetYaxis().SetTitle("MC sin^{2}(2#theta_{ee})")
    _h.SetTitle("All Error Bands -- \Delta m^2 = {}".format(template["name"][-2:]))
    canv = ROOT.gPad

    Length = array('d',[0.0,1,4,9,16,25,100])
    _h.SetContour(7,Length)
    _h.GetZaxis().SetTitle("#chi^{2}")
    #canv.GetCanvas().SetLogz()
    
    _h.Draw("colz")
    PlotTools.Print(AnalysisConfig.PlotPath(templates.plot_name,sideband,"sin2_chi2"))            
    canv.Clear()

    CanvasConfig = template.setdefault("canvastemplate",lambda x:True)
    unoscillated = mcTotalWithSys.Clone().GetCVHistoWithError()
    for i in range(gr.GetNbinsX()):
        if fakedata.GetBinContent(i) < 0:
            fakedata.SetBinContent(i,0)
        P_ee = 1 - 0.44*gr.GetBinContent(i)
        before = fakedata.GetBinContent(i)
        after = P_ee*before
        try:
            newerr = ((errors.GetCVHistoWithError(False).GetBinError(i)*after/before)**2 + (after**.5)**2)**.5
        except:
            newerr = ((errors.GetCVHistoWithError(False).GetBinError(i))**2 + (after**.5)**2)**.5
        fakedata.SetBinContent(i,after)
        fakedata.SetBinError(i,newerr)
    #fakedata.SetBinError(i,after**.5)
        unoscillated.SetBinError(i,unoscillated.GetBinError(i))
    #unoscillated.SetBinError(i,unoscillated.GetBinContent(i)**.5)

    bottomFraction = .2
    
    margin = .12
    overall = ROOT.TCanvas("Data/MC for " + str(template))
    top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()
    MNVPLOTTER.DrawDataMCWithErrorBand(fakedata,unoscillated,1,"TR")
    MNVPLOTTER.AddChi2Label(fakedata,unoscillated,1,"BC",.04,.45)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    ratio = fakedata.Clone()
    ratio.Divide(ratio, unoscillated)

    #TODO: I need GetCVHistoWithError() from mcRatio, but THStack doesn't keep a MnvH1D.  I have to Add() the histograms myself.

    #Now fill mcRatio with 1 for bin content and fractional error
    mcRatio = mcTotalWithSys.GetTotalError(False, True, False) #The second "true" makes this fractional error
    for whichBin in range(1, mcRatio.GetXaxis().GetNbins()+1): 
        mcRatio.SetBinError(whichBin, max(mcRatio.GetBinContent(whichBin), 1e-9))
        mcRatio.SetBinContent(whichBin, 1)

    ratio.SetTitle("")
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(3)
    ratio.SetTitleSize(0)

    ratio.GetYaxis().SetTitle("Data / MC")
    ratio.GetYaxis().SetLabelSize(.15)
    ratio.GetYaxis().SetTitleSize(0.16)
    ratio.GetYaxis().SetTitleOffset(0.4)
    ratio.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.

    ratio.GetXaxis().SetTitleSize(0.16)
    ratio.GetXaxis().SetTitleOffset(0.9)
    ratio.GetXaxis().SetLabelSize(.15)

    ratio.SetMinimum(.5)
    ratio.SetMaximum(1.5)
    ratio.Draw()

    #Error envelope for the MC
    mcRatio.SetLineColor(ROOT.kRed)
    mcRatio.SetLineWidth(3)
    mcRatio.SetFillColorAlpha(ROOT.kPink + 1, 0.4)
    mcRatio.Draw("E2 SAME")

    #Draw a flat line at 1 for ratio of MC to itself
    straightLine = mcRatio.Clone()
    straightLine.SetFillStyle(0)
    straightLine.Draw("HIST SAME")
    top.cd()
    title = ROOT.TPaveText(0.3, 0.91, 0.7, 1.0, "nbNDC") #no border and use Normalized Device Coordinates to place it
    title.SetFillStyle(0)
    title.SetLineColor(ROOT.kWhite)
    title.AddText("Pseudo Data Oscillation")
    title.Draw()

    MNVPLOTTER.WritePreliminary(0.3, 0.82, 5e-2, True)
    overall.Print(AnalysisConfig.PlotPath(templates.plot_name,sideband,"oscillated.png"))            

    MNVPLOTTER.DrawErrorSummary(errors)
    PlotTools.Print(AnalysisConfig.PlotPath(templates.plot_name,sideband,"data_errorsummary"))            

    MNVPLOTTER.DrawErrorSummary(mcsignal.GetHist())
    PlotTools.Print(AnalysisConfig.PlotPath(templates.plot_name,sideband,"mc_errorsummary"))            
    
    return True

if __name__ == "__main__":
    #input knobs
    playlist=AnalysisConfig.playlist
    #signalplaylist=AnalysisConfig.playlist[:-3]+"_Signal"+AnalysisConfig.playlist[-3:]

    bkg_file_path = AnalysisConfig.BackgroundFitPath(playlist,AnalysisConfig.bkgTune_tag,False)

    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    #signaltype_path_map = { t:AnalysisConfig.SelectionHistoPath(signalplaylist,t =="data",False) for t in AnalysisConfig.data_types}

    data_file,mc_file,pot_scale,data_pot,mc_pot = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag,True)
    #signaldata_file,signalmc_file,signalpot_scale,signaldata_pot,signalmc_pot = Utilities.getFilesAndPOTScale(signalplaylist,signaltype_path_map,AnalysisConfig.ntuple_tag,True)

    #bkg_file_path=bkg_file_path.replace("_Templates","")
    #bkg_file_path=bkg_file_path.replace("Templates","Eavail200MeV")
    #bkg_file_path=bkg_file_path.replace("Templates","dEdx")
    #bkg_file_path=bkg_file_path.replace("_Templates","")
    bkgdata = ROOT.TFile.Open(bkg_file_path)
    bkgmc = ROOT.TFile.Open(bkg_file_path)

    standPOT = data_pot if data_pot is not None else mc_pot 
    #signalstandPOT = signaldata_pot if signaldata_pot is not None else signalmc_pot 
    sideband_map = {}

    sinavg = []

    for template in PLOTS_TO_MAKE:
        sideband_group =  template.setdefault("sideband_group",["Signal"]+SELECTED_SIDEBANDS)
        if isinstance(sideband_group,list):
            for sideband in sideband_group:
                templates = HistHolder(template["name"] if "name" in template else template,mc_file,sideband,True,mc_pot,standPOT)
                mcsignal = HistHolder("Background Subbed MC",bkgmc,"Signal",True,mc_pot,standPOT)
                datasignal = HistHolder("Background Subbed Data",bkgdata,"Signal",False,mc_pot,standPOT)
                MakePlot(templates,mcsignal,datasignal)
                 
    ms = [i for i in range(1,16)] + [20,30,50,-1]
    deltam2 = 0
    legend = ROOT.TLegend(.45,.6,.9,.9)
    colors = [1,2,3,4,6]
    lines = []
    color = 0
    for gr in sinavg[1:]:
        line = gr.Clone()
        m2 = ms[deltam2]
        option = "SAME"
        if m2 in [1,6,11,20]:
            CanvasConfig = PLOTS_TO_MAKE[deltam2].setdefault("canvastemplate",lambda x:True)
            option = ""
            legend = ROOT.TLegend(.7,.6,.9,.9)
            gr.SetTitle("Average sin^{2}(1.27 #Delta m^{2} L/E")
            color = 0
        gr.GetXaxis().SetRangeUser(0, 20)
        gr.GetYaxis().SetRangeUser(0, 1) 
        gr.GetYaxis().SetTitle("average sin^{2}")
        gr.SetLineWidth(4)
        gr.SetLineStyle(1)
        #gr.SetMarkerColor(colors[color])
        gr.SetLineColor(colors[color])
        gr.SetFillColorAlpha(colors[color],.4)
        #gr.SetMarkerSize(4)
        gr.Draw(option+" E3")
        line.SetLineColor(colors[color])
        line.SetLineWidth(4)
        lines.append(line)
        legend.AddEntry(gr,"#Delta m^{2} = " + str(m2),"l")
        legend.Draw("SAME")
        deltam2 += 1
        m2 = ms[deltam2]
        color += 1
        if m2 in [6,11,20,-1]:
            for l in lines:
                l.Draw("SAME HIST L")
            PlotTools.Print(AnalysisConfig.PlotPath("SINAVG_TEST_{}".format(deltam2),sideband,"average"))
            lines = []
        
