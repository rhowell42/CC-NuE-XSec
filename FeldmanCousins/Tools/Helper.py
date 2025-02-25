import ROOT
import numpy as np
import ctypes

legend_text_size = .025

def TMatrix_to_Numpy(matrix):
    rows = matrix.GetNrows()
    cols = matrix.GetNcols()
    arr = np.empty((rows, cols), dtype=np.float64)
    for i in range(rows):
        for j in range(cols):
            arr[i, j] = matrix(i, j)
    return arr

def PlotWithRatio(MNVPLOTTER,plotName,h_cv,hists=[],titles=[],colors=[],styles=[]):
    canvas = ROOT.TCanvas()
    margin = .12
    bottomFraction = .2
    top = ROOT.TPad("DATAMC", "DATAMC", 0, bottomFraction, 1, 1)
    bottom = ROOT.TPad("Ratio", "Ratio", 0, 0, 1, bottomFraction+margin)

    top.Draw()
    bottom.Draw()

    top.cd()
    
    h_cv.Draw("hist")

    null = h_cv.GetCVHistoWithError()
    null.SetLineColor(ROOT.kBlack)
    null.SetLineWidth(2)
    null.SetMarkerStyle(0)
    null.SetFillColorAlpha(ROOT.kGray+1, 0.4)
    null.Draw("E2 SAME")

    leg = ROOT.TLegend(.4,.3)
    leg.SetTextSize(legend_text_size);
    leg.AddEntry(h_cv,"CV Prediction","l")

    for i,hist in enumerate(hists):
        if i < len(titles):
            hist.SetTitle(titles[i])
        if i < len(colors):
                hist.SetLineColor(colors[i])
        if i < len(styles):
            hist.SetLineStyle(styles[i])
        hist.Draw("hist same")
        leg.AddEntry(hist,hist.GetTitle(),"l")

    leg.Draw()

    nullRatio =  h_cv.Clone()
    nullRatio.Divide(nullRatio,h_cv)

    bottom.cd()
    bottom.SetTopMargin(0)
    bottom.SetBottomMargin(0.3)

    nullErrors = h_cv.GetTotalError(False, True, False) #The second "true" makes this fractional error, the third "true" makes this cov area normalized
    for whichBin in range(0, nullErrors.GetXaxis().GetNbins()+1): 
        nullErrors.SetBinError(whichBin, max(nullErrors.GetBinContent(whichBin), 1e-9))
        nullErrors.SetBinContent(whichBin, 1)

    nullRatio.SetLineColor(ROOT.kBlack)
    nullRatio.SetLineWidth(3)

    #Error envelope for the MC
    nullErrors.SetLineWidth(0)
    nullErrors.SetMarkerStyle(0)
    nullErrors.SetFillColorAlpha(ROOT.kGray+1, 0.4)
    nullErrors.GetYaxis().SetTitle("#splitline{Ratio to CV}{Prediction}")
    RatioAxis(nullErrors,MNVPLOTTER)
    nullErrors.SetMinimum(0.4)
    nullErrors.SetMaximum(1.6)
    nullErrors.Draw("E2")

    #Draw the data ratios
    nullRatio.Draw("same")

    for hist in hists:
        i_ratio = hist.Clone()
        i_ratio.Divide(i_ratio,h_cv)
        i_ratio.DrawClone("same hist")

    straightLine = nullErrors.Clone()
    straightLine.SetLineColor(ROOT.kBlack)
    straightLine.SetLineWidth(2)
    straightLine.SetFillColor(0)
    straightLine.Draw("HIST SAME")
    canvas.Print(plotName)

def RatioAxis(hist,MNVPLOTTER):
    #MNVPLOTTER.ApplyAxisStyle(hist)

    hist.SetTitle("")
    hist.GetYaxis().SetNdivisions(304) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
    hist.GetYaxis().SetTitleOffset(1.5)
    hist.GetYaxis().SetTitleFont(43)
    hist.GetYaxis().SetTitleSize(26)
    hist.GetYaxis().SetLabelFont(43)
    hist.GetYaxis().SetLabelSize(34)
    hist.GetYaxis().CenterTitle()
    hist.GetXaxis().SetNdivisions(510) #5 minor divisions between 9 major divisions.  I'm trying to match a specific paper here.
    hist.GetXaxis().SetTitleOffset(2.5)
    hist.GetXaxis().SetTitleFont(43)
    hist.GetXaxis().SetTitleSize(26)
    hist.GetXaxis().SetLabelFont(43)
    hist.GetXaxis().SetLabelSize(34)
    hist.GetXaxis().CenterTitle()

    #hist.GetYaxis().SetLabelSize(.13)
    #hist.GetYaxis().SetTitleSize(0.1)
    #hist.GetYaxis().SetTitleOffset(0.6)
    #hist.GetYaxis().SetNdivisions(505) #5 minor divisions between 5 major divisions.  I'm trying to match a specific paper here.
    #hist.GetXaxis().SetTitleSize(0.16)
    #hist.GetXaxis().SetTitleOffset(0.9)
    #hist.GetXaxis().SetLabelSize(.15)

def CanvasPartition(C, Nx, Ny, lMargin, rMargin, bMargin, tMargin):
    vSpacing = 0.0
    vStep = (1. - bMargin - tMargin - (Ny - 1) * vSpacing) / Ny

    hSpacing = 0.0
    hStep = (1. - lMargin - rMargin - (Nx - 1) * hSpacing) / Nx

    vposd, vposu, vmard, vmaru, vfactor = 0, 0, 0, 0, 0
    hposl, hposr, hmarl, hmarr, hfactor = 0, 0, 0, 0, 0

    for i in range(Nx):
        if (i == 0):
            hposl = 0.0
            hposr = lMargin + hStep
            hfactor = hposr - hposl
            hmarl = lMargin / hfactor
            hmarr = 0.0
        elif (i == Nx - 1):
            hposl = hposr + hSpacing
            hposr = hposl + hStep + rMargin
            hfactor = hposr - hposl
            hmarl = 0.0
            hmarr = rMargin / (hposr - hposl)
        else:
            hposl = hposr + hSpacing
            hposr = hposl + hStep
            hfactor = hposr - hposl
            hmarl = 0.0
            hmarr = 0.0

        for j in range(Ny):
            if (j == 0):
                vposd = 0.0
                vposu = bMargin + vStep
                vfactor = vposu - vposd
                vmard = bMargin / vfactor
                vmaru = 0.0
            elif (j == Ny - 1):
                vposd = vposu + vSpacing
                vposu = vposd + vStep + tMargin
                vfactor = vposu - vposd
                vmard = 0.0
                vmaru = tMargin / (vposu - vposd)
            else:
                vposd = vposu + vSpacing
                vposu = vposd + vStep
                vfactor = vposu - vposd
                vmard = 0.0
                vmaru = 0.0

            C.cd(0)

            name = "pad_{}_{}".format(i, j)
            pad = ROOT.TPad(name, "", hposl, vposd, hposr, vposu)
            pad.SetLeftMargin(hmarl)
            pad.SetRightMargin(hmarr)
            pad.SetBottomMargin(vmard)
            pad.SetTopMargin(vmaru)

            pad.SetFrameBorderMode(0)
            pad.SetBorderMode(0)
            pad.SetBorderSize(0)

            pad.Draw()
 
def XtoPad(x):
   xl = ctypes.c_double(0)
   yl = ctypes.c_double(0)
   xu = ctypes.c_double(0)
   yu = ctypes.c_double(0)
   ROOT.gPad.GetPadPar(xl, yl, xu, yu)
   pw = ctypes.c_double(xu.value - xl.value)
   lm = ROOT.gPad.GetLeftMargin()
   rm = ROOT.gPad.GetRightMargin()
   fw = pw.value - pw.value * lm - pw.value * rm
   return (x * fw + pw.value * lm) / pw.value
 
def YtoPad(y):
   xl = ctypes.c_double(0)
   yl = ctypes.c_double(0)
   xu = ctypes.c_double(0)
   yu = ctypes.c_double(0)
   ROOT.gPad.GetPadPar(xl, yl, xu, yu)
   ph = ctypes.c_double(yu.value - yl.value)
   tm = ROOT.gPad.GetTopMargin()
   bm = ROOT.gPad.GetBottomMargin()
   fh = ph.value - ph.value * bm - ph.value * tm
   return (y * fh + bm * ph.value) / ph.value
