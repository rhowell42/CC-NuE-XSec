import ROOT
import numpy as np
import ctypes
import json

legend_text_size = .025

def plot_side_by_side(hists,narrow_pads=None, narrow_factor=0.5):
    if not hists:
        return

    n = len(hists)
    if narrow_pads is None:
        narrow_pads = []

    # Compute global min/max for log scale
    global_max = max(h.GetMaximum() for h in hists)
    global_min = min(h.GetBinContent(b)
                     for h in hists
                     for b in range(1, h.GetNbinsX()+1)
                     if h.GetBinContent(b) > 0)

    # Desired widths before margin adjustments
    base_widths = []
    for i in range(1, n+1):
        w = narrow_factor if i in narrow_pads else 1.0
        base_widths.append(w)
    total_base_width = sum(base_widths)

    # Margins (in pad coordinates)
    left_margin = 0.14
    right_margin = 0.05
    bottom_margin = 0.14
    top_margin = 0.05

    # Convert to total drawable widths so first pad isn't visually smaller
    drawable_widths = []
    for i in range(1, n+1):
        # Pad width in canvas space
        pad_width = base_widths[i-1] / total_base_width
        # Drawable width = pad width * (1 - left_margin - right_margin)
        lm = left_margin if i == 1 else 0.0
        rm = right_margin if i == n else 0.0
        drawable_widths.append(pad_width * (1 - lm - rm))

    # Make all non-narrow pads have same drawable width
    max_drawable = max(drawable_widths[i-1] for i in range(1, n+1) if i not in narrow_pads)
    scale_factor = max_drawable / drawable_widths[0] if 1 not in narrow_pads else 1.0

    # Adjust base widths to equalize drawable space
    adjusted_widths = []
    for i in range(1, n+1):
        lm = left_margin if i == 1 else 0.0
        rm = right_margin if i == n else 0.0
        drawable_target = max_drawable if i not in narrow_pads else max_drawable * narrow_factor
        pad_width = drawable_target / (1 - lm - rm)
        adjusted_widths.append(pad_width)

    # Normalize adjusted widths to [0,1] total
    total_adjusted = sum(adjusted_widths)
    norm_widths = [w / total_adjusted for w in adjusted_widths]

    # Create canvas
    c = ROOT.TCanvas("c", "Custom Pads", 1200, 400)

    # Place pads with margins
    x_start = 0.0
    pads = []
    for i, wnorm in enumerate(norm_widths, start=1):
        x_end = x_start + wnorm
        pad = ROOT.TPad(f"pad{i}", "", x_start, 0, x_end, 1.0)
        pad.SetLogy(True)
        pad.SetBottomMargin(bottom_margin)
        pad.SetTopMargin(top_margin)
        pad.SetLeftMargin(left_margin if i == 1 else 0.0)
        pad.SetRightMargin(right_margin if i == n else 0.0)
        pad.Draw()
        pads.append(pad)
        x_start = x_end

    # Uniform label/tick label sizes in **canvas coordinates**
    # Scale label sizes by inverse pad size so they look identical
    target_label_size = 0.05  # fraction of canvas height
    target_title_size = 0.06  # fraction of canvas height

    for i, (pad, h) in enumerate(zip(pads, hists), start=1):
        pad.cd()
        h.SetMaximum(global_max * 1.2)
        h.SetMinimum(global_min / 2.0)

        pw = pad.GetWw()  # pad width in pixels
        ph = pad.GetWh()  # pad height in pixels
        cw = c.GetWw()    # canvas width in pixels
        ch = c.GetWh()    # canvas height in pixels

        # Adjust sizes to be consistent across pads
        x_scale = cw / pw
        y_scale = ch / ph
        h.GetXaxis().SetLabelSize(target_label_size * y_scale)
        h.GetXaxis().SetTitleSize(target_title_size * y_scale)
        h.GetYaxis().SetLabelSize(target_label_size * x_scale)
        h.GetYaxis().SetTitleSize(target_title_size * x_scale)

        if i in narrow_pads:
            h.GetXaxis().SetNdivisions(2,5,0,ROOT.kFALSE)
        else:
            h.GetXaxis().SetNdivisions(4,5,0,ROOT.kFALSE)

        if i != 1:
            h.GetYaxis().SetLabelSize(0)
            h.GetYaxis().SetTitleSize(0)

        h.Draw("HIST")

    c.Update()
    return c

def GetSliceIndices(fname,exclude,keys):
    bin_config = {}
    with open(fname, "r") as file:
        bin_config = json.load(file)

    sliceInds = []
    for h in keys:
        if h not in bin_config.keys():
            continue
        if not checkRemove(exclude,h):
            sliceInds.extend(list(range(bin_config[h]["start"],bin_config[h]["end"]+1)))
    return(sliceInds)

def slicer(arr,inds,axis=None):
    if arr.ndim == 1: # 1D array
        return arr[inds]
    elif arr.ndim == 2: # 2D array
        if axis == 1:
            return arr[:,inds]
        elif axis == 0:
            return arr[inds,:]
        else:
            ret = arr[:,inds]
            ret = ret[inds,:]
            return ret

def checkRemove(exclude,h):
    if "fhc" in exclude:
        if "fhc" in h and "selection" in h:
            return(True)
    if "rhc" in exclude:
        if "rhc" in h and "selection" in h:
            return(True)
    if "numu" in exclude and "numu" in h:
        return(True)
    if "nue" in exclude and "nue" in h:
        return(True)
    if "elastic" in exclude and "elastic" in h:
        return(True)
    if "imd" in exclude and "imd" in h:
        return(True)
    if "ratio" in exclude and "ratio" in h:
        return(True)

    return(False)

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
    nullErrors.GetYaxis().SetTitle("Ratio to CV")
    RatioAxis(nullErrors,MNVPLOTTER)
    nullErrors.SetMinimum(0.7)
    nullErrors.SetMaximum(1.3)
    nullErrors.GetXaxis().SetTitle("Neutrino Energy")
    nullErrors.Draw("E2")

    #Draw the data ratios
    nullRatio.Draw("same hist")

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
