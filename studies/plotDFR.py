import ROOT,PlotUtils

f =ROOT.TFile.Open("/minerva/data/users/hsu/nu_e/kin_dist_mcme_NCDIF_col10EelPt_Had.root")

hist_name = "inline_upstream_E_Psi"
energy_index = range(1,6)

c = ROOT.TCanvas("c1","c1",1920,1080)
for i in energy_index:
    h_excess1 = f.Get("{}{}_Excess_High_Inline".format(hist_name,i))
    h_excess2 = f.Get("{}{}_Excess_Low_Inline".format(hist_name,i))
    h_signal = f.Get("{}{}".format(hist_name,i))
    first = 2.5*min(i,4)+ 5*max(0,i-4)
    last = 2.5*min(i+1,4) + 5*max(0,i-3)
    h_signal.SetTitle("{}<Ee(GeV)<{}".format(first,last))
    h_signal.Draw("COLZ")
    c.Print("Signal{}.png".format(i))
    c.Clear()
    h_excess1.Add(h_excess2)
    h_excess1.SetTitle("{}<Ee(GeV)<{}".format(first,last))
    h_excess1.Draw("COLZ")
    c.Print("Excess{}.png".format(i))
    c.Clear()
