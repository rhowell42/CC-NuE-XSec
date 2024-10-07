#don't turn on display. It slow down the processing
import os
if 'DISPLAY' in os.environ:
    del os.environ['DISPLAY']

import ROOT
import math
from tools.EventClassification import EventClassifier
from tools.KinematicsCalculator import KinematicsCalculator
from tools.SystematicsUniverse import GetAllSystematicsUniverses

def MakeHist(chain,e,cate=None):
    hist_correlation = ROOT.TH2D("h1{}".format(e),";prev_dEdX;this dEdX;",50,0,10,50,0,10)
    hist_count = ROOT.TH2D("h2{}".format(e),";N-Plane DS(mm);N-Clusters;NEvents",40,0,40,20,0,10)
    hist_perplane = ROOT.TH3D("h3{}".format(e),";N-Plane DS; Energy(GeV);dE/dX(MeV/cm)",40,0,40,10,0,10,100,0,20)
    hist_LLR = ROOT.TH2D("h4{}".format(e),";Energy(Gev);LLR score; NEvents",10,0,10,40,-20,20)
    #hist_modelLLR = PlotUtils.MnvH2D("h4m{}".format(e),"",30,0,30,100,-10,10)
    hist_mindedx = ROOT.TH1D("h5{}".format(e),";dEdX(MeV/cm);N-hits",50,0,20)
    hist_spread = ROOT.TH2D("h8{}".format(e),"",40,0,80,100,0,1000)
    hist_Exuv = ROOT.TH2D("h6{}".format(e),"",100,-0.5,0.5,40,0,10000)
    hist_Euv = ROOT.TH2D("h7{}".format(e),"",100,-0.5,0.5,40,0,10000)
    #hist_nplane = ROOT.TH2D("h8{}".format(e),"",20,0,20,50,0,100)
    #hist_mindedx_proxy = ROOT.TH1D("h8{}".format(e),"",50,0,5)
    hist_strange0 = ROOT.TH1D("h9{}".format(e),"",50,0,10)
    hist_strange1 = ROOT.TH1D("h10{}".format(e),"",50,0,10)
    hist_strange2 = ROOT.TH1D("h11{}".format(e),"",50,0,10)

    count = 0
    for e in chain:
        if e.n_prongs == 0:
            continue
        
        # Exuv = (e.prong_XViewE[0]-e.prong_UViewE[0]-e.prong_VViewE[0])/(e.prong_XViewE[0]+e.prong_UViewE[0]+e.prong_VViewE[0])
        # Euv = (e.prong_UViewE[0]-e.prong_VViewE[0])/(e.prong_UViewE[0]+e.prong_VViewE[0])
        # Esum = e.prong_XViewE[0]+e.prong_UViewE[0]+e.prong_VViewE[0]

        # hist_Exuv.Fill(Exuv,e.prong_XViewE[0]+e.prong_UViewE[0]+e.prong_VViewE[0])
        # hist_Euv.Fill(Euv,e.prong_XViewE[0]+e.prong_UViewE[0]+e.prong_VViewE[0])

        

        if cate is not None:
            #chain.kin_cal.CalculateKinematics(e)
            eventClassifier.Classify(e)
            if not eventClassifier.is_reco_signal:
                continue
            if e.prong_part_E[0][3]<1500:
                continue

            if e.recoile_passive_idclus > 1500:
                continue

            if cate == "signal":
                #if e.mc_intType!=7:
                if abs(e.mc_incoming) != 12 or e.mc_current!=1 :
                    continue
            if cate == "NC":
                if e.mc_current!=2:# or e.mc_intType!=4:
                # or 211 in map(abs, e.mc_FSPartPDG) or 2212 in map(abs,e.mc_FSPartPDG):
                    continue
            # if abs(Exuv)>0.2 or abs(Euv)>0.3:
            #     continue

            #if e.recoile_passive >1200:
            #    continue
        energy = e.prong_part_E[0][3]/1e3
        #print energy
        dedx = e.prong_dEdXs[0]
        dedx_dz = e.prong_dEdXs_dz[0]
        dedx_dx = e.prong_dEdXs_projection[0]

        dEdX_sum = {}
        vertex_plane = Get_z_plane(e.mc_vtx[2])
        #print dedx.size()
        #print [e.prong_binned_energy_bin_contents[0][i] for i in range(e.prong_binned_energy_bin_contents[0].size())]
        for i in range(min(dedx.size(),dedx_dz.size())):
            try: 
                dEdX_sum[dedx_dz[i]] += dedx[i]
            except KeyError:
                dEdX_sum[dedx_dz[i]] = dedx[i]
        # print "count:",count
        print(e.mc_run, e.mc_subrun, e.mc_nthEvtInFile+1)
        print(e.prong_dEdXMeanFrontTracker[0])
        print(e.prong_TruePID[0], e.prong_part_E[0][3]-e.prong_TrueEnergy[0])
        print(e.prong_axis_vertex[0][2]-e.vtx[2] if e.vtx[2]>0 else -1)
        for key in sorted(dEdX_sum):
            print("("+ repr(key)+ ", dEdX:" + repr(dEdX_sum[key]) + ")")
        input("continue?")
        top_hit_counts = len(dEdX_sum)/3
        hit_count = 0
        LLR = [0,0,0]
        min_dedx = 99999
        model_LLR = 0
        strange_count = [0,0,0]
        #print len(dedx_sum)
        energy_index = min(11,int(math.ceil(energy)))
        for i in range(-10,84*2+1):
            if len(dedx_dz) == 0 or i < dedx_dz[0]:
                continue
            if hit_count>top_hit_counts :
                break
            hit_count += 1
            dplane = i-dedx_dz[0]
            dedx_step = dEdX_sum.setdefault(i,0)
            if dedx_step<0.5 and dedx_step>=0:
                strange_count[0] +=1
            if dedx_step<12.5 and dedx_step>=12:
                strange_count[1] +=1
            if dedx_step<9 and dedx_step>=8.5:
                strange_count[2] +=1
            hist_mindedx.Fill(dedx_step,1)
            # if (e.ev_run==10069 and  e.ev_subrun==44 and e.ev_gate==455):
            #     print "count:",count
            #     print e.prong_part_E[0][3]
            #     for key in sorted(dEdX_sum):
            #         print "("+ repr(key)+ ", dEdX:" + repr(dEdX_sum[key]) + ")"
            #     raw_input("continue?")
            #if dedx_sum[j] == 0:
            #    print "warning: zero dedx"
            #LLR *= ROOT.TMath.Landau(dedx_sum[j],3.75,1.46,True)/ROOT.TMath.Landau(dedx_sum[j],1.7,0.506,True)
            #dplane = plane - dedx_dz_trun[0]
            # prob_e = f.Get("elikelihood{}_{}".format(dplane,energy_index))
            # prob_pi = f.Get("glikelihood{}_{}".format(dplane,energy_index))
            # if not prob_e or not prob_pi:
            #     print dplane,energy_index
            # pe = max(prob_e.GetBinContent(prob_e.FindBin(dedx_sum[j])),0.0001)
            # ppi = max(prob_pi.GetBinContent(prob_pi.FindBin(dedx_sum[j])),0.0001)
            # k = 2 if dedx_dz_trun[j] % 2 == 1 else (dedx_dz_trun[j] % 4)/2
            # LLR[k]+=-math.log(pe)+math.log(ppi)

            # model_LLR += -math.log(ROOT.TMath.Landau(dedx_sum[j],1.7,0.506,True))+math.log(ROOT.TMath.Landau(dedx_sum[j],3.75,1.46,True))
            # hist_modelLLR.Fill(j,max(-10,min(10,model_LLR)),1)
            # if (j>=3):
            #     ave_dedx = sum(dedx_sum[j-3:j+1])/4
            #     min_dedx = min(ave_dedx,min_dedx)
            #     hist_mindedx.Fill(j-3,min_dedx,1)
            #hist_spread.Fill(j,min(999.9,dedx_spread[j]),1)
            hist_perplane.Fill(dplane,energy,dedx_step)
            # if j==15:
            #     if abs(dedx_sum[j] -2) < 0.5:
            #         print e.mc_run, e.mc_subrun, e.mc_nthEvtInFile+1
            #         print e.prong_axis_vertex[0][2]
            #         print [i for i in dedx_sum]
            #         print [e.prong_binned_energy_bin_contents[0][i] for i in range(e.prong_binned_energy_bin_contents[0].size())]
            #         print Exuv, Euv
            #         raw_input("prese ENTER to continue")
           # j+=1
        # if len(dedx_dz_trun) == 0:
        #      print e.mc_run, e.mc_subrun, e.mc_nthEvtInFile+1
        #      print e.prong_axis_vertex[0][2]
        #      print [i for i in dedx]
        #      print [i for i in dedx_dz]
        #      continue
        # if (dplane<15 and dedx_dz_trun[-1]< 84*2+1):
        #     for j in range(dplane+1, min(15, 84*2+1-dedx_dz_trun[0])):
        #         dplane = j
        #         hist_perplane.Fill(dplane,energy,0)
        #         prob_e = f.Get("elikelihood{}_{}".format(dplane,energy_index))
        #         prob_pi = f.Get("glikelihood{}_{}".format(dplane,energy_index))
        #         pe = max(prob_e.GetBinContent(1),0.001)
        #         ppi = max(prob_pi.GetBinContent(1),0.001)
        #         k = 2 if (dedx_dz_trun[0]+j) % 2 == 1 else ((dedx_dz_trun[0]+j) % 4)/2
        #         LLR[k]+=-math.log(pe)+math.log(ppi)

        # if len(dedx_dz_trun)>34:
        #     # that we can do something fancy:
        #     for k in range(0,16):
        #         prob_pi = f.Get("glikelihood{}_{}".format(k,energy_index))
        #         for l in range(len(startPilike)):
        #             ppi = max(prob_pi.GetBinContent(prob_pi.FindBin(dedx_sum[l+k])),0.001)
        #             startPilike[l]+=-math.log(ppi)
        #     # print "hey"
        #     # print dedx_sum
        #     # print startPilike
        #     min_value = min(startPilike)
        #     #print min_value
        #     hist_minPilike.Fill(min_value)
        #     hist_nplane.Fill(startPilike.index(min_value),min_value)
        #print sum(LLR)
        #score = sum(LLR)
        #hist_LLR.Fill(energy_index,max(-20,min(19.9,score)),1)
        #hist_KNN.Fill(min_LLR,max(min(dedx_spread[0::2]),min(dedx_spread[1::2])))
        #dtheta = e.prong_axis_vector[0][0]-e.prong_TrueTheta[0]
        #hist_dtheta.Fill(dtheta)
        #hist_KNN.Fill(e.prong_TrueTrajFracInProng[0],min(abs(e.prong_TruePID[0]),120))
        hist_strange0.Fill(min(strange_count[0],9))
        hist_strange1.Fill(min(strange_count[1],9))
        hist_strange2.Fill(min(strange_count[2],9))
              # print e.mc_run, e.mc_subrun, e.mc_nthEvtInFile+1
              # print dEdX_sum
              # print [e.prong_binned_energy_bin_contents[0][i] for i in range(e.prong_binned_energy_bin_contents[0].size())]
              # #print Exuv, Euv
              # print e.recoile_passive
              # raw_input("prese ENTER to continue")
        count+=1
        # print e.mc_run, e.mc_subrun, e.mc_nthEvtInFile+1
        # print [i for i in dedx]
        # print [i for i in dedx_dz]
        # raw_input("press enter to continue")
        #break

    print(count)
    #hist.Divide(hist,hist_count)
    hist_count.Scale(1.0/count)
    #hist_perplane.Scale(1.0/count)

    return hist_LLR,hist_mindedx,hist_perplane,hist_strange0,hist_strange1,hist_strange2

def MakeChain(path):
    chain = ROOT.TChain("NuECCQE");
    chain.Add(path)
    return chain

def calMCC(hist_true,hist_false,cut):
    TP = hist_true.Integral(1,cut)
    FN = hist_true.Integral(cut,101)
    FP = hist_false.Integral(1,cut)
    TN = hist_false.Integral(cut,101)
    try:
        MCC = (TP*TN-FP*FN) / (math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
        EFF = TP/(TP+FN)
        PUR = TP/(TP+FP)
    except ZeroDivisionError:
        print(TP,FN,FP,TN)
        return None,None,None
    print("MCC is {}".format( MCC))
    #print "efficiency is {}".format(EFF)
    #print "purity is {}".format(PUR)
    return MCC,EFF,PUR

def Get_z_plane(z):
    with open(os.environ["PLOTUTILSROOT"]+"/data/Minerva_Planes.txt") as f:
        for line in f.readlines():
            this_plane = line.split()
            if (abs(z -float(this_plane[5]))<25 or z<float(this_plane[5])):
                return int(this_plane[3])*2+int(this_plane[4])-1
    return None

if __name__ == "__main__":
    
    path_e = "/pnfs/minerva/persistent/users/hsu/electronPC-test/grid/central_value/minerva/ana/v21r1p1/00/00/00/01/*.root"
    path_g = "/pnfs/minerva/persistent/users/hsu/photonPCAna-oe2/grid/central_value/minerva/ana/v21r1p1/00/00/00/01/*.root"

    nue_IA=True
    
    if nue_IA:
        f = ROOT.TFile.Open("hist_IA_tmp.root","RECREATE")
        chain =  ROOT.TChain("NuECCQE");
        #for i in range(1):
        print(chain.Add("/pnfs/minerva/persistent/users/hsu/NuECCQE-v3-mc/me1B_merged/NuECCQE_mc_AnaTuple_*.root"))
        #print chain.Add("/pnfs/minerva/persistent/users/hsu/short_test/grid/central_value/minerva/ana/v21r1p1/00/11/10/00/*.root")
        #print chain.Add("/exp/minerva/data/users/hsu/test/nogrid/central_value/minerva/ana/v21r1p1/00/11/10/01/SIM_minerva_00111001_Subruns_0122_NuECCQE_Ana_Tuple_v21r1p1.root")
        kin_cal = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False, calc_true = True)
        eventClassifier = EventClassifier(classifiers=["Reco"], use_kin_cuts=False,use_sideband=[])
        chain.kin_cal = kin_cal
        chain.ShortName = lambda : "cv"
        hist_LLR_e,hist_dedx_e,hist_perplane_e,hist_count_e,hist_Exuv_e,hist_Euv_e = MakeHist(chain,"e","signal")
        hist_LLR_g,hist_dedx_g,hist_perplane_g,hist_count_g,hist_Exuv_g,hist_Euv_g= MakeHist(chain,"g","NC")
    else :
        input("will overwrite hist.root, are you sure?")
        f = ROOT.TFile.Open("hist_PC.root","RECREATE")
        chain_e = MakeChain(path_e)
        chain_pi = MakeChain(path_g)
        hist_LLR_e,hist_dedx_e,hist_perplane_e, hist_count_e,hist_Exuv_e,hist_Euv_e= MakeHist(chain_e,"e")
        #hist_LLR_g,hist_dedx_g,hist_perplane_g, hist_count_g, hist_Exuv_g,hist_Euv_g = MakeHist(chain_pi,"g")

    # f2 = ROOT.TFile.Open("plot_bal.root","RECREATE")
    hist_LLR_e.Write()
    hist_dedx_e.Write()
    hist_perplane_e.Write()
    hist_count_e.Write()
    hist_Exuv_e.Write()
    hist_Euv_e.Write()
    

    hist_LLR_g.Write()
    hist_dedx_g.Write()
    hist_perplane_g.Write()
    hist_count_g.Write()
    hist_Exuv_g.Write()
    hist_Euv_g.Write()
    # hist_perplane_e.Write()
    # hist_perplane_g.Write()

    # p1 = ROOT.TGraphErrors()
    # p2 = ROOT.TGraphErrors()
    # print "min dedx mcc"
    # for i in range(1,51):
    #     mcc,eff,pur = calMCC(hist_dedx_e,hist_dedx_g,i)
    #     if mcc is not None:
    #         p1.SetPoint(i,eff,pur)

    # print "LLR"
    # for i in range (1,41):
    #     mcc,eff,pur = calMCC(hist_LLR_e.ProjectionY("1"),hist_LLR_g.ProjectionY("2"),i)
    #     if mcc is not None:
    #         p2.SetPoint(i,eff,pur)
  

    # p1.Write("tg1")
    # p2.Write("tg2")

    for i in range(1,80):
        for j in range(1,12):
            hist_temp = hist_perplane_e.ProjectionZ("elikelihood{}_{}".format(i-1,j),i,i,j,j)
            if hist_temp.Integral(0,101)!= 0:
                hist_temp.Scale(1.0/hist_temp.Integral(0,101))
                hist_temp.Write()
            hist_temp = hist_perplane_g.ProjectionZ("glikelihood{}_{}".format(i-1,j),i,i,j,j)
            if hist_temp.Integral(0,101)!= 0:
                hist_temp.Scale(1.0/hist_temp.Integral(0,101))
                hist_temp.Write()

    for i in range(0,52):
        hist_temp = hist_LLR_e.ProjectionY("ecorrelation{}".format(i),i,i)
        if hist_temp.Integral(0,51)!= 0:
            hist_temp.Scale(1.0/hist_temp.Integral(0,51))
            hist_temp.Write()
        hist_temp = hist_LLR_g.ProjectionY("gcorrelation{}".format(i),i,i)
        if hist_temp.Integral(0,51)!= 0:
            hist_temp.Scale(1.0/hist_temp.Integral(0,51))
            hist_temp.Write()

    c= ROOT.TCanvas("c","c",1024,768)
    #for i in range(1,12):
    #hist_perplane_e.GetYaxis().SetRange(i,i);
    hist_perplane_e.Project3D("zx").Draw("COLZ")
    c.Print("e_profile{}.png".format(0))
    #hist_perplane_g.GetYaxis().SetRange(i,i)
    hist_perplane_g.Project3D("zx").Draw("COLZ")
    c.Print("g_profile{}.png".format(0))

    #calMCC(hist_LLR_e.ProjectionY("",i,i),hist_LLR_g.ProjectionY("",i,i),5)

    #print "total"
    #print "baseline"
    #calMCC(hist_dedx_e,hist_dedx_g,24)
    #calMCC(hist_dedx_e.ProjectionY("",20,20),hist_dedx_g.ProjectionY("",20,20),24)

    # for i in range(1,12):
    hist_LLR_e.SetLineColor(ROOT.kRed)
    hist_LLR_e.ProjectionY("").Draw("")
    hist_LLR_g.ProjectionY("").Draw("SAME")
    legend = ROOT.TLegend(0.3,0.7,0.7,0.9)
    legend.AddEntry(hist_LLR_e,"CCNuE")
    legend.AddEntry(hist_LLR_g,"NC")
    legend.Draw()
    c.Print("hist_LLR.png")
    #hist_perplane_e.Project3D("zx").Draw("COLZ")
    #c.Print("hist_nue_signal1.png")
    #hist_LLR_g.Draw("COLZ")
    #c.Print("hist_g.png")
    #hist_perplane_g.Project3D("zx").Draw("COLZ")
    #c.Print("hist_nue_NC1.png")
    #hist_LLR_g.Draw("COLZ")
    #c.Print("NCshowerwidth.png")
    #hist_perplane_g.Draw("COLZ")
    #c.Print("NCcount.png")
    #hist_dedx_e.Draw("COLZ")
    # c.Print("10deg_spread.png")
    f.Close()
