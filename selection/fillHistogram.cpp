#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>
#include "dk2nu/tree/dk2nu.h" 
#include "dk2nu/tree/dkmeta.h" 
#include "dk2nu/tree/readWeightLocations.h" 

void fillHistogram() {
    // Create a TChain and add multiple files
    TChain tree("dk2nuTree");
    tree.Add("merged.root");
    TChain meta("dkmetaTree");
    meta.Add("merged.root");

    // Create histograms for
    TH1D *hnue_EFrac = new TH1D("nue_energyfraction", "Neutrino Energy Fraction;2E^{CM}_{#nu}/M_{#mu};Entries", 30, 0, 1);
    TH1D *hnue_CosTheta = new TH1D("nue_pcostheta", "Angle Between Neutrino and Muon Polarization;cos(#theta);Entries", 15, -1, 1);
    TH1D *hnue_polarization = new TH1D("nue_distribution", "Muon Decay Neutrino x-sec;#frac{d^{2}N}{dxd#Omega};Entries", 30, 0,.35);

    TH2D *hnue_costheta_v_x = new TH2D("nue_costheta_vs_x", "#nu_{e} cos(#theta) vs x;cos(#theta);x", 30, -1, 1, 30, 0, 1);
    TH2D *hnumu_costheta_v_x = new TH2D("numu_costheta_vs_x", "#nu_{#mu} cos(#theta) vs x;cos(#theta);x", 30, -1, 1, 30, 0, 1);

    TH2D *hnue_x_v_d2n = new TH2D("nue_x_vs_d2n", "#nu_{e} x vs #frac{d^{2}N}{dxd#Omega};2E^{CM}_{#nu}/M_{#mu};#frac{d^{2}N}{dxd#Omega}", 30, 0, 1, 30, 0, .4);
    TH2D *hnumu_x_v_d2n = new TH2D("numu_x_vs_d2n", "#nu_{#mu} x vs #frac{d^{2}N}{dxd#Omega};2E^{CM}_{#nu}/M_{#mu};#frac{d^{2}N}{dxd#Omega}", 30, 0, 1, 30, 0, .4);

    TH1D *hnumu_EFrac = new TH1D("numu_energyfraction", "Neutrino Energy Fraction;2E^{CM}_{#nu}/M_{#mu};Entries", 30, 0, 1);
    TH1D *hnumu_CosTheta = new TH1D("numu_pcostheta", "Angle Between Neutrino and Muon Polarization;cos(#theta);Entries", 15, -1, 1);
    TH1D *hnumu_polarization = new TH1D("numu_distribution", "Muon Decay Neutrino x-sec;#frac{d^{2}N}{dxd#Omega};Entries", 30, 0, .35);

    TGraph *numu = new TGraph();
    TGraph *nue = new TGraph();

    // Variables to hold the branch data
    bsim::Dk2Nu*  dk2nu  = new bsim::Dk2Nu;
    tree.SetBranchAddress("dk2nu", &dk2nu);

    Long64_t nentries = tree.GetEntries();
    for (Long64_t i=0; i < nentries; ++i ) {
        tree.GetEntry(i);
        int ntype  = dk2nu->decay.ntype;
        int ptype  = dk2nu->decay.ptype;

        if (abs(ptype) == 13) {
            double pdpx = dk2nu->decay.pdpx;
            double pdpy = dk2nu->decay.pdpy;
            double pdpz = dk2nu->decay.pdpz;
            double necm = dk2nu->decay.necm;

            double nimpwt = dk2nu->decay.nimpwt;
            double wgt = dk2nu->nuray.at(0).wgt;

            double npE  = dk2nu->nuray.at(0).E;
            double npx  = dk2nu->nuray.at(0).px;
            double npy  = dk2nu->nuray.at(0).py;
            double npz  = dk2nu->nuray.at(0).pz;
            double mmu = 0.1057;

            double pd2 = pdpx*pdpx + pdpy*pdpy + pdpz*pdpz;
            double pd3 = TMath::Sqrt(pd2);
            double pdE2 = mmu*mmu + pdpx*pdpx + pdpy*pdpy + pdpz*pdpz;
            double pdE = TMath::Sqrt(pdE2);

            TLorentzVector nuvec;
            TLorentzVector muonvec;

            nuvec.SetPxPyPzE(npx,npy,npz,npE);
            muonvec.SetPxPyPzE(pdpx/pdE,pdpy/pdE,pdpz/pdE,pdE/pdE);

            TVector3 parent3vec = muonvec.Vect();
            nuvec.Boost(-parent3vec);

            double theta = 0;
            if (ptype < 0) { theta = nuvec.Angle(-muonvec.Vect()); }
            else if (ptype > 0) { theta = nuvec.Angle(muonvec.Vect()); }

            double costheta = TMath::Cos(theta);
            double x = 2*necm/mmu;

            double E_fraction = 2*necm/mmu;
            double pcostheta = E_fraction/(necm/mmu) - 1;

            if (abs(ntype) == 14) {
                int fact = 0;
                if (ntype > 0) { fact = -1; }
                else if (ntype < 0) { fact = 1; }
                hnumu_costheta_v_x->Fill(costheta,x,wgt);
                double d2n = 2*x*x/(4*TMath::Pi()) * ((3-2*x) + fact*(1-2*x)*costheta); 

                hnumu_x_v_d2n->Fill(x,d2n,wgt*nimpwt);
                hnumu_CosTheta->Fill(costheta,wgt*nimpwt); 
                hnumu_EFrac->Fill(E_fraction,wgt*nimpwt); 
                numu->SetPoint(numu->GetN(),2*necm/mmu,d2n);
                hnumu_polarization->Fill(2*x*x/(4*TMath::Pi()) * ((3-2*x) + fact*(1-2*x)*costheta));
            }
            else if (abs(ntype) == 12) {
                int fact = 0;
                if (ntype > 0) { fact = 1; }
                else if (ntype < 0) { fact = -1; }
                hnue_costheta_v_x->Fill(costheta,x,wgt*nimpwt);
                double d2n = 12*x*x/(4*TMath::Pi()) * ((1-x) + fact*(1-x)*costheta);

                hnue_x_v_d2n->Fill(x,d2n,wgt*nimpwt);
                hnue_CosTheta->Fill(TMath::Cos(theta),wgt); 
                hnue_EFrac->Fill(E_fraction,wgt*nimpwt); 
                nue->SetPoint(nue->GetN(),2*necm/mmu,d2n);
                hnue_polarization->Fill(12*x*x/(4*TMath::Pi()) * ((1-x) + fact*(1-x)*costheta));
            }
        }
    }
    // Create a canvas to draw the histogram
    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1", "efrac Histogram", 800, 600);
    hnumu_EFrac->SetLineColor(kBlue);
    hnue_EFrac->SetLineColor(kRed);
    hnumu_EFrac->SetLineWidth(3);
    hnue_EFrac->SetLineWidth(3);
    hnue_EFrac->Draw("hist");
    hnumu_EFrac->Draw("same hist");

    TLegend *leg1 = new TLegend(0.12, 0.7, 0.3, 0.88);
    leg1->AddEntry(hnumu_EFrac,"#nu_{#mu}","l");
    leg1->AddEntry(hnue_EFrac,"#nu_{e}","l");
    leg1->Draw();
    c1->SaveAs("energy_fraction.png");

    TCanvas *c2 = new TCanvas("c2", "TGraph", 800, 600);
    numu->SetMarkerColorAlpha(kBlue,0.4);
    numu->SetFillColor(kBlue);
    nue->SetMarkerColorAlpha(kRed,0.4);
    nue->SetFillColor(kRed);
    numu->SetMarkerSize(5);
    numu->SetMarkerStyle(7);
    nue->SetMarkerSize(5);
    nue->SetMarkerStyle(7);
    nue->SetTitle("x vs #frac{d^{2}N}{dxd#Omega}");
    nue->GetXaxis()->SetTitle("2E_{#nu}/m_{#mu}");
    nue->GetYaxis()->SetTitle("#frac{d^{2}N}{dxd#Omega}");
    numu->GetXaxis()->SetTitle("2E_{#nu}/m_{#mu}");
    numu->GetYaxis()->SetTitle("#frac{d^{2}N}{dxd#Omega}");
    nue->Draw("ap");
    numu->Draw("sp");

    TLegend *leg2 = new TLegend(0.12, 0.7, 0.3, 0.88);
    leg2->AddEntry(numu,"#nu_{#mu}","f");
    leg2->AddEntry(nue,"#nu_{e}","f");
    leg2->Draw();
    c2->SaveAs("costheta.png");

    TCanvas *c3 = new TCanvas("c3", "pcostheta Histogram", 800, 600);
    hnumu_CosTheta->SetLineColor(kBlue);
    hnue_CosTheta->SetLineColor(kRed);
    hnumu_CosTheta->SetLineWidth(3);
    hnue_CosTheta->SetLineWidth(3);
    hnumu_CosTheta->Draw("hist");
    hnue_CosTheta->Draw("same hist");

    TLegend *leg3 = new TLegend(0.6, 0.7, 0.78, 0.88);
    leg3->AddEntry(hnumu_CosTheta,"#nu_{#mu}","l");
    leg3->AddEntry(hnue_CosTheta,"#nu_{e}","l");
    leg3->Draw();
    c3->SaveAs("pcostheta.png");

    TCanvas *c4 = new TCanvas("c4", "xsec distribution", 800, 600);
    hnumu_polarization->SetLineColor(kBlue);
    hnue_polarization->SetLineColor(kRed);
    hnumu_polarization->SetLineWidth(3);
    hnue_polarization->SetLineWidth(3);
    hnumu_polarization->Draw("hist");
    hnue_polarization->Draw("same hist");

    TLegend *leg4 = new TLegend(0.6, 0.7, 0.78, 0.88);
    leg4->AddEntry(hnumu_polarization,"#nu_{#mu}","l");
    leg4->AddEntry(hnue_polarization,"#nu_{e}","l");
    leg4->Draw();
    c4->SaveAs("neutrino_xsec.png");

    TCanvas *c5 = new TCanvas("c5", "2d_numu", 800, 600);
    hnumu_costheta_v_x->Draw("colz");
    c5->SaveAs("2dhist_numu.png");
    
    TCanvas *c6 = new TCanvas("c6", "2d_nue", 800, 600);
    hnue_costheta_v_x->Draw("colz");
    c6->SaveAs("2dhist_nue.png");
    
    TCanvas *c7 = new TCanvas("c7", "2d_numud2n", 800, 600);
    hnumu_x_v_d2n->Draw("colz");
    c7->SaveAs("2dhist_numud2n.png");
    
    TCanvas *c8 = new TCanvas("c8", "2d_nued2n", 800, 600);
    hnue_x_v_d2n->Draw("colz");
    c8->SaveAs("2dhist_nued2n.png");
}
