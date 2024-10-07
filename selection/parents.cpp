#include <TFile.h>
#include <TTree.h>
#include <TEntryList.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <iostream>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>

void filterEntriesInDirectory(std::string playlist, TH1D* hist) {
    std::string combine = "/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me"+playlist+"_DualVertex_p4/";
    const char* dirPath = combine.c_str();
    TSystemDirectory dir("rootDir",dirPath);
    TList *files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "Error opening directory" << std::endl;
        return;
    }
    TIter next(files);
    TSystemFile *file;

    std::vector<Double_t> muonE;
    std::vector<Double_t> kaonE;
    std::vector<Double_t> pionE;
    std::vector<Double_t> otherE;

    std::vector<Double_t> muonL;
    std::vector<Double_t> kaonL;
    std::vector<Double_t> pionL;
    std::vector<Double_t> otherL;

    int filecount = 0;

    while ((file = (TSystemFile*)next())) {
        TString fileName = file->GetName();
        if (!file->IsDirectory() && fileName.EndsWith(".root")) {
            TString filePath = TString(dirPath) + "/" + fileName;
            TFile *rootFile = TFile::Open(filePath);
            // Open the ROOT file
            if (!file || file->IsZombie()) {
                std::cerr << "Error opening file" << std::endl;
                return ;
            }
            std::cout << "Filling from " << filePath << std::endl;
            filecount++;

            // Get the TTree from the file
            TTree *tree;
            rootFile->GetObject("MasterAnaDev", tree);
            if (!tree) {
                std::cerr << "Error getting TTree from file: " << filePath << std::endl;
                rootFile->Close();
                continue;
            }

            Double_t energy;
            tree->SetBranchAddress("mc_incomingE", &energy);
            Double_t vtx[4];
            tree->SetBranchAddress("mc_vtx", &vtx);
            Double_t parent_vtx[4];
            tree->SetBranchAddress("mc_fr_nuParentDecVtx", &parent_vtx);
            Int_t parent_pdg;
            tree->SetBranchAddress("mc_fr_nuParentID", &parent_pdg);

            // Create an entry list
            tree->Draw(">>entryList", "mc_incoming == 12 || mc_incoming == -12", "entrylist");
            TEntryList *entryList = (TEntryList*)gDirectory->Get("entryList");
            if (!entryList) {
                std::cerr << "Error creating entry list" << std::endl;
                return;
            }

            // Set the entry list for the tree
            tree->SetEntryList(entryList);

            // Loop over the entries in the entry list
            Long64_t nEntries = entryList->GetN();
            for (Long64_t i = 0; i < nEntries; ++i) {
                Long64_t entry = entryList->GetEntry(i);
                tree->GetEntry(entry);
                Double_t l = (.9825+vtx[2]/1e6 - parent_vtx[2]/1e6);
                Double_t e = energy/1000;
                if (abs(parent_pdg) == 211) {
                        pionE.push_back(e);
                        pionL.push_back(l);
                }
                else if (abs(parent_pdg) == 130 || abs(parent_pdg) == 310 || abs(parent_pdg) == 311 || abs(parent_pdg) == 321) {
                        kaonE.push_back(e);
                        kaonL.push_back(l);
                }
                else if (abs(parent_pdg) == 13) {
                        muonE.push_back(e);
                        muonL.push_back(l);
                }
                else {
                        std::cout<<parent_pdg<<std::endl;
                }
            }
            // Cleanup
            delete entryList;
            rootFile->Close();
        }
        if (filecount == 10) {
            break;
        }
    }

    TCanvas *c1 = new TCanvas("c1","c1");
    c1->SetRightMargin(0.2);
    TLegend *leg = new TLegend(0.81, 0.38, 0.99, 0.62);

    double me_arr[muonE.size()];
    std::copy(muonE.begin(), muonE.end(), me_arr);
    double ml_arr[muonE.size()];
    std::copy(muonL.begin(), muonL.end(), ml_arr);

    double ke_arr[kaonE.size()];
    std::copy(kaonE.begin(), kaonE.end(), ke_arr);
    double kl_arr[kaonL.size()];
    std::copy(kaonL.begin(), kaonL.end(), kl_arr);

    double pe_arr[pionE.size()];
    std::copy(pionE.begin(), pionE.end(), pe_arr);
    double pl_arr[pionL.size()];
    std::copy(pionL.begin(), pionL.end(), pl_arr);

    TGraph *gpion = new TGraph(pionE.size(),pe_arr,pl_arr);
    TGraph *gmuon = new TGraph(muonE.size(),me_arr,ml_arr);
    TGraph *gkaon = new TGraph(kaonE.size(),ke_arr,kl_arr);

    gmuon->SetMarkerColorAlpha(kRed,1);
    gpion->SetMarkerColorAlpha(kBlue,1);
    gkaon->SetMarkerColorAlpha(kGreen,1);

    gmuon->SetMarkerStyle(kPlus);
    gpion->SetMarkerStyle(kCircle);
    gkaon->SetMarkerStyle(kMultiply);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gpion);
    mg->Add(gmuon);
    mg->Add(gkaon);
    mg->GetXaxis()->SetLimits(0,20);
    mg->GetXaxis()->SetTitle("True Neutrino Energy (GeV)");
    mg->GetYaxis()->SetTitle("True Neutrino Length (km)");
    mg->GetHistogram()->SetTitle("Electron Neutrino L/E");
    mg->Draw("ap");
    leg->AddEntry(gpion,"Pion","P");
    leg->AddEntry(gmuon,"Muon","P");
    leg->AddEntry(gkaon,"Kaon","P");
    leg->Draw();
    c1->Print("rhc_muonneutrino.png");
    
}


int parents() {
    // Replace with the path to your directory containing ROOT files and your TTree name
    TH1D *hist = new TH1D("imd_scattering_template", "Histogram of Another Branch", 30, 0, 0.5); // Adjust binning and range as needed
    //filterEntriesInDirectory("5A", hist);
    filterEntriesInDirectory("6A", hist);
    //filterEntriesInDirectory("6B", hist);
    
    TFile outFile("RHC_nue_parentdist.root", "RECREATE");
    hist->Write();
    outFile.Close();
    delete hist;
    return 0;
}

