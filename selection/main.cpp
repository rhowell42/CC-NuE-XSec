#include <TFile.h>
#include <TTree.h>
#include <TEntryList.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <iostream>
#include <TH1D.h>

void filterEntriesInDirectory(std::string playlist, TH1D* hist) {
    std::string combine = "/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me"+playlist+"_DualVertex_p4/";
    const char* dirPath = combine.c_str();
    TSystemDirectory dir("rootDir",dirPath);
    TList *files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "Error opening directory" << std::endl;
        return 0;
    }
    TIter next(files);
    TSystemFile *file;
    while ((file = (TSystemFile*)next())) {
        TString fileName = file->GetName();
        if (!file->IsDirectory() && fileName.EndsWith(".root")) {
            TString filePath = TString(dirPath) + "/" + fileName;
            TFile *rootFile = TFile::Open(filePath);
            // Open the ROOT file
            if (!file || file->IsZombie()) {
                std::cerr << "Error opening file" << std::endl;
                return 0;
            }
            std::cout << "Filling from " << filePath << std::endl;


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
            Double_t vtx[3];
            tree->SetBranchAddress("mc_vtx", &vtx);
            Double_t parent_vtx[3];
            tree->SetBranchAddress("mc_fr_nuParentDecVtx", &parent_vtx);

            // Create an entry list
            tree->Draw(">>entryList", "mc_intType == 6", "entrylist");
            TEntryList *entryList = (TEntryList*)gDirectory->Get("entryList");
            if (!entryList) {
                std::cerr << "Error creating entry list" << std::endl;
                return 0;
            }

            // Set the entry list for the tree
            tree->SetEntryList(entryList);

            // Loop over the entries in the entry list
            Long64_t nEntries = entryList->GetN();
            for (Long64_t i = 0; i < nEntries; ++i) {
                Long64_t entry = entryList->GetEntry(i);
                tree->GetEntry(entry);
                Double_t l_over_e = (.9825+vtx[2]/1e6 - parent_vtx[2]/1e6)/(energy/1000);
                hist->Fill(l_over_e);
                // Add your processing code here
            }
            // Cleanup
            delete entryList;
            rootFile->Close();
        }
    }
    return 0;
}


int main() {
    // Replace with the path to your directory containing ROOT files and your TTree name
    TH1D *hist = new TH1D("imd_scattering_template", "Histogram of Another Branch", 30, 0, 0.5); // Adjust binning and range as needed
    filterEntriesInDirectory("5A", hist);
    filterEntriesInDirectory("6A", hist);
    filterEntriesInDirectory("6B", hist);
    
    TFile outFile("RHC_imd_scattering_template.root", "RECREATE");
    hist->Write();
    outFile.Close();
    delete hist;
    return 0;
}

