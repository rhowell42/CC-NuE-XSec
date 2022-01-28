#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include <string.h>


void mergeMultipleStategies() {
  TH1::AddDirectory(false);
  std::string file_path = "/minerva/data/users/hsu/nu_e/";
  std::vector<std::string> hist_names;
  hist_names.push_back("Eavail_Lepton_Pt_bkg_unfolding");
  hist_names.push_back("Eavail_q3_bkg_unfolding");
  std::vector<std::string> file_names;
  file_names.push_back("unfolded_me1D-P-tune1_col12.4_MAD_Global.root");
  file_names.push_back("unfolded_me1D-P-tune2_col12.4_MAD_Global.root");
  //file_names.push_back("bkgfit_me_pt_tune_FHC_Scale_Sarah.root");
  //file_names.push_back("bkgfit_me_q3_tune_FHC_Scale_Sarah.root");
  //file_names.push_back("bkgfit_me_pt_tune_FHC_Scale_Sarah.root");
  //file_names.push_back("bkgfit_me_regpt_tune_no_Scale_Sarah.root");
  //file_names.push_back("bkgfit_me_visE_tune_no_Scale_Sarah.root");
  char* cv_file = "unfolded_me1D-P_col12.4_MAD_Global.root";
  TH2D* cv_base = NULL;
  char* out_file = cv_file;

  //TCanvas c1("c1","c1",1920,1080);
  for (int i =0; i<hist_names.size();++i) {
    std::vector<TH2D*> errorband_vector;
    std::string path = file_path;
    path+=cv_file;
    TFile * f = TFile::Open(path.data());
    cv_base = dynamic_cast<TH2D*>(f->Get(hist_names[i].data())->Clone());
    //cv_base->Draw("COLZ");
    //c1.Print("cv.png");
    for (int j =0; j<file_names.size();++j) {
      std::string path = file_path;
      path+=file_names[j];
      TFile* f = TFile::Open(path.data());
      TH2D* hist = dynamic_cast<TH2D*>(f->Get(hist_names[i].data())->Clone());
      if (!hist) {
        std::cout << file_names[j] <<std::endl;
        continue;
      }
      errorband_vector.push_back(hist);
      f->Close();
    }
    PlotUtils::MnvVertErrorBand2D ErrorBand("bkgfit",cv_base,errorband_vector);
    path = file_path;
    path+=out_file;
    PlotUtils::MnvH2D * hist = dynamic_cast<PlotUtils::MnvH2D*>(f->Get(hist_names[i].data()));
    hist->PushErrorBand("bkgfit",&ErrorBand);
    TFile * f2 = TFile::Open(path.data(),"UPDATE");
    f2->cd();
    hist->Write(hist_names[i].data());
    //TH2D error_hist = ErrorBand.GetErrorBand(true);
    //error_hist.SetAxisRange(0,1,"Z");
    //error_hist.Draw("COLZ");
    //std::cout << error_hist.Integral();
    //c1.Print("err.png");
    //TH2D stat = hist->GetStatError(true);
    //stat.SetAxisRange(0,1,"Z");
    //stat.Draw("COLZ");
    //c1.Print("stat.png");
    f2->Close();
  }

}
