#include <string>
#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TH1.h"
#include <vector>
double count(TFile* const f, const char* const histname) {
  TH1D * hist = dynamic_cast<TH1D*>(f->Get(histname));
  if (hist) {
    return hist->Integral(0,-1);
  } //else {
  //  std::cout << "wrong hist name: " << histname << std::endl;
  //}
  //if (hist) return hist->Integral(0,100);
  return 0;
}

void confusionM(double selected,double selected_signal,double total_event,double total_signal,double* CM) {
  double TP = selected_signal;
  double FP = selected-selected_signal;
  double FN = total_signal-selected_signal;
  double TN = total_event-TP-FP-FN;
  double TPR = TP/(TP+FN);
  double TNR =  TN/(TN+FP);
  CM[0]=TPR+TNR-1;
  CM[1] = (TP*TN-FP*FN) / (sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)));
  //std::cout << "confusion Matrix:" << std::endl;
  //std::cout << TP << "," << FP <<std::endl;
  //std::cout << FN << "," << TN <<std::endl;
}

int printEfficiencyNPurity() {
  TFile* f = TFile::Open("/minerva/data/users/hsu/nu_e/cut_study_mcme1B_nx_wcut_real_Sarah.root");
  //TFile* fsig = TFile::Open("/minerva/data/users/hsu/nu_e/truth_dist_mcme1B_nx_collab1_LLR.root")
  const char* const var = "tEavail";
  std::vector<std::string> signal_list;
  //{"CCNuEQE","CCNuEDelta","CCNuEDIS","CCNuE2p2h","CCNuE"};
  //signal_list.push_back("CCNuEQE");
  //signal_list.push_back("CCNuEDelta");
  //signal_list.push_back("CCNuEDIS");
  //signal_list.push_back("CCNuE2p2h");
  signal_list.push_back("CCNuE");
  //signal_list.push_back("NCPi0");
  //signal_list.push_back("NCOther");
  double total_signal = -1;
  double total_event = -1;
  double prev_selected = -1;
  double ConM[2];
  int i = 0;
  while (1) {
    double selected =0;
    std::stringstream ss;
    ss << var << "_" << i << "cuts";
    selected += count(f,ss.str().c_str());
    ss.str("");
    ss << var << "_" << i << "cut";
    selected += count(f,ss.str().c_str());
    double selected_signal = 0;
    for (int j = 0;j<signal_list.size();++j) { 
      ss.str("");
      ss << var << "_" <<signal_list[j]<<"_" << i <<"cuts";
      selected_signal += count(f,ss.str().c_str());
    }
    //backward compatability
    ss.str("");
    ss << var << "_" << i <<"cut_signal" ;
    selected_signal += count(f,ss.str().c_str());
    if (selected == 0 || selected_signal == 0) {
      break;
    }
    if (i==0) {
      total_signal = selected_signal;
      total_event = selected;
      std::cout << "total event: " << total_event  << std::endl;
      std::cout << "total signal: " << total_signal  << std::endl;
    }
    confusionM(selected,selected_signal,total_event,total_signal,ConM);
    std::cout << "Cut number: " << i;
    std::cout << " &" << selected;
    std::cout << " & " << selected_signal/selected;
    std::cout << " & " << selected_signal/total_signal;
    std::cout << "& " << selected_signal / prev_selected;
    std::cout << ", MCC: " << ConM[1] <<std::endl;
    ++i;
    prev_selected = selected_signal;
  }
  return 0;
}
