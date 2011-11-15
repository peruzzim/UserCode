#include "corrections_E.C"

int check_E(){

  TFile *f = new TFile("./histCorrections_E.root","read");  
  TH1F *hEB = (TH1F*)f->Get("h_CBE_EB");
  TH1F *hEE = (TH1F*)f->Get("h_CBE_EE");
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);

//   TCanvas *cEB = new TCanvas("cEB","cEB",1200,600);
//   cEB->Divide(2,1);
//   cEB->cd(1);
//   hEB->Draw();
//  
//   TF1 *corr = new TF1("corr",F_applyScCorrectionsE_EB,-100,5000,0);
//   corr->SetLineColor(4);
//   corr->Draw("same");
//   hEB->Draw("same");
//
//   TH1F* res_EEB = new TH1F("res_EEB","res_EEB",100,-0.0005,0.0005);
//   for (Int_t i = 1; i< hEB->GetNbinsX()+1; ++i){
//     res_EEB->Fill((hEB->GetBinContent(i) - corr->Eval(hEB->GetBinCenter(i))));
//   }
//   cEB->cd(2);
//   res_EEB->Draw();

  TCanvas *cEE = new TCanvas("cEE","cEE",1200,600);
  cEE->Divide(2,1);
  cEE->cd(1);
  hEE->Draw();
  
  TF1 *corr = new TF1("corr",F_applyScCorrectionsE_EE,-100,5000,0);
  corr->SetLineColor(4);
  corr->Draw("same");
  hEE->Draw("same");

  TH1F* res_EEE = new TH1F("res_EEE","res_EEE",100,-0.0005,0.0005);
  for (Int_t i = 1; i< hEE->GetNbinsX()+1; ++i){
    res_EEE->Fill((hEE->GetBinContent(i) - corr->Eval(hEE->GetBinCenter(i))));
  }
  cEE->cd(2);
  res_EEE->Draw();
}
