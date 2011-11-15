#include "corrections_ET.C"

int check_ET(){

  TFile *f = new TFile("./histCorrections_ET.root","read");  
  TH1F *hEB = (TH1F*)f->Get("h_CBET_EB");
  TH1F *hEE = (TH1F*)f->Get("h_CBET_EE");
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);

  TCanvas *cEB = new TCanvas("cEB","cEB",1200,600);
  cEB->Divide(2,1);
  cEB->cd(1);
  hEB->Draw();
  
  TF1 *corr = new TF1("corr",F_applyScCorrectionsET_EB,-100,500,0);
  corr->SetLineColor(4);
  corr->Draw("same");
  hEB->Draw("same");

  TH1F* res_ETEB = new TH1F("res_ETEB","res_ETEB",100,-0.0005,0.0005);
  for (Int_t i = 1; i< hEB->GetNbinsX()+1; ++i){
    res_ETEB->Fill((hEB->GetBinContent(i) - corr->Eval(hEB->GetBinCenter(i))));
  }
  cEB->cd(2);
  res_ETEB->Draw();

  TCanvas *cEE = new TCanvas("cEE","cEE",1200,600);
  cEE->Divide(2,1);
  cEE->cd(1);
  hEE->Draw();
  
  TF1 *corr = new TF1("corr",F_applyScCorrectionsET_EE,-100,500,0);
  corr->SetLineColor(4);
  corr->Draw("same");
  hEE->Draw("same");

  TH1F* res_ETEE = new TH1F("res_ETEE","res_ETEE",100,-0.0005,0.0005);
  for (Int_t i = 1; i< hEE->GetNbinsX()+1; ++i){
    res_ETEE->Fill((hEE->GetBinContent(i) - corr->Eval(hEE->GetBinCenter(i))));
  }
  cEE->cd(2);
  res_ETEE->Draw();
}
