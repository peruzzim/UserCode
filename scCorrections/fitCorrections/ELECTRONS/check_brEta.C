#include "../../binning.h"
#include "corrections_brEta.C"
#include "TF2.h"
#include "TFile.h"

int check_brEta(){ 
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111111);

  TFile *f = new TFile("./histCorrections_brEta.root","read");  
  TF2 *f2 = new TF2("f2",F_applyScCorrectionsBrEta,0,3,0,10,2);

  for (Int_t i = 0 ; i< nBinsEta; ++i){
    
    TH1F *h = (TH1F*)f->Get(Form("h_corr_%d",i));    
    TH1F *h2 = new TH1F("h2","h2",nBinsBr*2-1,brbins);
    TH1F* res = new TH1F("res","res",100,-0.001,0.001);
    
    TCanvas *c = new TCanvas("c","c",1200,600);
    c->Divide(2);
    c->cd(1);
    h->GetYaxis()->SetRangeUser(0.90,1.01);
    if (i>6) h->GetYaxis()->SetRangeUser(0.85,1.01);
    h->GetXaxis()->SetRangeUser(0,5);
    h->Draw("EP");
    TLatex l(3.5,0.95,Form("h_corr_%d",i));
    l.Draw();
    
    for (Int_t j = 0; j< nBinsBr; ++j){
      h2->SetBinContent(2*j+1, f2->Eval( 0.5*(rightEta[i]-leftEta[i])+leftEta[i], h->GetBinCenter(2*j+1) ) );
      h2->SetBinError(j,0);
      res->Fill(h->GetBinContent(j) - h2->GetBinContent(j));
    }
    h2->Draw("same");
    
    c->cd(2);
    res->Draw();

    c->Update();
    getchar();
  }
  
  return 0;
}

