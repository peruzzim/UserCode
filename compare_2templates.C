#include "binsdef.h"
#include <assert>
#include <iostream>


compare_2templates(TString file1, TString file2, TString name1, TString name2, float maxrange=5, int rbin=1){

  TFile *f1 = new TFile(file1.Data(),"read");
  TFile *f2 = new TFile(file2.Data(),"read");

  TH1F *h[2];

  f1->GetObject(name1.Data(),h[0]);
  f2->GetObject(name2.Data(),h[1]);
  
  assert (h[0]!=NULL);
  assert (h[1]!=NULL);

  // Rebin
  for (int i=0; i<2; i++) h[i]->Rebin(rbin);

  for (int i=0; i<2; i++){
    h[i]->Scale(1.0/h[i]->Integral());
    h[i]->SetLineColor(i+1);
    h[i]->SetLineWidth(2);
  }

  // Kolmogorov test
  cout << "Kolmogorov test" << endl;
  cout << h[0]->KolmogorovTest(h[1]) << endl;

  TCanvas *c1 = new TCanvas();
  c1->cd();

  h[0]->GetXaxis()->SetRangeUser(0,maxrange);
  h[0]->Draw();
  h[1]->Draw("same");

  c1->Update();

}
