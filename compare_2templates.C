#include "binsdef.h"
#include <assert>
#include <iostream>


compare_2templates(TString pref1, TString pref2, TString dset1, TString dset2, TString temp1, TString temp2, TString reg="EB", TString colorstring="br", int rbin=1){

  TString file1=pref1;
  file1.Append("_");
  file1.Append(dset1);
  file1.Append("_");
  file1.Append(temp1);
  file1.Append(".root");

  TString file2=pref2;
  file2.Append("_");
  file2.Append(dset2);
  file2.Append("_");
  file2.Append(temp2);
  file2.Append(".root");


  TFile *f1 = new TFile(file1.Data(),"read");
  TFile *f2 = new TFile(file2.Data(),"read");

  TH1F *h[2];

  TString name1="";
  if (dset1=="data") name1.Append("data_Tree_"); else name1.Append("mc_Tree_");
  if (temp1=="bkg") name1.Append("background_template/template_background_");
  if (temp1=="sig") name1.Append("signal_template/template_signal_");
  if (temp1=="rcone") name1.Append("randomcone_signal_template/template_signal_");
  if (temp1=="impinging") name1.Append("impinging_track_template/template_background_");
  if (temp1=="sieiesideband") name1.Append("sieiesideband_sel/template_background_");
  if (temp1=="combisosideband") name1.Append("combisosideband_sel/template_background_");
  name1.Append(reg);
  name1.Append("_b9");

  TString name2="";
  if (dset2=="data") name2.Append("data_Tree_"); else name2.Append("mc_Tree_");
  if (temp2=="bkg") name2.Append("background_template/template_background_");
  if (temp2=="sig") name2.Append("signal_template/template_signal_");
  if (temp2=="rcone") name2.Append("randomcone_signal_template/template_signal_");
  if (temp2=="impinging") name2.Append("impinging_track_template/template_background_");
  if (temp2=="sieiesideband") name2.Append("sieiesideband_sel/template_background_");
  if (temp2=="combisosideband") name2.Append("combisosideband_sel/template_background_");
  name2.Append(reg);
  name2.Append("_b9");



  std::cout << name1.Data() << std::endl;
  f1->GetObject(name1.Data(),h[0]);
  assert (h[0]!=NULL);


  std::cout << name2.Data() << std::endl;
  f2->GetObject(name2.Data(),h[1]);
  assert (h[1]!=NULL);

  int colors[2];
  for (int i=0; i<2; i++){
    char c = colorstring[i];
    if (TString(c)=="b") colors[i]=1;
    else if (TString(c)=="r") colors[i]=2;
    else if (TString(c)=="p") colors[i]=6;
    else if (TString(c)=="u") colors[i]=4;
    else colors[i]=3;
  }


  // Rebin
  for (int i=0; i<2; i++) h[i]->Rebin(rbin);

  for (int i=0; i<2; i++){
    h[i]->Scale(1.0/h[i]->Integral());
    h[i]->SetLineColor(colors[i]);
    h[i]->SetLineWidth(2);
  }

//  // Kolmogorov test
//  cout << "Kolmogorov test" << endl;
//  cout << h[0]->KolmogorovTest(h[1]) << endl;

  TCanvas *c1 = new TCanvas();
  c1->cd();
  c1->Divide(2);

  c1->cd(1);
  h[0]->GetXaxis()->SetRangeUser(-5,5);
  h[0]->Draw();
  h[1]->Draw("same");

  c1->cd(2);
  c1->GetPad(2)->SetLogy();
  h[0]->GetXaxis()->SetRangeUser(-5,5);
  h[0]->Draw();
  h[1]->Draw("same");

  c1->Update();

}
