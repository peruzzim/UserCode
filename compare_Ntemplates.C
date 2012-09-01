#include "binsdef.h"
#include <assert>
#include <iostream>


void compare_Ntemplates(TString pref1, TString pref2, TString dset1, TString dset2, TString temp1, TString temp2, TString reg="EB", int rbin=1, int choosebin1=-1, int choosebin2=-1, float maxrange=5){

  const int nbins=10;

  int choosebin[2];
  choosebin[0]=choosebin1;
  choosebin[1]=choosebin2;

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

  TH1F *h[2][nbins];

  TString name1="";
  if (dset1=="data") name1.Append("data_Tree_"); else name1.Append("mc_Tree_");
  if (temp1=="bkg") name1.Append("background_template/template_background_");
  if (temp1=="sig") name1.Append("signal_template/template_signal_");
  if (temp1=="rcone") name1.Append("randomcone_signal_template/template_signal_");
  if (temp1=="impinging") name1.Append("impinging_track_template/template_background_");
  if (temp1=="sieiesideband") name1.Append("sieiesideband_sel/template_background_");
  name1.Append(reg);
  name1.Append("_b");

  TString name2="";
  if (dset2=="data") name2.Append("data_Tree_"); else name2.Append("mc_Tree_");
  if (temp2=="bkg") name2.Append("background_template/template_background_");
  if (temp2=="sig") name2.Append("signal_template/template_signal_");
  if (temp2=="rcone") name2.Append("randomcone_signal_template/template_signal_");
  if (temp2=="impinging") name2.Append("impinging_track_template/template_background_");
  if (temp2=="sieiesideband") name2.Append("sieiesideband_sel/template_background_");
  name2.Append(reg);
  name2.Append("_b");


  std::cout << name1.Data() << std::endl;
  std::cout << name2.Data() << std::endl;

  for (int i=0; i<nbins; i++){
    TString name = name1;
    name.Append(Form("%d",i));
    f1->GetObject(name,h[0][i]);
    assert (h[0][i]!=NULL);
  }

  for (int i=0; i<nbins; i++){
    TString name = name2;
    name.Append(Form("%d",i));
    f2->GetObject(name,h[1][i]);
    assert (h[1][i]!=NULL);
  }

  int colors[nbins];
  for (int i=0; i<nbins; i++) colors[i]=20+i;
  colors[nbins-1]=1;

  for (int j=0; j<2; j++){
    for (int i=0; i<nbins; i++){
      h[j][i]->SetLineColor(colors[i]);
      h[j][i]->SetMarkerColor(colors[i]);
      if (j==0) {
	h[j][i]->SetLineStyle(1);
	h[j][i]->SetMarkerStyle(20);
      }
      if (j==1) {
	h[j][i]->SetLineStyle(2);
	h[j][i]->SetMarkerStyle(22);
      }
      h[j][i]->SetLineWidth(2);
      h[j][i]->Rebin(rbin);
      h[j][i]->Scale(1.0/h[j][i]->Integral());
    }
  }
  

//  // Kolmogorov test
//  cout << "Kolmogorov test" << endl;
//  cout << h[0]->KolmogorovTest(h[1]) << endl;

  TCanvas *c1 = new TCanvas();
  c1->cd();
  c1->Divide(2);

  c1->cd(1);
  h[0][nbins-1]->GetXaxis()->SetRangeUser(0,maxrange);
  h[0][nbins-1]->Draw("axis");
  for (int j=0; j<2; j++){
    //    for (int i=0; i<nbins-1; i++){
    for (int i=0; i<nbins; i++){
      if (choosebin[j]>=0 && i!=choosebin[j]) continue;
      h[j][i]->Draw("same");
    }
  }
  //  h[0][nbins-1]->Draw("same");
  //  h[1][nbins-1]->Draw("same");

  c1->cd(2);
  c1->GetPad(2)->SetLogy();
  h[0][nbins-1]->GetXaxis()->SetRangeUser(0,maxrange);
  h[0][nbins-1]->Draw("axis");
  for (int j=0; j<2; j++){
    //    for (int i=0; i<nbins-1; i++){
    for (int i=0; i<nbins; i++){
      if (choosebin[j]>=0 && i!=choosebin[j]) continue;
      h[j][i]->Draw("same");
    }   
  }
  //  h[0][nbins-1]->Draw("same");
  //  h[1][nbins-1]->Draw("same");

  c1->Update();

}
