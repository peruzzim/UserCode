#include "binsdef.h"
#include <assert>
#include <iostream>


compare_4templates(TString pref, TString reg, int rbin=-1, int binchoice=9, float normmin=-9999, float normmax=9999, float minplotrange=-5, float maxplotrange=5){

  TString temp="sig";

  TString pref1;
  TString dset1("gjet");
  TString temp1;
  TString pref2;
  TString dset2("gjet");
  TString temp2;
  TString pref3;
  TString dset3("data");
  TString temp3;
  TString pref4;
  TString dset4("data");
  TString temp4;

  if (rbin==-1) rbin = (temp=="sig") ? 1 : 10;

  pref1=pref; pref2=pref; pref3=pref; pref4=pref;

  if (temp=="sig"){
    temp1=TString("sig");
    temp2=TString("rcone");
    temp3=TString("rcone");
    temp4=TString("muon");
  }

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

  TString file3=pref3;
  file3.Append("_");
  file3.Append(dset3);
  file3.Append("_");
  file3.Append(temp3);
  file3.Append(".root");

  TString file4=pref4;
  file4.Append("_");
  file4.Append(dset4);
  file4.Append("_");
  file4.Append(temp4);
  file4.Append(".root");


  TFile *f1 = new TFile(file1.Data(),"read");
  TFile *f2 = new TFile(file2.Data(),"read");
  TFile *f3 = new TFile(file3.Data(),"read");
  TFile *f4 = new TFile(file4.Data(),"read");

  TH1F *h[4];

  TString name1="";
  if (dset1=="data") name1.Append("data_Tree_"); else name1.Append("mc_Tree_");
  if (temp1=="bkg") name1.Append("background_template/template_background_");
  if (temp1=="sig") name1.Append("signal_template/template_signal_");
  if (temp1=="rcone") name1.Append("randomcone_signal_template/template_signal_");
  if (temp1=="impinging") name1.Append("impinging_track_template/template_background_");
  if (temp1=="sieiesideband") name1.Append("sieiesideband_sel/template_background_");
  if (temp1=="combisosideband") name1.Append("combisosideband_sel/template_background_");
  if (temp1=="muon") name1.Append("muoncone/template_signal_");
  name1.Append(reg);
  name1.Append(Form("_b%d",binchoice));

  TString name2="";
  if (dset2=="data") name2.Append("data_Tree_"); else name2.Append("mc_Tree_");
  if (temp2=="bkg") name2.Append("background_template/template_background_");
  if (temp2=="sig") name2.Append("signal_template/template_signal_");
  if (temp2=="rcone") name2.Append("randomcone_signal_template/template_signal_");
  if (temp2=="impinging") name2.Append("impinging_track_template/template_background_");
  if (temp2=="sieiesideband") name2.Append("sieiesideband_sel/template_background_");
  if (temp2=="combisosideband") name2.Append("combisosideband_sel/template_background_");
  if (temp2=="muon") name2.Append("muoncone/template_signal_");
  name2.Append(reg);
  name2.Append(Form("_b%d",binchoice));

  TString name3="";
  if (dset3=="data") name3.Append("data_Tree_"); else name3.Append("mc_Tree_");
  if (temp3=="bkg") name3.Append("background_template/template_background_");
  if (temp3=="sig") name3.Append("signal_template/template_signal_");
  if (temp3=="rcone") name3.Append("randomcone_signal_template/template_signal_");
  if (temp3=="impinging") name3.Append("impinging_track_template/template_background_");
  if (temp3=="sieiesideband") name3.Append("sieiesideband_sel/template_background_");
  if (temp3=="combisosideband") name3.Append("combisosideband_sel/template_background_");
  if (temp3=="muon") name3.Append("muoncone/template_signal_");
  name3.Append(reg);
  name3.Append(Form("_b%d",binchoice));

  TString name4="";
  if (dset4=="data") name4.Append("data_Tree_"); else name4.Append("mc_Tree_");
  if (temp4=="bkg") name4.Append("background_template/template_background_");
  if (temp4=="sig") name4.Append("signal_template/template_signal_");
  if (temp4=="rcone") name4.Append("randomcone_signal_template/template_signal_");
  if (temp4=="impinging") name4.Append("impinging_track_template/template_background_");
  if (temp4=="sieiesideband") name4.Append("sieiesideband_sel/template_background_");
  if (temp4=="combisosideband") name4.Append("combisosideband_sel/template_background_");
  if (temp4=="muon") name4.Append("muoncone/template_signal_");
  name4.Append(reg);
  name4.Append(Form("_b%d",binchoice));



  std::cout << name1.Data() << std::endl;
  f1->GetObject(name1.Data(),h[0]);
  assert (h[0]!=NULL);

  std::cout << name2.Data() << std::endl;
  f2->GetObject(name2.Data(),h[1]);
  assert (h[1]!=NULL);

  std::cout << name3.Data() << std::endl;
  f3->GetObject(name3.Data(),h[2]);
  assert (h[2]!=NULL);

  std::cout << name4.Data() << std::endl;
  f4->GetObject(name4.Data(),h[3]);
  assert (h[3]!=NULL);


  int colors[4];
  colors[0]=2; // rosso (MC driven)
  colors[1]=4; // blu (MC data-driven tech.)
  colors[2]=1; // nero (dati)
  colors[3]=3; // muons

  // Rebin
  for (int i=0; i<4; i++) h[i]->Rebin(rbin);

  for (int i=0; i<4; i++){
    h[i]->Scale(1.0/h[i]->Integral(h[i]->FindBin(normmin),h[i]->FindBin(normmax)));
    h[i]->SetAxisRange(minplotrange,maxplotrange);
    h[i]->SetLineColor(colors[i]);
    h[i]->SetLineWidth(2);
  }

//  // Kolmogorov test
//  cout << "Kolmogorov test" << endl;
//  cout << h[0]->KolmogorovTest(h[1]) << endl;

//  TLegend *leg = new TLegend(0.1,0.7,0.3,0.9);
//  leg->AddEntry(h[0],Form("%s_%s",dset1.Data(),temp1.Data()),"l");
//  leg->AddEntry(h[1],Form("%s_%s",dset2.Data(),temp2.Data()),"l");
//  leg->AddEntry(h[2],Form("%s_%s",dset3.Data(),temp3.Data()),"l");

  TCanvas *c1 = new TCanvas("iso","iso",1024,768);
  c1->cd();
  //  c1->Divide(2);

  h[2]->GetXaxis()->SetTitle(pref.Data());
  //  h[0]->GetXaxis()->SetRangeUser(0,maxrange);
  
  //  h[2]->GetYaxis()->SetRangeUser(0,0.09);

  float max=0;
  for (int i=0; i<4; i++){
    float thismax = h[i]->GetBinContent(h[i]->GetMaximumBin());
    max = (thismax>max) ? thismax : max;
  }
  h[2]->GetYaxis()->SetRangeUser(0,max*1.05);

  //  c1->cd(1);
  h[2]->Draw();
  h[1]->Draw("same");
  h[0]->Draw("same");
  h[3]->Draw("same");
  //  leg->Draw();	

//  c1->cd(2);
//  c1->GetPad(2)->SetLogy();
//  h[2]->Draw();
//  h[1]->Draw("same");
//  h[0]->Draw("same");
//  //  leg->Draw();	
//
  c1->Update();

  //  c1->SaveAs(Form("plots/comparison_%s_%s_%s_b%d.root",pref.Data(),temp.Data(),reg.Data(),binchoice));
  c1->SaveAs(Form("plots/comparison_%s_%s_%s_b%d.png",pref.Data(),temp.Data(),reg.Data(),binchoice));

}
