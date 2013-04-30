{
#include "binsdef.h"
#include <assert>
#include <iostream>

  TString reg="EB";
  bool dosignal=0;
  int rbin=1;


  TH1F* h[100];

  TString name;
  if (dosignal) name="signal";
  else name="background";

  int n_templates=0;
  if (reg=="EB") n_templates=n_templates_EB;
  if (reg=="EE") n_templates=n_templates_EE;

  assert(n_templates>0);

  for (int i=0; i<n_templates; i++) {

    gDirectory->GetObject(Form("mc_Tree_standard_sel/template_%s_%s_b%d",name.Data(),reg.Data(),i),h[i]);
    h[i]->Scale(1.0/h[i]->Integral());
    h[i]->SetLineColor(i+1);

  }

  gDirectory->GetObject(Form("mc_Tree_standard_sel/template_%s_%s_b%d",name.Data(),reg.Data(),n_bins),h[n_bins]);
  h[n_bins]->Scale(1.0/h[n_bins]->Integral());
  h[n_bins]->SetLineColor(kBlack);
  h[n_bins]->SetLineWidth(2);

  // Kolmogorov test
  cout << "Kolmogorov test" << endl;
  for (int i=0; i<n_templates; i++) cout << "b" << i << " " << h[n_bins]->KolmogorovTest(h[i]) << endl;

  // Rebin
  for (int i=0; i<n_templates; i++) h[i]->Rebin(rbin);
  h[n_bins]->Rebin(rbin);


  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  for (int i=0; i<n_templates; i++) leg->AddEntry(h[i],Form("b%d",i),"l");
  leg->AddEntry(h[n_bins],"all","l");



  TCanvas *c1 = new TCanvas();
  c1->cd();

  h[n_bins]->Draw();
  for (int i=0;i<n_templates; i++) h[i]->Draw("Esame");
  leg->Draw();


  TCanvas *c2 = new TCanvas();
  c2->cd();

  c2->Divide(3,3);

  for (int i=0; i<n_templates; i++) {
    c2->cd(i+1);
    h[n_bins]->Draw("axis");
    h[i]->Draw("same");
    h[n_bins]->Draw("same");
    h[n_bins]->SetAxisRange(0,.15,"y");
  }

  for (int i=0; i<n_templates; i++) {
    c2->cd(i+1);
    c2->Update();
  }





}
