#include <string>
#include <vector>
#include <stdlib>
#include <iostream>

using namespace std;

void Tree_compare(const char* plot, const char* cut="", const char* filename1="output_mc_ETETA_nophietacracks/MiniTree_Zee.root", const char* filename2="output_data_ETETA_nophietacracks/MiniTree_Zee.root"){

  TFile *file1 = TFile::Open(filename1,"read");
  TFile *file2 = TFile::Open(filename2,"read");

  TTree *t1;
  TTree *t2;

  file1->GetObject("Tree",t1);
  file2->GetObject("Tree",t2);

  TString plotstring(plot);
  plotstring.Append(">>htemp");
  t1->Draw(plotstring.Data(),cut);
  TH1F *hmc = new TH1F(*htemp);
  t2->Draw(plotstring.Data(),cut);
  TH1F *hdata = new TH1F(*htemp);

  hmc->SetFillStyle(0);
  hmc->SetLineColor(kRed);
  hdata->SetFillStyle(0);
  hdata->SetLineColor(kBlue);

  hmc->Scale(1.0/hmc->Integral());
  hdata->Scale(1.0/hdata->Integral());

  hmc->Draw();
  hdata->Draw("same");

}
