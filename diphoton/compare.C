void compare_w(TString file1, TString file2){

  TFile *f1 = new TFile(file1.Data());
  TFile *f2 = new TFile(file2.Data());

  std::vector<TString> histos;
  histos.push_back("purity_sigsig");
  histos.push_back("purity_sigbkg");
  histos.push_back("purity_bkgsig");
  histos.push_back("purity_bkgbkg");
  histos.push_back("purity_fsig1");
  histos.push_back("purity_fsig2");

  TH1F **h1 = new TH1F*[histos.size()];
  TH1F **h2 = new TH1F*[histos.size()];
  TCanvas *c = new TCanvas("canv","canv",1800,1000);
  c->Divide(6,2);

  for (int i=0; i<histos.size(); i++){
    f1->GetObject(histos.at(i).Data(),h1[i]);
    f2->GetObject(histos.at(i).Data(),h2[i]);
    if (!h1[i] || !h2[i]) continue;
    h1[i]->Print();
    h2[i]->Print();
    h1[i]->Sumw2();
    h2[i]->Sumw2();
    h1[i]->SetLineWidth(1);
    h2[i]->SetLineWidth(1);
    c->cd(i+1);
    h1[i]->Draw("PE1");
    h1[i]->SetTitle(histos.at(i).Data());
    h2[i]->SetMarkerStyle(21);
    h2[i]->Draw("samePE1");
    c->cd(i+6+1);
    TH1F *temp = (TH1F*)(h1[i]->Clone("temp"));
    temp->Divide(h2[i]);
    temp->SetAxisRange(0.5,1.5,"Y");
    temp->Draw();
  }




};

void compare(TString dir1, TString dir2, TString var, TString reg){
  TString file1(Form("%s/histo_purity_%s_%s_allbins.root",dir1.Data(),var.Data(),reg.Data()));
  TString file2(Form("%s/histo_purity_%s_%s_allbins.root",dir2.Data(),var.Data(),reg.Data()));
  compare_w(file1,file2);
}
