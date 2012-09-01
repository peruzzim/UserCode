TH1F* prepare_reweight_histo(TString dset1, TString dset2, TString temp1, TString temp2, TString reg="EB"){


  TString file1="out_";
  file1.Append(dset1);
  file1.Append("_");
  file1.Append(temp1);
  file1.Append(".root");

  TString file2="out_";
  file2.Append(dset2);
  file2.Append("_");
  file2.Append(temp2);
  file2.Append(".root");


  TFile *f1 = new TFile(file1.Data(),"read");
  TFile *f2 = new TFile(file2.Data(),"read");


  TString name1="";
  if (dset1=="data") name1.Append("data_Tree_"); else name1.Append("mc_Tree_");
  if (temp1=="bkg") name1.Append("background_template/");
  if (temp1=="sig") name1.Append("signal_template/");
  if (temp1=="rcone") name1.Append("randomcone_signal_template/");
  if (temp1=="impinging") name1.Append("impinging_track_template/");
  if (temp1=="sieiesideband") name1.Append("sieiesideband_sel/");
  name1.Append("histo_pt_");
  name1.Append(reg);

  TString name2="";
  if (dset2=="data") name2.Append("data_Tree_"); else name2.Append("mc_Tree_");
  if (temp2=="bkg") name2.Append("background_template/");
  if (temp2=="sig") name2.Append("signal_template/");
  if (temp2=="rcone") name2.Append("randomcone_signal_template/");
  if (temp2=="impinging") name2.Append("impinging_track_template/");
  if (temp2=="sieiesideband") name2.Append("sieiesideband_sel/");
  name2.Append("histo_pt_");
  name2.Append(reg);

  TH1F *h[2];
  f1->GetObject(name1,h[0]);
  f2->GetObject(name2,h[1]);
  assert(h[0]!=NULL);
  assert(h[1]!=NULL);

  h[0]->Print();
  h[1]->Print();

  TH1F *newhist = (TH1F*)(h[1]->Clone("reweight"));
  assert(newhist!=NULL);
  newhist->Print();

  newhist->Divide(h[0]);

  return newhist;

}
