


TH1F* convert_histo_float(TFile *a, const char* filename="pileup.root", const char* name="analyze/PileUpStats"){

  TFile *b = new TFile(filename,"recreate");

  TH1I *ist;
  a->GetObject(name,ist);

  TH1F *res = new TH1F("pileup","pileup",36,-0.5,35.5);

  for (int i=0;i<35;i++) res->Fill(i,ist->GetBinContent(ist->FindBin(i)));

  b->cd();

  res->Write();

  b->Close();

      






}
