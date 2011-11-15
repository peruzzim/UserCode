{
  TChain data("analyze/Analysis");

  Int_t n = 0;
  TString path = "/Users/peruzzi/work/ntuples/";
  n+=data.Add(path+"DiPhotonGun*.root",0);
    
  cout << "#files = " << n << endl;
}
