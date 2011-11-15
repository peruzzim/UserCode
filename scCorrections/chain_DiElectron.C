{
  TChain data("ggNtuplizer/EventTree");

  Int_t n = 0;
  TString path = "../DATA/NTUPLES/ElectronGun_Summer11_NPC/";
  n+=data.Add(path+"*.root",0);
  
  cout << "#files = " << n << endl;
}
