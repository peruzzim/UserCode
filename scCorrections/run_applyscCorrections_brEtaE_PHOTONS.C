{
  gROOT->Reset();
  gROOT->ProcessLine(".x ./chain_DiPhoton.C");
  gROOT->ProcessLine(".L scCorrections.C+");
  gROOT->ProcessLine("scCorrections a((TTree*) data)");
  gROOT->ProcessLine("a.run_apply_brEtaE()");

}
