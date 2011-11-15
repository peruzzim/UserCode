{
  gROOT->Reset();
  gROOT->ProcessLine(".x ./chain_DiElectron.C");
  gROOT->ProcessLine(".L scCorrections.C+");
  gROOT->ProcessLine("scCorrections a((TTree*) data)");
  gROOT->ProcessLine("a.run_brEta()");
}  
