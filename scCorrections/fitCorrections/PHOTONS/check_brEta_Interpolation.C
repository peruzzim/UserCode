#include "corrections_brEta.C"

int check_brEta_Interpolation(){

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  TF2 *f2 = new TF2("f2",F_applyScCorrectionsBrEta,-3.0,3.0,0,10,2);
  f2->SetNpx(20);
  f2->SetNpy(20);
  f2->Draw("lego2");
  return 0;
}
