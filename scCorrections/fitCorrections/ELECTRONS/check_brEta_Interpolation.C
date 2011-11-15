#include "corrections_brEta.C"

int check_brEta_Interpolation(){

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  TF2 *f2 = new TF2("f2",F_applyScCorrectionsBrEta,-3.0,3.0,-1,15,2);
  f2->SetNpx(50);
  f2->SetNpy(50);
  f2->Draw("lego2");
  return 0;
}
