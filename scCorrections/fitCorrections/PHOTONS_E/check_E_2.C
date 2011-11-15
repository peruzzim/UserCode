#include "corrections_E.C"

int check_E_2(){
  //  TF1 *corr = new TF1("corr",F_applyScCorrectionsE_EB,-100,5000,0);
  TF1 *corr = new TF1("corr",F_applyScCorrectionsE_EE,-100,5000,0);
  corr->SetLineColor(4);
  corr->Draw("");
  return 0;
}
