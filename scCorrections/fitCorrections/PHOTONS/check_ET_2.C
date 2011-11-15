#include "corrections_ET.C"

int check_ET_2(){
  //  TF1 *corr = new TF1("corr",F_applyScCorrectionsET_EB,-100,500,0);
  TF1 *corr = new TF1("corr",F_applyScCorrectionsET_EE,-100,500,0);
  corr->SetLineColor(4);
  corr->Draw("");
  return 0;
}
