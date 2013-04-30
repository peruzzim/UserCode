{
  gROOT->ProcessLine(".L template_production.C+O");

  TString file("gg_minitree_020616_mc_full2011purew_nopf02presel_noeffarea/DiPhotonJets_7TeV_madgraph_Fall11_PU_S6_START42_V14B_v1_AODSIM.root");
  //TString file("./input_data.root");

  get_eff_area(file.Data(),1,"photoniso");
//  get_eff_area(file.Data(),1,"chargediso");
//  get_eff_area(file.Data(),1,"neutraliso");
//  get_eff_area(file.Data(),1,"combiso");
  get_eff_area(file.Data(),0,"photoniso");
//  get_eff_area(file.Data(),0,"chargediso");
//  get_eff_area(file.Data(),0,"neutraliso");
//  get_eff_area(file.Data(),0,"combiso");

}
