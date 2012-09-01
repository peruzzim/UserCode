{

  TString dir("gg_minitree_020615_purew2011B_witheffarea/");

  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"signal",0,"outphoton_gjet_sig.root","photoniso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"background",0,"outphoton_gjet_bkg.root","photoniso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"randomcone",0,"outphoton_gjet_rcone.root","photoniso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"sieiesideband",0,"outphoton_gjet_sieiesideband.root","photoniso");
  //  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"impinging",0,"outphoton_gjet_impinging.root","photoniso");

  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"signal",0,"outcomb_gjet_sig.root","combiso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"background",0,"outcomb_gjet_bkg.root","combiso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"randomcone",0,"outcomb_gjet_rcone.root","combiso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"sieiesideband",0,"outcomb_gjet_sieiesideband.root","combiso");
  //  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"impinging",0,"outcomb_gjet_impinging.root","combiso");

  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"signal",0,"outcharged_gjet_sig.root","chargediso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"background",0,"outcharged_gjet_bkg.root","chargediso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"randomcone",0,"outcharged_gjet_rcone.root","chargediso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"sieiesideband",0,"outcharged_gjet_sieiesideband.root","chargediso");
  //  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"impinging",0,"outcharged_gjet_impinging.root","chargediso");

  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"signal",0,"outneutral_gjet_sig.root","neutraliso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"background",0,"outneutral_gjet_bkg.root","neutraliso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"randomcone",0,"outneutral_gjet_rcone.root","neutraliso");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"sieiesideband",0,"outneutral_gjet_sieiesideband.root","neutraliso");
  //  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"impinging",0,"outneutral_gjet_impinging.root","neutraliso");


  
 dir="gg_minitree_020615_dataB_witheffarea/";
 
 gen_templates(dir+TString("input_data.root"),"standard",1,"outphoton_data_standard.root","photoniso");
 gen_templates(dir+TString("input_data.root"),"randomcone",1,"outphoton_data_rcone.root","photoniso");
 gen_templates(dir+TString("input_data.root"),"sieiesideband",1,"outphoton_data_sieiesideband.root","photoniso");
 // gen_templates(dir+TString("input_data.root"),"impinging",1,"outphoton_data_impinging.root","photoniso");

 gen_templates(dir+TString("input_data.root"),"standard",1,"outcomb_data_standard.root","combiso");
 gen_templates(dir+TString("input_data.root"),"randomcone",1,"outcomb_data_rcone.root","combiso");
 gen_templates(dir+TString("input_data.root"),"sieiesideband",1,"outcomb_data_sieiesideband.root","combiso");
 // gen_templates(dir+TString("input_data.root"),"impinging",1,"outcomb_data_impinging.root","combiso");

 gen_templates(dir+TString("input_data.root"),"standard",1,"outcharged_data_standard.root","chargediso");
 gen_templates(dir+TString("input_data.root"),"randomcone",1,"outcharged_data_rcone.root","chargediso");
 gen_templates(dir+TString("input_data.root"),"sieiesideband",1,"outcharged_data_sieiesideband.root","chargediso");
 // gen_templates(dir+TString("input_data.root"),"impinging",1,"outcharged_data_impinging.root","chargediso");

 gen_templates(dir+TString("input_data.root"),"standard",1,"outneutral_data_standard.root","neutraliso");
 gen_templates(dir+TString("input_data.root"),"randomcone",1,"outneutral_data_rcone.root","neutraliso");
 gen_templates(dir+TString("input_data.root"),"sieiesideband",1,"outneutral_data_sieiesideband.root","neutraliso");
 // gen_templates(dir+TString("input_data.root"),"impinging",1,"outneutral_data_impinging.root","neutraliso");
 


}
