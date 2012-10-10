{
  gROOT->ProcessLine(".L template_production.C+O");
  TString dir("gg_minitree_020615_full2011_data_freeze/");
  gen_templates(dir+TString("input_data.root"),"randomcone",1,"outphoton_data_rcone.root","photoniso");
  gROOT->ProcessLine(".! cp outphoton_data_rcone.root forreweight_data_rcone.root");
  gen_templates(dir+TString("input_data.root"),"standard",1,"outphoton_data_standard.root","photoniso");
  gen_templates(dir+TString("input_data.root"),"sieiesideband",1,"outphoton_data_sieiesideband.root","photoniso");
  gROOT->ProcessLine(".! cp outphoton_data_sieiesideband.root forreweight_data_sieiesideband.root");

  dir = TString("gg_minitree_020615_purewfull2011_freeze/");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"randomcone",0,"outphoton_gjet_rcone.root","photoniso");
  gROOT->ProcessLine(".! mv outphoton_gjet_rcone.root forreweight_gjet_rcone.root");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"randomcone",0,"outphoton_gjet_rcone.root","photoniso",1,"gjet","data","rcone","rcone");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"background",0,"outphoton_gjet_bkg.root","photoniso");
  gROOT->ProcessLine(".! mv outphoton_gjet_bkg.root forreweight_gjet_bkg.root");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"background",0,"outphoton_gjet_bkg.root","photoniso",1,"gjet","data","bkg","sieiesideband");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"signal",0,"outphoton_gjet_sig.root","photoniso");
  gROOT->ProcessLine(".! mv outphoton_gjet_sig.root forreweight_gjet_sig.root");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"signal",0,"outphoton_gjet_sig.root","photoniso",1,"gjet","data","sig","rcone");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"sieiesideband",0,"outphoton_gjet_sieiesideband.root","photoniso");
  gROOT->ProcessLine(".! mv outphoton_gjet_sieiesideband.root forreweight_gjet_sieiesideband.root");
  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"sieiesideband",0,"outphoton_gjet_sieiesideband.root","photoniso",1,"gjet","data","sieiesideband","sieiesideband");

  gen_templates(dir+TString("GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Fall11_PU_S6_START42_V14B_v1_AODSIM.root"),"standard",0,"outphoton_gjet_standard.root","photoniso");
  gen_templates(dir+TString("DiPhotonJets_7TeV_madgraph_Fall11_PU_S6_START42_V14B_v1_AODSIM_triggerskim.root"),"standard",0,"outphoton_ggjets_standard.root","photoniso");


  //  dir="gg_minitree_020615_dataB_muons_muonveto/";
  //  gen_templates(dir+TString("DoubleMu_Run2011B_16Jan2012_v1_AOD_dimuonskim.root"),"muon",1,"outphoton_data_muon.root","photoniso");

  //  TFile *f1 = new TFile("outphoton_data_rcone.root");
  //  TFile *f2 = new TFile("outphoton_gjet_rcone.root");
  //  TBrowser *tb = new TBrowser();

  gROOT->ProcessLine(".L compare_3templates.C+");
  for (int i=0; i<7; i++) compare_3templates("outphoton","sig","EB",1,i,0.1,3,0.1,3);
  for (int i=0; i<5; i++) compare_3templates("outphoton","sig","EE",2,i,0.1,3,0.1,3);
  compare_3templates("outphoton","sig","EB",1,13,0.1,3,0.1,3);
  compare_3templates("outphoton","sig","EE",2,13,0.1,3,0.1,3);
  //  for (int i=0; i<7; i++) compare_3templates("outphoton","bkg","EB",2,i,0.1,3,0.1,3);
  //  for (int i=0; i<5; i++) compare_3templates("outphoton","bkg","EE",5,i,0.1,3,0.1,3);
  compare_3templates("outphoton","bkg","EB",5,13,0.1,3,0.1,3);
  compare_3templates("outphoton","bkg","EE",5,13,0.1,3,0.1,3);

//  gROOT->ProcessLine(".L compare_4templates.C");
//  for (int i=0; i<7; i++) compare_4templates("outphoton","EB",1,i,0.1,3,0.1,3);
//  for (int i=0; i<5; i++) compare_4templates("outphoton","EE",2,i,0.1,3,0.1,3);
//  compare_4templates("outphoton","EB",1,9,0.1,3,0.1,3);
//  compare_4templates("outphoton","EE",2,9,0.1,3,0.1,3);






  gROOT->ProcessLine(".! open plots/*.png");
}


