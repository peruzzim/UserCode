{

  
  gROOT->ProcessLine(".L template_production.C+O");
  TString dir;


 dir=TString("./");
 gen_templates(dir+TString("input_data.root"),"sigsig",1,"outphoton_data_sigsig.root","photoniso",5e5);
 gen_templates(dir+TString("input_data.root"),"sigbkg",1,"outphoton_data_sigbkg.root","photoniso",-1);
 gen_templates(dir+TString("input_data.root"),"bkgbkg",1,"outphoton_data_bkgbkg.root","photoniso",-1);
 gen_templates(dir+TString("input_data.root"),"standard",1,"outphoton_data_standard.root","photoniso");



 gen_templates(dir+TString("input_data.root"),"zee",1,"outphoton_data_zee.root","photoniso",-1);
 gen_templates(dir+TString("input_data.root"),"randomcone",1,"outphoton_data_rcone.root","photoniso",5e5);
 gen_templates(dir+TString("input_data.root"),"sieiesideband",1,"outphoton_data_sieiesideband.root","photoniso",-1);



}


