#define efficiency_measure_cxx
#include "efficiency_measure.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void efficiency_measure::Loop(){
  TFile *outf = new TFile(outname.Data(),"recreate");
  //  diffvariables_list.push_back(TString("dosingle"));
  for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
    LoopOne(*diffvariable,outf);
  }
  outf->Close();
};

void efficiency_measure::LoopOne(TString diffvariable, TFile *outf)
{

  assert (outf->GetName()=="eff_ggjets_num.root" || outf->GetName()=="eff_ggjets_den.root");
  bool isnum = (outf->GetName()=="eff_ggjets_num.root");

//   In a ROOT session, you can do:
//      Root > .L efficiency_measure.C
//      Root > efficiency_measure t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1::SetDefaultSumw2(kTRUE);

   TString sg_name[2]={"EB","EE"};
   TString gg_name[3]={"EBEB","EBEE","EEEE"};

   int n_templates_sg[2];
   float *binsdef_sg[2];

   int n_templates_gg[3];
   float *binsdef_gg[3];

   bool dosingle=0;

  if (diffvariable=="invmass"){
    n_templates_gg[0]=n_templates_invmass_EBEB;
    n_templates_gg[1]=n_templates_invmass_EBEE;
    n_templates_gg[2]=n_templates_invmass_EEEE; 
    binsdef_gg[0]=binsdef_diphoton_invmass_EBEB;
    binsdef_gg[1]=binsdef_diphoton_invmass_EBEE;
    binsdef_gg[2]=binsdef_diphoton_invmass_EEEE;
  }
  if (diffvariable=="diphotonpt"){
    n_templates_gg[0]=n_templates_diphotonpt_EBEB;
    n_templates_gg[1]=n_templates_diphotonpt_EBEE;
    n_templates_gg[2]=n_templates_diphotonpt_EEEE; 
    binsdef_gg[0]=binsdef_diphoton_diphotonpt_EBEB;
    binsdef_gg[1]=binsdef_diphoton_diphotonpt_EBEE;
    binsdef_gg[2]=binsdef_diphoton_diphotonpt_EEEE;
  }
  if (diffvariable=="costhetastar"){
    n_templates_gg[0]=n_templates_costhetastar_EBEB;
    n_templates_gg[1]=n_templates_costhetastar_EBEE;
    n_templates_gg[2]=n_templates_costhetastar_EEEE; 
    binsdef_gg[0]=binsdef_diphoton_costhetastar_EBEB;
    binsdef_gg[1]=binsdef_diphoton_costhetastar_EBEE;
    binsdef_gg[2]=binsdef_diphoton_costhetastar_EEEE;
  }
  if (diffvariable=="dphi"){
    n_templates_gg[0]=n_templates_dphi_EBEB;
    n_templates_gg[1]=n_templates_dphi_EBEE;
    n_templates_gg[2]=n_templates_dphi_EEEE; 
    binsdef_gg[0]=binsdef_diphoton_dphi_EBEB;
    binsdef_gg[1]=binsdef_diphoton_dphi_EBEE;
    binsdef_gg[2]=binsdef_diphoton_dphi_EEEE;
  }
  if (diffvariable=="dosingle"){
    dosingle=1;
    n_templates_sg[0]=n_templates_EB;
    n_templates_sg[1]=n_templates_EE;
    binsdef_sg[0]=binsdef_single_gamma_EB_eta;
    binsdef_sg[1]=binsdef_single_gamma_EE_eta;
    diffvariable="singlegamma_eta";
  }



  if (dosingle){
    for (int i=0; i<2; i++) w_tot_sg[i] = new TH1F(Form("w_tot_sg_%s",sg_name[i].Data()),Form("w_tot_sg_%s",sg_name[i].Data()),n_templates_sg[i],binsdef_sg[i]);
    for (int i=0; i<2; i++) w_tot_sg[i]->Sumw2();
  }
  else {
    for (int i=0; i<3; i++) w_tot_gg[i] = new TH1F(Form("w_tot_gg_%s_%s",gg_name[i].Data(),diffvariable.Data()),Form("w_tot_gg_%s_%s",gg_name[i].Data(),diffvariable.Data()),n_templates_gg[i],binsdef_gg[i]);
    for (int i=0; i<3; i++) w_tot_gg[i]->Sumw2();
  }



   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;



      Float_t weight=event_luminormfactor*event_Kfactor*event_weight;



      Int_t event_ok_for_dataset=-1;

      Int_t reg_lead;
      Int_t reg_trail;

      // categorization:
      // templates 0:EB 1:EE (only EBEB and EEEE events used)
      // dataset 0:EBEB 3/4->1:EBEE 2:EEEE

      if (fabs(pholead_SCeta)<1.4442 && fabs(photrail_SCeta)<1.4442) {
	event_ok_for_dataset=0;
	reg_lead=0;
	reg_trail=0;
      }
      else if (fabs(pholead_SCeta)>1.56 && fabs(photrail_SCeta)>1.56) {
	event_ok_for_dataset=2;
	reg_lead=1;
	reg_trail=1;
      }
      else if (fabs(pholead_SCeta)<1.4442 && fabs(photrail_SCeta)>1.56) {
	event_ok_for_dataset=3;
	reg_lead=0;
	reg_trail=1;
      }
      else if (fabs(pholead_SCeta)>1.56 && fabs(photrail_SCeta)<1.4442) {
	event_ok_for_dataset=4;
	reg_lead=1;
	reg_trail=0;
      }

      if (event_ok_for_dataset==3 || event_ok_for_dataset==4) event_ok_for_dataset=1;

      //      if (pholead_pt<40 || photrail_pt<25 || dipho_mgg_photon<80) continue;
      if (pholead_pt<40 || photrail_pt<25) continue;


    bool recalc_lead =  true;
    bool recalc_trail = true;

    float pholead_outvar = -999;
    float photrail_outvar = -999;


    if (recalc_lead){
      float et_recalc = 0;
      float e_recalc = 0;
      int number_recalc = 0;

      bool printout=false;

      if (printout) std::cout << "---" << std::endl;

      for (int i=0; i<pholead_Npfcandphotonincone; i++) {

	float et=pholead_photonpfcandets[i];
	float e=pholead_photonpfcandenergies[i];
	float deta=pholead_photonpfcanddetas[i];
	float dphi=pholead_photonpfcanddphis[i];
	float dR=sqrt(deta*deta+dphi*dphi);
	float eta=fabs(TMath::ACosH(e/et));
	if (printout) {
	  std::cout << et << " " << e << " " << deta << " " << dphi << " " << dR << " " << eta << std::endl;
	}
	if (eta>1.4442 && eta<1.56) continue;
	if (eta>2.5) continue;

#include "cleaning.cc"

	if (fabs(pholead_SCeta)<1.4442 && eta>1.4442) continue;
	if (fabs(pholead_SCeta)>1.56 && eta<1.56) continue;
	et_recalc+=et;
	e_recalc+=e;
	number_recalc++;

	//	hist2d_singlecandet->Fill(et,eta,weight*ptweight_lead);
	//	hist2d_singlecandenergy->Fill(e,eta,weight*ptweight_lead);
	//	hist2d_singlecandet->Fill(et/pholead_pt,eta,weight*ptweight_lead);
	//	hist2d_singlecandenergy->Fill(e/pholead_energy,eta,weight*ptweight_lead);
	//	hist2d_singlecanddeta->Fill(deta,eta,weight*ptweight_lead);
	//	hist2d_singlecanddphi->Fill(dphi,eta,weight*ptweight_lead);
	//	hist2d_singlecanddR->Fill(dR,eta,weight*ptweight_lead);
      }

//      hist2d_coneet->Fill(et_recalc,fabs(pholead_SCeta),weight*ptweight_lead);
//      hist2d_coneenergy->Fill(e_recalc,fabs(pholead_SCeta),weight*ptweight_lead);
//      hist2d_iso_ncand[reg_lead][bin_lead]->Fill(et_recalc,number_recalc,weight*ptweight_lead);

      pholead_outvar=et_recalc;

      if (printout) std::cout << "---" << std::endl;

    }

    if (recalc_trail){
      float et_recalc = 0;
      float e_recalc = 0;
      int number_recalc = 0;

      bool printout=false;

      if (printout) std::cout << "---" << std::endl;

      for (int i=0; i<photrail_Npfcandphotonincone; i++) {

	float et=photrail_photonpfcandets[i];
	float e=photrail_photonpfcandenergies[i];
	float deta=photrail_photonpfcanddetas[i];
	float dphi=photrail_photonpfcanddphis[i];
	float dR=sqrt(deta*deta+dphi*dphi);
	float eta=fabs(TMath::ACosH(e/et));
	if (printout) {
	  std::cout << et << " " << e << " " << deta << " " << dphi << " " << dR << " " << eta << std::endl;
	}
	if (eta>1.4442 && eta<1.56) continue;
	if (eta>2.5) continue;

#include "cleaning.cc"

	if (fabs(photrail_SCeta)<1.4442 && eta>1.4442) continue;
	if (fabs(photrail_SCeta)>1.56 && eta<1.56) continue;
	et_recalc+=et;
	e_recalc+=e;
	number_recalc++;

      }


      photrail_outvar=et_recalc;

      if (printout) std::cout << "---" << std::endl;

    }
     
    pholead_outvar-=getpuenergy(reg_lead,pholead_SCeta);
    photrail_outvar-=getpuenergy(reg_trail,photrail_SCeta);

    assert (pholead_outvar>-100 && photrail_outvar>-100);

      bool pass1=true;
      if (!(pholead_PhoMCmatchexitcode==1 || pholead_PhoMCmatchexitcode==2)) pass1=false;
      if (!(photrail_PhoMCmatchexitcode==1 || photrail_PhoMCmatchexitcode==2)) pass1=false;
      if (pholead_GenPhotonIsoDR04>5 || photrail_GenPhotonIsoDR04>5) pass1=false;
//      if (pholead_outvar<leftrange)   pass1=false;
//      if (photrail_outvar<leftrange)  pass1=false;

      if (isnum){
	if (pholead_outvar>=rightrange) pass1=false;
	if (photrail_outvar>=rightrange)pass1=false;
      }

      if (!dosingle){
      float fillvar=0;
      if (diffvariable==TString("invmass")) fillvar=dipho_mgg_photon;
	if (diffvariable==TString("diphotonpt")){
	  float px = pholead_px+photrail_px;
	  float py = pholead_py+photrail_py;
	  float pt = sqrt(px*px+py*py);
	  fillvar = pt;
	}
	if (diffvariable==TString("costhetastar")){
	  TLorentzVector pho1(pholead_px,pholead_py,pholead_pz,pholead_energy);
	  TLorentzVector pho2(photrail_px,photrail_py,photrail_pz,photrail_energy);
	  TVector3 boost = (pho1+pho2).BoostVector();
	  TLorentzVector boostedpho1 = pho1;
	  boostedpho1.Boost(-boost);
	  float thetastar1 = boostedpho1.Angle(boost);
	  if (thetastar1>TMath::Pi()/2) thetastar1 = TMath::Pi()-thetastar1;
	  fillvar=TMath::Cos(thetastar1);
	}
	if (diffvariable==TString("dphi")){
	  float phi1 = pholead_SCphi;
	  float phi2 = photrail_SCphi;
	  float dphi = AbsDeltaPhi(phi1,phi2);
	  fillvar=dphi;
	}
	

      if (pass1) w_tot_gg[event_ok_for_dataset]->Fill(fillvar,weight);
      }

      if (dosingle){
	if (pass1) w_tot_sg[reg_lead]->Fill(fabs(pholead_SCeta),weight);
	if (pass1) w_tot_sg[reg_trail]->Fill(fabs(photrail_SCeta),weight);
      }




   } // end event loop


   outf->cd();
   if (!dosingle) for (int i=0; i<3; i++) w_tot_gg[i]->Write();   
   else for (int i=0; i<2; i++) w_tot_sg[i]->Write();   



}

void divide_eff_histosOne(TString numerator, TString denominator, TString diffvariable, TFile *outf){

  TString gg_name[3]={"EBEB","EBEE","EEEE"};
  TString sg_name[2]={"EB","EE"};

  TFile *nfile = new TFile(numerator.Data(),"read");
  TFile *dfile = new TFile(denominator.Data(),"read");

  TH1F w_eff_gg[3];
  TH1F *w_num_gg[3];
  TH1F *w_den_gg[3];
  TH1F w_eff_sg[2];
  TH1F *w_num_sg[2];
  TH1F *w_den_sg[2];

  bool dosingle=0;
  if (diffvariable=="dosingle"){
    dosingle=1;
    diffvariable="singlegamma_eta";    
  }

  if (dosingle){
  for (int i=0; i<2; i++) nfile->GetObject(Form("w_tot_sg_%s",sg_name[i].Data()),w_num_sg[i]);
  for (int i=0; i<2; i++) dfile->GetObject(Form("w_tot_sg_%s",sg_name[i].Data()),w_den_sg[i]);
  for (int i=0; i<2; i++) w_num_sg[i]->Copy(w_eff_sg[i]);
  for (int i=0; i<2; i++) w_eff_sg[i].Divide(w_num_sg[i],w_den_sg[i],1,1,"B");  
  for (int i=0; i<2; i++) { w_eff_sg[i].SetName(Form("w_eff_sg_%s",sg_name[i].Data())); w_eff_sg[i].SetTitle(Form("w_eff_sg_%s",sg_name[i].Data()));}
  }
  else {
  for (int i=0; i<3; i++) nfile->GetObject(Form("w_tot_gg_%s_%s",gg_name[i].Data(),diffvariable.Data()),w_num_gg[i]);
  for (int i=0; i<3; i++) dfile->GetObject(Form("w_tot_gg_%s_%s",gg_name[i].Data(),diffvariable.Data()),w_den_gg[i]);
  for (int i=0; i<3; i++) w_num_gg[i]->Copy(w_eff_gg[i]);
  for (int i=0; i<3; i++) w_eff_gg[i].Divide(w_num_gg[i],w_den_gg[i],1,1,"B");  
  for (int i=0; i<3; i++) { w_eff_gg[i].SetName(Form("w_eff_gg_%s_%s",gg_name[i].Data(),diffvariable.Data())); w_eff_gg[i].SetTitle(Form("w_eff_gg_%s_%s",gg_name[i].Data(),diffvariable.Data()));}
  }


  outf->cd();
  if (!dosingle) for (int i=0; i<3; i++) w_eff_gg[i].Write();  
  else for (int i=0; i<2; i++) w_eff_sg[i].Write();  

};

void divide_eff_histos(TString numerator, TString denominator){
  TFile *outf = new TFile("efficiencies.root","recreate");
  //  diffvariables_list.push_back(TString("dosingle"));
  for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
    divide_eff_histosOne(numerator,denominator,*diffvariable,outf);
  }  
  outf->Close();
};

float efficiency_measure::getpuenergy(int reg, float eta){

  int bin = Choose_bin_eta(fabs(eta),reg);
  float eff_area;
  // MC eff areas
  eff_area = (reg==0) ? eff_areas_EB_mc[bin] : eff_areas_EE_mc[bin];

  return 0.4*0.4*3.14*event_rho*eff_area;

};
