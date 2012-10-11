#ifndef template_production_cxx
#define template_production_cxx
#include "template_production.h"

using namespace std;


void template_production::Loop()
{
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
  if (!initialized){
    std::cout << "Not initialized! Call Setup() first." << std::endl;
    return;
  }


  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100000==0) std::cout << "Processing entry " << jentry << std::endl;

    // initial kinematic selection
    //    if (dodistribution) if (pholead_pt<40 || photrail_pt<30 || dipho_mgg_photon<80) continue;
    if (dodistribution) if (pholead_pt<40 || photrail_pt<30) continue;
    float cutpt = (mode=="muon") ? 10 : 30;
    if (dosignaltemplate || dobackgroundtemplate) if (pholead_pt<cutpt) continue;

    if (isdata && event_CSCTightHaloID>0) continue;
    //if (mode!="muon" && event_NMuons>0) continue;
    //if (mode!="muon" && event_NMuonsTot>0) continue;
    //    if (!isdata && (event_PUOOTnumInteractionsEarly<20 && event_PUOOTnumInteractionsLate<20)) continue;
    //    if (!isdata && event_PUOOTnumInteractionsEarly>3) continue;

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
    else std::cout << "We have a problem here!!!" << std::endl;





//    if (differentialvariable=="photoniso"){
//      pholead_outvar=pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;
//      photrail_outvar=photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;
//      candcounter=pholead_Npfcandphotonincone;
//    }
//    else if (differentialvariable=="combiso"){
//      pholead_outvar=pholead_pho_Cone04PFCombinedIso;
//      photrail_outvar=photrail_pho_Cone04PFCombinedIso;
//      candcounter=pholead_Npfcandphotonincone+pholead_Npfcandchargedincone+pholead_Npfcandneutralincone;
//    }
//    else if (differentialvariable=="chargediso"){
//      pholead_outvar=pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
//      photrail_outvar=photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
//      candcounter=pholead_Npfcandchargedincone;
//    }
//    else if (differentialvariable=="neutraliso"){
//      pholead_outvar=pholead_pho_Cone04NeutralHadronIso_mvVtx;
//      photrail_outvar=photrail_pho_Cone04NeutralHadronIso_mvVtx;
//      candcounter=pholead_Npfcandneutralincone;
//    }


    pholead_outvar = -999;
    photrail_outvar = -999;
    candcounter = -999;

    bool recalc_lead = false;
    bool recalc_trail = false;

    if (dosignaltemplate||dobackgroundtemplate) recalc_lead=true;
    if (dodistribution) {recalc_lead=true; recalc_trail=true;}


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


    if (pholead_outvar<leftrange) pholead_outvar=leftrange;
    if (photrail_outvar<leftrange) photrail_outvar=leftrange;
    if (pholead_outvar>=rightrange) pholead_outvar=rightrange-1e-5; // overflow in last bin
    if (photrail_outvar>=rightrange) photrail_outvar=rightrange-1e-5; // overflow in last bin

    // OCCHIO ALL'ORDINE DI QUANDO DI METTE IL TAGLIO LEFT/RIGHTRANGE QUI!!!!


//    float eff_area[2] = {0,0};
//    if (mode!="muon"){
//      // use this to UNDO effarea corrections
//      if (differentialvariable=="photoniso") {eff_area[0] = 0.221; eff_area[1] = 0.130;}
//      if (differentialvariable=="chargediso") {eff_area[0] = 0.016; eff_area[1] = 0.017;}
//      if (differentialvariable=="neutraliso") {eff_area[0] = 0.097; eff_area[1] = 0.132;}
//    }
//    float puincone = 0.4*0.4*3.14*event_rho;
//    pholead_outvar+=puincone*eff_area[reg_lead];
//    photrail_outvar+=puincone*eff_area[reg_trail];




//    float scale_lead = 1;
//    float scale_trail = 1;
//    if (mode=="muon") {scale_lead = 1.066;} // not using trail here
//    else if (differentialvariable=="chargediso") {scale_lead = 1+2.5e-3; scale_trail=1+2.5e-3;}
//    else if (differentialvariable=="photoniso") {scale_lead = pholead_scareaSF; scale_trail = photrail_scareaSF;}
//    pholead_outvar/=scale_lead;
//    photrail_outvar/=scale_trail;


    Float_t weight=event_luminormfactor*event_Kfactor*event_weight;
    float ptweight_lead = 1;
    float ptweight_trail = 1;

    if (do_pt_reweighting) ptweight_lead*=FindPtWeight(pholead_pt,pholead_SCeta);
    if (do_eta_reweighting) ptweight_lead*=FindEtaWeight(pholead_SCeta);
    if (do_pt_reweighting) ptweight_trail*=FindPtWeight(photrail_pt,photrail_SCeta);
    if (do_eta_reweighting) ptweight_trail*=FindEtaWeight(photrail_SCeta);

    if (do_pt_eta_reweighting) {
      ptweight_lead*=FindPtEtaWeight(pholead_pt,pholead_SCeta);
      ptweight_trail*=FindPtEtaWeight(photrail_pt,photrail_SCeta);
    }


    if (dosignaltemplate||dobackgroundtemplate){

      Int_t bin_lead = Choose_bin_eta(pholead_SCeta,reg_lead);
   
      if (dosignaltemplate){
	template_signal[reg_lead][bin_lead]->Fill(pholead_outvar,weight*ptweight_lead);
	histo_pt[reg_lead]->Fill(pholead_pt,weight*ptweight_lead);
	histo_eta->Fill(fabs(pholead_SCeta),weight*ptweight_lead);
	histo_pt_eta->Fill(pholead_pt,fabs(pholead_SCeta),weight*ptweight_lead);
	//      std::cout << weight << " " << weight*ptweight_lead << std::endl;
      }
      
      if (dobackgroundtemplate){
	if (mode=="sieiesideband"){ 
	  if (fabs(pholead_SCeta)<1.4442 && pholead_sieie>0.011 && pholead_sieie<0.014) template_background[reg_lead][bin_lead]->Fill(pholead_outvar,weight*ptweight_lead);
	  if (fabs(pholead_SCeta)>1.56 && pholead_sieie>0.030 && pholead_sieie<0.031) template_background[reg_lead][bin_lead]->Fill(pholead_outvar,weight*ptweight_lead);
	  if ((fabs(pholead_SCeta)<1.4442 && pholead_sieie>0.011 && pholead_sieie<0.014) || (fabs(pholead_SCeta)>1.56 && pholead_sieie>0.030 && pholead_sieie<0.031)){
	    histo_pt[reg_lead]->Fill(pholead_pt,weight*ptweight_lead);
	    histo_eta->Fill(fabs(pholead_SCeta),weight*ptweight_lead);
	    histo_pt_eta->Fill(pholead_pt,fabs(pholead_SCeta),weight*ptweight_lead);
	  }
	}
	else {
	  template_background[reg_lead][bin_lead]->Fill(pholead_outvar,weight*ptweight_lead);
	  histo_pt[reg_lead]->Fill(pholead_pt,weight*ptweight_lead);
	  histo_eta->Fill(fabs(pholead_SCeta),weight*ptweight_lead);
	  histo_pt_eta->Fill(pholead_pt,fabs(pholead_SCeta),weight*ptweight_lead);
	}
      }
      
    }

    if (dodistribution && event_ok_for_dataset>-1){
    
      for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){

	Int_t bin_lead = Choose_bin_eta(pholead_SCeta,reg_lead);
	Int_t bin_trail = Choose_bin_eta(photrail_SCeta,reg_trail);
	Int_t bin_couple = -999;
	
	if (*diffvariable==TString("invmass")) bin_couple = Choose_bin_invmass(dipho_mgg_photon,event_ok_for_dataset);
	if (*diffvariable==TString("diphotonpt")){
	  float px = pholead_px+photrail_px;
	  float py = pholead_py+photrail_py;
	  float pt = sqrt(px*px+py*py);
	  bin_couple = Choose_bin_diphotonpt(pt,event_ok_for_dataset);
	}
	if (*diffvariable==TString("costhetastar")){
	  TLorentzVector pho1(pholead_px,pholead_py,pholead_pz,pholead_energy);
	  TLorentzVector pho2(photrail_px,photrail_py,photrail_pz,photrail_energy);
	  TVector3 boost = (pho1+pho2).BoostVector();
	  TLorentzVector boostedpho1 = pho1;
	  boostedpho1.Boost(-boost);
	  float thetastar1 = boostedpho1.Angle(boost);
	  if (thetastar1>TMath::Pi()/2) thetastar1 = TMath::Pi()-thetastar1;
	  bin_couple = Choose_bin_costhetastar(TMath::Cos(thetastar1),event_ok_for_dataset);
	}
	if (*diffvariable==TString("dphi")){
	  float phi1 = pholead_SCphi;
	  float phi2 = photrail_SCphi;
	  float dphi = AbsDeltaPhi(phi1,phi2);
	  bin_couple = Choose_bin_dphi(dphi,event_ok_for_dataset);
	}
      
	obs_hist_single[get_name_obs_single(reg_lead,bin_lead)]->Fill(pholead_outvar,weight*ptweight_lead);
	obs_hist_single[get_name_obs_single(reg_trail,bin_trail)]->Fill(photrail_outvar,weight*ptweight_trail);
	
	float in1=pholead_outvar;
	float in2=photrail_outvar;

	bool doswap=false;

	if ((event_ok_for_dataset==0 || event_ok_for_dataset==2) && (randomgen->Uniform()>0.5)) doswap=true;
	
	if (event_ok_for_dataset==4) doswap=true;

	if (event_ok_for_dataset==3 || event_ok_for_dataset==4) event_ok_for_dataset=1;

	if (doswap){
	  in2=pholead_outvar;
	  in1=photrail_outvar;
	}
      
	obs_hist[get_name_obs(event_ok_for_dataset,*diffvariable,bin_couple)]->Fill(in1,in2,weight);
		 
		 
      }
      
    }
    



  } // end event loop
  std::cout << "ended event loop" << std::endl;
};


#endif

void gen_templates(TString filename="input.root", TString mode="", bool isdata=1, const char* outfile="out.root", TString differentialvariable="photoniso", bool doreweight=false, TString dset1="", TString dset2="", TString temp1="", TString temp2=""){
  
  TFile *outF = TFile::Open(outfile,"recreate");
  outF->Close();

  TFile *file;
  file = TFile::Open(filename.Data(),"read");

  TTree *t;

  TString treename[11];
  treename[0] = TString("Tree_standard_sel");
  treename[1] = TString("Tree_signal_template");
  treename[2] = TString("Tree_background_template");
  treename[3] = TString("Tree_DY_sel");
  treename[4] = TString("Tree_randomcone_signal_template");
  treename[5] = TString("Tree_impinging_track_template");
  treename[6] = TString("Tree_singlephoton_sel_noimpingingcut");
  treename[7] = TString("Tree_onlypreselection");
  treename[8] = TString("Tree_sieiesideband_sel");
  treename[9] = TString("Tree_combisosideband_sel");
  treename[10] = TString("Tree_muoncone");

  TString treename_chosen="";
  if (mode=="standard") treename_chosen=treename[0];
  if (mode=="signal") treename_chosen=treename[1];
  if (mode=="background") treename_chosen=treename[2];
  if (mode=="randomcone") treename_chosen=treename[4];
  if (mode=="impinging") treename_chosen=treename[5];
  if (mode=="sieiesideband") treename_chosen=treename[8];
  if (mode=="combisosideband") treename_chosen=treename[9];
  if (mode=="muon") treename_chosen=treename[10];

  file->GetObject(treename_chosen.Data(),t);

  std::cout << "Processing selection " << treename_chosen.Data() << std::endl;
  template_production *temp = new template_production(t);
  temp->Setup(isdata,mode,differentialvariable);
  if (!doreweight) {
    temp->SetNoPtReweighting();
    temp->SetNoEtaReweighting();
    temp->SetNoPtEtaReweighting();
  }
  else {
    temp->SetNoPtReweighting();
    temp->SetNoEtaReweighting();
    //    temp->Initialize_Pt_Reweighting(dset1,dset2,temp1,temp2);
    //    temp->Initialize_Eta_Reweighting(dset1,dset2,temp1,temp2);
    temp->Initialize_Pt_Eta_Reweighting(dset1,dset2,temp1,temp2);
  }
  temp->Loop();
  std::cout << "exited from event loop" << std::endl;
  temp->WriteOutput(outfile,treename_chosen.Data());
  std::cout << "written output" << std::endl;

  file->Close();

};



void get_eff_area(bool doEB, TString comp){

  const char* outfile=Form("puscaling_%s.root",(doEB) ? "EB" : "EE");
  
  TFile *outF = TFile::Open(outfile,"recreate");
  
  TString dir("./");

  TString filenameMC=dir;
  filenameMC.Append("outfile.root");

  TFile *file[1];
  file[0] = TFile::Open(filenameMC.Data());

  TTree *t;

  TString treename;
  //  treename = TString("Tree_randomcone_signal_template");
  //  treename[3] = TString("Tree_DYnocombisowithinvmasscut_sel");
  treename = TString("Tree_randomcone_nocombisocut");

  TProfile** output;

  file[0]->GetObject(treename.Data(),t);
  template_production *temp = new template_production(t);
  output=temp->GetPUScaling(doEB,comp);

  output[0]->Print();
  output[1]->Print();
  output[2]->Print();

  TF1 *f_iso = new TF1("f_iso","pol1(0)",5,20);
  TF1 *f_rho = new TF1("f_rho","pol1(0)",5,20);
  TF1 *f_iso_pu = new TF1("f_iso_pu","pol1(0)",5,20);

  output[0]->Fit(f_iso,"R");
  output[1]->Fit(f_rho,"R");
  output[2]->Fit(f_iso_pu,"R");  
   
  std::cout << "RESULT (eff. area)" << std::endl << f_iso->GetParameter(1)/f_rho->GetParameter(1) << std::endl;

  outF->cd();
  output[0]->Write();
  output[1]->Write();
  output[2]->Write();
  f_iso->Write();
  f_rho->Write();
  f_iso_pu->Write();

  file[0]->Close();

  outF->Close();

};

TProfile** template_production::GetPUScaling(bool doEB, TString diffvar){


  TProfile *prof_iso = new TProfile("prof_iso","prof_iso",30,0,30);
  TProfile *prof_rho = new TProfile("prof_rho","prof_rho",30,0,30);
  TProfile *prof_iso_pu = new TProfile("prof_iso_pu","prof_iso_pu",30,0,30);

  prof_rho->SetLineColor(kRed);
  prof_iso_pu->SetLineColor(kBlue);

  Init();

  if (fChain == 0){
    std::cout << "No chain!" << std::endl;
    return NULL;
  } 

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100000==0) { 
      std::cout << "Processing entry " << jentry << std::endl;
    }

    float puincone = 0.4*0.4*3.14*event_rho;

    //    bool isbarrel = (fabs(pholead_SCeta)<1.4442);

    float eff_area = 0;

    // use this to UNDO effarea corrections
//    if (diffvar=="photoniso") eff_area = isbarrel ? 0.221 : 0.130;
//    if (diffvar=="chargediso") eff_area = isbarrel ? 0.016 : 0.017;
//    if (diffvar=="neutraliso") eff_area = isbarrel ? 0.097 : 0.132;

    if (diffvar=="photoniso"){
      pholead_outvar=pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;
    }
    else if (diffvar=="combiso"){
      pholead_outvar=pholead_pho_Cone04PFCombinedIso;
    }
    else if (diffvar=="chargediso"){
      pholead_outvar=pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
    }
    else if (diffvar=="neutraliso"){
      pholead_outvar=pholead_pho_Cone04NeutralHadronIso_mvVtx;
    }
    pholead_outvar+=puincone*eff_area;

    if (pholead_outvar==-999) continue;

    //    if (pholead_PhoMCmatchexitcode!=1 && pholead_PhoMCmatchexitcode!=2) continue;
    //    if (dipho_mgg_photon>95 || dipho_mgg_photon<85) continue;

    if (pholead_pt<30) continue;

    Float_t weight=event_luminormfactor*event_Kfactor*event_weight;

    if (!doEB){
      if (fabs(pholead_SCeta)<1.4442) continue;
    }
    else {
      if (fabs(pholead_SCeta)>1.56) continue;
    }

    prof_iso->Fill(event_nRecVtx,pholead_outvar,weight/2);

     float eff_area_fraction_EB=0;
     float eff_area_fraction_EE=0;
     const float dR=0.4;
     pholead_outvar-=event_rho*3.14*dR*dR*((fabs(pholead_SCeta)<1.4442) ? eff_area_fraction_EB : eff_area_fraction_EE);
    
    prof_iso_pu->Fill(event_nRecVtx,pholead_outvar,weight/2);


    prof_rho->Fill(event_nRecVtx,event_rho*3.14*0.4*0.4,weight);

    

  }

  TProfile** out = new TProfile*[3];
  out[0]=prof_iso;
  out[1]=prof_rho;
  out[2]=prof_iso_pu;

  prof_iso->Print();
  prof_rho->Print();
  prof_iso_pu->Print();

  return out;

};

