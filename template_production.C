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
  int limit_entries = 1e+6;
  //int limit_entries = -1;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100000==0) std::cout << "Processing entry " << jentry << std::endl;

    if (limit_entries>0){
      if (randomgen->Uniform(0,1) > float(limit_entries)/float(nentries)) continue;
    }

    //    if (event_nRecVtx>5) continue;

    //    if (jentry==1e+4) break;

    // initial kinematic selection
    //    if (dodistribution) if (pholead_pt<40 || photrail_pt<30 || dipho_mgg_photon<80) continue;
    if (dodistribution) if (pholead_pt<40 || photrail_pt<25) continue;
    if (do2dtemplate) if (pholead_pt<40 || photrail_pt<25) continue;
    float cutpt = (mode=="muon") ? 10 : 25;
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
    if (do2ptemplate || do1p1ftemplate || do2ftemplate) {recalc_lead=true; recalc_trail=true;} 
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

    roorho->setVal(event_rho);
    roosigma->setVal(event_sigma);

    Int_t bin_lead = Choose_bin_eta(pholead_SCeta,reg_lead);
    Int_t bin_trail = dodistribution ? Choose_bin_eta(photrail_SCeta,reg_trail) : -999;
    pholead_outvar-=getpuenergy(reg_lead,pholead_SCeta);
    if (dodistribution) photrail_outvar-=getpuenergy(reg_trail,photrail_SCeta);
    if (do2ptemplate || do1p1ftemplate || do2ftemplate) photrail_outvar-=getpuenergy(reg_trail,photrail_SCeta);


//    float scale_lead = 1;
//    float scale_trail = 1;
//    if (mode=="muon") {scale_lead = 1.066;} // not using trail here
//    else if (differentialvariable=="chargediso") {scale_lead = 1+2.5e-3; scale_trail=1+2.5e-3;}
//    else if (differentialvariable=="photoniso") {scale_lead = pholead_scareaSF; scale_trail = photrail_scareaSF;}
//    pholead_outvar/=scale_lead;
//    photrail_outvar/=scale_trail;


    if (pholead_outvar<leftrange) pholead_outvar=leftrange;
    if (photrail_outvar<leftrange) photrail_outvar=leftrange;
    if (pholead_outvar>=rightrange) pholead_outvar=rightrange-1e-5; // overflow in last bin
    if (photrail_outvar>=rightrange) photrail_outvar=rightrange-1e-5; // overflow in last bin
//    if (pholead_outvar==0) pholead_outvar=randomgen->Uniform(0,0.1);
//    if (photrail_outvar==0) photrail_outvar=randomgen->Uniform(0,0.1);

    if (purew_initialized) event_weight=FindNewPUWeight(event_nPU);
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

    histo_pu_nvtx->Fill(event_nPU,event_nRecVtx,weight);


    if (dosignaltemplate||dobackgroundtemplate){

      if (dosignaltemplate){
	template_signal[reg_lead][bin_lead]->Fill(pholead_outvar,weight*ptweight_lead);
	roovar1->setVal(pholead_outvar);
	roovar2->setVal(pholead_outvar);
	roopt1->setVal(pholead_pt);
	roopt2->setVal(pholead_pt);
	rooeta1->setVal(fabs(pholead_SCeta));
	rooeta2->setVal(fabs(pholead_SCeta));
	rooweight->setVal(weight*ptweight_lead);
	roodset_signal[reg_lead][bin_lead][0]->add(RooArgList(*roovar1,*roopt1,*rooeta1,*roorho,*roosigma),weight*ptweight_lead);
	roodset_signal[reg_lead][bin_lead][1]->add(RooArgList(*roovar2,*roopt2,*rooeta2,*roorho,*roosigma),weight*ptweight_lead);
	histo_pt[reg_lead]->Fill(pholead_pt,weight*ptweight_lead);
	histo_eta->Fill(fabs(pholead_SCeta),weight*ptweight_lead);
	histo_pt_eta->Fill(pholead_pt,fabs(pholead_SCeta),weight*ptweight_lead);
	histo_rho_sigma->Fill(event_rho,event_sigma,weight*ptweight_lead);
	//      std::cout << weight << " " << weight*ptweight_lead << std::endl;
      }
      
      if (dobackgroundtemplate){
	if (mode=="sieiesideband"){ 
	  if (fabs(pholead_SCeta)<1.4442 && pholead_sieie>0.011 && pholead_sieie<0.014){
	    template_background[reg_lead][bin_lead]->Fill(pholead_outvar,weight*ptweight_lead);
	    roovar1->setVal(pholead_outvar);
	    roovar2->setVal(pholead_outvar);
	    roopt1->setVal(pholead_pt);
	    roopt2->setVal(pholead_pt);
	    rooeta1->setVal(fabs(pholead_SCeta));
	    rooeta2->setVal(fabs(pholead_SCeta));
	    rooweight->setVal(weight*ptweight_lead);
	    roodset_background[reg_lead][bin_lead][0]->add(RooArgList(*roovar1,*roopt1,*rooeta1,*roorho,*roosigma),weight*ptweight_lead);
	    roodset_background[reg_lead][bin_lead][1]->add(RooArgList(*roovar2,*roopt2,*rooeta2,*roorho,*roosigma),weight*ptweight_lead);
	  }
	  if (fabs(pholead_SCeta)>1.56 && pholead_sieie>0.030 && pholead_sieie<0.031){
	    template_background[reg_lead][bin_lead]->Fill(pholead_outvar,weight*ptweight_lead);
	    roovar1->setVal(pholead_outvar);
	    roovar2->setVal(pholead_outvar);
	    roopt1->setVal(pholead_pt);
	    roopt2->setVal(pholead_pt);
	    rooeta1->setVal(fabs(pholead_SCeta));
	    rooeta2->setVal(fabs(pholead_SCeta));
	    rooweight->setVal(weight*ptweight_lead);
	    roodset_background[reg_lead][bin_lead][0]->add(RooArgList(*roovar1,*roopt1,*rooeta1,*roorho,*roosigma),weight*ptweight_lead);
	    roodset_background[reg_lead][bin_lead][1]->add(RooArgList(*roovar2,*roopt2,*rooeta2,*roorho,*roosigma),weight*ptweight_lead);
	  }
	  if ((fabs(pholead_SCeta)<1.4442 && pholead_sieie>0.011 && pholead_sieie<0.014) || (fabs(pholead_SCeta)>1.56 && pholead_sieie>0.030 && pholead_sieie<0.031)){
	    histo_pt[reg_lead]->Fill(pholead_pt,weight*ptweight_lead);
	    histo_eta->Fill(fabs(pholead_SCeta),weight*ptweight_lead);
	    histo_pt_eta->Fill(pholead_pt,fabs(pholead_SCeta),weight*ptweight_lead);
	    histo_rho_sigma->Fill(event_rho,event_sigma,weight*ptweight_lead);
	  }
	}
	else {
	  template_background[reg_lead][bin_lead]->Fill(pholead_outvar,weight*ptweight_lead);
	  roovar1->setVal(pholead_outvar);
	  roovar2->setVal(pholead_outvar); 
	  roopt1->setVal(pholead_pt);
	  roopt2->setVal(pholead_pt);
	  rooeta1->setVal(fabs(pholead_SCeta));
	  rooeta2->setVal(fabs(pholead_SCeta));
	  rooweight->setVal(weight*ptweight_lead);
	  roodset_background[reg_lead][bin_lead][0]->add(RooArgList(*roovar1,*roopt1,*rooeta1,*roorho,*roosigma),weight*ptweight_lead);
	  roodset_background[reg_lead][bin_lead][1]->add(RooArgList(*roovar2,*roopt2,*rooeta2,*roorho,*roosigma),weight*ptweight_lead);
	  histo_pt[reg_lead]->Fill(pholead_pt,weight*ptweight_lead);
	  histo_eta->Fill(fabs(pholead_SCeta),weight*ptweight_lead);
	  histo_pt_eta->Fill(pholead_pt,fabs(pholead_SCeta),weight*ptweight_lead);
	  histo_rho_sigma->Fill(event_rho,event_sigma,weight*ptweight_lead);
	}
      }
      
    }


    if (do2dtemplate){

	float in1=pholead_outvar;
	float in2=photrail_outvar;
	float ptin1=pholead_pt;
	float ptin2=photrail_pt;
	float etain1=fabs(pholead_SCeta);
	float etain2=fabs(photrail_SCeta);

	bool doswap = false;
	if ((event_ok_for_dataset==0 || event_ok_for_dataset==2) && (randomgen->Uniform()>0.5)) doswap=true;
	if (event_ok_for_dataset==4) doswap=true;
	if (event_ok_for_dataset==3 || event_ok_for_dataset==4) event_ok_for_dataset=1;

	if (doswap){
	  float temp;
	  temp=in1; in1=in2; in2=temp;
	  temp=ptin1; ptin1=ptin2; ptin2=temp;
	  temp=etain1; etain1=etain2; etain2=temp;
	  event_pass12whoisrcone=!event_pass12whoisrcone;
	}

	TString sigorbkg;
	if (do2ptemplate) sigorbkg=TString("sigsig");
	if (do2ftemplate) sigorbkg=TString("bkgbkg");
	if (do1p1ftemplate) sigorbkg=TString("sigbkg");
	if (do1p1ftemplate && event_ok_for_dataset==1 && event_pass12whoisrcone==1) sigorbkg=TString("bkgsig");
	
	roovar1->setVal(in1);
	roovar2->setVal(in2);
	roopt1->setVal(ptin1);
	roopt2->setVal(ptin2);
	rooeta1->setVal(etain1);
	rooeta2->setVal(etain2);
	rooweight->setVal(weight);
	template2d_roodset[get_name_template2d_roodset(event_ok_for_dataset,sigorbkg)]->add(RooArgSet(*roovar1,*roovar2,*roopt1,*roopt2,*rooeta1,*rooeta2,*roorho,*roosigma),weight);

    }


    if (dodistribution && event_ok_for_dataset>-1){

      obs_hist_single[get_name_obs_single(reg_lead,bin_lead)]->Fill(pholead_outvar,weight*ptweight_lead);
      obs_hist_single[get_name_obs_single(reg_trail,bin_trail)]->Fill(photrail_outvar,weight*ptweight_trail);

    
      for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){

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
      
	
	float in1=pholead_outvar;
	float in2=photrail_outvar;
	float ptin1=pholead_pt;
	float ptin2=photrail_pt;
	float etain1=fabs(pholead_SCeta);
	float etain2=fabs(photrail_SCeta);

	bool doswap=false;

	if ((event_ok_for_dataset==0 || event_ok_for_dataset==2) && (randomgen->Uniform()>0.5)) doswap=true;
	
	if (event_ok_for_dataset==4) doswap=true;

	if (event_ok_for_dataset==3 || event_ok_for_dataset==4) event_ok_for_dataset=1;

	if (doswap){
	  float temp;
	  temp=in1; in1=in2; in2=temp;
	  temp=ptin1; ptin1=ptin2; ptin2=temp;
	  temp=etain1; etain1=etain2; etain2=temp;
	}
      
	obs_hist[get_name_obs(event_ok_for_dataset,*diffvariable,bin_couple)]->Fill(in1,in2,weight);
	roovar1->setVal(in1);
	roovar2->setVal(in2);
	roopt1->setVal(ptin1);
	roopt2->setVal(ptin2);
	rooeta1->setVal(etain1);
	rooeta2->setVal(etain2);
	rooweight->setVal(weight);
	obs_roodset[get_name_obs_roodset(event_ok_for_dataset,*diffvariable,bin_couple)]->add(RooArgSet(*roovar1,*roovar2,*roopt1,*roopt2,*rooeta1,*rooeta2,*roorho,*roosigma),weight);
		 
      }
      
    }
    



  } // end event loop
  std::cout << "ended event loop" << std::endl;
};


#endif

void gen_templates(TString filename="input.root", TString mode="", bool isdata=1, const char* outfile="out.root", TString differentialvariable="photoniso", TString targetpufile="", bool doreweight=false, TString dset1="", TString dset2="", TString temp1="", TString temp2=""){
  
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
  treename[5] = TString("Tree_doublerandomcone_sel");
  treename[6] = TString("Tree_singlephoton_sel_noimpingingcut");
  treename[7] = TString("Tree_onlypreselection");
  treename[8] = TString("Tree_sieiesideband_sel");
  treename[9] = TString("Tree_combisosideband_sel");
  treename[10] = TString("Tree_muoncone");

  treename[0] =  TString("Tree_standard_sel");
  treename[1] =  TString("Tree_signal_template");
  treename[2] =  TString("Tree_background_template");
  treename[3] =  TString("Tree_DY_sel");
  treename[4] =  TString("Tree_randomcone_signal_template");
  treename[5] =  TString("Tree_doublerandomcone_sel");
  treename[6] =  TString("Tree_DYnocombisowithinvmasscut_sel");
  treename[7] =  TString("Tree_onlypreselection");
  treename[8] =  TString("Tree_sieiesideband_sel");
  treename[9] =  TString("Tree_nocombisocut_sel");
  treename[10] = TString("Tree_randomcone_nocombisocut");
  treename[11] = TString("Tree_muoncone");
  treename[12] = TString("Tree_randomconesideband_sel");
  treename[13] = TString("Tree_doublesieiesideband_sel");


  TString treename_chosen="";
  if (mode=="standard") treename_chosen=treename[0];
  if (mode=="signal") treename_chosen=treename[1];
  if (mode=="background") treename_chosen=treename[2];
  if (mode=="randomcone") treename_chosen=treename[4];
  if (mode=="sieiesideband") treename_chosen=treename[8];
  if (mode=="muon") treename_chosen=treename[11];
  if (mode=="sigsig") treename_chosen=treename[5];
  if (mode=="sigbkg") treename_chosen=treename[12];
  if (mode=="bkgbkg") treename_chosen=treename[13];

  file->GetObject(treename_chosen.Data(),t);

  std::cout << "Processing selection " << treename_chosen.Data() << std::endl;
  template_production *temp = new template_production(t);
  temp->Setup(isdata,mode,differentialvariable);
  if (targetpufile!="") temp->InitializeNewPUReweighting(filename,targetpufile);
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



void get_eff_area(TString filename, bool doEB, TString comp){

  const char* outfile=Form("puscaling_%s.root",(doEB) ? "EB" : "EE");
  
  TFile *outF = TFile::Open(outfile,"recreate");
  

  TFile *file[1];
  file[0] = TFile::Open(filename.Data());

  TTree *t;

  TString treename;
  //  treename = TString("Tree_randomcone_signal_template");
  //  treename[3] = TString("Tree_DYnocombisowithinvmasscut_sel");
  treename = TString("Tree_randomcone_nocombisocut");


  std::vector<std::vector<TProfile*> > output;

  file[0]->GetObject(treename.Data(),t);
  template_production *temp = new template_production(t);
  output=temp->GetPUScaling(doEB,comp);

//  output[0]->Print();
//  output[1]->Print();
//  output[2]->Print();

  TF1 *f_iso[n_bins];
  TF1 *f_rho[n_bins];
  TF1 *f_iso_pu[n_bins];


  for (int i=0; i<n_bins; i++) f_iso[i] = new TF1(Form("f_iso_%d",i),"pol1(0)",5,20);
  for (int i=0; i<n_bins; i++) f_rho[i] = new TF1(Form("f_rho_%d",i),"pol1(0)",5,20);
  for (int i=0; i<n_bins; i++) f_iso_pu[i] = new TF1(Form("f_iso_pu_%d",i),"pol1(0)",5,20);

  for (int i=0; i<n_bins; i++){
    output[0][i]->Fit(f_iso[i],"R");
    output[1][i]->Fit(f_rho[i],"R");
    output[2][i]->Fit(f_iso_pu[i],"R");  
    std::cout << "RESULT (eff. area)" << std::endl << f_iso[i]->GetParameter(1)/f_rho[i]->GetParameter(1) << std::endl;
  }


  for (int i=0; i<n_bins; i++) gROOT->ProcessLine(Form(".! echo %s %s bin%d %e >> puscaling.txt",comp.Data(),(doEB) ? "EB" : "EE",i,f_iso[i]->GetParameter(1)/f_rho[i]->GetParameter(1)));

  outF->cd();

  for (int i=0; i<n_bins; i++){
    output[0][i]->Write();
    output[1][i]->Write();
    output[2][i]->Write();
    f_iso[i]->Write();
    f_rho[i]->Write();
    f_iso_pu[i]->Write();
  }

  file[0]->Close();

  outF->Close();

};

std::vector<std::vector<TProfile*> > template_production::GetPUScaling(bool doEB, TString diffvar){

  TProfile *prof_iso[n_bins];
  TProfile *prof_rho[n_bins];
  TProfile *prof_iso_pu[n_bins];

  for (int i=0; i<n_bins; i++) prof_iso[i] = new TProfile(Form("prof_iso_%d",i),Form("prof_iso_%d",i),30,0,30);
  for (int i=0; i<n_bins; i++) prof_rho[i] = new TProfile(Form("prof_rho_%d",i),Form("prof_rho_%d",i),30,0,30);
  for (int i=0; i<n_bins; i++) prof_iso_pu[i] = new TProfile(Form("prof_iso_pu_%d",i),Form("prof_iso_pu_%d",i),30,0,30);

  //  prof_rho->SetLineColor(kRed);
  //  prof_iso_pu->SetLineColor(kBlue);

  Init();

  if (fChain == 0){
    std::cout << "No chain!" << std::endl;
    return std::vector<std::vector<TProfile*> >();
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

    //    bool isbarrel = (fabs(pholead_SCeta)<1.4442);

    event_weight = FindNewPUWeight(event_nPU);


    float pholead_phoiso=-999;
    { // recalc phoiso w/cleaning
      float et_recalc = 0;
      float e_recalc = 0;
      int number_recalc = 0;


      for (int i=0; i<pholead_Npfcandphotonincone; i++) {

	float et=pholead_photonpfcandets[i];
	float e=pholead_photonpfcandenergies[i];
	float deta=pholead_photonpfcanddetas[i];
	float dphi=pholead_photonpfcanddphis[i];
	//	float dR=sqrt(deta*deta+dphi*dphi);
	float eta=fabs(TMath::ACosH(e/et));

	if (eta>1.4442 && eta<1.56) continue;
	if (eta>2.5) continue;

#include "cleaning.cc"

	if (fabs(pholead_SCeta)<1.4442 && eta>1.4442) continue;
	if (fabs(pholead_SCeta)>1.56 && eta<1.56) continue;
	et_recalc+=et;
	e_recalc+=e;
	number_recalc++;

      }


      pholead_phoiso=et_recalc;
    }



    if (diffvar=="photoniso"){
      //      pholead_outvar=pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;
      pholead_outvar=pholead_phoiso;
    }
    else if (diffvar=="combiso"){
      //      pholead_outvar=pholead_pho_Cone04PFCombinedIso;
      pholead_outvar=pholead_phoiso+pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01+pholead_pho_Cone04NeutralHadronIso_mvVtx;
    }
    else if (diffvar=="chargediso"){
      pholead_outvar=pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
    }
    else if (diffvar=="neutraliso"){
      pholead_outvar=pholead_pho_Cone04NeutralHadronIso_mvVtx;
    }


    if (pholead_outvar==-999) continue;

    //    if (pholead_PhoMCmatchexitcode!=1 && pholead_PhoMCmatchexitcode!=2) continue;
    //    if (dipho_mgg_photon>95 || dipho_mgg_photon<85) continue;

    if (pholead_pt<25) continue;

    Float_t weight=event_luminormfactor*event_Kfactor*event_weight;

    if (!doEB){
      if (fabs(pholead_SCeta)<1.56) continue;
    }
    else {
      if (fabs(pholead_SCeta)>1.4442) continue;
    }

    Int_t bin_lead = Choose_bin_eta(pholead_SCeta,(doEB) ? 0 : 1);

    prof_iso[bin_lead]->Fill(event_nRecVtx,pholead_outvar,weight);

    
    prof_iso_pu[bin_lead]->Fill(event_nRecVtx,pholead_outvar,weight);


    prof_rho[bin_lead]->Fill(event_nRecVtx,event_rho*3.14*0.4*0.4,weight);

    

  }

  std::vector<std::vector<TProfile*> > out;

  out.resize(3);
  out[0].resize(n_bins);
  out[1].resize(n_bins);
  out[2].resize(n_bins);

  for (int i=0; i<n_bins; i++) out[0][i]=prof_iso[i];
  for (int i=0; i<n_bins; i++) out[1][i]=prof_rho[i];
  for (int i=0; i<n_bins; i++) out[2][i]=prof_iso_pu[i];

//  prof_iso->Print();
//  prof_rho->Print();
//  prof_iso_pu->Print();

  return out;

};

float template_production::getpuenergy(int reg, float eta){

  //  return 0;

  int bin = Choose_bin_eta(fabs(eta),reg);
  float eff_area = (reg==0) ? eff_areas_EB[bin] : eff_areas_EE[bin];

  return 0.4*0.4*3.14*event_rho*eff_area;

};
