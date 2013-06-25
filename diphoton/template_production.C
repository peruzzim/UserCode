#ifndef template_production_cxx
#define template_production_cxx
#include "template_production.h"

using namespace std;


void template_production::Loop(int maxevents)
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
  int limit_entries = maxevents;
  //int limit_entries = 1e+4;
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
    if (mode=="zee") {
      if (fabs(dipho_mgg_photon-91.2)>10) continue;
    }
    //    if (mode=="background") if (pholead_PhoMCmatchexitcode!=3) continue;

    //    if (mode=="background") if (pholead_PhoMCmatchexitcode!=1 && pholead_PhoMCmatchexitcode!=2) continue;

    //    if (mode=="background") if (pholead_PhoMCmatchexitcode!=0) continue;  

    //    if (isdata && event_CSCTightHaloID>0) continue;
    //if (mode!="muon" && event_NMuons>0) continue;
    //if (mode!="muon" && event_NMuonsTot>0) continue;
    //    if (!isdata && (event_PUOOTnumInteractionsEarly<20 && event_PUOOTnumInteractionsLate<20)) continue;
    //    if (!isdata && event_PUOOTnumInteractionsEarly>3) continue;

    if (mode=="fragmentation") if (pholead_PhoMCmatchexitcode!=1) continue;
    if (mode=="nofragmentation") if (pholead_PhoMCmatchexitcode!=2) continue;

    if (mode=="cutPFchargediso_signal" || mode=="cutPFchargediso_background" || mode=="cutPFchargediso_randomcone" || mode=="cutPFchargediso_sieiesideband") if (pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01>0.1) continue;

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

    //    std::cout << event_ok_for_dataset << " " << fabs(pholead_SCeta) << " " << fabs(photrail_SCeta) << std::endl;


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

    //    std::cout << pholead_outvar << " " << photrail_outvar << std::endl;
    
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

    if (recalc_lead){
      if (pholead_outvar<-100) std::cout << "PROBLEM WITH ISOLATION CALCULATION!!!" << std::endl;
      assert (pholead_outvar>=-100);
      if (pholead_outvar<leftrange) {/*std::cout << "Warning: fixing underflow " << pholead_outvar << std::endl;*/ pholead_outvar=leftrange+1e-5;}
      //      if (pholead_outvar>=rightrange) continue;
      if (pholead_outvar>=rightrange) pholead_outvar=rightrange-1e-5; // overflow in last bin 
    }
    if (recalc_trail){
      if (photrail_outvar<-100) std::cout << "PROBLEM WITH ISOLATION CALCULATION!!!" << std::endl;
      assert (photrail_outvar>=-100);
      if (photrail_outvar<leftrange) {/*std::cout << "Warning: fixing underflow " << photrail_outvar << std::endl;*/ photrail_outvar=leftrange+1e-5;}
      //      if (photrail_outvar>=rightrange) continue;
      if (photrail_outvar>=rightrange) photrail_outvar=rightrange-1e-5; // overflow in last bin 
    }

//    if (pholead_outvar<-10 || photrail_outvar<-10) std::cout << "PROBLEM WITH ISOLATION CALCULATION!!!" << std::endl;
//    if (pholead_outvar<leftrange) pholead_outvar=leftrange;
//    if (photrail_outvar<leftrange) photrail_outvar=leftrange;
//    if (pholead_outvar>=rightrange) pholead_outvar=rightrange-1e-5; // overflow in last bin
//    if (photrail_outvar>=rightrange) photrail_outvar=rightrange-1e-5; // overflow in last bin
//    if (pholead_outvar<leftrange)   continue;
//    if (pholead_outvar>=rightrange) continue;
//    if (dodistribution||do2dtemplate){
//      if (photrail_outvar<leftrange)  continue;
//      if (photrail_outvar>=rightrange)continue;
//    }


//    if (purew_initialized) event_weight=FindNewPUWeight(event_nPU);
    Float_t weight=event_luminormfactor*event_Kfactor*event_weight;
    float ptweight_lead = 1;
    float ptweight_trail = 1;

    if (mode=="standard_2frag" || mode=="2pgen_2frag" || mode=="1p1fbothgen_2frag" || mode=="1pgen1fside_2frag") {
      if (pholead_PhoMCmatchexitcode==1 && pholead_GenPhotonIsoDR04<5) weight*=2;
      if (photrail_PhoMCmatchexitcode==1 && photrail_GenPhotonIsoDR04<5) weight*=2;
    }
    if (mode=="signal_2frag"){
      if (pholead_PhoMCmatchexitcode==1 && pholead_GenPhotonIsoDR04<5) weight*=2;
    }

//    if (do_pt_reweighting) ptweight_lead*=FindPtWeight(pholead_pt,pholead_SCeta);
//    if (do_eta_reweighting) ptweight_lead*=FindEtaWeight(pholead_SCeta);
//    if (do_pt_reweighting) ptweight_trail*=FindPtWeight(photrail_pt,photrail_SCeta);
//    if (do_eta_reweighting) ptweight_trail*=FindEtaWeight(photrail_SCeta);
//
//    if (do_pt_eta_reweighting) {
//      ptweight_lead*=FindPtEtaWeight(pholead_pt,pholead_SCeta);
//      ptweight_trail*=FindPtEtaWeight(photrail_pt,photrail_SCeta);
//    }
//

    histo_pu_nvtx->Fill(event_nPU,event_nRecVtx,weight);


    if (dosignaltemplate||dobackgroundtemplate){

      if (dosignaltemplate){
	template_signal[reg_lead][bin_lead]->Fill(pholead_outvar,weight*ptweight_lead);
	roovar1->setVal(pholead_outvar);
	roovar2->setVal(pholead_outvar);
	roopt1->setVal(pholead_pt);
	roopt2->setVal(pholead_pt);
	roosieie1->setVal(pholead_sieie);
	roosieie2->setVal(pholead_sieie);
	rooeta1->setVal(fabs(pholead_SCeta));
	rooeta2->setVal(fabs(pholead_SCeta));
	rooweight->setVal(weight*ptweight_lead);
	roodset_signal[reg_lead][bin_lead][0]->add(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma),weight*ptweight_lead);
	roodset_signal[reg_lead][bin_lead][1]->add(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma),weight*ptweight_lead);
	histo_pt[reg_lead]->Fill(pholead_pt,weight*ptweight_lead);
	histo_eta->Fill(fabs(pholead_SCeta),weight*ptweight_lead);
	histo_pt_eta->Fill(pholead_pt,fabs(pholead_SCeta),weight*ptweight_lead);
	histo_rho_sigma->Fill(event_rho,event_sigma,weight*ptweight_lead);
	//	std::cout << weight << " " << weight*ptweight_lead << std::endl;
      }
      
      if (dobackgroundtemplate){
	  template_background[reg_lead][bin_lead]->Fill(pholead_outvar,weight*ptweight_lead);
	  roovar1->setVal(pholead_outvar);
	  roovar2->setVal(pholead_outvar); 
	  roopt1->setVal(pholead_pt);
	  roopt2->setVal(pholead_pt);
	  roosieie1->setVal(pholead_sieie);
	  roosieie2->setVal(pholead_sieie);
	  rooeta1->setVal(fabs(pholead_SCeta));
	  rooeta2->setVal(fabs(pholead_SCeta));
	  rooweight->setVal(weight*ptweight_lead);
	  roodset_background[reg_lead][bin_lead][0]->add(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma),weight*ptweight_lead);
	  roodset_background[reg_lead][bin_lead][1]->add(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma),weight*ptweight_lead);
	  histo_pt[reg_lead]->Fill(pholead_pt,weight*ptweight_lead);
	  histo_eta->Fill(fabs(pholead_SCeta),weight*ptweight_lead);
	  histo_pt_eta->Fill(pholead_pt,fabs(pholead_SCeta),weight*ptweight_lead);
	  histo_rho_sigma->Fill(event_rho,event_sigma,weight*ptweight_lead);
      }
      
    }

    if (do2dtemplate){

	int event_ok_for_dataset_local = event_ok_for_dataset;

	float in1=pholead_outvar;
	float in2=photrail_outvar;
	float ptin1=pholead_pt;
	float ptin2=photrail_pt;
	float sieiein1=pholead_sieie;
	float sieiein2=photrail_sieie;	
	float etain1=fabs(pholead_SCeta);
	float etain2=fabs(photrail_SCeta);

	bool doswap = false;
	if ((event_ok_for_dataset_local==0 || event_ok_for_dataset_local==2) && (randomgen->Uniform()>0.5)) doswap=true;
	if (event_ok_for_dataset_local==4) doswap=true;
	if (event_ok_for_dataset_local==3 || event_ok_for_dataset_local==4) event_ok_for_dataset_local=1;

	if (doswap){
	  float temp;
	  temp=in1; in1=in2; in2=temp;
	  temp=ptin1; ptin1=ptin2; ptin2=temp;
	  temp=sieiein1; sieiein1=sieiein2; sieiein2=temp;
	  temp=etain1; etain1=etain2; etain2=temp;
	  event_pass12whoissiglike=!event_pass12whoissiglike;
	}

	TString sigorbkg;
	if (do2ptemplate) sigorbkg=TString("sigsig");
	if (do2ftemplate) sigorbkg=TString("bkgbkg");
	if (do1p1ftemplate && event_pass12whoissiglike==0) sigorbkg=TString("sigbkg");
	if (do1p1ftemplate && event_pass12whoissiglike==1) sigorbkg=TString("bkgsig");
	
	roovar1->setVal(in1);
	roovar2->setVal(in2);
	roopt1->setVal(ptin1);
	roopt2->setVal(ptin2);
	roosieie1->setVal(sieiein1);
	roosieie2->setVal(sieiein2);
	rooeta1->setVal(etain1);
	rooeta2->setVal(etain2);
	rooweight->setVal(weight);
	RooArgSet args(*roovar1,*roovar2,*roopt1,*roopt2,*roosieie1,*roosieie2,*rooeta1,*rooeta2);
	args.add(RooArgSet(*roorho,*roosigma));
	template2d_roodset[get_name_template2d_roodset(event_ok_for_dataset_local,sigorbkg)]->add(args,weight);

    }

    if (dodistribution && event_ok_for_dataset>-1){

      obs_hist_single[get_name_obs_single(reg_lead,bin_lead)]->Fill(pholead_outvar,weight*ptweight_lead);
      obs_hist_single[get_name_obs_single(reg_trail,bin_trail)]->Fill(photrail_outvar,weight*ptweight_trail);

    
      for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){

	Int_t bin_couple = -999;
	
	float value_diffvariable;
	int event_ok_for_dataset_local = event_ok_for_dataset;

	if (*diffvariable==TString("invmass")) {
	  bin_couple = Choose_bin_invmass(dipho_mgg_photon,event_ok_for_dataset_local);
	  value_diffvariable=dipho_mgg_photon;
	  invmass_vector.push_back(value_diffvariable);
	}
	if (*diffvariable==TString("diphotonpt")){
	  float px = pholead_px+photrail_px;
	  float py = pholead_py+photrail_py;
	  float pt = sqrt(px*px+py*py);
	  bin_couple = Choose_bin_diphotonpt(pt,event_ok_for_dataset_local);
	  value_diffvariable=pt;
	  diphotonpt_vector.push_back(value_diffvariable);
	}
	if (*diffvariable==TString("costhetastar")){
	  TLorentzVector pho1(pholead_px,pholead_py,pholead_pz,pholead_energy);
	  TLorentzVector pho2(photrail_px,photrail_py,photrail_pz,photrail_energy);

// COS THETASTAR HX
//	  TVector3 boost = (pho1+pho2).BoostVector();
//	  TLorentzVector boostedpho1 = pho1;
//	  boostedpho1.Boost(-boost);
//	  float thetastar1 = boostedpho1.Angle(boost);
//	  bin_couple = Choose_bin_costhetastar(fabs(TMath::Cos(thetastar1)),event_ok_for_dataset_local);
//	  value_diffvariable=fabs(TMath::Cos(thetastar1));

// DELTA ETA
//	  value_diffvariable = fabs(TMath::TanH((pho1.Rapidity()-pho2.Rapidity())/2));
//	  bin_couple = Choose_bin_costhetastar(value_diffvariable,event_ok_for_dataset_local);

// COS THETASTAR CS
	  TLorentzVector b1,b2,diphoton;
	  b1.SetPx(0); b1.SetPy(0); b1.SetPz( 3500); b1.SetE(3500);
	  b2.SetPx(0); b2.SetPy(0); b2.SetPz(-3500); b2.SetE(3500);
	  TLorentzVector boostedpho1 = pho1; 
	  TLorentzVector boostedpho2 = pho2; 
	  TLorentzVector boostedb1 = b1; 
	  TLorentzVector boostedb2 = b2; 
	  TVector3 boost = (pho1+pho2).BoostVector();
	  boostedpho1.Boost(-boost);
	  boostedpho2.Boost(-boost);
	  boostedb1.Boost(-boost);
	  boostedb2.Boost(-boost);
	  TVector3 direction_cs = (boostedb1.Vect().Unit()-boostedb2.Vect().Unit()).Unit();
	  value_diffvariable = fabs(TMath::Cos(direction_cs.Angle(boostedpho1.Vect())));
	  bin_couple = Choose_bin_costhetastar(value_diffvariable,event_ok_for_dataset_local); 

	}
	if (*diffvariable==TString("dphi")){
	  float phi1 = pholead_SCphi;
	  float phi2 = photrail_SCphi;
	  float dphi = AbsDeltaPhi(phi1,phi2);
	  bin_couple = Choose_bin_dphi(dphi,event_ok_for_dataset_local);
	  value_diffvariable=dphi;
	}
      
	if (bin_couple<0) continue;
	
	float in1=pholead_outvar;
	float in2=photrail_outvar;
	float ptin1=pholead_pt;
	float ptin2=photrail_pt;
	float sieiein1=pholead_sieie;
	float sieiein2=photrail_sieie;
	float etain1=fabs(pholead_SCeta);
	float etain2=fabs(photrail_SCeta);

	bool doswap=false;


	if ((event_ok_for_dataset_local==0 || event_ok_for_dataset_local==2) && (randomgen->Uniform()>0.5)) doswap=true;
	
	if (event_ok_for_dataset_local==4) doswap=true;

	if (event_ok_for_dataset_local==3 || event_ok_for_dataset_local==4) event_ok_for_dataset_local=1;

	if (doswap){
	  float temp;
	  temp=in1; in1=in2; in2=temp;
	  temp=ptin1; ptin1=ptin2; ptin2=temp;
	  temp=sieiein1; sieiein1=sieiein2; sieiein2=temp;
	  temp=etain1; etain1=etain2; etain2=temp;
	}
      

	obs_hist_distribution[get_name_obs_distribution(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,weight);
	obs_hist[get_name_obs(event_ok_for_dataset_local,*diffvariable,bin_couple)]->Fill(in1,in2,weight);
	roovar1->setVal(in1);
	roovar2->setVal(in2);
	roopt1->setVal(ptin1);
	roopt2->setVal(ptin2);
	roosieie1->setVal(sieiein1);
	roosieie2->setVal(sieiein2);
	rooeta1->setVal(etain1);
	rooeta2->setVal(etain2);
	rooweight->setVal(weight);
	RooArgSet args(*roovar1,*roovar2,*roopt1,*roopt2,*roosieie1,*roosieie2,*rooeta1,*rooeta2);
	args.add(RooArgSet(*roorho,*roosigma));
	obs_roodset[get_name_obs_roodset(event_ok_for_dataset_local,*diffvariable,bin_couple)]->add(args,weight);

	if (!isdata){
	  int isppevent = 0;
	  if ( ((pholead_PhoMCmatchexitcode==1 || pholead_PhoMCmatchexitcode==2) && pholead_GenPhotonIsoDR04<5) && ((photrail_PhoMCmatchexitcode==1 || photrail_PhoMCmatchexitcode==2) && photrail_GenPhotonIsoDR04<5) ) isppevent=1;
	  true_purity[get_name_true_purity(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,isppevent,weight);
	}

//	if (!isdata) { SISTEMARE;
//	  if (leadistruesig && trailistruesig) weights_2p[event_ok_for_dataset_local][*diffvariable][bin_couple]+=weight;
//	  else if (!leadistruesig && !trailistruesig) weights_2f[event_ok_for_dataset_local][*diffvariable][bin_couple]+=weight;
//	  else weights_1p1f[event_ok_for_dataset_local][*diffvariable][bin_couple]+=weight;
//	}

      }
      
    }
    



  } // end event loop
  std::cout << "ended event loop" << std::endl;

//  // gen-level purities printout
//  for (int k=0; k<3; k++){
//    if (k==0) std::cout << "EBEB" << std::endl;
//    if (k==1) std::cout << "EBEE" << std::endl;
//    if (k==2) std::cout << "EEEE" << std::endl;
//    for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
//      std::cout << diffvariable->Data() << std::endl;
//      for (int i=0; i<n_bins; i++) {
//	if (weights_2p[k][*diffvariable][i]==0) continue;
//	std::cout << "bin " << i << " " << weights_2p[k][*diffvariable][i] << " " << weights_1p1f[k][*diffvariable][i] << " " << weights_2f[k][*diffvariable][i] << std::endl;
//	std::cout << "bin " << i << " " << weights_2p[k][*diffvariable][i]/(weights_2p[k][*diffvariable][i]+weights_1p1f[k][*diffvariable][i]+weights_2f[k][*diffvariable][i]) << std::endl;
//      }
//    }
//  }

  if (invmass_vector.size()>0){
    std::sort(invmass_vector.begin(),invmass_vector.end());
    std::sort(diphotonpt_vector.begin(),diphotonpt_vector.end());
    for (int i=1; i<10; i++) std::cout << "invmass" << invmass_vector.at(invmass_vector.size()-i) << std::endl;
    for (int i=1; i<10; i++) std::cout << "diphotonpt" << diphotonpt_vector.at(diphotonpt_vector.size()-i) << std::endl;
  }

};


#endif

void gen_templates(TString filename="input.root", TString mode="", bool isdata=1, const char* outfile="out.root", TString differentialvariable="photoniso", int maxevents=-1){
  
  TFile *outF = TFile::Open(outfile,"recreate");
  outF->Close();

  TFile *file;
  file = TFile::Open(filename.Data(),"read");

  TTree *t;

  TString treename[18];

  treename[0] = TString("Tree_2Dstandard_selection");
  treename[1] = TString("Tree_1Drandomcone_template");
  treename[2] = TString("Tree_1Dsideband_template");
  treename[3] = TString("Tree_2DZee_pixelvetoreversed_selection");
  treename[4] = TString("Tree_1Dpreselection");
  treename[5] = TString("Tree_1Dselection");
  treename[6] = TString("Tree_2Drandomcone_template");
  treename[7] = TString("Tree_2Drandomconesideband_template");
  treename[8] = TString("Tree_2Dsideband_template");
  treename[9] = TString("Tree_2Dstandard_preselection");
  treename[10] = TString("Tree_2DZmumu_selection");
  treename[11] = TString("Tree_1Dsignal_template");
  treename[12] = TString("Tree_1Dbackground_template");
  treename[13] = TString("Tree_2Dtruesigsig_template");
  treename[14] = TString("Tree_2Dtruesigbkg_template");
  treename[15] = TString("Tree_2Dtruebkgbkg_template");
  treename[16] = TString("Tree_2Drconeplusgenfake_template");
  treename[17] = TString("Tree_2Dgenpromptplussideband_template");


  TString treename_chosen="";
  if (mode=="standard") treename_chosen=treename[0];
  if (mode=="standard_2frag") treename_chosen=treename[0];
  if (mode=="signal") treename_chosen=treename[11];
  if (mode=="signal_2frag") treename_chosen=treename[11];
  if (mode=="fragmentation") treename_chosen=treename[11];
  if (mode=="nofragmentation") treename_chosen=treename[11];
  if (mode=="cutPFchargediso_signal") treename_chosen=treename[11];
  if (mode=="background") treename_chosen=treename[12];
  if (mode=="cutPFchargediso_background") treename_chosen=treename[12];
  if (mode=="randomcone") treename_chosen=treename[1];
  if (mode=="cutPFchargediso_randomcone") treename_chosen=treename[1];
  if (mode=="sieiesideband") treename_chosen=treename[2];
  if (mode=="cutPFchargediso_sieiesideband") treename_chosen=treename[2];
  if (mode=="zmumu") treename_chosen=treename[10];
  if (mode=="sigsig") treename_chosen=treename[6];
  if (mode=="sigbkg") treename_chosen=treename[7];
  if (mode=="bkgbkg") treename_chosen=treename[8];
  if (mode=="zee") treename_chosen=treename[3];
  if (mode=="2pgen") treename_chosen=treename[13];
  if (mode=="2pgen_2frag") treename_chosen=treename[13];
  if (mode=="1p1fbothgen") treename_chosen=treename[14];
  if (mode=="1p1fbothgen_2frag") treename_chosen=treename[14];
  if (mode=="2fgen") treename_chosen=treename[15];
  if (mode=="1prcone1fgen") treename_chosen=treename[16];
  if (mode=="1pgen1fside") treename_chosen=treename[17];
  if (mode=="1pgen1fside_2frag") treename_chosen=treename[17];


  file->GetObject(treename_chosen.Data(),t);

  std::cout << "Processing selection " << treename_chosen.Data() << std::endl;
  
  template_production *temp = new template_production(t);
  temp->Setup(isdata,mode,differentialvariable);

  temp->SetNoPtReweighting();
  temp->SetNoEtaReweighting();
  temp->SetNoPtEtaReweighting();

  if (maxevents>0) temp->Loop(maxevents); else temp->Loop();
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
  treename = TString("Tree_1Drandomcone_template");


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

    //    event_weight = FindNewPUWeight(event_nPU);


    float pholead_phoiso=-999;
    { // recalc phoiso w/cleaning
      float et_recalc = 0;
      float e_recalc = 0;
      int number_recalc = 0;


      for (int i=0; i<pholead_Npfcandphotonincone; i++) {

	float et=pholead_photonpfcandets[i];
	float e=pholead_photonpfcandenergies[i];
//	float deta=pholead_photonpfcanddetas[i];
//	float dphi=pholead_photonpfcanddphis[i];
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
      pholead_outvar=pholead_phoiso;
      assert (pholead_outvar>-100);
    }
    else if (diffvar=="combiso"){
      pholead_outvar=pholead_phoiso+pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01+pholead_pho_Cone04NeutralHadronIso_mvVtx;
    }
    else if (diffvar=="chargediso"){
      pholead_outvar=pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
    }
    else if (diffvar=="neutraliso"){
      pholead_outvar=pholead_pho_Cone04NeutralHadronIso_mvVtx;
    }

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

  int bin = Choose_bin_eta(fabs(eta),reg);
  float eff_area;
  if (isdata) eff_area = (reg==0) ? eff_areas_EB_data[bin] : eff_areas_EE_data[bin];
  else eff_area = (reg==0) ? eff_areas_EB_mc[bin] : eff_areas_EE_mc[bin];

  return 0.4*0.4*3.14*event_rho*eff_area;

};
