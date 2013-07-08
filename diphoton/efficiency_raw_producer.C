#define efficiency_raw_producer_cxx
#include "efficiency_raw_producer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void efficiency_raw_producer::Loop()
{

  const bool apply_scale_factors = true;
  const bool do_energy_smearing = true;
  const bool apply_trigger_efficiency = true;


  const float trig_eff_EBEB_highr9 = 1.;
  const float trig_eff_EBEB_lowr9 = 0.993;
  const float trig_eff_notEBEB_highr9 = 1.;
  const float trig_eff_notEBEB_lowr9 = 0.988;


  // monitoring
  TH1F *true_reco = (TH1F*)(histo_pass[get_name_histo_pass(0,"invmass")]->Clone("reco"));
  TH1F *true_gen = (TH1F*)(histo_pass[get_name_histo_pass(0,"invmass")]->Clone("gen"));


  TRandom3 *rand = new TRandom3(0);

  TFile *file_scalefactor_zee = new TFile("histo_scalefactor_Zee_totaluncertainty.root","read");
  TH2F *h_zee=NULL; file_scalefactor_zee->GetObject("histo_withsyst",h_zee);
  TFile *file_scalefactor_zuug = new TFile("histo_scalefactor_Zuug_totaluncertainty.root","read");
  TH1F *h_zuug=NULL; file_scalefactor_zuug->GetObject("h_zuug",h_zuug);
  assert (h_zee);
  assert (h_zuug);


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // BUGFIX for productions up to jul7 included 
      gen_in_acc_has_no_matched_reco = gen_in_acc && !gen_in_acc_has_matched_reco;

      Int_t event_ok_for_dataset=-1;

      // dataset 0:EBEB 3/4->1:EBEE 2:EEEE

      if (reco_has_matched_gen_no_acceptance){
	if (fabs(pholead_SCeta)<1.4442 && fabs(photrail_SCeta)<1.4442) {
	  event_ok_for_dataset=0;
	}
	else if (fabs(pholead_SCeta)>1.56 && fabs(photrail_SCeta)>1.56) {
	  event_ok_for_dataset=2;
	}
	else if (fabs(pholead_SCeta)<1.4442 && fabs(photrail_SCeta)>1.56) {
	  event_ok_for_dataset=1;
	}
	else if (fabs(pholead_SCeta)>1.56 && fabs(photrail_SCeta)<1.4442) {
	  event_ok_for_dataset=1;
	}
	else std::cout << "We have a problem here!!!" << std::endl;
      }
      else if (gen_in_acc){
	if (fabs(pholead_GEN_eta)<1.4442 && fabs(photrail_GEN_eta)<1.4442) {
	  event_ok_for_dataset=0;
	}
	else if (fabs(pholead_GEN_eta)>1.566 && fabs(photrail_GEN_eta)>1.566) {
	  event_ok_for_dataset=2;
	}
	else if (fabs(pholead_GEN_eta)<1.4442 && fabs(photrail_GEN_eta)>1.566) {
	  event_ok_for_dataset=1;
	}
	else if (fabs(pholead_GEN_eta)>1.566 && fabs(photrail_GEN_eta)<1.4442) {
	  event_ok_for_dataset=1;
	}
	else std::cout << "We have a problem here!!!" << std::endl;
      }
      
      Float_t weight=event_luminormfactor*event_Kfactor*event_weight;

      if (event_ok_for_dataset<=-1) continue;
	
      if (reco_has_matched_gen_no_acceptance && do_energy_smearing){
	pholead_pt*=rand->Gaus(1,Smearing(pholead_SCeta,pholead_r9));
	photrail_pt*=rand->Gaus(1,Smearing(photrail_SCeta,photrail_r9));
      }

	for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){

	  Int_t bin_couple = -999;	 
 	  Int_t bin_coupleGEN = -999;

	  float value_diffvariable;
	  float value_diffvariableGEN;
	  int event_ok_for_dataset_local = event_ok_for_dataset;

	  if (gen_in_acc || reco_has_matched_gen_no_acceptance){

	  FillDiffVariablesGEN(); // WARNING: WORKS ONLY IF DIFF VARIABLES ARE NOT SENSITIVE TO SWAPPING 1 WITH 2

	  if (*diffvariable==TString("invmass")) {
	    value_diffvariableGEN=localvarGEN_invmass;
	    bin_coupleGEN = Choose_bin_invmass(value_diffvariableGEN,event_ok_for_dataset_local);
	  }
	  if (*diffvariable==TString("diphotonpt")){
	    value_diffvariableGEN=localvarGEN_diphotonpt;
	    bin_coupleGEN = Choose_bin_diphotonpt(value_diffvariableGEN,event_ok_for_dataset_local);
	  }
	  if (*diffvariable==TString("costhetastar")){
	    value_diffvariableGEN = localvarGEN_costhetastar;
	    bin_coupleGEN = Choose_bin_costhetastar(value_diffvariableGEN,event_ok_for_dataset_local); 
	  }
	  if (*diffvariable==TString("dphi")){
	    value_diffvariableGEN=localvarGEN_dphi;
	    bin_coupleGEN = Choose_bin_dphi(value_diffvariableGEN,event_ok_for_dataset_local);
	  }
	  if (*diffvariable==TString("dR")){
	    value_diffvariableGEN=localvarGEN_dR;
	    bin_coupleGEN = Choose_bin_dR(value_diffvariableGEN,event_ok_for_dataset_local);
	  }
      
	  }

	  if (gen_in_acc && bin_coupleGEN>=0){

	    float addweight=1;

	  if (apply_scale_factors && gen_in_acc_has_matched_reco) {
	    addweight *= h_zee->GetBinContent(h_zee->FindBin(pholead_pt<h_zee->GetXaxis()->GetXmax() ? pholead_pt : h_zee->GetXaxis()->GetXmax()-1e-5 ,fabs(pholead_SCeta)));
	    addweight *= h_zee->GetBinContent(h_zee->FindBin(photrail_pt<h_zee->GetXaxis()->GetXmax() ? photrail_pt : h_zee->GetXaxis()->GetXmax()-1e-5 ,fabs(photrail_SCeta)));
	    addweight *= h_zuug->GetBinContent(h_zuug->FindBin(fabs(pholead_SCeta)));
	    addweight *= h_zuug->GetBinContent(h_zuug->FindBin(fabs(photrail_SCeta)));
	  }
	  if (apply_trigger_efficiency && gen_in_acc_has_matched_reco){
	    if (event_ok_for_dataset_local==0) addweight *= (pholead_r9>0.94 && photrail_r9>0.94) ? trig_eff_EBEB_highr9 : trig_eff_EBEB_lowr9;
	    else addweight *= (pholead_r9>0.94 && photrail_r9>0.94) ? trig_eff_notEBEB_highr9 : trig_eff_notEBEB_lowr9;
	  }

	  if (gen_in_acc_has_matched_reco) histo_pass[get_name_histo_pass(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariableGEN,weight*addweight);
	  else histo_fail[get_name_histo_fail(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariableGEN,weight*addweight);

	  }


	  if (reco_has_matched_gen_no_acceptance){

	  FillDiffVariables(); // WARNING: WORKS ONLY IF DIFF VARIABLES ARE NOT SENSITIVE TO SWAPPING 1 WITH 2

	  if (*diffvariable==TString("invmass")) {
	    value_diffvariable=localvar_invmass;
	    bin_couple = Choose_bin_invmass(value_diffvariable,event_ok_for_dataset_local,true);
	  }
	  if (*diffvariable==TString("diphotonpt")){
	    value_diffvariable=localvar_diphotonpt;
	    bin_couple = Choose_bin_diphotonpt(value_diffvariable,event_ok_for_dataset_local,true);
	  }
	  if (*diffvariable==TString("costhetastar")){
	    value_diffvariable = localvar_costhetastar;
	    bin_couple = Choose_bin_costhetastar(value_diffvariable,event_ok_for_dataset_local,true); 
	  }
	  if (*diffvariable==TString("dphi")){
	    value_diffvariable=localvar_dphi;
	    bin_couple = Choose_bin_dphi(value_diffvariable,event_ok_for_dataset_local,true);
	  }
	  if (*diffvariable==TString("dR")){
	    value_diffvariable=localvar_dR;
	    bin_couple = Choose_bin_dR(value_diffvariable,event_ok_for_dataset_local,true);
	  }

	  if (bin_couple<0) {reco_has_matched_gen_no_acceptance=false; reco_has_matched_gen_within_acceptance=false; reco_has_matched_gen_outside_acceptance=false;}

	  }

	  if (reco_has_matched_gen_no_acceptance){
	    if (reco_has_matched_gen_within_acceptance) response[get_name_response(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,value_diffvariableGEN,weight);
	    else if (reco_has_matched_gen_outside_acceptance) response[get_name_response(event_ok_for_dataset_local,*diffvariable)]->Fake(value_diffvariable,weight);
	  }

	  if (reco_has_matched_gen_no_acceptance || gen_in_acc){
	    if (reco_has_matched_gen_within_acceptance || gen_in_acc_has_matched_reco) responsewitheff[get_name_responsewitheff(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,value_diffvariableGEN,weight);
	    else if (reco_has_matched_gen_outside_acceptance) responsewitheff[get_name_responsewitheff(event_ok_for_dataset_local,*diffvariable)]->Fake(value_diffvariable,weight);
	    else if (gen_in_acc && !gen_in_acc_has_matched_reco) responsewitheff[get_name_responsewitheff(event_ok_for_dataset_local,*diffvariable)]->Miss(value_diffvariableGEN,weight);
	  }

	  // monitoring
	  if (reco_has_matched_gen_no_acceptance && event_ok_for_dataset_local==0 && (*diffvariable==TString("invmass"))) true_reco->Fill(value_diffvariable,weight);
	  if (gen_in_acc && event_ok_for_dataset_local==0 &&(*diffvariable==TString("invmass"))) true_gen->Fill(value_diffvariableGEN,weight);

	}


   } // end event loop


   TFile *fu = new TFile("unfolding.root","recreate");
   fu->cd();
   for (int i=0; i<3; i++)
     for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
       response[get_name_response(i,*diffvariable)]->Write(get_name_response(i,*diffvariable).Data());
       responsewitheff[get_name_responsewitheff(i,*diffvariable)]->Write(get_name_responsewitheff(i,*diffvariable).Data());
     }
   fu->Close();

   for (int i=0; i<3; i++)
     for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
       histo_eff[get_name_histo_eff(i,*diffvariable)] = (TH1F*)(histo_pass[get_name_histo_pass(i,*diffvariable)]->Clone(get_name_histo_eff(i,*diffvariable).Data()));
       histo_eff[get_name_histo_eff(i,*diffvariable)]->SetTitle(get_name_histo_eff(i,*diffvariable).Data());
       TH1F *histo_tot = (TH1F*)(histo_pass[get_name_histo_pass(i,*diffvariable)]->Clone("tot"));
       histo_tot->Add(histo_fail[get_name_histo_fail(i,*diffvariable)]);
       histo_eff[get_name_histo_eff(i,*diffvariable)]->Divide(histo_eff[get_name_histo_eff(i,*diffvariable)],histo_tot,1,1,"B");
       delete histo_tot;
     }

   TFile *f = new TFile ("efficiency.root","recreate");
   f->cd();

   for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++)
     for (int i=0; i<3; i++){
       histo_eff[get_name_histo_eff(i,*diffvariable)]->SetAxisRange(0,1,"Y");
       histo_eff[get_name_histo_eff(i,*diffvariable)]->Write();
     }

   f->mkdir("other_histos");
   f->cd("other_histos");
   for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++)
     for (int i=0; i<3; i++){
       histo_pass[get_name_histo_pass(i,*diffvariable)]->Write();
       histo_fail[get_name_histo_fail(i,*diffvariable)]->Write();
     }

   f->Close();




   // monitoring

   TH1F *true_reco_effonly = (TH1F*)(true_reco->Clone("true_reco_effonly"));
   TH1F *true_reco_effunf = (TH1F*)(true_reco->Clone("true_reco_effunf"));

   true_gen->SetLineColor(kRed);
   true_gen->SetMarkerColor(kRed);
   true_gen->SetMarkerStyle(20);

   true_reco->SetLineColor(kBlack);
   true_reco->SetMarkerColor(kBlack);
   true_reco->SetMarkerStyle(20);

   // NO UNFOLDING
   true_reco_effonly->SetLineColor(kBlack);
   true_reco_effonly->SetMarkerColor(kBlack);
   true_reco_effonly->SetMarkerStyle(21);
   true_reco_effonly->Divide(histo_eff[get_name_histo_eff(0,"invmass")]);

   // ONE SHOT, EFF+UNFOLDING INSIEME, NO SCALE FACTORS / HLT
   RooUnfoldBayes unf2(responsewitheff[get_name_responsewitheff(0,"invmass")],true_reco,4);
   TH1D *u2 = (TH1D*)(unf2.Hreco());
   u2->SetLineColor(kGreen);
   u2->SetMarkerColor(kGreen);
   u2->SetMarkerStyle(20);

   // EFFICIENZA E UNFOLDING IN SEQUENZA (ATTUALE)
   true_reco_effunf->Divide(histo_eff[get_name_histo_eff(0,"invmass")]);
   RooUnfoldBayes unf3(response[get_name_response(0,"invmass")],true_reco_effunf,4);
   TH1D *u3 = (TH1D*)(unf3.Hreco());
   u3->SetLineColor(kBlue);
   u3->SetMarkerColor(kBlue);
   u3->SetMarkerStyle(20);
   RooUnfoldBinByBin unf3b(response[get_name_response(0,"invmass")],true_reco_effunf);
   TH1D *u3b = (TH1D*)(unf3b.Hreco());
   u3b->SetLineColor(kYellow);
   u3b->SetMarkerColor(kYellow);
   u3b->SetMarkerStyle(20);

   DivideByBinSize(true_gen);
   DivideByBinSize(true_reco);
   DivideByBinSize(true_reco_effonly);
   DivideByBinSize(u2);
   DivideByBinSize(u3);
   DivideByBinSize(u3b);

   true_gen->Draw("P");
   true_reco->Draw("sameP");
   true_reco_effonly->Draw("sameP");
   u2->Draw("sameP");
   u3->Draw("sameP");
   u3b->Draw("sameP");



   



}


float efficiency_raw_producer::Smearing(float eta, float r9){

  // Smearings from Shervin for regression energy, 16Jan rereco

  eta = fabs(eta);

  float r=0;

  if (r9<0.94){
    if (eta<1) r= 0.96;
    else if (eta<1.4442) r= 1.96;
    else if (eta<2) r= 2.79;
    else r= 3.01;
  }
  else {
    if (eta<1) r= 0.74;
    else if (eta<1.4442) r= 1.43;
    else if (eta<2) r= 2.68;
    else r= 2.93;
  }

  return r/100.;

}

void efficiency_raw_producer::DivideByBinSize(TH1* h){

  for (int i=1; i<=h->GetNbinsX(); i++) h->SetBinContent(i,h->GetBinContent(i)/h->GetBinWidth(i));
  for (int i=1; i<=h->GetNbinsX(); i++) h->SetBinError(i,h->GetBinError(i)/h->GetBinWidth(i));

}
