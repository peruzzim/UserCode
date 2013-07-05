#define efficiency_raw_producer_cxx
#include "efficiency_raw_producer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void efficiency_raw_producer::Loop()
{

  TH1F *true_reco = (TH1F*)(histo_pass[get_name_histo_pass(0,"invmass")]->Clone("reco"));
  TH1F *true_reco_nofakes = (TH1F*)(histo_pass[get_name_histo_pass(0,"invmass")]->Clone("reco_nofakes"));
  TH1F *true_gen = (TH1F*)(histo_pass[get_name_histo_pass(0,"invmass")]->Clone("gen"));

  const bool apply_scale_factors = true;

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
      
      if (!tree_found_gen && !tree_found_reco) continue;
      
      Int_t event_ok_for_dataset=-1;

      // dataset 0:EBEB 3/4->1:EBEE 2:EEEE

      if (tree_found_gen){
	if (fabs(pholead_GEN_eta)<1.4442 && fabs(photrail_GEN_eta)<1.4442) {
	  event_ok_for_dataset=0;
	}
	else if (fabs(pholead_GEN_eta)>1.56 && fabs(photrail_GEN_eta)>1.56) {
	  event_ok_for_dataset=2;
	}
	else if (fabs(pholead_GEN_eta)<1.4442 && fabs(photrail_GEN_eta)>1.56) {
	  event_ok_for_dataset=1;
	}
	else if (fabs(pholead_GEN_eta)>1.56 && fabs(photrail_GEN_eta)<1.4442) {
	  event_ok_for_dataset=1;
	}
	else std::cout << "We have a problem here!!!" << std::endl;
      }
      else if (tree_found_reco){
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
      
      Float_t weight=event_luminormfactor*event_Kfactor*event_weight;

      if (event_ok_for_dataset<=-1) continue;
	
	for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){

	  Int_t bin_couple = -999;	 
 	  Int_t bin_coupleGEN = -999;

	  float value_diffvariable;
	  float value_diffvariableGEN;
	  int event_ok_for_dataset_local = event_ok_for_dataset;

	  if (tree_found_gen){

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

	  if (tree_found_gen && bin_coupleGEN>=0){

	    float addweight=1;

	  if (apply_scale_factors && tree_found_match) {
	    addweight *= h_zee->GetBinContent(h_zee->FindBin(pholead_GEN_pt<h_zee->GetXaxis()->GetXmax() ? pholead_GEN_pt : h_zee->GetXaxis()->GetXmax()-1e-5 ,fabs(pholead_GEN_eta)));
	    addweight *= h_zee->GetBinContent(h_zee->FindBin(photrail_GEN_pt<h_zee->GetXaxis()->GetXmax() ? photrail_GEN_pt : h_zee->GetXaxis()->GetXmax()-1e-5 ,fabs(photrail_GEN_eta)));
	    addweight *= h_zuug->GetBinContent(h_zuug->FindBin(fabs(pholead_GEN_eta)));
	    addweight *= h_zuug->GetBinContent(h_zuug->FindBin(fabs(photrail_GEN_eta)));
	  }

	  if (tree_found_match) histo_pass[get_name_histo_pass(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariableGEN,weight*addweight);
	  else histo_fail[get_name_histo_fail(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariableGEN,weight*addweight);

	  }


	  if (tree_found_reco){

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

	  if (bin_couple<0) {tree_found_reco=false; tree_found_match=false;}

	  }



	  if (tree_found_match) responsewithmatch[get_name_responsewithmatch(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,value_diffvariableGEN,weight);

	  if (tree_found_reco){
	    if (tree_found_match) responsewithfakes[get_name_responsewithfakes(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,value_diffvariableGEN,weight);
	    else responsewithfakes[get_name_responsewithfakes(event_ok_for_dataset_local,*diffvariable)]->Fake(value_diffvariable,weight);
	  }

	  if (tree_found_reco || tree_found_gen){
	    if (tree_found_match) responsewitheff[get_name_responsewitheff(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,value_diffvariableGEN,weight);
	    else if (tree_found_gen) responsewitheff[get_name_responsewitheff(event_ok_for_dataset_local,*diffvariable)]->Miss(value_diffvariableGEN,weight);
	    else if (tree_found_reco) responsewitheff[get_name_responsewitheff(event_ok_for_dataset_local,*diffvariable)]->Fake(value_diffvariable,weight);
	  }

	  if (tree_found_reco && event_ok_for_dataset_local==0 && (*diffvariable==TString("invmass"))) true_reco->Fill(value_diffvariable,weight);
	  if (tree_found_match && event_ok_for_dataset_local==0 && (*diffvariable==TString("invmass"))) true_reco_nofakes->Fill(value_diffvariable,weight);
	  if (tree_found_gen && event_ok_for_dataset_local==0 &&(*diffvariable==TString("invmass"))) true_gen->Fill(value_diffvariableGEN,weight);
	  
	}


   } // end event loop


   TFile *fu = new TFile("unfolding.root","recreate");
   fu->cd();
   for (int i=0; i<3; i++)
     for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
       responsewithmatch[get_name_responsewithmatch(i,*diffvariable)]->Write(get_name_responsewithmatch(i,*diffvariable).Data());
       responsewithfakes[get_name_responsewithfakes(i,*diffvariable)]->Write(get_name_responsewithfakes(i,*diffvariable).Data());
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





   TH1F *true_reco_effnounf = (TH1F*)(true_reco_nofakes->Clone("true_reco_effnounf"));
   TH1F *true_reco_effunf = (TH1F*)(true_reco_nofakes->Clone("true_reco_effunf"));

   true_gen->SetLineColor(kRed);
   true_gen->SetMarkerColor(kRed);
   true_gen->SetMarkerStyle(20);
   true_gen->Draw("P");

   true_reco->SetLineColor(kBlack);
   true_reco->SetMarkerColor(kBlack);
   true_reco->SetMarkerStyle(20);
   true_reco->Draw("sameP");

   true_reco_effnounf->SetLineColor(kMagenta);
   true_reco_effnounf->Divide(histo_eff[get_name_histo_eff(0,"invmass")]);
   true_reco_effnounf->Draw("same");


   RooUnfoldBayes unf1(responsewithfakes[get_name_responsewithfakes(0,"invmass")],true_reco,4);
   TH1D *u1 = (TH1D*)(unf1.Hreco());
   u1->Divide(histo_eff[get_name_histo_eff(0,"invmass")]);
   u1->SetLineColor(kBlue);
   u1->Draw("same");
  
   RooUnfoldBayes unf2(responsewitheff[get_name_responsewitheff(0,"invmass")],true_reco,4);
   TH1D *u2 = (TH1D*)(unf2.Hreco());
   u2->SetLineColor(kGreen);
   u2->Draw("same");

   true_reco_effunf->Divide(histo_eff[get_name_histo_eff(0,"invmass")]);
   RooUnfoldBayes unf3(responsewithmatch[get_name_responsewithmatch(0,"invmass")],true_reco_effunf,4);
   TH1D *u3 = (TH1D*)(unf3.Hreco());
   u3->SetLineColor(kYellow);
   u3->Draw("same");




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
