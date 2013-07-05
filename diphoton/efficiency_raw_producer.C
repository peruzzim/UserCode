#define efficiency_raw_producer_cxx
#include "efficiency_raw_producer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void efficiency_raw_producer::Loop()
{


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
      
	  if (bin_coupleGEN<0) break; // this should never happen
	    
	  if (apply_scale_factors && tree_found_match) {
	    weight *= h_zee->GetBinContent(h_zee->FindBin(pholead_GEN_pt<h_zee->GetXaxis()->GetXmax() ? pholead_GEN_pt : h_zee->GetXaxis()->GetXmax()-1e-5 ,fabs(pholead_GEN_eta)));
	    weight *= h_zee->GetBinContent(h_zee->FindBin(photrail_GEN_pt<h_zee->GetXaxis()->GetXmax() ? photrail_GEN_pt : h_zee->GetXaxis()->GetXmax()-1e-5 ,fabs(photrail_GEN_eta)));
	    weight *= h_zuug->GetBinContent(h_zuug->FindBin(fabs(pholead_GEN_eta)));
	    weight *= h_zuug->GetBinContent(h_zuug->FindBin(fabs(photrail_GEN_eta)));
	  }

	  if (tree_found_match) histo_pass[get_name_histo_pass(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariableGEN,weight);
	  else histo_fail[get_name_histo_fail(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariableGEN,weight);

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

	  if (bin_couple<0) break;

	  if (tree_found_match) response[get_name_response(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,value_diffvariableGEN,weight);
	  else response[get_name_response(event_ok_for_dataset_local,*diffvariable)]->Fake(value_diffvariable,weight);

	  }

	}


   } // end event loop


   TFile *fu = new TFile("unfolding.root","recreate");
   fu->cd();
   for (int i=0; i<3; i++)
     for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
       response[get_name_response(i,*diffvariable)]->Write(get_name_response(i,*diffvariable).Data());
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

}
