#define efficiency_raw_producer_cxx
#include "efficiency_raw_producer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void efficiency_raw_producer::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if (!tree_found_gen) continue;

      
      Int_t event_ok_for_dataset=-1;

      // dataset 0:EBEB 3/4->1:EBEE 2:EEEE

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

      Float_t weight=event_luminormfactor*event_Kfactor*event_weight;

      if (event_ok_for_dataset>-1){
	
	for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
	  
	  Int_t bin_couple = -999;
	    
	  float value_diffvariable;
	  int event_ok_for_dataset_local = event_ok_for_dataset;

	  FillDiffVariablesGEN(); // WARNING: WORKS ONLY IF DIFF VARIABLES ARE NOT SENSITIVE TO SWAPPING 1 WITH 2

	  if (*diffvariable==TString("invmass")) {
	    value_diffvariable=localvar_invmass;
	    bin_couple = Choose_bin_invmass(value_diffvariable,event_ok_for_dataset_local);
	    //  invmass_vector.push_back(value_diffvariable);
	  }
	  if (*diffvariable==TString("diphotonpt")){
	    value_diffvariable=localvar_diphotonpt;
	    bin_couple = Choose_bin_diphotonpt(value_diffvariable,event_ok_for_dataset_local);
	    //  diphotonpt_vector.push_back(value_diffvariable);
	  }
	  if (*diffvariable==TString("costhetastar")){
	    value_diffvariable = localvar_costhetastar;
	    bin_couple = Choose_bin_costhetastar(value_diffvariable,event_ok_for_dataset_local); 
	  }
	  if (*diffvariable==TString("dphi")){
	    value_diffvariable=localvar_dphi;
	    bin_couple = Choose_bin_dphi(value_diffvariable,event_ok_for_dataset_local);
	  }
	  if (*diffvariable==TString("dR")){
	    value_diffvariable=localvar_dR;
	    bin_couple = Choose_bin_dR(value_diffvariable,event_ok_for_dataset_local);
	  }
      
	  if (bin_couple<0) continue;
	    
	  if (tree_found_match) histo_pass[get_name_histo_pass(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,weight);
	  else histo_fail[get_name_histo_fail(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,weight);

	}
      
      }


   } // end event loop


   for (int i=0; i<3; i++)
     for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
       histo_raweff[get_name_histo_raweff(i,*diffvariable)] = (TH1F*)(histo_pass[get_name_histo_pass(i,*diffvariable)]->Clone(get_name_histo_raweff(i,*diffvariable).Data()));
       histo_raweff[get_name_histo_raweff(i,*diffvariable)]->SetTitle(get_name_histo_raweff(i,*diffvariable).Data());
       TH1F *histo_tot = (TH1F*)(histo_pass[get_name_histo_pass(i,*diffvariable)]->Clone("tot"));
       histo_tot->Add(histo_fail[get_name_histo_fail(i,*diffvariable)]);
       histo_raweff[get_name_histo_raweff(i,*diffvariable)]->Divide(histo_tot);
       delete histo_tot;
     }

   TFile *f = new TFile ("raw_efficiency.root","recreate");
   f->cd();

   for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++)
     for (int i=0; i<3; i++){
       histo_raweff[get_name_histo_raweff(i,*diffvariable)]->SetAxisRange(0,1,"Y");
       histo_raweff[get_name_histo_raweff(i,*diffvariable)]->Write();
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
