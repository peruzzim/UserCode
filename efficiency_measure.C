#define efficiency_measure_cxx
#include "efficiency_measure.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void efficiency_measure::Loop()
{
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

   static const int n_bins=15;

   float binsdef_single_gamma[n_bins+1]={30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180};
   float binsdef_diphoton[n_bins+1]={80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230};

   TH1::SetDefaultSumw2(kTRUE);

   TString gg_name[3]={"EBEB","EBEE","EEEE"};
   TString sg_name[2]={"EB","EE"};

   for (int i=0; i<3; i++) w_eff_gg[i] = new TH1F(Form("w_eff_gg_%s",gg_name[i].Data()),Form("w_eff_gg_%s",gg_name[i].Data()),n_bins,binsdef_diphoton);
   for (int i=0; i<2; i++) w_eff_1g[i] = new TH1F(Form("w_eff_1g_%s",sg_name[i].Data()),Form("w_eff_1g_%s",sg_name[i].Data()),n_bins,binsdef_single_gamma);
   for (int i=0; i<3; i++) w_tot_gg[i] = new TH1F(Form("w_tot_gg_%s",gg_name[i].Data()),Form("w_tot_gg_%s",gg_name[i].Data()),n_bins,binsdef_diphoton);
   for (int i=0; i<2; i++) w_tot_1g[i] = new TH1F(Form("w_tot_1g_%s",sg_name[i].Data()),Form("w_tot_1g_%s",sg_name[i].Data()),n_bins,binsdef_single_gamma);
   for (int i=0; i<3; i++) w_passing_gg[i] = new TH1F(Form("w_passing_gg_%s",gg_name[i].Data()),Form("w_passing_gg_%s",gg_name[i].Data()),n_bins,binsdef_diphoton);
   for (int i=0; i<2; i++) w_passing_1g[i] = new TH1F(Form("w_passing_1g_%s",sg_name[i].Data()),Form("w_passing_1g_%s",sg_name[i].Data()),n_bins,binsdef_single_gamma);

   for (int i=0; i<3; i++) w_eff_gg[i]->Sumw2();
   for (int i=0; i<2; i++) w_eff_1g[i]->Sumw2();
   for (int i=0; i<3; i++) w_tot_gg[i]->Sumw2();
   for (int i=0; i<2; i++) w_tot_1g[i]->Sumw2();
   for (int i=0; i<3; i++) w_passing_gg[i]->Sumw2();
   for (int i=0; i<2; i++) w_passing_1g[i]->Sumw2();

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

      if (pholead_pt<40 || photrail_pt<30 || dipho_mgg_photon<80) continue;

      { // 2 prompt photon efficiency

      bool pass1=true;

      if (!(pholead_PhoMCmatchexitcode==1 || pholead_PhoMCmatchexitcode==2)) pass1=false;
      if (!(photrail_PhoMCmatchexitcode==1 || photrail_PhoMCmatchexitcode==2)) pass1=false;

      if (pass1) w_tot_gg[event_ok_for_dataset]->Fill(dipho_mgg_photon,weight);
      
      bool pass2;
      pass2=true;

      float cutUP, cutLOW;
      if (fabs(pholead_SCeta)<1.4442) {cutLOW=0; cutUP=0.011;} // EB
      else {cutLOW=0; cutUP=0.028;} // EE

      if (pholead_hoe>0.05) pass2=false;
      if (pholead_sieie>cutUP) pass2=false;
      if (pholead_sieie<cutLOW) pass2=false;

      if (fabs(photrail_SCeta)<1.4442) {cutLOW=0; cutUP=0.011;} // EB
      else {cutLOW=0; cutUP=0.028;} // EE

      if (photrail_hoe>0.05) pass2=false;
      if (photrail_sieie>cutUP) pass2=false;
      if (photrail_sieie<cutLOW) pass2=false;

      if (pass1 && pass2) w_passing_gg[event_ok_for_dataset]->Fill(dipho_mgg_photon,weight);

      } // end 2 prompt efficiency    

      { // single prompt photon efficiency

	bool pass1;
	bool pass2;
      float cutUP, cutLOW;

	pass1=true;
	pass2=true;
	
	if (!(pholead_PhoMCmatchexitcode==1 || pholead_PhoMCmatchexitcode==2)) pass1=false;
	if (pass1) w_tot_1g[reg_lead]->Fill(pholead_pt,weight);


      if (fabs(pholead_SCeta)<1.4442) {cutLOW=0; cutUP=0.011;} // EB
      else {cutLOW=0; cutUP=0.028;} // EE

      if (pholead_hoe>0.05) pass2=false;
      if (pholead_sieie>cutUP) pass2=false;
      if (pholead_sieie<cutLOW) pass2=false;

      if (pass1 && pass2) w_passing_1g[reg_lead]->Fill(pholead_pt,weight);

      pass1=true;
      pass2=true;

      if (!(photrail_PhoMCmatchexitcode==1 || photrail_PhoMCmatchexitcode==2)) pass1=false;
      if (pass1) w_tot_1g[reg_trail]->Fill(photrail_pt,weight);

      if (fabs(photrail_SCeta)<1.4442) {cutLOW=0; cutUP=0.011;} // EB
      else {cutLOW=0; cutUP=0.028;} // EE

      if (photrail_hoe>0.05) pass2=false;
      if (photrail_sieie>cutUP) pass2=false;
      if (photrail_sieie<cutLOW) pass2=false;

      if (pass1 && pass2) w_passing_1g[reg_trail]->Fill(photrail_pt,weight);

      } // end 2 prompt efficiency    

   } // end event loop

   for (int i=0; i<3; i++) w_eff_gg[i]->Divide(w_passing_gg[i],w_tot_gg[i],1,1,"B");
   for (int i=0; i<2; i++) w_eff_1g[i]->Divide(w_passing_1g[i],w_tot_1g[i],1,1,"B");

   if (true){
   TCanvas *c = new TCanvas();
   c->Divide(3,2);
   for (int i=0; i<3; i++) {c->cd(i+1); w_eff_gg[i]->Draw("e1");}
   for (int i=0; i<2; i++) {c->cd(i+4); w_eff_1g[i]->Draw("e1");}
   }

   if (true){
     TFile *outf = new TFile("efficiencies.root","recreate");
     outf->cd();
     for (int i=0; i<3; i++) w_eff_gg[i]->Write();
     for (int i=0; i<2; i++) w_eff_1g[i]->Write();
     outf->Close();
   }

}
