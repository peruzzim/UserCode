//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul  5 00:44:42 2013 by ROOT version 5.34/08
// from TTree LightTreeGenReco/LightTreeGenReco
// found on file: DiPhotonJets_M0_7TeV_madgraph_Fall11_PU_S6_START42_V14B_v1_AODSIM.root
//////////////////////////////////////////////////////////

#ifndef efficiency_raw_producer_h
#define efficiency_raw_producer_h

#include "RooUnfold-1.1.1/src/RooUnfold.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"

#include "binsdef.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TRandom3.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TF1.h"
#include "TString.h"
#include "TH1.h"
#include "TProfile.h"
#include "TMath.h"
#include "TRandom3.h"
#include <map>
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>
#include <algorithm> 

using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class efficiency_raw_producer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   // Declaration of leaf types
   Int_t event_run;
   Int_t event_lumi;
   Int_t event_number;
   Float_t         event_luminormfactor;
   Float_t         event_Kfactor;
   Float_t         event_weight;
   Float_t         pholead_pt;
   Float_t         photrail_pt;
   Float_t         pholead_eta;
   Float_t         photrail_eta;
   Float_t         pholead_phi;
   Float_t         photrail_phi;
   Float_t         pholead_SCeta;
   Float_t         photrail_SCeta;
   Float_t         pholead_SCphi;
   Float_t         photrail_SCphi;
   Float_t         pholead_r9;
   Float_t         photrail_r9;
   Float_t         pholead_GEN_pt;
   Float_t         photrail_GEN_pt;
   Float_t         pholead_GEN_eta;
   Float_t         photrail_GEN_eta;
   Float_t         pholead_GEN_phi;
   Float_t         photrail_GEN_phi;
   Bool_t          reco_has_matched_gen_no_acceptance;
   Bool_t          reco_has_matched_gen_within_acceptance;
   Bool_t          reco_has_matched_gen_outside_acceptance;
   Bool_t          gen_in_acc;
   Bool_t          gen_in_acc_has_matched_reco;
   Bool_t          gen_in_acc_has_no_matched_reco;

   // List of branches
   TBranch *b_event_run;
   TBranch *b_event_lumi;
   TBranch *b_event_number;
   TBranch        *b_event_luminormfactor;   //!
   TBranch        *b_event_Kfactor;   //!
   TBranch        *b_event_weight;   //!
   TBranch        *b_pholead_pt;   //!
   TBranch        *b_photrail_pt;   //!
   TBranch        *b_pholead_eta;   //!
   TBranch        *b_photrail_eta;   //!
   TBranch        *b_pholead_phi;   //!
   TBranch        *b_photrail_phi;   //!
   TBranch        *b_pholead_SCeta;   //!
   TBranch        *b_photrail_SCeta;   //!
   TBranch        *b_pholead_SCphi;   //!
   TBranch        *b_photrail_SCphi;   //!
   TBranch        *b_pholead_r9;   //!
   TBranch        *b_photrail_r9;   //!
   TBranch        *b_pholead_GEN_pt;   //!
   TBranch        *b_photrail_GEN_pt;   //!
   TBranch        *b_pholead_GEN_eta;   //!
   TBranch        *b_photrail_GEN_eta;   //!
   TBranch        *b_pholead_GEN_phi;   //!
   TBranch        *b_photrail_GEN_phi;   //!
   TBranch        *b_reco_has_matched_gen_no_acceptance;   //!
   TBranch        *b_reco_has_matched_gen_within_acceptance;   //!
   TBranch        *b_reco_has_matched_gen_outside_acceptance;   //!
   TBranch        *b_gen_in_acc;   //!
   TBranch        *b_gen_in_acc_has_matched_reco;   //!
   TBranch        *b_gen_in_acc_has_no_matched_reco;   //!

   efficiency_raw_producer(TTree *tree=0);
   virtual ~efficiency_raw_producer();
   //   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void FillDiffVariables();
   void FillDiffVariablesGEN();

   TString get_name_histo_den(int region, TString diffvariable);
   TString get_name_histo_num(int region, TString diffvariable);
   TString get_name_histo_eff(int region, TString diffvariable);
   TString get_name_response(int region, TString diffvariable);

   std::map<TString, TH1F*> histo_eff;
   std::map<TString, TH1F*> histo_den;
   std::map<TString, TH1F*> histo_num;

   std::map<TString, RooUnfoldResponse*> response;

   Float_t   localvar_invmass;
   Float_t   localvar_diphotonpt;
   Float_t   localvar_costhetastar;
   Float_t   localvar_dphi;
   Float_t   localvar_dR;

   Float_t   localvarGEN_invmass;
   Float_t   localvarGEN_diphotonpt;
   Float_t   localvarGEN_costhetastar;
   Float_t   localvarGEN_dphi;
   Float_t   localvarGEN_dR;

   Int_t Choose_bin_invmass(float invmass, int region, bool *is_overflow=NULL);
   Int_t Choose_bin_diphotonpt(float diphotonpt, int region, bool *is_overflow=NULL);
   Int_t Choose_bin_costhetastar(float costhetastar, int region, bool *is_overflow=NULL);
   Int_t Choose_bin_dphi(float dphi, int region, bool *is_overflow=NULL);
   Int_t Choose_bin_dR(float dR, int region, bool *is_overflow=NULL);

   float Smearing(float eta, float r9);

   void DivideByBinSize(TH1* h);

};

#endif

#ifdef efficiency_raw_producer_cxx
efficiency_raw_producer::efficiency_raw_producer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  assert(tree!=NULL);
  Init(tree);


   for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){

     for (int i=0; i<3; i++) {
       TString reg;
       if (i==0) reg="EBEB"; else if (i==1) reg="EBEE"; else if (i==2) reg="EEEE"; else if (i==3) reg="EEEB";

       int bins_to_run=-1; 
       float *binsdef=NULL;

       TString splitting = reg;
       if (splitting=="EEEB") splitting="EBEE";
       if (*diffvariable=="invmass"){
	 if (splitting=="EBEB")      bins_to_run+=n_templates_invmass_EBEB;
	 else if (splitting=="EBEE") bins_to_run+=n_templates_invmass_EBEE;
	 else if (splitting=="EEEE") bins_to_run+=n_templates_invmass_EEEE; 
	 if (splitting=="EBEB")      binsdef=binsdef_diphoton_invmass_EBEB;
	 else if (splitting=="EBEE") binsdef=binsdef_diphoton_invmass_EBEE;
	 else if (splitting=="EEEE") binsdef=binsdef_diphoton_invmass_EEEE;
       }
       if (*diffvariable=="diphotonpt"){
	 if (splitting=="EBEB")      bins_to_run+=n_templates_diphotonpt_EBEB;
	 else if (splitting=="EBEE") bins_to_run+=n_templates_diphotonpt_EBEE;
	 else if (splitting=="EEEE") bins_to_run+=n_templates_diphotonpt_EEEE; 
	 if (splitting=="EBEB")      binsdef=binsdef_diphoton_diphotonpt_EBEB;
	 else if (splitting=="EBEE") binsdef=binsdef_diphoton_diphotonpt_EBEE;
	 else if (splitting=="EEEE") binsdef=binsdef_diphoton_diphotonpt_EEEE;
       }
       if (*diffvariable=="costhetastar"){
	 if (splitting=="EBEB")      bins_to_run+=n_templates_costhetastar_EBEB;
	 else if (splitting=="EBEE") bins_to_run+=n_templates_costhetastar_EBEE;
	 else if (splitting=="EEEE") bins_to_run+=n_templates_costhetastar_EEEE; 
	 if (splitting=="EBEB")      binsdef=binsdef_diphoton_costhetastar_EBEB;
	 else if (splitting=="EBEE") binsdef=binsdef_diphoton_costhetastar_EBEE;
	 else if (splitting=="EEEE") binsdef=binsdef_diphoton_costhetastar_EEEE;
       }
       if (*diffvariable=="dphi"){
	 if (splitting=="EBEB")      bins_to_run+=n_templates_dphi_EBEB;
	 else if (splitting=="EBEE") bins_to_run+=n_templates_dphi_EBEE;
	 else if (splitting=="EEEE") bins_to_run+=n_templates_dphi_EEEE; 
	 if (splitting=="EBEB")      binsdef=binsdef_diphoton_dphi_EBEB;
	 else if (splitting=="EBEE") binsdef=binsdef_diphoton_dphi_EBEE;
	 else if (splitting=="EEEE") binsdef=binsdef_diphoton_dphi_EEEE;
       }
       if (*diffvariable=="dR"){
	 if (splitting=="EBEB")      bins_to_run+=n_templates_dR_EBEB;
	 else if (splitting=="EBEE") bins_to_run+=n_templates_dR_EBEE;
	 else if (splitting=="EEEE") bins_to_run+=n_templates_dR_EEEE; 
	 if (splitting=="EBEB")      binsdef=binsdef_diphoton_dR_EBEB;
	 else if (splitting=="EBEE") binsdef=binsdef_diphoton_dR_EBEE;
	 else if (splitting=="EEEE") binsdef=binsdef_diphoton_dR_EEEE;
       }



//       TString t3=get_name_histo_eff(i,diffvariable->Data());
//       histo_eff[t3] = new TH1F(t3.Data(),t3.Data(),bins_to_run,binsdef);
       TString t3_pass=get_name_histo_den(i,diffvariable->Data());
       histo_den[t3_pass] = new TH1F(t3_pass.Data(),t3_pass.Data(),bins_to_run,binsdef);
       histo_den[t3_pass]->Sumw2();
       TString t3_fail=get_name_histo_num(i,diffvariable->Data());
       histo_num[t3_fail] = new TH1F(t3_fail.Data(),t3_fail.Data(),bins_to_run,binsdef);
       histo_num[t3_fail]->Sumw2();

       TH1F *proto_h_gen = (TH1F*)(histo_den[t3_pass]->Clone("GEN"));
       proto_h_gen->SetTitle("GEN");
       TH1F *proto_h_reco = (TH1F*)(histo_den[t3_pass]->Clone("RECO"));
       proto_h_reco->SetTitle("RECO");

       response[get_name_response(i,*diffvariable)] = new RooUnfoldResponse(proto_h_gen,proto_h_reco); // histo_den is used as template for binning here, contents are not used

     }


   }

}

efficiency_raw_producer::~efficiency_raw_producer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t efficiency_raw_producer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t efficiency_raw_producer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void efficiency_raw_producer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_run",&event_run,&b_event_run);
   fChain->SetBranchAddress("event_lumi",&event_lumi,&b_event_lumi);
   fChain->SetBranchAddress("event_number",&event_number,&b_event_number);
   fChain->SetBranchAddress("event_luminormfactor", &event_luminormfactor, &b_event_luminormfactor);
   fChain->SetBranchAddress("event_Kfactor", &event_Kfactor, &b_event_Kfactor);
   fChain->SetBranchAddress("event_weight", &event_weight, &b_event_weight);
   fChain->SetBranchAddress("pholead_pt", &pholead_pt, &b_pholead_pt);
   fChain->SetBranchAddress("photrail_pt", &photrail_pt, &b_photrail_pt);
   fChain->SetBranchAddress("pholead_eta", &pholead_eta, &b_pholead_eta);
   fChain->SetBranchAddress("photrail_eta", &photrail_eta, &b_photrail_eta);
   fChain->SetBranchAddress("pholead_phi", &pholead_phi, &b_pholead_phi);
   fChain->SetBranchAddress("photrail_phi", &photrail_phi, &b_photrail_phi);
   fChain->SetBranchAddress("pholead_SCeta", &pholead_SCeta, &b_pholead_SCeta);
   fChain->SetBranchAddress("photrail_SCeta", &photrail_SCeta, &b_photrail_SCeta);
   fChain->SetBranchAddress("pholead_SCphi", &pholead_SCphi, &b_pholead_SCphi);
   fChain->SetBranchAddress("photrail_SCphi", &photrail_SCphi, &b_photrail_SCphi);
   fChain->SetBranchAddress("pholead_r9", &pholead_r9, &b_pholead_r9);
   fChain->SetBranchAddress("photrail_r9", &photrail_r9, &b_photrail_r9);
   fChain->SetBranchAddress("pholead_GEN_pt", &pholead_GEN_pt, &b_pholead_GEN_pt);
   fChain->SetBranchAddress("photrail_GEN_pt", &photrail_GEN_pt, &b_photrail_GEN_pt);
   fChain->SetBranchAddress("pholead_GEN_eta", &pholead_GEN_eta, &b_pholead_GEN_eta);
   fChain->SetBranchAddress("photrail_GEN_eta", &photrail_GEN_eta, &b_photrail_GEN_eta);
   fChain->SetBranchAddress("pholead_GEN_phi", &pholead_GEN_phi, &b_pholead_GEN_phi);
   fChain->SetBranchAddress("photrail_GEN_phi", &photrail_GEN_phi, &b_photrail_GEN_phi);
   fChain->SetBranchAddress("reco_has_matched_gen_no_acceptance", &reco_has_matched_gen_no_acceptance, &b_reco_has_matched_gen_no_acceptance);
   fChain->SetBranchAddress("reco_has_matched_gen_within_acceptance", &reco_has_matched_gen_within_acceptance, &b_reco_has_matched_gen_within_acceptance);
   fChain->SetBranchAddress("reco_has_matched_gen_outside_acceptance", &reco_has_matched_gen_outside_acceptance, &b_reco_has_matched_gen_outside_acceptance);
   fChain->SetBranchAddress("gen_in_acc", &gen_in_acc, &b_gen_in_acc);
   fChain->SetBranchAddress("gen_in_acc_has_matched_reco", &gen_in_acc_has_matched_reco, &b_gen_in_acc_has_matched_reco);
   fChain->SetBranchAddress("gen_in_acc_has_no_matched_reco", &gen_in_acc_has_no_matched_reco, &b_gen_in_acc_has_no_matched_reco);
   Notify();
}

Bool_t efficiency_raw_producer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void efficiency_raw_producer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
//Int_t efficiency_raw_producer::Cut(Long64_t entry)
//{
//// This function may be called from Loop.
//// returns  1 if entry is accepted.
//// returns -1 otherwise.
//   return 1;
//}


void efficiency_raw_producer::FillDiffVariables(){

  TLorentzVector pho1; pho1.SetPtEtaPhiM(pholead_pt,pholead_eta,pholead_phi,0); 
  TLorentzVector pho2; pho2.SetPtEtaPhiM(photrail_pt,photrail_eta,photrail_phi,0);

  localvar_invmass = (pho1+pho2).M();
  localvar_diphotonpt = (pho1+pho2).Pt();

  {
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
    localvar_costhetastar=( fabs(TMath::Cos(direction_cs.Angle(boostedpho1.Vect()))) );
  }
  {
    float phi1 = pholead_phi;
    float phi2 = photrail_phi;
    float dphi = AbsDeltaPhi(phi1,phi2);
    localvar_dphi=(dphi);
  }
  {
    float phi1 = pholead_phi;
    float phi2 = photrail_phi;
    float dphi = AbsDeltaPhi(phi1,phi2);
    float deta = pholead_eta-photrail_eta;
    float dR = sqrt(deta*deta+dphi*dphi);
    localvar_dR=(dR);
  }



  return;

}

void efficiency_raw_producer::FillDiffVariablesGEN(){

  TLorentzVector pho1; pho1.SetPtEtaPhiM(pholead_GEN_pt,pholead_GEN_eta,pholead_GEN_phi,0);
  TLorentzVector pho2; pho2.SetPtEtaPhiM(photrail_GEN_pt,photrail_GEN_eta,photrail_GEN_phi,0);

  localvarGEN_invmass = (pho1+pho2).M();
  localvarGEN_diphotonpt = (pho1+pho2).Pt();
 
  {
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
    localvarGEN_costhetastar = fabs(TMath::Cos(direction_cs.Angle(boostedpho1.Vect()))) ;
  }
  {
    float phi1 = pholead_GEN_phi;
    float phi2 = photrail_GEN_phi;
    float dphi = AbsDeltaPhi(phi1,phi2);
    localvarGEN_dphi = dphi;
  }
  {
    float phi1 = pholead_GEN_phi;
    float phi2 = photrail_GEN_phi;
    float dphi = AbsDeltaPhi(phi1,phi2);
    float deta = pholead_GEN_eta-photrail_GEN_eta;
    float dR = sqrt(deta*deta+dphi*dphi);
    if (!(dR>0 && dR<6)) {
      std::cout << dR << " " << phi1 << " " << phi2 << " " << deta << std::endl;
      std::cout << gen_in_acc << " " << reco_has_matched_gen_no_acceptance << std::endl;
      std::cout << pholead_GEN_pt << " " << pholead_GEN_eta << " " << pholead_GEN_phi << std::endl;
      std::cout << pholead_pt << " " << pholead_eta << " " << pholead_phi  << std::endl;
      std::cout << photrail_GEN_pt << " " << photrail_GEN_eta << " " << photrail_GEN_phi << std::endl;
      std::cout << photrail_pt << " " << photrail_eta << " " << photrail_phi  << std::endl;
    }
    localvarGEN_dR = dR;
  }
  
  return;

}

TString efficiency_raw_producer::get_name_histo_den(int region, TString diffvariable){
  TString name_signal="histo_den";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s",name_signal.Data(),reg.Data(),diffvariable.Data());
  return t;
}

TString efficiency_raw_producer::get_name_histo_num(int region, TString diffvariable){
  TString name_signal="histo_num";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s",name_signal.Data(),reg.Data(),diffvariable.Data());
  return t;
}

TString efficiency_raw_producer::get_name_histo_eff(int region, TString diffvariable){
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  std::map<TString,TString> translation;
  translation.insert(std::pair<TString,TString>(TString("invmass"),TString("mgg")));
  translation.insert(std::pair<TString,TString>(TString("diphotonpt"),TString("qtgg")));
  translation.insert(std::pair<TString,TString>(TString("costhetastar"),TString("costhetastar")));
  translation.insert(std::pair<TString,TString>(TString("dphi"),TString("deltaphi")));
  translation.insert(std::pair<TString,TString>(TString("dR"),TString("dR")));
  TString t=Form("h_%s_%s_WithTotErr",translation[diffvariable].Data(),reg.Data());
  return t;
}

TString efficiency_raw_producer::get_name_response(int region, TString diffvariable){
  TString name_signal="response";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s",name_signal.Data(),reg.Data(),diffvariable.Data());
  return t;
}

Int_t efficiency_raw_producer::Choose_bin_invmass(float invmass, int region, bool *is_overflow){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_diphoton_invmass_EBEB; index=n_templates_invmass_EBEB;}
  if (region==2) {cuts=binsdef_diphoton_invmass_EEEE; index=n_templates_invmass_EEEE;}
  if (region==3) {cuts=binsdef_diphoton_invmass_EBEE; index=n_templates_invmass_EBEE;}
  if (region==4) {cuts=binsdef_diphoton_invmass_EBEE; index=n_templates_invmass_EBEE;}
  if (region==1) {cuts=binsdef_diphoton_invmass_EBEE; index=n_templates_invmass_EBEE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]= 9999;

  if (invmass<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value mass " << invmass << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((invmass>=cuts[i]) && (invmass<cuts[i+1])) {if (is_overflow) *is_overflow = (i==index-1); return i;}
  
   std::cout << "WARNING: called bin choice for out-of-range value mass " << invmass << std::endl;
  return -999;


};

Int_t efficiency_raw_producer::Choose_bin_diphotonpt(float diphotonpt, int region, bool *is_overflow){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_diphoton_diphotonpt_EBEB; index=n_templates_diphotonpt_EBEB;}
  if (region==2) {cuts=binsdef_diphoton_diphotonpt_EEEE; index=n_templates_diphotonpt_EEEE;}
  if (region==3) {cuts=binsdef_diphoton_diphotonpt_EBEE; index=n_templates_diphotonpt_EBEE;}
  if (region==4) {cuts=binsdef_diphoton_diphotonpt_EBEE; index=n_templates_diphotonpt_EBEE;}
  if (region==1) {cuts=binsdef_diphoton_diphotonpt_EBEE; index=n_templates_diphotonpt_EBEE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]= 9999;

  if (diphotonpt<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value diphotonpt " << diphotonpt << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((diphotonpt>=cuts[i]) && (diphotonpt<cuts[i+1])) {if (is_overflow) *is_overflow = (i==index-1); return i;}
  
   std::cout << "WARNING: called bin choice for out-of-range value diphotonpt " << diphotonpt << std::endl;
  return -999;


};

Int_t efficiency_raw_producer::Choose_bin_costhetastar(float costhetastar, int region, bool *is_overflow){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_diphoton_costhetastar_EBEB; index=n_templates_costhetastar_EBEB;}
  if (region==2) {cuts=binsdef_diphoton_costhetastar_EEEE; index=n_templates_costhetastar_EEEE;}
  if (region==3) {cuts=binsdef_diphoton_costhetastar_EBEE; index=n_templates_costhetastar_EBEE;}
  if (region==4) {cuts=binsdef_diphoton_costhetastar_EBEE; index=n_templates_costhetastar_EBEE;}
  if (region==1) {cuts=binsdef_diphoton_costhetastar_EBEE; index=n_templates_costhetastar_EBEE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]= 9999;

  if (costhetastar<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << costhetastar << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((costhetastar>=cuts[i]) && (costhetastar<cuts[i+1])) {if (is_overflow) *is_overflow = (i==index-1); return i;}
  
   std::cout << "WARNING: called bin choice for out-of-range value " << costhetastar << std::endl;
  return -999;


};

Int_t efficiency_raw_producer::Choose_bin_dphi(float dphi, int region, bool *is_overflow){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_diphoton_dphi_EBEB; index=n_templates_dphi_EBEB;}
  if (region==2) {cuts=binsdef_diphoton_dphi_EEEE; index=n_templates_dphi_EEEE;}
  if (region==3) {cuts=binsdef_diphoton_dphi_EBEE; index=n_templates_dphi_EBEE;}
  if (region==4) {cuts=binsdef_diphoton_dphi_EBEE; index=n_templates_dphi_EBEE;}
  if (region==1) {cuts=binsdef_diphoton_dphi_EBEE; index=n_templates_dphi_EBEE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]= 9999;

  if (dphi<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << dphi << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((dphi>=cuts[i]) && (dphi<cuts[i+1])) {if (is_overflow) *is_overflow = (i==index-1); return i;}
  
   std::cout << "WARNING: called bin choice for out-of-range value " << dphi << std::endl;
  return -999;


};

Int_t efficiency_raw_producer::Choose_bin_dR(float dR, int region, bool *is_overflow){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_diphoton_dR_EBEB; index=n_templates_dR_EBEB;}
  if (region==2) {cuts=binsdef_diphoton_dR_EEEE; index=n_templates_dR_EEEE;}
  if (region==3) {cuts=binsdef_diphoton_dR_EBEE; index=n_templates_dR_EBEE;}
  if (region==4) {cuts=binsdef_diphoton_dR_EBEE; index=n_templates_dR_EBEE;}
  if (region==1) {cuts=binsdef_diphoton_dR_EBEE; index=n_templates_dR_EBEE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]= 9999;

  if (dR<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << dR << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((dR>=cuts[i]) && (dR<cuts[i+1])) {if (is_overflow) *is_overflow = (i==index-1); return i;}
  
   std::cout << "WARNING: called bin choice for out-of-range value " << dR << std::endl;
  return -999;


};





#endif // #ifdef efficiency_raw_producer_cxx
