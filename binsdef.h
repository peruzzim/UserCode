#ifndef __BINSDEF__
#define __BINSDEF__

#include <vector>
#include "TString.h"
#include "TMath.h"
#include <iostream>

TString __variables__[] = {TString("invmass"),TString("diphotonpt"),TString("costhetastar"),TString("dphi")};
std::vector<TString> diffvariables_list (__variables__, __variables__ + sizeof(__variables__) / sizeof(TString) );

const Int_t n_histobins = 96;
const Float_t leftrange = -3;
const Float_t rightrange = 9;

const int n_templatebins_full = 10; 
Double_t templatebinsboundaries_full[n_templatebins_full+1] = {-3,-0.75,-0.5,-0.25,0,0.5,1,2,4,6,9};
const int n_templatebins_reduced = 5; 
Double_t templatebinsboundaries_reduced[n_templatebins_reduced+1] = {-3,-0.5,0,1,4,9};
int n_templatebins = 0;
Double_t *templatebinsboundaries = NULL;

static const int n_bins=26;

int n_templates_EB=7;
int n_templates_EE=5;
float binsdef_single_gamma_EB[n_bins+1]={30,40,50,60,70,80,90,110,140,150};
float binsdef_single_gamma_EE[n_bins+1]={30,40,50,60,70,80,90,110,120};
float binsdef_single_gamma_EB_eta[n_bins+1]={0,0.2,0.4,0.6,0.8,1,1.2,1.4442};
float binsdef_single_gamma_EE_eta[n_bins+1]={1.56,1.653,1.8,2,2.2,2.5};

int n_templates_invmass_EBEB=13;
int n_templates_invmass_EBEE=13;
int n_templates_invmass_EEEE=13;
float binsdef_diphoton_invmass_EBEB[n_bins+1]={0,60,75,80,85,90,95,100,110,120,150,250,400,401};
float binsdef_diphoton_invmass_EBEE[n_bins+1]={0,60,75,80,85,90,95,100,110,120,150,250,400,401};
float binsdef_diphoton_invmass_EEEE[n_bins+1]={0,60,75,80,85,90,95,100,110,120,150,250,400,401};

int n_templates_diphotonpt_EBEB=17;
int n_templates_diphotonpt_EBEE=17;
int n_templates_diphotonpt_EEEE=17;
float binsdef_diphoton_diphotonpt_EBEB[n_bins+1]={0,6,10,12,14,16,18,20,22,24,28,34,40,50,60,80,120,121};
float binsdef_diphoton_diphotonpt_EBEE[n_bins+1]={0,6,10,12,14,16,18,20,22,24,28,34,40,50,60,80,120,121};
float binsdef_diphoton_diphotonpt_EEEE[n_bins+1]={0,6,10,12,14,16,18,20,22,24,28,34,40,50,60,80,120,121};

int n_templates_costhetastar_EBEB=13;
int n_templates_costhetastar_EBEE=13;
int n_templates_costhetastar_EEEE=13;
float binsdef_diphoton_costhetastar_EBEB[n_bins+1]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.8,1};
float binsdef_diphoton_costhetastar_EBEE[n_bins+1]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.8,1};
float binsdef_diphoton_costhetastar_EEEE[n_bins+1]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.8,1};


int n_templates_dphi_EBEB=13;
int n_templates_dphi_EBEE=13;
int n_templates_dphi_EEEE=13;
const float Pi = TMath::Pi();
float binsdef_diphoton_dphi_EBEB[n_bins+1]={0,0.2*Pi,0.4*Pi,0.6*Pi,0.7*Pi,0.8*Pi,0.84*Pi,0.88*Pi,0.90*Pi,0.92*Pi,0.94*Pi,0.96*Pi,0.98*Pi,1.0*Pi};
float binsdef_diphoton_dphi_EBEE[n_bins+1]={0,0.2*Pi,0.4*Pi,0.6*Pi,0.7*Pi,0.8*Pi,0.84*Pi,0.88*Pi,0.90*Pi,0.92*Pi,0.94*Pi,0.96*Pi,0.98*Pi,1.0*Pi};
float binsdef_diphoton_dphi_EEEE[n_bins+1]={0,0.2*Pi,0.4*Pi,0.6*Pi,0.7*Pi,0.8*Pi,0.84*Pi,0.88*Pi,0.90*Pi,0.92*Pi,0.94*Pi,0.96*Pi,0.98*Pi,1.0*Pi};


float AbsDeltaPhi(double phi1, double phi2){
  // From cmssw reco::deltaPhi()
  double result = phi1 - phi2;
  while( result >   TMath::Pi() ) result -= TMath::TwoPi();
  while( result <= -TMath::Pi() ) result += TMath::TwoPi();
  return TMath::Abs(result);
};

//const int n_rho_cats=12; // = (rhobins-1)
//const int n_sigma_cats=9; // = (sigmabins-1)
//int n_rhosigma_cats=n_rho_cats*n_sigma_cats; // = (rhobins-1)*(sigmabins-1)
//float rhobins[n_rho_cats+1]={0,2,4,6,8,10,12,14,16,18,20,25,100};
//float sigmabins[n_sigma_cats+1]={0,1,2,3,4,5,6,7,8,100};
const int n_rho_cats=1; // = (rhobins-1)
//const int n_sigma_cats=7; // = (sigmabins-1)
const int n_sigma_cats=1; // = (sigmabins-1)
int n_rhosigma_cats=n_rho_cats*n_sigma_cats; // = (rhobins-1)*(sigmabins-1)
float rhobins[n_rho_cats+1]={0,100};
//float sigmabins[n_sigma_cats+1]={0,1,2,3,4,5,6,100};
float sigmabins[n_sigma_cats+1]={0,100};

const int n_eta_cats = n_templates_EB;
int n_eta1eta2_cats = n_eta_cats*n_eta_cats;
float *etabins = binsdef_single_gamma_EB_eta+0;

// FOR PHOTON COMPONENT
// 020616 from data, no cleaning, no pf charged cut in presel
float eff_areas_EB_data[n_bins] = {2.615381e-01,2.587228e-01,2.641340e-01,2.640072e-01,2.569331e-01,2.395977e-01,1.651901e-01};
float eff_areas_EE_data[n_bins] = {5.848781e-02,8.411012e-02,1.189097e-01,1.446802e-01,1.148867e-01};
float eff_areas_EB_mc[n_bins] = {2.719313e-01,2.748353e-01,2.746818e-01,2.704172e-01,2.662482e-01,2.427560e-01,1.708246e-01};
float eff_areas_EE_mc[n_bins] = {5.107922e-02,8.485798e-02,1.354249e-01,1.671490e-01,1.454153e-01};



const int n_ptbins_forreweighting = 1;
Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={0,300};

//const int n_ptbins_forreweighting = 6;
//Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={0,35,45,55,65,85,999};


#endif
