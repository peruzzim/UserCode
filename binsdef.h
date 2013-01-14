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

const int n_templatebins = 10; 
Double_t templatebinsboundaries[n_templatebins+1] = {-3,-0.75,-0.5,-0.25,0,0.5,1,2,4,6,9};


static const int n_bins=13;

int n_templates_EB=7;
int n_templates_EE=5;
float binsdef_single_gamma_EB[n_bins+1]={30,40,50,60,70,80,90,110,140,150};
float binsdef_single_gamma_EE[n_bins+1]={30,40,50,60,70,80,90,110,120};
float binsdef_single_gamma_EB_eta[n_bins+1]={0,0.2,0.4,0.6,0.8,1,1.2,1.4442};
float binsdef_single_gamma_EE_eta[n_bins+1]={1.56,1.653,1.8,2,2.2,2.5};

int n_templates_invmass_EBEB=13;
int n_templates_invmass_EBEE=6;
int n_templates_invmass_EEEE=6;
//float binsdef_diphoton_invmass_EBEB[n_bins+1]={0,80,90,100,120,160,200};
float binsdef_diphoton_invmass_EBEB[n_bins+1]={0,50,60,70,80,90,100,110,120,130,140,150,170,200};
float binsdef_diphoton_invmass_EBEE[n_bins+1]={0,80,90,100,120,160,200};
float binsdef_diphoton_invmass_EEEE[n_bins+1]={0,80,90,100,120,160,200};
//float binsdef_diphoton_invmass_EEEE[n_bins+1]={0,20,40,60,80,90,100,110,120,130,140,160,190,200};
//float binsdef_diphoton_invmass_EBEB[n_bins+1]={80,90,100,110,120,130,140,160,190,200};
//float binsdef_diphoton_invmass_EBEE[n_bins+1]={80,90,100,110,120,130,140,160,190,200};
//float binsdef_diphoton_invmass_EEEE[n_bins+1]={80,90,100,110,120,130,140,160,190,200};

int n_templates_diphotonpt_EBEB=5;
int n_templates_diphotonpt_EBEE=5;
int n_templates_diphotonpt_EEEE=5;
float binsdef_diphoton_diphotonpt_EBEB[n_bins+1]={0,15,30,60,120,200};
float binsdef_diphoton_diphotonpt_EBEE[n_bins+1]={0,15,30,60,120,200};
float binsdef_diphoton_diphotonpt_EEEE[n_bins+1]={0,15,30,60,120,200};
//float binsdef_diphoton_diphotonpt_EEEE[n_bins+1]={0,20,40,60,80,90,100,110,120,130,140,160,190,200};

int n_templates_costhetastar_EBEB=7;
int n_templates_costhetastar_EBEE=7;
int n_templates_costhetastar_EEEE=7;
float binsdef_diphoton_costhetastar_EBEB[n_bins+1]={0,0.1,0.2,0.3,0.4,0.6,0.7,1};
float binsdef_diphoton_costhetastar_EBEE[n_bins+1]={0,0.1,0.2,0.3,0.4,0.6,0.7,1};
float binsdef_diphoton_costhetastar_EEEE[n_bins+1]={0,0.1,0.2,0.3,0.4,0.6,0.7,1};
//float binsdef_diphoton_costhetastar_EEEE[n_bins+1]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};

int n_templates_dphi_EBEB=5;
int n_templates_dphi_EBEE=5;
int n_templates_dphi_EEEE=5;
float binsdef_diphoton_dphi_EBEB[n_bins+1]={0,1.5,2.1,2.5,3.0,3.14};
float binsdef_diphoton_dphi_EBEE[n_bins+1]={0,1.5,2.1,2.5,3.0,3.14};
float binsdef_diphoton_dphi_EEEE[n_bins+1]={0,1.5,2.1,2.5,3.0,3.14};


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
// 020616 from data, no cleaning
float eff_areas_EB[n_bins] = {2.601118e-01,2.584915e-01,2.640072e-01,2.656851e-01,2.564615e-01,2.396511e-01,1.645776e-01};
float eff_areas_EE[n_bins] = {5.783452e-02,8.321881e-02,1.177009e-01,1.422445e-01,1.139434e-01};

// 020615rho with cleaning numbers
//float eff_areas_EB[n_bins] = {2.004668e-01,2.000222e-01,2.083325e-01,2.129163e-01,2.082317e-01,1.982015e-01,1.383834e-01};
//float eff_areas_EE[n_bins] = {3.727486e-02,5.494237e-02,7.876623e-02,1.006998e-01,8.432818e-02};

const int n_ptbins_forreweighting = 1;
//Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={0,40,60,300};
//Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={0,20,40,60,80,100,120,140,160,180,200,300};
Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={0,300};
//Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={0,30,300};

#endif
