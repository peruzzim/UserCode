#ifndef __BINSDEF__
#define __BINSDEF__

#include <vector>
#include "TString.h"
#include "TMath.h"
#include <iostream>

TString __variables__[] = {TString("invmass"),TString("diphotonpt"),TString("costhetastar"),TString("dphi")};
std::vector<TString> diffvariables_list (__variables__, __variables__ + sizeof(__variables__) / sizeof(TString) );

const Int_t n_histobins = 20;
const Float_t leftrange = 0.0;
const Float_t rightrange = 10.0;





static const int n_bins=13;

int n_templates_EB=7;
int n_templates_EE=5;
float binsdef_single_gamma_EB[n_bins+1]={30,40,50,60,70,80,90,110,140,150};
float binsdef_single_gamma_EE[n_bins+1]={30,40,50,60,70,80,90,110,120};
float binsdef_single_gamma_EB_eta[n_bins+1]={0,0.2,0.4,0.6,0.8,1,1.2,1.442};
float binsdef_single_gamma_EE_eta[n_bins+1]={1.56,1.653,1.8,2,2.2,2.5};

int n_templates_invmass_EBEB=6;
int n_templates_invmass_EBEE=6;
int n_templates_invmass_EEEE=6;
float binsdef_diphoton_invmass_EBEB[n_bins+1]={0,80,90,100,120,160,200};
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






#endif
