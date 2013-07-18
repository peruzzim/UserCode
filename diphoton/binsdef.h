#ifndef __BINSDEF__
#define __BINSDEF__

#include <map>
#include <vector>
#include "TString.h"
#include "TMath.h"
#include <iostream>

TString __variables__[] = {TString("invmass"),TString("diphotonpt"),TString("costhetastar"),TString("dphi"),TString("dR")};
std::vector<TString> diffvariables_list (__variables__, __variables__ + sizeof(__variables__) / sizeof(TString) );

std::map<TString,TString> diffvariables_names_list_;
TString diffvariables_names_list(TString diffvariable){
  if (diffvariables_names_list_.size()==0) {
    diffvariables_names_list_[TString("invmass")]=TString("m_{#gamma #gamma}");
    diffvariables_names_list_[TString("diphotonpt")]=TString("p_{T, #gamma #gamma}");
    diffvariables_names_list_[TString("costhetastar")]=TString("|cos #theta^{*}|");
    diffvariables_names_list_[TString("dphi")]=TString("#Delta #phi_{#gamma #gamma}");
    diffvariables_names_list_[TString("dR")]=TString("#Delta R_{#gamma #gamma}");
  }
  return diffvariables_names_list_[diffvariable];
};

std::map<TString,TString> diffvariables_units_list_;
TString diffvariables_units_list(TString diffvariable){
  if (diffvariables_units_list_.size()==0) {
    diffvariables_units_list_[TString("invmass")]=TString("GeV");
    diffvariables_units_list_[TString("diphotonpt")]=TString("GeV");
    diffvariables_units_list_[TString("costhetastar")]=TString("");
    diffvariables_units_list_[TString("dphi")]=TString("");
    diffvariables_units_list_[TString("dR")]=TString("");
  }
  return diffvariables_units_list_[diffvariable];
};




const Int_t n_histobins = 96;
const Float_t leftrange = -3;
const Float_t rightrange = 9;

const float default_threshold_adaptive_binning = -999;

const int n_templatebins_max = 1000; 
int n_templatebins = 0;
Double_t templatebinsboundaries[n_templatebins_max+1];

static const int n_bins=26;

int n_templates_pt=4;
float binsdef_single_gamma_pt[n_bins+1]={20,35,50,80,150};
int n_templates_EB_eta=7;
int n_templates_EE_eta=5;
float binsdef_single_gamma_EB_eta[n_bins+1]={0,0.2,0.4,0.6,0.8,1,1.2,1.4442};
float binsdef_single_gamma_EE_eta[n_bins+1]={1.56,1.653,1.8,2,2.2,2.5};

int n_templates_invmass_EBEB=16;
int n_templates_invmass_EBEE=16;
int n_templates_invmass_EEEE=16;
float binsdef_diphoton_invmass_EBEB[n_bins+1]={0,40,60,70,75,80,85,90,95,100,110,120,150,250,400,800,801};
float binsdef_diphoton_invmass_EBEE[n_bins+1]={0,40,60,70,75,80,85,90,95,100,110,120,150,250,400,800,801};
float binsdef_diphoton_invmass_EEEE[n_bins+1]={0,40,60,70,75,80,85,90,95,100,110,120,150,250,400,800,801};

int n_templates_diphotonpt_EBEB=21;
int n_templates_diphotonpt_EBEE=21;
int n_templates_diphotonpt_EEEE=21;
float binsdef_diphoton_diphotonpt_EBEB[n_bins+1]={0,6,10,12,14,16,18,20,22,24,28,34,40,50,60,70,80,90,100,120,200,201};
float binsdef_diphoton_diphotonpt_EBEE[n_bins+1]={0,6,10,12,14,16,18,20,22,24,28,34,40,50,60,70,80,90,100,120,200,201};
float binsdef_diphoton_diphotonpt_EEEE[n_bins+1]={0,6,10,12,14,16,18,20,22,24,28,34,40,50,60,70,80,90,100,120,200,201};

int n_templates_costhetastar_EBEB=8;
int n_templates_costhetastar_EBEE=8;
int n_templates_costhetastar_EEEE=8;
float binsdef_diphoton_costhetastar_EBEB[n_bins+1]={0,0.20,0.28,0.36,0.44,0.60,0.90,1.00,1.01};
float binsdef_diphoton_costhetastar_EBEE[n_bins+1]={0,0.20,0.28,0.36,0.44,0.60,0.90,1.00,1.01};
float binsdef_diphoton_costhetastar_EEEE[n_bins+1]={0,0.20,0.28,0.36,0.44,0.60,0.90,1.00,1.01};


int n_templates_dphi_EBEB=14;
int n_templates_dphi_EBEE=14;
int n_templates_dphi_EEEE=14;
const float Pi = TMath::Pi();
float binsdef_diphoton_dphi_EBEB[n_bins+1]={0,0.2*Pi,0.4*Pi,0.6*Pi,0.7*Pi,0.8*Pi,0.84*Pi,0.88*Pi,0.90*Pi,0.92*Pi,0.94*Pi,0.96*Pi,0.98*Pi,1.0*Pi,1.01*Pi};
float binsdef_diphoton_dphi_EBEE[n_bins+1]={0,0.2*Pi,0.4*Pi,0.6*Pi,0.7*Pi,0.8*Pi,0.84*Pi,0.88*Pi,0.90*Pi,0.92*Pi,0.94*Pi,0.96*Pi,0.98*Pi,1.0*Pi,1.01*Pi};
float binsdef_diphoton_dphi_EEEE[n_bins+1]={0,0.2*Pi,0.4*Pi,0.6*Pi,0.7*Pi,0.8*Pi,0.84*Pi,0.88*Pi,0.90*Pi,0.92*Pi,0.94*Pi,0.96*Pi,0.98*Pi,1.0*Pi,1.01*Pi};

int n_templates_dR_EBEB=8;
int n_templates_dR_EBEE=8;
int n_templates_dR_EEEE=8;
const float MaxDrExperiment = sqrt(pow(5.,2)+pow(Pi,2)); // = 5.90
float binsdef_diphoton_dR_EBEB[n_bins+1]={0.45,1.0,2.0,2.5,3.0,3.5,4.5,MaxDrExperiment,MaxDrExperiment+0.01};
float binsdef_diphoton_dR_EBEE[n_bins+1]={0.45,1.0,2.0,2.5,3.0,3.5,4.5,MaxDrExperiment,MaxDrExperiment+0.01};
float binsdef_diphoton_dR_EEEE[n_bins+1]={0.45,1.0,2.0,2.5,3.0,3.5,4.5,MaxDrExperiment,MaxDrExperiment+0.01};


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

const int n_eta_cats = n_templates_EB_eta;
int n_eta1eta2_cats = n_eta_cats*n_eta_cats;
float *etabins = binsdef_single_gamma_EB_eta+0;

// FOR PHOTON COMPONENT
// 020616 from data, no cleaning, no pf charged cut in presel
float eff_areas_EB_data[n_bins] = {2.615381e-01,2.587228e-01,2.641340e-01,2.640072e-01,2.569331e-01,2.395977e-01,1.651901e-01};
float eff_areas_EE_data[n_bins] = {5.848781e-02,8.411012e-02,1.189097e-01,1.446802e-01,1.148867e-01};
float eff_areas_EB_mc[n_bins] = {2.719313e-01,2.748353e-01,2.746818e-01,2.704172e-01,2.662482e-01,2.427560e-01,1.708246e-01};
float eff_areas_EE_mc[n_bins] = {5.107922e-02,8.485798e-02,1.354249e-01,1.671490e-01,1.454153e-01};



//const int n_ptbins_forreweighting = 1;
//Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={0,300};

const int n_ptbins_forreweighting = 3;
Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={0,40,60,999};


#endif
