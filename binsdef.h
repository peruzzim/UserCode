#ifndef __BINSDEF__
#define __BINSDEF__

static const int n_bins=9;

float binsdef_single_gamma_EB[n_bins+1]={30,40,50,60,70,80,90,110,140,150};
float binsdef_single_gamma_EE[n_bins+1]={30,40,50,60,70,80,90,110,120};

float binsdef_single_gamma_EB_eta[n_bins+1]={0,0.2,0.4,0.6,0.8,1,1.2,1.442};
float binsdef_single_gamma_EE_eta[n_bins+1]={1.56,1.653,1.8,2,2.2,2.5};

float binsdef_diphoton_EBEB[n_bins+1]={80,90,100,110,120,130,140,160,190,200};
//float binsdef_diphoton_EBEB[n_bins+1]={80,90,100,120,130};
float binsdef_diphoton_EBEE[n_bins+1]={80,90,100,120,130};
float binsdef_diphoton_EEEE[n_bins+1]={80,90,100,120,130};

int n_templates_EB=9;
int n_templates_EE=8;

//int n_templates_EB=7;
//int n_templates_EE=5;

int n_templates_EBEB=9;
//int n_templates_EBEB=4;
int n_templates_EBEE=4;
int n_templates_EEEE=4;

#endif
