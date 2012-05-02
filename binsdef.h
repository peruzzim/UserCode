#ifndef __BINSDEF__
#define __BINSDEF__

static const int n_bins=9;
float binsdef_single_gamma_EB[n_bins+1]={30,40,50,60,70,80,90,110,140,150};
float binsdef_single_gamma_EE[n_bins+1]={30,40,50,60,70,80,90,110,120};
float binsdef_diphoton_EBEB[n_bins+1]={80,90,100,110,120,130,140,160,190,200};
float binsdef_diphoton_EBEE[n_bins+1]={80,90};
float binsdef_diphoton_EEEE[n_bins+1]={80,90,100,120,130};

int n_templates_EB=9;
int n_templates_EE=8;
int n_templates_EBEB=9;
int n_templates_EBEE=1;
int n_templates_EEEE=4;

#endif
