bool global_doplots = true;
bool doxcheckstemplates = false;
bool dolightcomparisonwithstandardselsig = false;
bool dolightcomparisonwithstandardselbkg = false;

#include <assert.h>

#include "tdrstyle.C"

#include "binsdef.h"
#include "RooFitResult.h"
#include "TLatex.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooMCStudy.h"
#include "RooBinning.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include <stdio.h>
#include "RooNLLVar.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooMinimizer.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooClassFactory.h"
#include "TCanvas.h"
#include "RooConstraintSum.h"
#include "RooAddition.h"
#include "RooAbsDataStore.h"
#include "RooCachedPdf.h"
#include "RooThresholdCategory.h"
#include "TF1.h"
#include "TF2.h"
#include "TLegend.h"
#include "TSystem.h"

TRandom3 *_random_generator = new TRandom3(0);

using namespace std;
using namespace RooFit;

typedef struct {
  RooDataSet *dset;
  TString legend;
  Int_t color;
} plot_dataset_struct;

typedef struct {
  RooFitResult *fr_pass1;
  RooFitResult *fr_pass2constraint;
  RooFitResult *fr;
  float tot_events;
  float pp;
  float pp_err;
  float pf;
  float pf_err;
  float fp;
  float fp_err;
  float ff;
  float ff_err;
  float eff_overflow_removal_pp;
  RooAbsPdf *pdf_forgen[8];
  float fsig1_firstpass;
  float fsig2_firstpass;
  float chi2;
  int ndof;
  float probchi2;
  float fitpulls[100];
  bool lowstatbin[100];
} fit_output; 

const int numcpu=1;

ProcInfo_t procinfo;

void reweight_pt_2d(RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_pt_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void reweight_eta_2d(RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_eta_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void reweight_rho(RooDataSet **dset, RooDataSet *dsetdestination, RooPlot *plot);
void reweight_sigma(RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_rhosigma(RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold = kTRUE);
void validate_reweighting(RooDataSet *dset, RooDataSet *dsetdestination, int numvar);
void plot_datasets_axis1(std::vector<plot_dataset_struct> dsets, TString outname, TString legtitle, bool legendup=true, bool dolin=false);
void plot_template_dependency_axis1(RooDataSet *dset, TString variable, float min, float max, int bins, bool dobinned=0);
void produce_category_binning(RooDataSet **dset, bool deleteold=kTRUE);
void randomize_dataset_statistically_binned(RooDataSet **dset);
void create_histo_from_dataset_binned(RooDataSet *dset, TH1F **h1out, TH2F **h2out);
void create_histo_from_dataset_variablebins(RooDataSet *dset, TH1F **h1out, TH2F **h2out);
//void generate_toy_dataset_1d(RooDataSet **target, RooAbsPdf *sigpdf, RooAbsPdf *bkgpdf, float fsig1toy);
void generate_toy_dataset_2d(RooDataSet **target, RooAbsPdf *sigsigpdf, RooAbsPdf *sigbkgpdf, RooAbsPdf *bkgsigpdf, RooAbsPdf *bkgbkgpdf, float pptoy, float pftoy, float fptoy);
void print_mem();
void find_adaptive_binning(RooDataSet *dset, int *n_found_bins, Double_t *array_bounds, int axis=-1, float threshold = default_threshold_adaptive_binning);

TH1F* AddTHInQuadrature(std::vector<TH1F*> vector, TString name);

RooRealVar *roovar1=NULL;
RooRealVar *roovar2=NULL;
RooRealVar *roopt1=NULL;
RooRealVar *roopt2=NULL;
RooRealVar *roosieie1=NULL;
RooRealVar *roosieie2=NULL;
RooRealVar *rooeta1=NULL;
RooRealVar *rooeta2=NULL;
RooRealVar *roorho=NULL;
RooRealVar *roosigma=NULL;
RooRealVar *rooweight=NULL;
RooThresholdCategory *binning_roovar1_threshold=NULL;
RooThresholdCategory *binning_roovar2_threshold=NULL;
RooRealVar *binning_roovar1=NULL;
RooRealVar *binning_roovar2=NULL;


TFile *inputfile_t2p  = NULL;
TFile *inputfile_t1p1f = NULL;
TFile *inputfile_t2f   = NULL;
TFile *inputfile_d = NULL;
TDirectoryFile *dir_t2p=NULL;
TDirectoryFile *dir_t1p1f=NULL;
TDirectoryFile *dir_t2f=NULL;
TDirectoryFile *dir_d=NULL;

fit_output* fit_dataset(TString diffvariable, TString splitting, int bin, const TString do_syst_string=TString("")){

  std::cout << "Calling fit_dataset " << diffvariable.Data() << " " << splitting.Data() << " " << bin << " " << do_syst_string.Data() << std::endl;

  bool doplots = global_doplots;

  if (do_syst_string=="savepdfMCtrue2D")   {std::cout << "RUNNING FOR SAVE PDF 2D" << std::endl; doplots=false;}
  if (do_syst_string=="savepdfMCtrue1D") {std::cout << "RUNNING FOR SAVE PDF 1D, IGNORE THE RESULTS OF THE FIT!!!" << std::endl; doplots=false;}

  if (do_syst_string=="data_donotwriteoutpurity") doplots=false;
  if (do_syst_string=="subtractionZee") doplots=false;

  bool doplots_b = doplots;
  bool doplots_ub = doplots;
//  bool doplots_b = true;
//  bool doplots_ub = false;

  TH1F::SetDefaultSumw2(kTRUE);

  int bins_to_run=-1; 
  float *binsdef=NULL;

  if (diffvariable=="invmass"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_invmass_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_invmass_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_invmass_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_invmass_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_invmass_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_invmass_EEEE;
  }
  if (diffvariable=="diphotonpt"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_diphotonpt_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_diphotonpt_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_diphotonpt_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_diphotonpt_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_diphotonpt_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_diphotonpt_EEEE;
  }
  if (diffvariable=="costhetastar"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_costhetastar_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_costhetastar_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_costhetastar_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_costhetastar_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_costhetastar_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_costhetastar_EEEE;
  }
  if (diffvariable=="dphi"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_dphi_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_dphi_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_dphi_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_dphi_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_dphi_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_dphi_EEEE;
  }
  if (diffvariable=="dR"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_dR_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_dR_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_dR_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_dR_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_dR_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_dR_EEEE;
  }

  if (bins_to_run==bin) bins_to_run=-1;

  fit_output *out=NULL;

  for (int k=0; k<2; k++) std::cout << std::endl;
  std::cout << "Process " << diffvariable.Data() << " bin " << bin << std::endl;
  for (int k=0; k<2; k++) std::cout << std::endl;

  TString inputfilename_t2p;
  TString inputfilename_t1p1f;
  TString inputfilename_t2f;
  TString inputfilename_d;

  if (do_syst_string=="savepdfMCtrue2D" ||  do_syst_string=="doMCtrue"){
    inputfilename_t2p   = "outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton_allmc_1p1fbothgen.root";
    inputfilename_t2f   = "outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton_allmc_standard.root";
  }
  else if (do_syst_string=="doMCtrue_2frag"){
    inputfilename_t2p   = "outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton_allmc_1p1fbothgen.root";
    inputfilename_t2f   = "outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton_allmc_standard_2frag.root";
  }
  else if (do_syst_string=="savepdfMCtrue1D" || do_syst_string=="templateshape2frag"){
    inputfilename_t2p   = "outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton_allmc_1p1fbothgen.root";
    inputfilename_t2f   = "outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton_data_standard.root";
  }
  else if (do_syst_string=="doMCpromptdriven"){
    inputfilename_t2p   = "outphoton_allmc_sigsig.root";
    inputfilename_t1p1f = "outphoton_allmc_1prcone1fgen.root";
    inputfilename_t2f   = "outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton_allmc_standard.root";
  }
  else if (do_syst_string=="templateshapeMCpromptdrivenEB" || do_syst_string=="templateshapeMCpromptdrivenEE"){
    inputfilename_t2p   = "outphoton_allmc_sigsig.root";
    inputfilename_t1p1f = "outphoton_allmc_1prcone1fgen.root";
    inputfilename_t2f   = "outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton_data_standard.root";
  }
  else if (do_syst_string=="doMCfakedriven"){
    inputfilename_t2p   = "outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton_allmc_1pgen1fside.root";
    inputfilename_t2f   = "outphoton_allmc_bkgbkg.root";
    inputfilename_d     = "outphoton_allmc_standard.root";
  }
  else if (do_syst_string=="templateshapeMCfakedrivenEB" || do_syst_string=="templateshapeMCfakedrivenEE"){
    inputfilename_t2p   = "outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton_allmc_1pgen1fside.root";
    inputfilename_t2f   = "outphoton_allmc_bkgbkg.root";
    inputfilename_d     = "outphoton_data_standard.root";
  }
  else if (do_syst_string=="subtractionZee"){
    inputfilename_t2p   = "outphoton_allmc_2pgen.root";
    inputfilename_t1p1f = "outphoton_allmc_1p1fbothgen.root";
    inputfilename_t2f   = "outphoton_allmc_2fgen.root";
    inputfilename_d     = "outphoton_allmc_standard.root";
  }  
  else if (do_syst_string=="doMCfulldriven") {
    inputfilename_t2p   = "outphoton_allmc_sigsig.root";
    inputfilename_t1p1f = "outphoton_allmc_sigbkg.root";
    inputfilename_t2f   = "outphoton_allmc_bkgbkg.root";
    inputfilename_d     = "outphoton_allmc_standard_1p1fbothgen_fakeaxis1.root";
  }  
  else if (do_syst_string=="newtemplates") {
    inputfilename_t2p   = "outphoton_data_sigsig_step2_1event.root";
    inputfilename_t1p1f = "outphoton_data_sigbkg_step2_1event.root";
    inputfilename_t2f   = "outphoton_data_bkgbkg_step2_1event.root";
//    inputfilename_t2p   = "outphoton_data_sigsig.root";
//    inputfilename_t1p1f = "outphoton_data_sigbkg.root";
//    inputfilename_t2f   = "outphoton_data_bkgbkg.root";
    inputfilename_d     = "outphoton_data_standard.root";
  }  
  else {
    inputfilename_t2p   = "outphoton_data_sigsig.root";
    inputfilename_t1p1f = "outphoton_data_sigbkg.root";
    inputfilename_t2f   = "outphoton_data_bkgbkg.root";
    inputfilename_d     = "outphoton_data_standard.root";
  }  

  if ((!inputfile_t2p)   ||  (TString(inputfile_t2p->GetName())   != TString(inputfilename_t2p)  )) {inputfile_t2p = TFile::Open(inputfilename_t2p);     dir_t2p=NULL;  }
  if ((!inputfile_t1p1f) ||  (TString(inputfile_t1p1f->GetName()) != TString(inputfilename_t1p1f))) {inputfile_t1p1f = TFile::Open(inputfilename_t1p1f); dir_t1p1f=NULL;}
  if ((!inputfile_t2f)   ||  (TString(inputfile_t2f->GetName())   != TString(inputfilename_t2f)  )) {inputfile_t2f = TFile::Open(inputfilename_t2f);     dir_t2f=NULL;  }
  if ((!inputfile_d)     ||  (TString(inputfile_d->GetName())     != TString(inputfilename_d)    )) {inputfile_d = TFile::Open(inputfilename_d);         dir_d=NULL;    }

  if (splitting=="EEEB") splitting="EBEE";

  TH1::SetDefaultSumw2(kTRUE);
  
  TString s1; TString s2;
  if (splitting=="EBEB") {s1="EB"; s2="EB";}
  else if (splitting=="EEEE") {s1="EE"; s2="EE";}
  else if (splitting=="EBEE") {s1="EB"; s2="EE";}
  bool sym  = (s1==s2);
  
  if(!dir_t2p)   inputfile_t2p->GetObject("roofit",dir_t2p);
  if(!dir_t1p1f) inputfile_t1p1f->GetObject("roofit",dir_t1p1f);
  if(!dir_t2f)   inputfile_t2f->GetObject("roofit",dir_t2f);
  if(!dir_d)     inputfile_d->GetObject("roofit",dir_d);

  assert(dir_t2p);
  assert(dir_t1p1f);
  assert(dir_t2f);
  assert(dir_d);

  dir_d->GetObject("roovar1",roovar1);
  dir_d->GetObject("roovar2",roovar2);
  roovar1->setRange(leftrange,rightrange);
  roovar2->setRange(leftrange,rightrange);
  roovar1->setBins(n_histobins);
  roovar2->setBins(n_histobins);
  roovar1->SetTitle("Iso_{1}");
  roovar2->SetTitle("Iso_{2}");
  {
    TString unit = diffvariables_units_list(diffvariable);
    if(unit!=TString("")) {
      roovar1->setUnit(unit.Data());
      roovar2->setUnit(unit.Data());
    }
  }

  dir_d->GetObject("roopt1",roopt1); 
  dir_d->GetObject("roosieie1",roosieie1); 
  dir_d->GetObject("rooeta1",rooeta1); 
  dir_d->GetObject("roopt2",roopt2); 
  dir_d->GetObject("roosieie2",roosieie2); 
  dir_d->GetObject("rooeta2",rooeta2); 
  dir_d->GetObject("roorho",roorho); 
  dir_d->GetObject("roosigma",roosigma); 
  rooweight = new RooRealVar("rooweight","rooweight",0,100);
  assert (roovar1);
  assert (roovar2);
  assert (roopt1);
  assert (roosieie1);
  assert (rooeta1);
  assert (roopt2);
  assert (roosieie2);
  assert (rooeta2);
  assert (roorho);
  assert (roosigma);
  assert (rooweight);

  const float pp_init = 0.2;
  const float pf_init = 0.23;
  const float fp_init = 0.23;
    
  RooDataSet *dataset_sigsig_orig = NULL;
  RooDataSet *dataset_sigbkg_orig = NULL;
  RooDataSet *dataset_bkgsig_orig = NULL;
  RooDataSet *dataset_bkgbkg_orig = NULL;
  RooDataSet *dataset_orig =        NULL;

  if (do_syst_string!="newtemplates"){
    dir_t2p->GetObject(Form("template_roodset_%s_sigsig",splitting.Data()),dataset_sigsig_orig);
    dir_t1p1f->GetObject(Form("template_roodset_%s_sigbkg",splitting.Data()),dataset_sigbkg_orig);
    dir_t1p1f->GetObject(Form("template_roodset_%s_bkgsig",splitting.Data()),dataset_bkgsig_orig);
    dir_t2f->GetObject(Form("template_roodset_%s_bkgbkg",splitting.Data()),dataset_bkgbkg_orig);
  }
  else {
    dir_t2p->GetObject(Form("newtempl_roodset_%s_%s_b%d_sigsig",splitting.Data(),diffvariable.Data(),bin),dataset_sigsig_orig);
    dir_t1p1f->GetObject(Form("newtempl_roodset_%s_%s_b%d_sigbkg",splitting.Data(),diffvariable.Data(),bin),dataset_sigbkg_orig);
    dir_t1p1f->GetObject(Form("newtempl_roodset_%s_%s_b%d_bkgsig",splitting.Data(),diffvariable.Data(),bin),dataset_bkgsig_orig);
    dir_t2f->GetObject(Form("newtempl_roodset_%s_%s_b%d_bkgbkg",splitting.Data(),diffvariable.Data(),bin),dataset_bkgbkg_orig);
//    dir_t2p->GetObject(Form("template_roodset_%s_sigsig",splitting.Data()),dataset_sigsig_orig);
//    dir_t1p1f->GetObject(Form("template_roodset_%s_sigbkg",splitting.Data()),dataset_sigbkg_orig);
//    dir_t1p1f->GetObject(Form("template_roodset_%s_bkgsig",splitting.Data()),dataset_bkgsig_orig);
//    dir_t2f->GetObject(Form("template_roodset_%s_bkgbkg",splitting.Data()),dataset_bkgbkg_orig);
  }
  dir_d->GetObject(Form("obs_roodset_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin),dataset_orig);
  assert(dataset_sigsig_orig);
  assert(dataset_sigbkg_orig);
  assert(dataset_bkgsig_orig);
  assert(dataset_bkgbkg_orig);
  assert(dataset_orig);

  RooDataSet *dataset_sigsig = (RooDataSet*)(dataset_sigsig_orig->reduce(Name("dataset_sigsig"),Cut(Form("roovar1<%f && roovar2<%f",rightrange-1e-5,rightrange-1e-5))));
  RooDataSet *dataset_sigbkg = (RooDataSet*)(dataset_sigbkg_orig->reduce(Name("dataset_sigbkg"),Cut(Form("roovar1<%f && roovar2<%f",rightrange-1e-5,rightrange-1e-5))));
  RooDataSet *dataset_bkgsig = (RooDataSet*)(dataset_bkgsig_orig->reduce(Name("dataset_bkgsig"),Cut(Form("roovar1<%f && roovar2<%f",rightrange-1e-5,rightrange-1e-5))));
  RooDataSet *dataset_bkgbkg = (RooDataSet*)(dataset_bkgbkg_orig->reduce(Name("dataset_bkgbkg"),Cut(Form("roovar1<%f && roovar2<%f",rightrange-1e-5,rightrange-1e-5))));
  RooDataSet *dataset =        (RooDataSet*)(dataset_orig->reduce(Name("dataset"),Cut(Form("roovar1<%f && roovar2<%f",rightrange-1e-5,rightrange-1e-5))));

  std::cout << "2D datasets" << std::endl;
  dataset_sigsig->Print();
  dataset_sigbkg->Print();
  dataset_bkgsig->Print();
  dataset_bkgbkg->Print();
  dataset_orig->Print();
  dataset->Print();
  const float eff_overflow_removal = dataset_sigsig->sumEntries()/dataset_sigsig_orig->sumEntries();
  std::cout << "TO BE DONE BETTER AFTER REWEIGHTING: Efficiency of overflow removal: " << eff_overflow_removal << std::endl;
  //  assert (eff_overflow_removal>0.995);


  RooDataSet *dataset_sig_axis1 = (RooDataSet*)(dataset_sigsig->reduce(Name("dataset_sig_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma))));
  RooDataSet *dataset_bkg_axis1 = (RooDataSet*)(dataset_bkgsig->reduce(Name("dataset_bkg_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma))));
  RooDataSet *dataset_sig_axis2 = (RooDataSet*)(dataset_sigsig->reduce(Name("dataset_sig_axis2"),SelectVars(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma))));
  RooDataSet *dataset_bkg_axis2 = (RooDataSet*)(dataset_sigbkg->reduce(Name("dataset_bkg_axis2"),SelectVars(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma))));

  RooDataSet *dataset_axis1 = (RooDataSet*)(dataset->reduce(Name("dataset_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma))));
  RooDataSet *dataset_axis2 = (RooDataSet*)(dataset->reduce(Name("dataset_axis2"),SelectVars(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma))));

  std::cout << "1D datasets" << std::endl;
  dataset_sig_axis1->Print();
  dataset_sig_axis2->Print();
  dataset_bkg_axis1->Print();
  dataset_bkg_axis2->Print();
  dataset_axis1->Print();
  dataset_axis2->Print();
  
  if (do_syst_string!="newtemplates"){
    { // rhosigma reweighting
      reweight_rhosigma(&dataset_sigsig,dataset);
      reweight_rhosigma(&dataset_sigbkg,dataset);
      reweight_rhosigma(&dataset_bkgsig,dataset);
      reweight_rhosigma(&dataset_bkgbkg,dataset);
      reweight_rhosigma(&dataset_sig_axis1,dataset_axis1);
      reweight_rhosigma(&dataset_bkg_axis1,dataset_axis1);
      reweight_rhosigma(&dataset_sig_axis2,dataset_axis2);
      reweight_rhosigma(&dataset_bkg_axis2,dataset_axis2);
    }
    
    { // eta reweighting
      reweight_eta_2d(&dataset_sigsig,dataset);
      reweight_eta_2d(&dataset_sigbkg,dataset);
      reweight_eta_2d(&dataset_bkgsig,dataset);
      reweight_eta_2d(&dataset_bkgbkg,dataset);
      reweight_eta_1d(&dataset_sig_axis1,dataset_axis1,1);
      reweight_eta_1d(&dataset_bkg_axis1,dataset_axis1,1);
      reweight_eta_1d(&dataset_sig_axis2,dataset_axis2,2);
      reweight_eta_1d(&dataset_bkg_axis2,dataset_axis2,2);
    }
    
    if (!(diffvariable=="invmass" && splitting=="EEEE" && bin==14)) { // pt reweighting
      //      reweight_pt_2d(&dataset_sigsig,dataset);
      reweight_pt_2d(&dataset_sigbkg,dataset);
      reweight_pt_2d(&dataset_bkgsig,dataset);
      reweight_pt_2d(&dataset_bkgbkg,dataset);
      //      reweight_pt_1d(&dataset_sig_axis1,dataset_axis1,1);
      reweight_pt_1d(&dataset_bkg_axis1,dataset_axis1,1);
      //      reweight_pt_1d(&dataset_sig_axis2,dataset_axis2,2);
      reweight_pt_1d(&dataset_bkg_axis2,dataset_axis2,2);
    }
  }

  /*
    { // validate reweighting
    validate_reweighting(dataset_sigsig,dataset,1);
    validate_reweighting(dataset_sigbkg,dataset,1);
    validate_reweighting(dataset_bkgsig,dataset,1);
    validate_reweighting(dataset_bkgbkg,dataset,1);
    validate_reweighting(dataset_sigsig,dataset,2);
    validate_reweighting(dataset_sigbkg,dataset,2);
    validate_reweighting(dataset_bkgsig,dataset,2);
    validate_reweighting(dataset_bkgbkg,dataset,2);
    validate_reweighting(dataset_sig_axis1,dataset_axis1,1);
    validate_reweighting(dataset_bkg_axis1,dataset_axis1,1);
    validate_reweighting(dataset_sig_axis2,dataset_axis2,2);
    validate_reweighting(dataset_bkg_axis2,dataset_axis2,2);
    }
  */




  RooDataSet *dset_mctrue_s = NULL;
  RooDataSet *dset_mcfrag_s = NULL;
  RooDataSet *dset_mcnofrag_s = NULL;
  RooDataSet *dset_mcrcone_s = NULL;
  RooDataSet *dset_zee_s = NULL;
  RooDataSet *dset_mctrue_b = NULL;
  RooDataSet *dset_mcrcone_b = NULL;
  RooDataSet *dset_mctrue_noEM = NULL;
  RooDataSet *dset_mcrcone_noEM = NULL;

  if (doxcheckstemplates || dolightcomparisonwithstandardselsig || dolightcomparisonwithstandardselbkg) {

    if (dolightcomparisonwithstandardselsig || dolightcomparisonwithstandardselbkg) if (splitting=="EBEE") return NULL;

    TFile *fmctrue_s = new TFile("outphoton_allmc_sig.root","read");
    fmctrue_s->GetObject(Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mctrue_s);
    assert(dset_mctrue_s);

    TFile *fmcfrag_s = new TFile("outphoton_allmc_frag.root","read");
    fmcfrag_s->GetObject(Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mcfrag_s);
    assert(dset_mcfrag_s);
  
    TFile *fmcnofrag_s = new TFile("outphoton_allmc_nofrag.root","read");
    fmcnofrag_s->GetObject(Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mcnofrag_s);
    assert(dset_mcnofrag_s);
  
    TFile *fmcrcone_s = new TFile("outphoton_allmc_rcone.root","read");
    fmcrcone_s->GetObject(Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mcrcone_s);
    assert(dset_mcrcone_s);
  
    TFile *fzee_s = new TFile("outphoton_data_zee.root","read");
    RooDataSet *dset_zee_s_2d = NULL;
    fzee_s->GetObject(Form("roofit/template_roodset_%s_sigsig",splitting.Data()),dset_zee_s_2d);
    assert(dset_zee_s_2d);
    dset_zee_s = (RooDataSet*)(dset_zee_s_2d->reduce(Name("dset_zee_s"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma))));
    assert(dset_zee_s);

    TFile *fmctrue_b = new TFile("outphoton_allmc_bkg.root","read");
    fmctrue_b->GetObject(Form("roofit/roodset_background_%s_rv1",s1.Data()),dset_mctrue_b);
    assert(dset_mctrue_b);

    TFile *fmcrcone_b = new TFile("outphoton_allmc_sieiesideband.root","read");
    fmcrcone_b->GetObject(Form("roofit/roodset_background_%s_rv1",s1.Data()),dset_mcrcone_b);
    assert(dset_mcrcone_b);

    TFile *fmctrue_noEM = new TFile("outphoton_allmc_bkg_noEMenr.root","read");
    fmctrue_noEM->GetObject(Form("roofit/roodset_background_%s_rv1",s1.Data()),dset_mctrue_noEM);
    assert(dset_mctrue_noEM);

    TFile *fmcrcone_noEM = new TFile("outphoton_allmc_sieiesideband_noEMenr.root","read");
    fmcrcone_noEM->GetObject(Form("roofit/roodset_background_%s_rv1",s1.Data()),dset_mcrcone_noEM);
    assert(dset_mcrcone_noEM);
  
    dset_mctrue_s = (RooDataSet*)(dset_mctrue_s->reduce(Name("dset_mctrue_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mcfrag_s = (RooDataSet*)(dset_mcfrag_s->reduce(Name("dset_mcfrag_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mcnofrag_s = (RooDataSet*)(dset_mcnofrag_s->reduce(Name("dset_mcnofrag_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mcrcone_s = (RooDataSet*)(dset_mcrcone_s->reduce(Name("dset_mcrcone_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mctrue_b = (RooDataSet*)(dset_mctrue_b->reduce(Name("dset_mctrue_b"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mcrcone_b = (RooDataSet*)(dset_mcrcone_b->reduce(Name("dset_mcrcone_b"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mctrue_noEM = (RooDataSet*)(dset_mctrue_noEM->reduce(Name("dset_mctrue_noEM"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mcrcone_noEM = (RooDataSet*)(dset_mcrcone_noEM->reduce(Name("dset_mcrcone_noEM"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_zee_s = (RooDataSet*)(dset_zee_s->reduce(Name("dset_zee_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));

    std::cout << "MC datasets" << std::endl;
    dset_mctrue_s->Print();
    dset_mcfrag_s->Print();
    dset_mcnofrag_s->Print();
    dset_mcrcone_s->Print();
    dset_zee_s->Print();
    dset_mctrue_b->Print();
    dset_mcrcone_b->Print();
    dset_mctrue_noEM->Print();
    dset_mcrcone_noEM->Print();

    reweight_rhosigma(&dset_mctrue_s,dataset_axis1);
    reweight_rhosigma(&dset_mcfrag_s,dataset_axis1);
    reweight_rhosigma(&dset_mcnofrag_s,dataset_axis1);
    reweight_rhosigma(&dset_mcrcone_s,dataset_axis1);
    reweight_rhosigma(&dset_zee_s,dataset_axis1);
    reweight_eta_1d(&dset_mctrue_s,dataset_axis1,1);
    reweight_eta_1d(&dset_mcfrag_s,dataset_axis1,1);
    reweight_eta_1d(&dset_mcnofrag_s,dataset_axis1,1);
    reweight_eta_1d(&dset_mcrcone_s,dataset_axis1,1);
    reweight_eta_1d(&dset_zee_s,dataset_axis1,1);
    reweight_pt_1d(&dset_mctrue_s,dataset_axis1,1);
    reweight_pt_1d(&dset_mcfrag_s,dataset_axis1,1);
    reweight_pt_1d(&dset_mcnofrag_s,dataset_axis1,1);
    //reweight_pt_1d(&dset_mcrcone_s,dataset_axis1,1);
    reweight_pt_1d(&dset_zee_s,dataset_axis1,1);
    reweight_rhosigma(&dset_mctrue_b,dataset_axis1);
    reweight_rhosigma(&dset_mcrcone_b,dataset_axis1);
    reweight_rhosigma(&dset_mctrue_noEM,dataset_axis1);
    reweight_rhosigma(&dset_mcrcone_noEM,dataset_axis1);
    reweight_eta_1d(&dset_mctrue_b,dataset_axis1,1);
    reweight_eta_1d(&dset_mcrcone_b,dataset_axis1,1);
    reweight_eta_1d(&dset_mctrue_noEM,dataset_axis1,1);
    reweight_eta_1d(&dset_mcrcone_noEM,dataset_axis1,1);
    reweight_pt_1d(&dset_mctrue_b,dataset_axis1,1);
    reweight_pt_1d(&dset_mcrcone_b,dataset_axis1,1);
    reweight_pt_1d(&dset_mctrue_noEM,dataset_axis1,1);
    reweight_pt_1d(&dset_mcrcone_noEM,dataset_axis1,1);

    RooDataSet *dset_mcrcone_b1 = NULL;
    RooDataSet *dset_mcrcone_b2 = NULL;
    RooDataSet *dataset_bkg_axis1_1 = NULL;
    RooDataSet *dataset_bkg_axis1_2 = NULL;
    if (s1=="EB") {
      dset_mcrcone_b1 = (RooDataSet*)(dset_mcrcone_b->reduce(Name("dset_mcrcone_b1"),Cut(Form("roosieie1<%f",0.0125))));
      dset_mcrcone_b2 = (RooDataSet*)(dset_mcrcone_b->reduce(Name("dset_mcrcone_b2"),Cut(Form("roosieie1>%f",0.0125))));
      dataset_bkg_axis1_1 = (RooDataSet*)(dataset_bkg_axis1->reduce(Name("dataset_bkg_axis1_1"),Cut(Form("roosieie1<%f",0.0125))));
      dataset_bkg_axis1_2 = (RooDataSet*)(dataset_bkg_axis1->reduce(Name("dataset_bkg_axis1_2"),Cut(Form("roosieie1>%f",0.0125))));
    }
    if (s1=="EE") {
      dset_mcrcone_b1 = (RooDataSet*)(dset_mcrcone_b->reduce(Name("dset_mcrcone_b1"),Cut(Form("roosieie1<%f",0.032))));
      dset_mcrcone_b2 = (RooDataSet*)(dset_mcrcone_b->reduce(Name("dset_mcrcone_b2"),Cut(Form("roosieie1>%f",0.032))));
      dataset_bkg_axis1_1 = (RooDataSet*)(dataset_bkg_axis1->reduce(Name("dataset_bkg_axis1_1"),Cut(Form("roosieie1<%f",0.032))));
      dataset_bkg_axis1_2 = (RooDataSet*)(dataset_bkg_axis1->reduce(Name("dataset_bkg_axis1_2"),Cut(Form("roosieie1>%f",0.032))));
    }


    plot_dataset_struct str_dataset_axis1;
    str_dataset_axis1.dset = dataset_axis1;
    str_dataset_axis1.legend = "Photon Iso in selection";
    str_dataset_axis1.color = kGreen;
    plot_dataset_struct str_dataset_sig_axis1;
    str_dataset_sig_axis1.dset = dataset_sig_axis1;
    str_dataset_sig_axis1.legend = "Rand. cone in data";
    str_dataset_sig_axis1.color = kBlack;
    plot_dataset_struct str_dset_mctrue_s;
    str_dset_mctrue_s.dset = dset_mctrue_s;
    str_dset_mctrue_s.legend = "Photon Iso in MC";
    str_dset_mctrue_s.color = kRed;
    plot_dataset_struct str_dset_mcrcone_s;
    str_dset_mcrcone_s.dset = dset_mcrcone_s;
    str_dset_mcrcone_s.legend = "Rand. cone in MC";
    str_dset_mcrcone_s.color = kBlue;
    plot_dataset_struct str_dset_zee_s;
    str_dset_zee_s.dset = dset_zee_s;
    str_dset_zee_s.legend = "Zee in data";
    str_dset_zee_s.color = kGreen+2;
    plot_dataset_struct str_dset_mcnofrag_s;
    str_dset_mcnofrag_s.dset = dset_mcnofrag_s;
    str_dset_mcnofrag_s.legend = "Direct photon Iso in MC";
    str_dset_mcnofrag_s.color = kCyan;
    plot_dataset_struct str_dset_mcfrag_s;
    str_dset_mcfrag_s.dset = dset_mcfrag_s;
    str_dset_mcfrag_s.legend = "Frag. photon Iso in MC";
    str_dset_mcfrag_s.color = kOrange;
    plot_dataset_struct str_dataset_bkg_axis1;
    str_dataset_bkg_axis1.dset = dataset_bkg_axis1;
    str_dataset_bkg_axis1.legend = "Sieie sideband in data";
    str_dataset_bkg_axis1.color = kBlack;
    plot_dataset_struct str_dataset_bkg_axis1_1;
    str_dataset_bkg_axis1_1.dset = dataset_bkg_axis1_1;
    str_dataset_bkg_axis1_1.legend = "Sieie sideband in data / left";
    str_dataset_bkg_axis1_1.color = kBlack;
    plot_dataset_struct str_dataset_bkg_axis1_2;
    str_dataset_bkg_axis1_2.dset = dataset_bkg_axis1_2;
    str_dataset_bkg_axis1_2.legend = "Sieie sideband in data / right";
    str_dataset_bkg_axis1_2.color = kBlack;
    plot_dataset_struct str_dset_mctrue_b;
    str_dset_mctrue_b.dset = dset_mctrue_b;
    str_dset_mctrue_b.legend = "Photon Iso in MC fakes";
    str_dset_mctrue_b.color = kRed;
    plot_dataset_struct str_dset_mcrcone_b;
    str_dset_mcrcone_b.dset = dset_mcrcone_b;
    str_dset_mcrcone_b.legend = "Sieie sideband in MC";
    str_dset_mcrcone_b.color = kBlue;
    plot_dataset_struct str_dset_mctrue_noEM;
    str_dset_mctrue_noEM.dset = dset_mctrue_noEM;
    str_dset_mctrue_noEM.legend = "Photon Iso in MC fakes, no EM enr.";
    str_dset_mctrue_noEM.color = kOrange;
    plot_dataset_struct str_dset_mcrcone_noEM;
    str_dset_mcrcone_noEM.dset = dset_mcrcone_noEM;
    str_dset_mcrcone_noEM.legend = "Sieie sideband in MC, no EM enr.";
    str_dset_mcrcone_noEM.color = kMagenta;
    plot_dataset_struct str_dset_mcrcone_b1;
    str_dset_mcrcone_b1.dset = dset_mcrcone_b1;
    str_dset_mcrcone_b1.legend = "Sieie sideband in MC / left";
    str_dset_mcrcone_b1.color = kBlue;
    plot_dataset_struct str_dset_mcrcone_b2;
    str_dset_mcrcone_b2.dset = dset_mcrcone_b2;
    str_dset_mcrcone_b2.legend = "Sieie sideband in MC / right";
    str_dset_mcrcone_b2.color = kBlue;

    if (dolightcomparisonwithstandardselsig){
      {
	std::vector<plot_dataset_struct> vec;
	vec.push_back(str_dataset_axis1);
	vec.push_back(str_dataset_sig_axis1);
	vec.push_back(str_dset_mctrue_s);
	vec.push_back(str_dset_mcrcone_s);
	plot_datasets_axis1(vec,Form("plots/histo_template_sig_compwithsel_%s_log_%s_%s_b%d",s1.Data(),diffvariable.Data(),splitting.Data(),bin),Form("Signal template %s",s1.Data()));
	plot_datasets_axis1(vec,Form("plots/histo_template_sig_compwithsel_%s_lin_%s_%s_b%d",s1.Data(),diffvariable.Data(),splitting.Data(),bin),Form("Signal template %s",s1.Data()),true,true);
      }
      return NULL;
    }
    if (dolightcomparisonwithstandardselbkg){
      {
	std::vector<plot_dataset_struct> vec;
	vec.push_back(str_dataset_axis1);
	vec.push_back(str_dataset_bkg_axis1);
	vec.push_back(str_dset_mctrue_b);
	vec.push_back(str_dset_mcrcone_b);
	plot_datasets_axis1(vec,Form("plots/histo_template_bkg_compwithsel_%s_log_%s_%s_b%d",s1.Data(),diffvariable.Data(),splitting.Data(),bin),Form("Background template %s",s1.Data()),false);
	plot_datasets_axis1(vec,Form("plots/histo_template_bkg_compwithsel_%s_lin_%s_%s_b%d",s1.Data(),diffvariable.Data(),splitting.Data(),bin),Form("Background template %s",s1.Data()),true,true);
      }
      return NULL;
    }

    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dataset_sig_axis1);
    vec.push_back(str_dset_mctrue_s);
    vec.push_back(str_dset_mcrcone_s);
    vec.push_back(str_dset_zee_s);
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_withZ_%s_log",s1.Data()),Form("Signal template %s",s1.Data()));
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_withZ_%s_lin",s1.Data()),Form("Signal template %s",s1.Data()),true,true);
    }
    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dataset_sig_axis1);
    vec.push_back(str_dset_mctrue_s);
    vec.push_back(str_dset_mcrcone_s);
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_%s_log",s1.Data()),Form("Signal template %s",s1.Data()));
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_%s_lin",s1.Data()),Form("Signal template %s",s1.Data()),true,true);
    }
    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dset_mctrue_s);
    vec.push_back(str_dset_mcrcone_s);
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_onlyMC_%s_log",s1.Data()),Form("Signal template %s",s1.Data()));
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_onlyMC_%s_lin",s1.Data()),Form("Signal template %s",s1.Data()),true,true);
    }
    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dset_mctrue_s);
    vec.push_back(str_dset_mcnofrag_s);
    vec.push_back(str_dset_mcfrag_s);
    vec.push_back(str_dset_mcrcone_s);
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_onlyMCwithfrag_%s_log",s1.Data()),Form("Signal template %s",s1.Data()));
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_onlyMCwithfrag_%s_lin",s1.Data()),Form("Signal template %s",s1.Data()),true,true);
    }
    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dataset_sig_axis1);
    vec.push_back(str_dset_mctrue_s);
    vec.push_back(str_dset_mcrcone_s);
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_%s_log",s1.Data()),Form("Signal template %s",s1.Data()));
    plot_datasets_axis1(vec,Form("plots/histo_template_sig_%s_lin",s1.Data()),Form("Signal template %s",s1.Data()),true,true);
    }

    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dataset_bkg_axis1);
    vec.push_back(str_dset_mctrue_b);
    vec.push_back(str_dset_mcrcone_b);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_%s_log",s1.Data()),Form("Background template %s",s1.Data()),false);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_%s_lin",s1.Data()),Form("Background template %s",s1.Data()),true,true);
    }
    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dset_mctrue_b);
    vec.push_back(str_dset_mcrcone_b);
    vec.push_back(str_dset_mctrue_noEM);
    vec.push_back(str_dset_mcrcone_noEM);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_onlyMC_%s_log",s1.Data()),Form("Background template %s",s1.Data()),false);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_onlyMC_%s_lin",s1.Data()),Form("Background template %s",s1.Data()),true,true);
    }
    {
    std::vector<plot_dataset_struct> vec;
    vec.push_back(str_dataset_bkg_axis1);
    vec.push_back(str_dataset_bkg_axis1_1);
    vec.push_back(str_dataset_bkg_axis1_2);
    vec.push_back(str_dset_mctrue_b);
    vec.push_back(str_dset_mcrcone_b);
    vec.push_back(str_dset_mcrcone_b1);
    vec.push_back(str_dset_mcrcone_b2);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_slicesieie_%s_log",s1.Data()),Form("Background template %s",s1.Data()),false);
    plot_datasets_axis1(vec,Form("plots/histo_template_bkg_slicesieie_%s_lin",s1.Data()),Form("Background template %s",s1.Data()),true,true);
    }

    return NULL;
  
  }



  RooDataSet *dset_mctrue_s_rv1 = NULL;
  RooDataSet *dset_mcrcone_s_rv1 = NULL;
  RooDataSet *dset_mctrue_s_rv2 = NULL;
  RooDataSet *dset_mcrcone_s_rv2 = NULL;
  RooDataSet *dset_mctrue_b_rv1 = NULL;
  RooDataSet *dset_mcrcone_b_rv1 = NULL;
  RooDataSet *dset_mctrue_b_rv2 = NULL;
  RooDataSet *dset_mcrcone_b_rv2 = NULL;

  if (do_syst_string==TString("savepdfMCtrue1D") || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")) {

    TFile *fmctrue_s = new TFile("outphoton_allmc_sig.root","read");
    fmctrue_s->GetObject(Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mctrue_s_rv1);
    fmctrue_s->GetObject(Form("roofit/roodset_signal_%s_rv2",s2.Data()),dset_mctrue_s_rv2);
    assert(dset_mctrue_s_rv1);
    assert(dset_mctrue_s_rv2);
    TFile *fmcrcone_s;
    if (do_syst_string==TString("templateshape2frag")){
      fmcrcone_s = new TFile("outphoton_allmc_sig_2frag.root","read");
      fmcrcone_s->GetObject(Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mcrcone_s_rv1);
      fmcrcone_s->GetObject(Form("roofit/roodset_signal_%s_rv2",s2.Data()),dset_mcrcone_s_rv2);
    }
    else {
      fmcrcone_s = new TFile("outphoton_allmc_rcone.root","read");
      fmcrcone_s->GetObject(Form("roofit/roodset_signal_%s_rv1",s1.Data()),dset_mcrcone_s_rv1);
      fmcrcone_s->GetObject(Form("roofit/roodset_signal_%s_rv2",s2.Data()),dset_mcrcone_s_rv2);
    }
    assert(dset_mcrcone_s_rv1);
    assert(dset_mcrcone_s_rv2);
    TFile *fmctrue_b = new TFile("outphoton_allmc_bkg.root","read");
    fmctrue_b->GetObject(Form("roofit/roodset_background_%s_rv1",s1.Data()),dset_mctrue_b_rv1);
    fmctrue_b->GetObject(Form("roofit/roodset_background_%s_rv2",s2.Data()),dset_mctrue_b_rv2);
    assert(dset_mctrue_b_rv1);
    assert(dset_mctrue_b_rv2);
    TFile *fmcrcone_b = new TFile("outphoton_allmc_sieiesideband.root","read");
    fmcrcone_b->GetObject(Form("roofit/roodset_background_%s_rv1",s1.Data()),dset_mcrcone_b_rv1);
    fmcrcone_b->GetObject(Form("roofit/roodset_background_%s_rv2",s2.Data()),dset_mcrcone_b_rv2);
    assert(dset_mcrcone_b_rv1);
    assert(dset_mcrcone_b_rv2);

    dset_mctrue_s_rv1 = (RooDataSet*)(dset_mctrue_s_rv1->reduce(Name("dset_mctrue_rv1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mcrcone_s_rv1 = (RooDataSet*)(dset_mcrcone_s_rv1->reduce(Name("dset_mcrcone_rv1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mctrue_s_rv2 = (RooDataSet*)(dset_mctrue_s_rv2->reduce(Name("dset_mctrue_rv2"),Cut(Form("roovar2<%f",rightrange-1e-5))));
    dset_mcrcone_s_rv2 = (RooDataSet*)(dset_mcrcone_s_rv2->reduce(Name("dset_mcrcone_rv2"),Cut(Form("roovar2<%f",rightrange-1e-5))));
    dset_mctrue_b_rv1 = (RooDataSet*)(dset_mctrue_b_rv1->reduce(Name("dset_mctrue_b_rv1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mcrcone_b_rv1 = (RooDataSet*)(dset_mcrcone_b_rv1->reduce(Name("dset_mcrcone_b_rv1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
    dset_mctrue_b_rv2 = (RooDataSet*)(dset_mctrue_b_rv2->reduce(Name("dset_mctrue_b_rv2"),Cut(Form("roovar2<%f",rightrange-1e-5))));
    dset_mcrcone_b_rv2 = (RooDataSet*)(dset_mcrcone_b_rv2->reduce(Name("dset_mcrcone_b_rv2"),Cut(Form("roovar2<%f",rightrange-1e-5))));

    std::cout << "MC true/rcone datasets" << std::endl;
    dset_mctrue_s_rv1->Print();
    dset_mcrcone_s_rv1->Print();
    dset_mctrue_s_rv2->Print();
    dset_mcrcone_s_rv2->Print();
    dset_mctrue_b_rv1->Print();
    dset_mcrcone_b_rv1->Print();
    dset_mctrue_b_rv2->Print();
    dset_mcrcone_b_rv2->Print();

    reweight_rhosigma(&dset_mctrue_s_rv1,dataset_axis1);
    reweight_rhosigma(&dset_mcrcone_s_rv1,dataset_axis1);
    reweight_eta_1d(&dset_mctrue_s_rv1,dataset_axis1,1);
    reweight_eta_1d(&dset_mcrcone_s_rv1,dataset_axis1,1);
    reweight_pt_1d(&dset_mctrue_s_rv1,dataset_axis1,1);
    //    reweight_pt_1d(&dset_mcrcone_s_rv1,dataset_axis1,1);
    reweight_rhosigma(&dset_mctrue_s_rv2,dataset_axis2);
    reweight_rhosigma(&dset_mcrcone_s_rv2,dataset_axis2);
    reweight_eta_1d(&dset_mctrue_s_rv2,dataset_axis2,2);
    reweight_eta_1d(&dset_mcrcone_s_rv2,dataset_axis2,2);
    reweight_pt_1d(&dset_mctrue_s_rv2,dataset_axis2,2);
    //    reweight_pt_1d(&dset_mcrcone_s_rv2,dataset_axis2,2);
    reweight_rhosigma(&dset_mctrue_b_rv1,dataset_axis1);
    reweight_rhosigma(&dset_mcrcone_b_rv1,dataset_axis1);
    reweight_eta_1d(&dset_mctrue_b_rv1,dataset_axis1,1);
    reweight_eta_1d(&dset_mcrcone_b_rv1,dataset_axis1,1);
    reweight_pt_1d(&dset_mctrue_b_rv1,dataset_axis1,1);
    reweight_pt_1d(&dset_mcrcone_b_rv1,dataset_axis1,1);
    reweight_rhosigma(&dset_mctrue_b_rv2,dataset_axis2);
    reweight_rhosigma(&dset_mcrcone_b_rv2,dataset_axis2);
    reweight_eta_1d(&dset_mctrue_b_rv2,dataset_axis2,2);
    reweight_eta_1d(&dset_mcrcone_b_rv2,dataset_axis2,2);
    reweight_pt_1d(&dset_mctrue_b_rv2,dataset_axis2,2);
    reweight_pt_1d(&dset_mcrcone_b_rv2,dataset_axis2,2);

  }

  bool islowstatcat = false;
  /*
  if (diffvariable=="costhetastar" && splitting=="EEEE" && bin==1) islowstatcat=true;
  if (diffvariable=="costhetastar" && splitting=="EEEE" && bin==5) islowstatcat=true;
  if (diffvariable=="costhetastar" && splitting=="EEEE" && bin==6) islowstatcat=true;
  if (diffvariable=="invmass" && bin>=13) islowstatcat=true;
  if (diffvariable=="diphotonpt" && splitting=="EEEE" && bin>=16) islowstatcat=true;
  if (diffvariable=="diphotonpt" && splitting=="EEEE" && bin<=4) islowstatcat=true;
  if (diffvariable=="dR" && splitting=="EBEE" && bin==6) islowstatcat=true;
  */

  find_adaptive_binning(dataset,&n_templatebins,templatebinsboundaries+0,1,islowstatcat ? -999 : -1);

  if (binning_roovar1==NULL || binning_roovar2==NULL){
  Double_t templatebinsboundaries_diagonal[n_templatebins_max+1];
  for (int i=0; i<n_templatebins+1; i++) templatebinsboundaries_diagonal[i]=templatebinsboundaries[i]*sqrt(2);
  binning_roovar1_threshold = new RooThresholdCategory("binning_roovar1_threshold","binning_roovar1_threshold",*roovar1);
  for (int i=1; i<n_templatebins+1; i++) binning_roovar1_threshold->addThreshold(templatebinsboundaries[i],Form("rv1_templatebin_thr_%d",i));
  binning_roovar2_threshold = new RooThresholdCategory("binning_roovar2_threshold","binning_roovar2_threshold",*roovar2);
  for (int i=1; i<n_templatebins+1; i++) binning_roovar2_threshold->addThreshold(templatebinsboundaries[i],Form("rv2_templatebin_thr_%d",i));
  binning_roovar1 = new RooRealVar("binning_roovar1","Binned Iso_{1}",0.5,n_templatebins+0.5); binning_roovar1->setBins(n_templatebins);
  binning_roovar2 = new RooRealVar("binning_roovar2","Binned Iso_{2}",0.5,n_templatebins+0.5); binning_roovar2->setBins(n_templatebins);
  }

  produce_category_binning(&dataset_sigsig);
  produce_category_binning(&dataset_sigbkg);
  produce_category_binning(&dataset_bkgsig);
  produce_category_binning(&dataset_bkgbkg);

  produce_category_binning(&dataset_sig_axis1);
  produce_category_binning(&dataset_bkg_axis1);
  produce_category_binning(&dataset_sig_axis2);
  produce_category_binning(&dataset_bkg_axis2);
  produce_category_binning(&dataset);
  produce_category_binning(&dataset_axis1);
  produce_category_binning(&dataset_axis2);

  if (do_syst_string==TString("savepdfMCtrue1D") || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")) {
    produce_category_binning(&dset_mctrue_s_rv1);
    produce_category_binning(&dset_mcrcone_s_rv1);
    produce_category_binning(&dset_mctrue_s_rv2);
    produce_category_binning(&dset_mcrcone_s_rv2);
    produce_category_binning(&dset_mctrue_b_rv1);
    produce_category_binning(&dset_mcrcone_b_rv1);
    produce_category_binning(&dset_mctrue_b_rv2);
    produce_category_binning(&dset_mcrcone_b_rv2);
  }




  int times_to_run = 1;
  const int ntoys = 10;

  std::vector<fit_output*> do_syst_templatestatistics_outputvector;
  std::vector<fit_output*> do_syst_purefitbias_outputvector;
  std::vector<fit_output*> do_syst_MCpromptdriven_outputvector;
  std::vector<fit_output*> do_syst_MCfakedriven_outputvector;
  std::vector<fit_output*> do_syst_2frag_outputvector;

  if (do_syst_string==TString("templatestatistics") || do_syst_string==TString("purefitbias") || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")) times_to_run = ntoys;

  RooAbsPdf *sigsigpdf_forgen = NULL;
  RooAbsPdf *sigbkgpdf_forgen = NULL;
  RooAbsPdf *bkgsigpdf_forgen = NULL;
  RooAbsPdf *bkgbkgpdf_forgen = NULL;
  RooAbsPdf *sigpdf_axis1_forgen = NULL;
  RooAbsPdf *bkgpdf_axis1_forgen = NULL;
  RooAbsPdf *sigpdf_axis2_forgen = NULL;
  RooAbsPdf *bkgpdf_axis2_forgen = NULL;
  fit_output *mctruthfr = NULL;
  fit_output *datafr = NULL;

  for (int runcount=0; runcount<times_to_run; runcount++){

      RooDataSet* original_dataset_sigsig   =NULL;
      RooDataSet* original_dataset_sigbkg   =NULL;
      RooDataSet* original_dataset_bkgsig   =NULL;
      RooDataSet* original_dataset_bkgbkg   =NULL;
      RooDataSet* original_dataset_sig_axis1=NULL;
      RooDataSet* original_dataset_bkg_axis1=NULL;
      RooDataSet* original_dataset_sig_axis2=NULL;
      RooDataSet* original_dataset_bkg_axis2=NULL;
      RooDataSet* original_dataset_axis1=NULL;
      RooDataSet* original_dataset_axis2=NULL;
      RooDataSet* original_dataset=NULL;


      if (do_syst_string==TString("templatestatistics") || (do_syst_string==TString("purefitbias") && runcount>0) || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")){
	original_dataset_sigsig   =dataset_sigsig   ;
        original_dataset_sigbkg   =dataset_sigbkg   ;
        original_dataset_bkgsig   =dataset_bkgsig   ;
        original_dataset_bkgbkg   =dataset_bkgbkg   ;
        original_dataset_sig_axis1=dataset_sig_axis1;
        original_dataset_bkg_axis1=dataset_bkg_axis1;
        original_dataset_sig_axis2=dataset_sig_axis2;
        original_dataset_bkg_axis2=dataset_bkg_axis2;
	dataset_sigsig   =(RooDataSet*)(original_dataset_sigsig->Clone("dataset_sigsig_forsyst"));
        dataset_sigbkg   =(RooDataSet*)(original_dataset_sigbkg->Clone("dataset_sigbkg_forsyst"));
        dataset_bkgsig   =(RooDataSet*)(original_dataset_bkgsig->Clone("dataset_bkgsig_forsyst"));
        dataset_bkgbkg   =(RooDataSet*)(original_dataset_bkgbkg->Clone("dataset_bkgbkg_forsyst"));
        dataset_sig_axis1=(RooDataSet*)(original_dataset_sig_axis1->Clone("dataset_sig_axis1_forsyst"));
        dataset_bkg_axis1=(RooDataSet*)(original_dataset_bkg_axis1->Clone("dataset_bkg_axis1_forsyst"));
        dataset_sig_axis2=(RooDataSet*)(original_dataset_sig_axis2->Clone("dataset_sig_axis2_forsyst"));
        dataset_bkg_axis2=(RooDataSet*)(original_dataset_bkg_axis2->Clone("dataset_bkg_axis2_forsyst"));
	original_dataset_axis1    =dataset_axis1;
	original_dataset_axis2    =dataset_axis2;
	original_dataset          =dataset;
	dataset_axis1    =(RooDataSet*)(original_dataset_axis1->Clone("dataset_axis1_forsyst"));
	dataset_axis2    =(RooDataSet*)(original_dataset_axis2->Clone("dataset_axis2_forsyst"));
	dataset          =(RooDataSet*)(original_dataset->Clone("dataset_forsyst"));
	}
	
      if (do_syst_string==TString("templatestatistics")){

	if (runcount==0) datafr = fit_dataset(diffvariable.Data(),splitting.Data(),bin,"data_donotwriteoutpurity");
	randomize_dataset_statistically_binned(&dataset_sigsig);
	randomize_dataset_statistically_binned(&dataset_sigbkg);
	randomize_dataset_statistically_binned(&dataset_bkgsig);
	randomize_dataset_statistically_binned(&dataset_bkgbkg);
	randomize_dataset_statistically_binned(&dataset_sig_axis1);
	randomize_dataset_statistically_binned(&dataset_bkg_axis1);
	randomize_dataset_statistically_binned(&dataset_sig_axis2);
	randomize_dataset_statistically_binned(&dataset_bkg_axis2);	

      }
	
      if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")) {
      
      if (runcount==0){
	datafr = fit_dataset(diffvariable.Data(),splitting.Data(),bin,"data_donotwriteoutpurity");
	mctruthfr = fit_dataset(diffvariable.Data(),splitting.Data(),bin,"savepdfMCtrue1D");
	sigsigpdf_forgen = mctruthfr->pdf_forgen[0];
	sigbkgpdf_forgen = mctruthfr->pdf_forgen[1];
	bkgsigpdf_forgen = mctruthfr->pdf_forgen[2];
	bkgbkgpdf_forgen = mctruthfr->pdf_forgen[3];
	sigpdf_axis1_forgen = mctruthfr->pdf_forgen[4];
	bkgpdf_axis1_forgen = mctruthfr->pdf_forgen[5];
	sigpdf_axis2_forgen = mctruthfr->pdf_forgen[6];
	bkgpdf_axis2_forgen = mctruthfr->pdf_forgen[7];
      }
      
      generate_toy_dataset_2d(&dataset,sigsigpdf_forgen,sigbkgpdf_forgen,bkgsigpdf_forgen,bkgbkgpdf_forgen,datafr->pp,datafr->pf,datafr->fp);

//      generate_toy_dataset_1d(&dataset_axis1,sigpdf_axis1_forgen,bkgpdf_axis1_forgen,mctruthfr->fsig1_firstpass);
//      generate_toy_dataset_1d(&dataset_axis2,sigpdf_axis2_forgen,bkgpdf_axis2_forgen,mctruthfr->fsig2_firstpass);
      delete dataset_axis1;
      delete dataset_axis2;
      dataset_axis1 = (RooDataSet*)(dataset->reduce(Name("dataset_axis1"),SelectVars(RooArgList(*binning_roovar1))));
      dataset_axis2 = (RooDataSet*)(dataset->reduce(Name("dataset_axis2"),SelectVars(RooArgList(*binning_roovar2))));
      }


    RooDataHist *sigsigdhist = new RooDataHist("sigsigdhist","sigsigdhist",RooArgList(*binning_roovar1,*binning_roovar2),*dataset_sigsig);
    RooDataHist *sigbkgdhist = new RooDataHist("sigbkgdhist","sigbkgdhist",RooArgList(*binning_roovar1,*binning_roovar2),*dataset_sigbkg);
    RooDataHist *bkgsigdhist = new RooDataHist("bkgsigdhist","bkgsigdhist",RooArgList(*binning_roovar1,*binning_roovar2),*dataset_bkgsig);
    RooDataHist *bkgbkgdhist = new RooDataHist("bkgbkgdhist","bkgbkgdhist",RooArgList(*binning_roovar1,*binning_roovar2),*dataset_bkgbkg);


    RooDataHist *dhist_mc_s_rv1 = NULL;
    RooDataHist *dhist_mc_s_rv2 = NULL;
    RooHistPdf *pdf_mc_s_rv1 = NULL;
    RooHistPdf *pdf_mc_s_rv2 = NULL;
    RooDataHist *dhist_mc_b_rv1 = NULL;
    RooDataHist *dhist_mc_b_rv2 = NULL;
    RooHistPdf *pdf_mc_b_rv1 = NULL;
    RooHistPdf *pdf_mc_b_rv2 = NULL;
    RooAbsPdf *sigsigpdf = NULL;
    RooAbsPdf *sigbkgpdf = NULL;
    RooAbsPdf *bkgsigpdf = NULL;
    RooAbsPdf *bkgbkgpdf = NULL;

    if (do_syst_string==TString("savepdfMCtrue1D") || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")){

//      This is for debug: if you uncomment this and comment the next block, you should get NO BIAS.
//      dhist_mc_s_rv1 = new RooDataHist("dhist_mc_s_rv1","dset_mc_s_rv1",RooArgList(*binning_roovar1),*dset_mctrue_s_rv1);
//      dhist_mc_s_rv2 = new RooDataHist("dhist_mc_s_rv2","dset_mc_s_rv2",RooArgList(*binning_roovar2),*dset_mctrue_s_rv2);
//      dhist_mc_b_rv1 = new RooDataHist("dhist_mc_b_rv1","dset_mc_b_rv1",RooArgList(*binning_roovar1),*dset_mctrue_b_rv1);
//      dhist_mc_b_rv2 = new RooDataHist("dhist_mc_b_rv2","dset_mc_b_rv2",RooArgList(*binning_roovar2),*dset_mctrue_b_rv2);
      {
	if (splitting=="EBEB"){        
	  dhist_mc_s_rv1 = new RooDataHist("dhist_mc_s_rv1","dset_mc_s_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEB")) ? *dset_mctrue_s_rv1 : *dset_mcrcone_s_rv1);
	  dhist_mc_s_rv2 = new RooDataHist("dhist_mc_s_rv2","dset_mc_s_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEB")) ? *dset_mctrue_s_rv2 : *dset_mcrcone_s_rv2);
	  dhist_mc_b_rv1 = new RooDataHist("dhist_mc_b_rv1","dset_mc_b_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshapeMCfakedrivenEB")) ? *dset_mctrue_b_rv1 : *dset_mcrcone_b_rv1);
	  dhist_mc_b_rv2 = new RooDataHist("dhist_mc_b_rv2","dset_mc_b_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshapeMCfakedrivenEB")) ? *dset_mctrue_b_rv2 : *dset_mcrcone_b_rv2);
	}
	else if (splitting=="EBEE"){        
	  dhist_mc_s_rv1 = new RooDataHist("dhist_mc_s_rv1","dset_mc_s_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEB")) ? *dset_mctrue_s_rv1 : *dset_mcrcone_s_rv1);
	  dhist_mc_s_rv2 = new RooDataHist("dhist_mc_s_rv2","dset_mc_s_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEE")) ? *dset_mctrue_s_rv2 : *dset_mcrcone_s_rv2);
	  dhist_mc_b_rv1 = new RooDataHist("dhist_mc_b_rv1","dset_mc_b_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshapeMCfakedrivenEB")) ? *dset_mctrue_b_rv1 : *dset_mcrcone_b_rv1);
	  dhist_mc_b_rv2 = new RooDataHist("dhist_mc_b_rv2","dset_mc_b_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshapeMCfakedrivenEE")) ? *dset_mctrue_b_rv2 : *dset_mcrcone_b_rv2);
	}
	else if (splitting=="EEEE"){        
	  dhist_mc_s_rv1 = new RooDataHist("dhist_mc_s_rv1","dset_mc_s_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEE")) ? *dset_mctrue_s_rv1 : *dset_mcrcone_s_rv1);
	  dhist_mc_s_rv2 = new RooDataHist("dhist_mc_s_rv2","dset_mc_s_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshape2frag") && do_syst_string!=TString("templateshapeMCpromptdrivenEE")) ? *dset_mctrue_s_rv2 : *dset_mcrcone_s_rv2);
	  dhist_mc_b_rv1 = new RooDataHist("dhist_mc_b_rv1","dset_mc_b_rv1",RooArgList(*binning_roovar1),(do_syst_string!=TString("templateshapeMCfakedrivenEE")) ? *dset_mctrue_b_rv1 : *dset_mcrcone_b_rv1);
	  dhist_mc_b_rv2 = new RooDataHist("dhist_mc_b_rv2","dset_mc_b_rv2",RooArgList(*binning_roovar2),(do_syst_string!=TString("templateshapeMCfakedrivenEE")) ? *dset_mctrue_b_rv2 : *dset_mcrcone_b_rv2);
	}
      }

      pdf_mc_s_rv1 = new RooHistPdf("pdf_mc_s_rv1","pdf_mc_s_rv1",RooArgList(*binning_roovar1),*dhist_mc_s_rv1);
      pdf_mc_s_rv2 = new RooHistPdf("pdf_mc_s_rv2","pdf_mc_s_rv2",RooArgList(*binning_roovar2),*dhist_mc_s_rv2);
      pdf_mc_b_rv1 = new RooHistPdf("pdf_mc_b_rv1","pdf_mc_b_rv1",RooArgList(*binning_roovar1),*dhist_mc_b_rv1);
      pdf_mc_b_rv2 = new RooHistPdf("pdf_mc_b_rv2","pdf_mc_b_rv2",RooArgList(*binning_roovar2),*dhist_mc_b_rv2);
      sigsigpdf = (RooAbsPdf*)(new RooProdPdf("sigsigpdf","sigsigpdf",*pdf_mc_s_rv1,*pdf_mc_s_rv2));
      sigbkgpdf = (RooAbsPdf*)(new RooProdPdf("sigbkgpdf","sigbkgpdf",*pdf_mc_s_rv1,*pdf_mc_b_rv2));
      bkgsigpdf = (RooAbsPdf*)(new RooProdPdf("bkgsigpdf","bkgsigpdf",*pdf_mc_b_rv1,*pdf_mc_s_rv2));
      bkgbkgpdf = (RooAbsPdf*)(new RooProdPdf("bkgbkgpdf","bkgbkgpdf",*pdf_mc_b_rv1,*pdf_mc_b_rv2));
    }
    else {
      sigsigpdf = (RooAbsPdf*)(new RooHistPdf("sigsigpdf","sigsigpdf",RooArgList(*binning_roovar1,*binning_roovar2),*sigsigdhist));
      sigbkgpdf = (RooAbsPdf*)(new RooHistPdf("sigbkgpdf","sigbkgpdf",RooArgList(*binning_roovar1,*binning_roovar2),*sigbkgdhist));
      bkgsigpdf = (RooAbsPdf*)(new RooHistPdf("bkgsigpdf","bkgsigpdf",RooArgList(*binning_roovar1,*binning_roovar2),*bkgsigdhist));
      bkgbkgpdf = (RooAbsPdf*)(new RooHistPdf("bkgbkgpdf","bkgbkgpdf",RooArgList(*binning_roovar1,*binning_roovar2),*bkgbkgdhist));
    }

    RooDataHist *sigdhist_axis1 = new RooDataHist("sigdhist_axis1","sigdhist_axis1",RooArgList(*binning_roovar1),*dataset_sig_axis1);
    RooDataHist *bkgdhist_axis1 = new RooDataHist("bkgdhist_axis1","bkgdhist_axis1",RooArgList(*binning_roovar1),*dataset_bkg_axis1);
    RooHistPdf *sigpdf_axis1 = new RooHistPdf("sigpdf_axis1","sigpdf_axis1",RooArgList(*binning_roovar1),*sigdhist_axis1);
    RooHistPdf *bkgpdf_axis1 = new RooHistPdf("bkgpdf_axis1","bkgpdf_axis1",RooArgList(*binning_roovar1),*bkgdhist_axis1);
    RooDataHist *sigdhist_axis2 = new RooDataHist("sigdhist_axis2","sigdhist_axis2",RooArgList(*binning_roovar2),*dataset_sig_axis2);
    RooDataHist *bkgdhist_axis2 = new RooDataHist("bkgdhist_axis2","bkgdhist_axis2",RooArgList(*binning_roovar2),*dataset_bkg_axis2);
    RooHistPdf *sigpdf_axis2 = new RooHistPdf("sigpdf_axis2","sigpdf_axis2",RooArgList(*binning_roovar2),*sigdhist_axis2);
    RooHistPdf *bkgpdf_axis2 = new RooHistPdf("bkgpdf_axis2","bkgpdf_axis2",RooArgList(*binning_roovar2),*bkgdhist_axis2);
  
    RooDataHist *sigsigdhist_unbinned = new RooDataHist("sigsigdhist_unbinned","sigsigdhist_unbinned",RooArgList(*roovar1,*roovar2),*dataset_sigsig);
    RooDataHist *sigbkgdhist_unbinned = new RooDataHist("sigbkgdhist_unbinned","sigbkgdhist_unbinned",RooArgList(*roovar1,*roovar2),*dataset_sigbkg);
    RooDataHist *bkgsigdhist_unbinned = new RooDataHist("bkgsigdhist_unbinned","bkgsigdhist_unbinned",RooArgList(*roovar1,*roovar2),*dataset_bkgsig);
    RooDataHist *bkgbkgdhist_unbinned = new RooDataHist("bkgbkgdhist_unbinned","bkgbkgdhist_unbinned",RooArgList(*roovar1,*roovar2),*dataset_bkgbkg);
    RooHistPdf *sigsigpdf_unbinned = new RooHistPdf("sigsigpdf_unbinned","sigsigpdf_unbinned",RooArgList(*roovar1,*roovar2),*sigsigdhist_unbinned);
    RooHistPdf *sigbkgpdf_unbinned = new RooHistPdf("sigbkgpdf_unbinned","sigbkgpdf_unbinned",RooArgList(*roovar1,*roovar2),*sigbkgdhist_unbinned);
    RooHistPdf *bkgsigpdf_unbinned = new RooHistPdf("bkgsigpdf_unbinned","bkgsigpdf_unbinned",RooArgList(*roovar1,*roovar2),*bkgsigdhist_unbinned);
    RooHistPdf *bkgbkgpdf_unbinned = new RooHistPdf("bkgbkgpdf_unbinned","bkgbkgpdf_unbinned",RooArgList(*roovar1,*roovar2),*bkgbkgdhist_unbinned);
    RooDataHist *sigdhist_axis1_unbinned = new RooDataHist("sigdhist_axis1_unbinned","sigdhist_axis1_unbinned",RooArgList(*roovar1),*dataset_sig_axis1);
    RooDataHist *bkgdhist_axis1_unbinned = new RooDataHist("bkgdhist_axis1_unbinned","bkgdhist_axis1_unbinned",RooArgList(*roovar1),*dataset_bkg_axis1);
    RooHistPdf *sigpdf_axis1_unbinned = new RooHistPdf("sigpdf_axis1_unbinned","sigpdf_axis1_unbinned",RooArgList(*roovar1),*sigdhist_axis1_unbinned);
    RooHistPdf *bkgpdf_axis1_unbinned = new RooHistPdf("bkgpdf_axis1_unbinned","bkgpdf_axis1_unbinned",RooArgList(*roovar1),*bkgdhist_axis1_unbinned);
    RooDataHist *sigdhist_axis2_unbinned = new RooDataHist("sigdhist_axis2_unbinned","sigdhist_axis2_unbinned",RooArgList(*roovar2),*dataset_sig_axis2);
    RooDataHist *bkgdhist_axis2_unbinned = new RooDataHist("bkgdhist_axis2_unbinned","bkgdhist_axis2_unbinned",RooArgList(*roovar2),*dataset_bkg_axis2);
    RooHistPdf *sigpdf_axis2_unbinned = new RooHistPdf("sigpdf_axis2_unbinned","sigpdf_axis2_unbinned",RooArgList(*roovar2),*sigdhist_axis2_unbinned);
    RooHistPdf *bkgpdf_axis2_unbinned = new RooHistPdf("bkgpdf_axis2_unbinned","bkgpdf_axis2_unbinned",RooArgList(*roovar2),*bkgdhist_axis2_unbinned);


    if (do_syst_string==TString("purefitbias") && (runcount>0)) {
      //	mctruthfr = fit_dataset(diffvariable.Data(),splitting.Data(),bin,"savepdfMCtrue2D");
      //      datafr = fit_dataset(diffvariable.Data(),splitting.Data(),bin,"data_donotwriteoutpurity");
      assert(datafr!=NULL);
      generate_toy_dataset_2d(&dataset,sigsigpdf,sigbkgpdf,bkgsigpdf,bkgbkgpdf,datafr->pp,datafr->pf,datafr->fp);
      delete dataset_axis1;
      delete dataset_axis2;
      dataset_axis1 = (RooDataSet*)(dataset->reduce(Name("dataset_axis1"),SelectVars(RooArgList(*binning_roovar1))));
      dataset_axis2 = (RooDataSet*)(dataset->reduce(Name("dataset_axis2"),SelectVars(RooArgList(*binning_roovar2))));
    }
	

    
    if (doplots_b) {
      TCanvas *c0 = new TCanvas(Form("c0"),Form("c0"),1200,800);
      c0->Divide(2,2);

      c0->cd(1);
      RooPlot *frame01sig = binning_roovar1->frame(Title("Signal template axis 1 - binned"));
      dataset_sigsig->plotOn(frame01sig,Name("proj"));
      //      sigsigpdf->plotOn(frame01sig);
      sigpdf_axis1->plotOn(frame01sig,LineColor(kRed),LineStyle(kDashed),Name("1d"));
      frame01sig->Draw();
      //    c0->GetPad(1)->SetLogy(1);
      TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
      leg->AddEntry("proj","proj. of 2-D template","lp");
      leg->AddEntry("1d","1-D template","l");
      leg->SetFillColor(kWhite);
      leg->Draw();

      c0->cd(2);
      RooPlot *frame02sig = binning_roovar2->frame(Title("Signal template axis 2 - binned"));
      dataset_sigsig->plotOn(frame02sig);
      //      sigsigpdf->plotOn(frame02sig);
      sigpdf_axis2->plotOn(frame02sig,LineColor(kRed),LineStyle(kDashed));
      frame02sig->Draw();
      //    c0->GetPad(2)->SetLogy(1);
      leg->Draw();

      c0->cd(3);
      RooPlot *frame01bkg = binning_roovar1->frame(Title("Background template axis 1 - binned"));
      dataset_bkgbkg->plotOn(frame01bkg,Name("proj"));
      //      bkgbkgpdf->plotOn(frame01bkg);
      bkgpdf_axis1->plotOn(frame01bkg,LineColor(kRed),LineStyle(kDashed),Name("1d"));
      frame01bkg->Draw();
      //    c0->GetPad(1)->SetLogy(1);
      TLegend *leg2 = new TLegend(0.3,0.7,0.5,0.9);
      leg2->AddEntry("proj","proj. of 2-D template","lp");
      leg2->AddEntry("1d","1-D template","l");
      leg2->SetFillColor(kWhite);
      leg2->Draw();
      leg2->Draw();

      c0->cd(4);
      RooPlot *frame02bkg = binning_roovar2->frame(Title("Background template axis 2 - binned"));
      dataset_bkgbkg->plotOn(frame02bkg);
      //      bkgbkgpdf->plotOn(frame02bkg);
      bkgpdf_axis2->plotOn(frame02bkg,LineColor(kRed),LineStyle(kDashed));
      frame02bkg->Draw();
      //    c0->GetPad(2)->SetLogy(1);
      leg2->Draw();

      c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
      c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));
      c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));
      c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));

    } // c0 binned

    if (doplots_ub) {
      TCanvas *c0_ub = new TCanvas(Form("c0_ub"),Form("c0_ub"),1200,800);
      c0_ub->Divide(2,2);

      
      c0_ub->cd(1);
      RooPlot *frame01sig = roovar1->frame(Title("Signal template axis 1"));
      dataset_sigsig->plotOn(frame01sig,Name("proj"));
      //      sigsigpdf_unbinned->plotOn(frame01sig);
      sigpdf_axis1_unbinned->plotOn(frame01sig,LineColor(kRed),LineStyle(kDashed),Name("1d"));
      frame01sig->Draw();
      //    c0_ub->GetPad(1)->SetLogy(1);
      TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
      leg->AddEntry("proj","proj. of 2-D template","lp");
      leg->AddEntry("1d","1-D template","l");
      leg->SetFillColor(kWhite);
      leg->Draw();


      c0_ub->cd(2);
      RooPlot *frame02sig = roovar2->frame(Title("Signal template axis 2"));
      dataset_sigsig->plotOn(frame02sig);
      //      sigsigpdf_unbinned->plotOn(frame02sig);
      sigpdf_axis2_unbinned->plotOn(frame02sig,LineColor(kRed),LineStyle(kDashed));
      frame02sig->Draw();
      //    c0_ub->GetPad(2)->SetLogy(1);
      leg->Draw();

      c0_ub->cd(3);
      RooPlot *frame01bkg = roovar1->frame(Title("Background template axis 1"));
      dataset_bkgbkg->plotOn(frame01bkg);
      //      bkgbkgpdf_unbinned->plotOn(frame01bkg);
      bkgpdf_axis1_unbinned->plotOn(frame01bkg,LineColor(kRed),LineStyle(kDashed));
      frame01bkg->Draw();
      //    c0_ub->GetPad(1)->SetLogy(1);
      leg->Draw();

      c0_ub->cd(4);
      RooPlot *frame02bkg = roovar2->frame(Title("Background template axis 2"));
      dataset_bkgbkg->plotOn(frame02bkg);
      //      bkgbkgpdf_unbinned->plotOn(frame02bkg);
      bkgpdf_axis2_unbinned->plotOn(frame02bkg,LineColor(kRed),LineStyle(kDashed));
      frame02bkg->Draw();
      //    c0_ub->GetPad(2)->SetLogy(1);
      leg->Draw();

      c0_ub->SaveAs(Form("plots/fittingplot0unbinned_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
      c0_ub->SaveAs(Form("plots/fittingplot0unbinned_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));
      c0_ub->SaveAs(Form("plots/fittingplot0unbinned_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));
      c0_ub->SaveAs(Form("plots/fittingplot0unbinned_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));


    } // c0_ub

    if (doplots_b) {
      //      plot_template_dependency_axis1(dataset_bkg_axis1,TString("pt"),20,70,2,kTRUE);
      //      plot_template_dependency_axis1(dataset_bkg_axis1,TString("sieie"),0.008,0.014,6,kTRUE);
    } // template dependency

    // inversione:
    //  pp = pp;
    //  pf = j1 - pp;
    //  fp = j2 - pp;
    //  ff = (1-j1-j2)+pp
   
    out = new fit_output();
    out->fr_pass1=NULL;
    out->fr_pass2constraint=NULL;
    out->fr=NULL;
    out->tot_events=0;
    out->pp=0;
    out->pp_err=0;
    out->pf=0;
    out->pf_err=0;
    out->fp=0;
    out->fp_err=0;
    out->ff=0;
    out->ff_err=0;
    out->eff_overflow_removal_pp=eff_overflow_removal;
    out->chi2=0;
    out->ndof=0;
    out->probchi2=0;
    for (int l=0; l<4; l++) out->pdf_forgen[l]=NULL;

    RooRealVar *pp = new RooRealVar("pp","pp",pp_init,0,1);
    RooRealVar *j1 = new RooRealVar("j1","j1",pp_init+pf_init,0,1);
    RooRealVar *j2 = new RooRealVar("j2","j2",pp_init+fp_init,0,1);

    RooFormulaVar *fsig1 = new RooFormulaVar("fsig1","fsig1","j1",RooArgList(*j1));
    RooFormulaVar *fsig2 = new RooFormulaVar("fsig2","fsig2","@0", (sym) ? RooArgList(*j1) : RooArgList(*j2) );

    sigpdf_axis1->Print();
    bkgpdf_axis1->Print();
    sigpdf_axis2->Print();
    bkgpdf_axis2->Print();
    sigpdf_axis1_unbinned->Print();
    bkgpdf_axis1_unbinned->Print();
    sigpdf_axis2_unbinned->Print();
    bkgpdf_axis2_unbinned->Print();
    

    RooAddPdf *model_axis1 = new RooAddPdf("model_axis1","model_axis1",RooArgList(*sigpdf_axis1,*bkgpdf_axis1),RooArgList(*fsig1));
    RooAddPdf *model_axis2 = new RooAddPdf("model_axis2","model_axis2",RooArgList(*sigpdf_axis2,*bkgpdf_axis2),RooArgList(*fsig2));
    RooAddPdf *model_axis1_unbinned = new RooAddPdf("model_axis1_unbinned","model_axis1_unbinned",RooArgList(*sigpdf_axis1_unbinned,*bkgpdf_axis1_unbinned),RooArgList(*fsig1));
    RooAddPdf *model_axis2_unbinned = new RooAddPdf("model_axis2_unbinned","model_axis2_unbinned",RooArgList(*sigpdf_axis2_unbinned,*bkgpdf_axis2_unbinned),RooArgList(*fsig2));

    RooFitResult *firstpass;

    RooNLLVar *model_axis1_noextended_nll = new RooNLLVar("model_axis1_noextended_nll","model_axis1_noextended_nll",*model_axis1,*dataset_axis1,NumCPU(numcpu>1 ? numcpu/2 : 1));
    RooNLLVar *model_axis1_nll = model_axis1_noextended_nll;
    RooNLLVar *model_axis2_noextended_nll = new RooNLLVar("model_axis2_noextended_nll","model_axis2_noextended_nll",*model_axis2,*dataset_axis2,NumCPU(numcpu>1 ? numcpu/2 : 1));
    RooNLLVar *model_axis2_nll = model_axis2_noextended_nll;
    RooAddition *model_2axes_nll = new RooAddition("model_2axes_nll","model_2axes_nll",RooArgSet(*model_axis1_nll,*model_axis2_nll));

    RooMinimizer *minuit_firstpass = new RooMinimizer(*model_2axes_nll);
    minuit_firstpass->migrad();
    minuit_firstpass->hesse();
    firstpass = minuit_firstpass->save("firstpass","firstpass");
    firstpass->Print();


    /*
      ofstream myfile;
      myfile.open(Form("plots/fitresults_%d.txt",bin));
      myfile << Form("bin %d",bin) << std::endl;
      myfile << "fsig1 " << fsig1->getVal() << " " << fsig1->getPropagatedError(*firstpass) << std::endl; 
      if (!sym) myfile << "fsig2 " << fsig2->getVal() << " " << fsig2->getPropagatedError(*firstpass) << std::endl; 
    */
  

    if (doplots_b) {
      TCanvas *c1 = new TCanvas("c1","c1",1200,800);
      c1->Divide(2);
      c1->cd(1);
      RooPlot *frame1bla = binning_roovar1->frame(Title("Fit axis 1 - first pass - binned"));
      dataset_axis1->plotOn(frame1bla,Name("data"));
      model_axis1->plotOn(frame1bla,Name("fit"));
      model_axis1->plotOn(frame1bla,Components("sigpdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
      model_axis1->plotOn(frame1bla,Components("bkgpdf_axis1"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
      frame1bla->Draw();
      //    c1->GetPad(1)->SetLogy(1);
      TLegend *leg = new TLegend(0.2,0.8,0.45,0.9);
      leg->AddEntry("data","data","lp");
      leg->AddEntry("fit","fit","l");
      leg->AddEntry("signal","signal comp.","l");
      leg->AddEntry("background","background comp.","l");
      leg->SetFillColor(kWhite);
      leg->Draw();

      c1->cd(2);
      RooPlot *frame2bla = binning_roovar2->frame(Title("Fit axis 2 - first pass - binned"));
      dataset_axis2->plotOn(frame2bla,Name("data"));
      model_axis2->plotOn(frame2bla,Name("fit"));
      model_axis2->plotOn(frame2bla,Components("sigpdf_axis2"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
      model_axis2->plotOn(frame2bla,Components("bkgpdf_axis2"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
      frame2bla->Draw();
      //    c1->GetPad(2)->SetLogy(1);
      leg->Draw();

      model_axis1->Print();
      model_axis2->Print();
      dataset_axis1->Print();
      dataset_axis2->Print();

      c1->SaveAs(Form("plots/fittingplot1_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
      c1->SaveAs(Form("plots/fittingplot1_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));
      c1->SaveAs(Form("plots/fittingplot1_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));
    } // c1

    if (doplots_ub) {
      TCanvas *c1_ub = new TCanvas("c1_ub","c1_ub",1200,800);
      c1_ub->Divide(2);
      c1_ub->cd(1);
      RooPlot *frame1bla = roovar1->frame(Title(""));
      frame1bla->GetXaxis()->SetLabelSize(0.03);
      frame1bla->GetYaxis()->SetLabelSize(0.03);
      frame1bla->GetYaxis()->SetTitleOffset(1.45);
      dataset_axis1->plotOn(frame1bla,Name("data"));
      model_axis1_unbinned->plotOn(frame1bla,Name("fit"));
      model_axis1_unbinned->plotOn(frame1bla,Components("sigpdf_axis1_unbinned"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
      model_axis1_unbinned->plotOn(frame1bla,Components("bkgpdf_axis1_unbinned"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
      frame1bla->Draw();
      //    c1_ub->GetPad(1)->SetLogy(1);
      TLegend *leg = new TLegend(0.55,0.7,0.9,0.9);
      leg->AddEntry("data","data","lp");
      leg->AddEntry("fit","fit","l");
      leg->AddEntry("signal","signal comp.","l");
      leg->AddEntry("background","background comp.","l");
      leg->SetFillColor(kWhite);
      leg->Draw();
      TLatex a;
      a.SetNDC();
      a.SetTextSize(0.03);
      a.DrawLatex(0.58,0.65,"#splitline{CMS Preliminary}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");

      c1_ub->cd(2);
      RooPlot *frame2bla = roovar2->frame(Title(""));
      frame2bla->GetXaxis()->SetLabelSize(0.03);
      frame2bla->GetYaxis()->SetLabelSize(0.03);
      frame2bla->GetYaxis()->SetTitleOffset(1.45);
      dataset_axis2->plotOn(frame2bla,Name("data"));
      model_axis2_unbinned->plotOn(frame2bla,Name("fit"));
      model_axis2_unbinned->plotOn(frame2bla,Components("sigpdf_axis2_unbinned"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
      model_axis2_unbinned->plotOn(frame2bla,Components("bkgpdf_axis2_unbinned"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
      frame2bla->Draw();
      //    c1_ub->GetPad(2)->SetLogy(1);
      leg->Draw();
      a.DrawLatex(0.58,0.65,"#splitline{CMS Preliminary}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");

      model_axis1_unbinned->Print();
      model_axis2_unbinned->Print();
      dataset_axis1->Print();
      dataset_axis2->Print();

      c1_ub->SaveAs(Form("plots/fittingplot1unbinned_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
      c1_ub->SaveAs(Form("plots/fittingplot1unbinned_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));
      c1_ub->SaveAs(Form("plots/fittingplot1unbinned_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));
    } // c1_ub
  
    return NULL;  
  
    RooFormulaVar *fsigsig = new RooFormulaVar("fsigsig","fsigsig","pp",RooArgList(*pp));
    RooFormulaVar *fsigbkg = new RooFormulaVar("fsigbkg","fsigbkg","fsig1-pp",RooArgList(*fsig1,*pp));  
    RooFormulaVar *fbkgsig = new RooFormulaVar("fbkgsig","fbkgsig","fsig2-pp",RooArgList(*fsig2,*pp));
    RooFormulaVar *fbkgbkg = new RooFormulaVar("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg-fbkgsig",RooArgList(*fsigsig,*fsigbkg,*fbkgsig));

  
    const float nsigma_tolerance = 0;
      
    float f1p = fsig1->getVal()+nsigma_tolerance*fsig1->getPropagatedError(*firstpass);
    float f2p = fsig2->getVal()+nsigma_tolerance*fsig2->getPropagatedError(*firstpass);
    float f1l = fsig1->getVal()-nsigma_tolerance*fsig1->getPropagatedError(*firstpass);
    float f2l = fsig2->getVal()-nsigma_tolerance*fsig2->getPropagatedError(*firstpass);

    float lowerbounds[4]={0,f1p-1,f2p-1,f1p+f2p-1};
    float upperbounds[4]={1,f1l,f2l,f1l+f2l};      

    float minpp = TMath::MaxElement(4,lowerbounds);
    float maxpp = TMath::MinElement(4,upperbounds);
    pp->setVal((minpp+maxpp)/2);
      
//    std::cout << "setting constrain pp val at " << pp->getVal() << " between " << minpp << " and " << maxpp << std::endl;
//    pp->setRange(minpp,maxpp);



    RooGaussian *constraint_gaussian_j1 = new RooGaussian("constraint_gaussian_j1","constraint_gaussian_j1",*j1,RooRealConstant::value(j1->getVal()),RooRealConstant::value(j1->getPropagatedError(*firstpass)));
    RooArgSet *constraint_pdf_set = new RooArgSet(*constraint_gaussian_j1);
    RooArgSet *constraint_parameters_set = new RooArgSet(*j1);
    RooGaussian *constraint_gaussian_j2=NULL;
    if (!sym){
      constraint_gaussian_j2 = new RooGaussian("constraint_gaussian_j2","constraint_gaussian_j2",*j2,RooRealConstant::value(j2->getVal()),RooRealConstant::value(j2->getPropagatedError(*firstpass)));
      constraint_pdf_set->add(*constraint_gaussian_j2);
      constraint_parameters_set->add(*j2);
    }
    RooConstraintSum *constraint_gaussian_nll = new RooConstraintSum("constraint_gaussian_nll","constraint_gaussian_nll",*constraint_pdf_set,*constraint_parameters_set);


    RooAddPdf *model_2D_uncorrelated = new RooAddPdf("model_2D_uncorrelated","model_2D_uncorrelated",RooArgList(*sigsigpdf,*sigbkgpdf,*bkgsigpdf,*bkgbkgpdf),RooArgList(*fsigsig,*fsigbkg,*fbkgsig,*fbkgbkg),kFALSE);
    RooAddPdf *model_2D_uncorrelated_unbinned = new RooAddPdf("model_2D_uncorrelated_unbinned","model_2D_uncorrelated_unbinned",RooArgList(*sigsigpdf_unbinned,*sigbkgpdf_unbinned,*bkgsigpdf_unbinned,*bkgbkgpdf_unbinned),RooArgList(*fsigsig,*fsigbkg,*fbkgsig,*fbkgbkg),kFALSE);
    RooNLLVar *model_2D_uncorrelated_noextended_nll = new RooNLLVar("model_2D_uncorrelated_noextended_nll","model_2D_uncorrelated_noextended_nll",*model_2D_uncorrelated,*dataset,NumCPU(numcpu));
    RooAddition *model_2D_uncorrelated_noextended_nll_constraint = new RooAddition("model_2D_uncorrelated_noextended_nll_constraint","model_2D_uncorrelated_noextended_nll_constraint",RooArgSet(*model_2D_uncorrelated_noextended_nll,*constraint_gaussian_nll));


    RooMinimizer *minuit_secondpass_constraint = new RooMinimizer(*model_2D_uncorrelated_noextended_nll_constraint);
    minuit_secondpass_constraint->migrad();
    minuit_secondpass_constraint->hesse();
    RooFitResult *secondpass_constraint;
    secondpass_constraint = minuit_secondpass_constraint->save("secondpass_constraint","secondpass_constraint");
    secondpass_constraint->Print();    

    /* 100% free per il secondpass finale
      std::cout << "setting constrain j1 val at " << j1->getVal() << " between " << j1->getVal()-nsigma_tolerance*j1->getPropagatedError(*firstpass) << " and " << j1->getVal()+nsigma_tolerance*j1->getPropagatedError(*firstpass) << std::endl;
      j1->setRange(j1->getVal()-nsigma_tolerance*j1->getPropagatedError(*firstpass),j1->getVal()+nsigma_tolerance*j1->getPropagatedError(*firstpass));
      
      if (!sym){
      std::cout << "setting constrain j2 val at " << j2->getVal() << " between " << j2->getVal()-nsigma_tolerance*j2->getPropagatedError(*firstpass) << " and " << j2->getVal()+nsigma_tolerance*j2->getPropagatedError(*firstpass) << std::endl;
      j2->setRange(j2->getVal()-nsigma_tolerance*j2->getPropagatedError(*firstpass),j2->getVal()+nsigma_tolerance*j2->getPropagatedError(*firstpass));
      }
    */

    RooMinimizer *minuit_secondpass = new RooMinimizer(*model_2D_uncorrelated_noextended_nll);
    minuit_secondpass->migrad();
    minuit_secondpass->hesse();
    RooFitResult *secondpass;
    secondpass = minuit_secondpass->save("secondpass","secondpass");
    secondpass->Print();

    /*
      myfile << "pp " << fsigsig->getVal() << " " << fsigsig->getPropagatedError(*secondpass) << std::endl; 
      myfile << "pf " << fsigbkg->getVal() << " " << fsigbkg->getPropagatedError(*secondpass) << std::endl; 
      myfile << "fp " << fbkgsig->getVal() << " " << fbkgsig->getPropagatedError(*secondpass) << std::endl; 
      myfile << "ff " << fbkgbkg->getVal() << " " << fbkgbkg->getPropagatedError(*secondpass) << std::endl; 
    */

    out->fr_pass1=firstpass;
    out->fr_pass2constraint=secondpass_constraint;
    out->fr=secondpass;
    out->tot_events=dataset->sumEntries();
    out->pp=fsigsig->getVal();
    out->pp_err=fsigsig->getPropagatedError(*secondpass);
    out->pf=fsigbkg->getVal();
    out->pf_err=fsigbkg->getPropagatedError(*secondpass);
    out->fp=fbkgsig->getVal();
    out->fp_err=fbkgsig->getPropagatedError(*secondpass);
    out->ff=fbkgbkg->getVal();
    out->ff_err=fbkgbkg->getPropagatedError(*secondpass);
    out->fsig1_firstpass=fsig1->getVal();
    out->fsig2_firstpass=fsig2->getVal();

    std::cout << "RAW YIELD " << out->pp*out->tot_events/out->eff_overflow_removal_pp << std::endl;

    {
      RooDataHist *dataset_dhist = new RooDataHist("dataset_dhist","dataset_dhist",RooArgList(*binning_roovar1,*binning_roovar2),*dataset);
      RooAbsReal *chi2var = model_2D_uncorrelated->createChi2(*dataset_dhist);
      float chi2 = chi2var->getVal();
      int nparams = (sym) ? 2 : 3;
      int ndof = n_templatebins*n_templatebins-nparams-1;
      std::cout << "OFFICIALCHI2 " << chi2 << " " << ndof << " " << chi2/ndof << " " << TMath::Prob(chi2,ndof) << std::endl; 

      float oldchi2 = 0;
      float newchi2 = 0;

      ndof=0;

      for (int i=0; i<dataset_dhist->numEntries(); i++){
        float r1 = dataset_dhist->get(i)->getRealValue("binning_roovar1");
        float r2 = dataset_dhist->get(i)->getRealValue("binning_roovar2");
        float w = dataset_dhist->store()->weight(i);
	float data_werr = dataset_dhist->weightError(RooAbsData::SumW2);

	if (w==0) std::cout << "WARNING: EMPTY DATA BIN " << r1 << " " << r2 << std::endl;
	if (w<10) std::cout << "WARNING: LOW STAT DATA BIN " << r1 << " " << r2 << " events:" << w << std::endl;
	//	if (data_werr==0) std::cout << "WARNING: ZERO UNCERTAINTY ON DATA " << r1 << " " << r2 << std::endl;

	binning_roovar1->setVal(r1);
	binning_roovar2->setVal(r2);
	float fitvalue = dataset_dhist->sumEntries()*model_2D_uncorrelated->getVal(RooArgSet(*binning_roovar1,*binning_roovar2));

	assert (sigsigdhist->get(i)->getRealValue("binning_roovar1")==r1);
	assert (sigsigdhist->get(i)->getRealValue("binning_roovar2")==r2);
	sigsigdhist->store()->get(i);
	float sigma_sigsig = sigsigdhist->store()->weightError(RooAbsData::SumW2)/sigsigdhist->sumEntries()*dataset_dhist->sumEntries()*fsigsig->getVal();
	assert (sigbkgdhist->get(i)->getRealValue("binning_roovar1")==r1);
	assert (sigbkgdhist->get(i)->getRealValue("binning_roovar2")==r2);
	sigbkgdhist->store()->get(i);
	float sigma_sigbkg = sigbkgdhist->store()->weightError(RooAbsData::SumW2)/sigbkgdhist->sumEntries()*dataset_dhist->sumEntries()*fsigbkg->getVal();
	assert (bkgsigdhist->get(i)->getRealValue("binning_roovar1")==r1);
	assert (bkgsigdhist->get(i)->getRealValue("binning_roovar2")==r2);
	bkgsigdhist->store()->get(i);
	float sigma_bkgsig = bkgsigdhist->store()->weightError(RooAbsData::SumW2)/bkgsigdhist->sumEntries()*dataset_dhist->sumEntries()*fbkgsig->getVal();
	assert (bkgbkgdhist->get(i)->getRealValue("binning_roovar1")==r1);
	assert (bkgbkgdhist->get(i)->getRealValue("binning_roovar2")==r2);
	bkgbkgdhist->store()->get(i);
	float sigma_bkgbkg = bkgbkgdhist->store()->weightError(RooAbsData::SumW2)/bkgbkgdhist->sumEntries()*dataset_dhist->sumEntries()*fbkgbkg->getVal();
   
//	if (sigsigdhist->store()->weightError(RooAbsData::SumW2)/sigsigdhist->store()->weight()>0.3) std::cout << "WARNING: bin known at " << sigsigdhist->store()->weightError(RooAbsData::SumW2)/sigsigdhist->store()->weight() << " sigsig " << r1 << " " << r2 << " " << sigsigdhist->store()->weight() << " +/- " << sigsigdhist->store()->weightError(RooAbsData::SumW2) << std::endl;
//	if (sigbkgdhist->store()->weightError(RooAbsData::SumW2)/sigbkgdhist->store()->weight()>0.3) std::cout << "WARNING: bin known at " << sigbkgdhist->store()->weightError(RooAbsData::SumW2)/sigbkgdhist->store()->weight() << " sigbkg " << r1 << " " << r2 << " " << sigbkgdhist->store()->weight() << " +/- " << sigbkgdhist->store()->weightError(RooAbsData::SumW2) << std::endl;
//	if (bkgsigdhist->store()->weightError(RooAbsData::SumW2)/bkgsigdhist->store()->weight()>0.3) std::cout << "WARNING: bin known at " << bkgsigdhist->store()->weightError(RooAbsData::SumW2)/bkgsigdhist->store()->weight() << " bkgsig " << r1 << " " << r2 << " " << bkgsigdhist->store()->weight() << " +/- " << bkgsigdhist->store()->weightError(RooAbsData::SumW2) << std::endl;
//	if (bkgbkgdhist->store()->weightError(RooAbsData::SumW2)/bkgbkgdhist->store()->weight()>0.3) std::cout << "WARNING: bin known at " << bkgbkgdhist->store()->weightError(RooAbsData::SumW2)/bkgbkgdhist->store()->weight() << " bkgbkg " << r1 << " " << r2 << " " << bkgbkgdhist->store()->weight() << " +/- " << bkgbkgdhist->store()->weightError(RooAbsData::SumW2) << std::endl;

	//	std::cout << " " <<  sigma_sigsig << " " <<  sigma_sigbkg << " " <<  sigma_bkgsig << " " <<  sigma_bkgbkg << " " <<  std::endl;

	float den = sqrt(pow(data_werr,2) + pow(sigma_sigsig,2) + pow(sigma_sigbkg,2) + pow(sigma_bkgsig,2) + pow(sigma_bkgbkg,2));

	if (den==0) {std::cout << "WARNING: 0 denominator" << std::endl; out->fitpulls[i]=0; continue;}

//	model_2D_uncorrelated->Print();
//	model_2D_uncorrelated->Print("v");
//	std::cout << dataset_dhist->sumEntries()*model_2D_uncorrelated->getVal(RooArgSet(*binning_roovar1,*binning_roovar2)) << std::endl;
	//	std::cout << (w-model_2D_uncorrelated_noextended_nll->evaluate())/sqrt(w) << std::endl;
	out->fitpulls[i]=(w-fitvalue)/den;

//	std::cout << r1 << " " << r2 << " " << w << " " << fitvalue << std::endl;
//	std::cout << (w-fitvalue)/data_werr << " " << (w-fitvalue)/den << std::endl;

	out->lowstatbin[i]=false;
	if (w<20) out->lowstatbin[i]=true;
	if (sigsigdhist->store()->weight()<20) out->lowstatbin[i]=true;
	if (sigbkgdhist->store()->weight()<20) out->lowstatbin[i]=true;
	if (bkgsigdhist->store()->weight()<20) out->lowstatbin[i]=true;
	if (bkgbkgdhist->store()->weight()<20) out->lowstatbin[i]=true;

	if (!(out->lowstatbin[i])){
	  ndof++;
	  oldchi2+=pow((w-fitvalue)/data_werr,2);
	  newchi2+=pow((w-fitvalue)/den,2);
	}

      }

      ndof -= (nparams+1);

      if (ndof<1) {
	std::cout << "NOT ENOUGH BINS TO DO THE CHI2" << std::endl;
	out->chi2=0;
	out->ndof=1;
	out->probchi2=0;
      }
      else {
	std::cout << "OLDCHI2 " << oldchi2 << std::endl;
	std::cout << "NEWCHI2 " << newchi2 << " " << ndof << " " << newchi2/ndof << " " << TMath::Prob(newchi2,ndof) << std::endl; 
	out->chi2=newchi2;
	out->ndof=ndof;
	out->probchi2=TMath::Prob(newchi2,ndof);
      }

      delete chi2var;
      delete dataset_dhist;
    }

    if (doplots_ub){
      TCanvas *c2_ub = new TCanvas("c2_ub","c2_ub",1200,800);
      c2_ub->Divide(2);   
       
      c2_ub->cd(1);
      RooPlot *frame1bla = roovar1->frame(Title(""));
      frame1bla->GetXaxis()->SetLabelSize(0.03);
      frame1bla->GetYaxis()->SetLabelSize(0.03);
      frame1bla->GetYaxis()->SetTitleOffset(1.45);
      dataset->plotOn(frame1bla,Name("data"));
      model_axis1_unbinned->plotOn(frame1bla,Name("fit"));
      model_axis1_unbinned->plotOn(frame1bla,Components("sigpdf_axis1_unbinned"),Name("plot_sigsig_axis1_unbinned"),Normalization(fsigsig->getVal()/fsig1->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kRed));
      model_axis1_unbinned->plotOn(frame1bla,Components("sigpdf_axis1_unbinned"),Name("plot_sigbkg_axis1_unbinned"),Normalization(fsigbkg->getVal()/fsig1->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kGreen));
      model_axis1_unbinned->plotOn(frame1bla,Components("bkgpdf_axis1_unbinned"),Name("plot_bkgsig_axis1_unbinned"),Normalization(fbkgsig->getVal()/(1-fsig1->getVal()),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kCyan));
      model_axis1_unbinned->plotOn(frame1bla,Components("bkgpdf_axis1_unbinned"),Name("plot_bkgbkg_axis1_unbinned"),Normalization(fbkgbkg->getVal()/(1-fsig1->getVal()),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kBlack));
      frame1bla->Draw();
      TLegend *leg = new TLegend(0.55,0.6,0.9,0.9);
      leg->AddEntry("data","data","lp");
      leg->AddEntry("fit","fit","l");
      leg->AddEntry("plot_sigsig_axis1_unbinned","prompt-prompt","l");
      leg->AddEntry("plot_sigbkg_axis1_unbinned","prompt-fake","l");
      leg->AddEntry("plot_bkgsig_axis1_unbinned","fake-prompt","l");
      leg->AddEntry("plot_bkgbkg_axis1_unbinned","fake-fake","l");
      leg->SetFillColor(kWhite);
      leg->Draw();
    TLatex a;
    a.SetNDC();
    a.SetTextSize(0.03);
    a.DrawLatex(0.58,0.55,"#splitline{CMS Preliminary}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");

      c2_ub->cd(2);
      RooPlot *frame2bla = roovar2->frame(Title(""));
      frame2bla->GetXaxis()->SetLabelSize(0.03);
      frame2bla->GetYaxis()->SetLabelSize(0.03);
      frame2bla->GetYaxis()->SetTitleOffset(1.45);
      dataset->plotOn(frame2bla);
      model_axis2_unbinned->plotOn(frame2bla);
      model_axis2_unbinned->plotOn(frame2bla,Components("sigpdf_axis2_unbinned"),Normalization(fsigsig->getVal()/fsig2->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kRed));
      model_axis2_unbinned->plotOn(frame2bla,Components("sigpdf_axis2_unbinned"),Normalization(fbkgsig->getVal()/fsig2->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kCyan));
      model_axis2_unbinned->plotOn(frame2bla,Components("bkgpdf_axis2_unbinned"),Normalization(fsigbkg->getVal()/(1-fsig2->getVal()),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kGreen));
      model_axis2_unbinned->plotOn(frame2bla,Components("bkgpdf_axis2_unbinned"),Normalization(fbkgbkg->getVal()/(1-fsig2->getVal()),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kBlack));
      frame2bla->Draw();
      leg->Draw();
    a.DrawLatex(0.58,0.55,"#splitline{CMS Preliminary}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");

      c2_ub->SaveAs(Form("plots/fittingplot2unbinned_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));   
      c2_ub->SaveAs(Form("plots/fittingplot2unbinned_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));   
      c2_ub->SaveAs(Form("plots/fittingplot2unbinned_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));   
      c2_ub->SaveAs(Form("plots/fittingplot2unbinned_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));   


    } // c2_ub

    if (doplots_b) {
      TCanvas *c2 = new TCanvas("c2","c2",1200,800);
      c2->Divide(2,2);   

      {       
	c2->cd(1);
	RooPlot *frame1final = binning_roovar1->frame(Title("Fit axis 1 - binned"));
	dataset->plotOn(frame1final,Name("data"));
	model_2D_uncorrelated->plotOn(frame1final,Name("fit"));
	sigsigpdf->plotOn(frame1final,Normalization(fsigsig->getVal(),RooAbsPdf::Relative),Name("plot_sigsig_axis1_unbinned"),LineStyle(kDashed),LineColor(kRed));	  
	sigbkgpdf->plotOn(frame1final,Normalization(fsigbkg->getVal(),RooAbsPdf::Relative),Name("plot_sigbkg_axis1_unbinned"),LineStyle(kDashed),LineColor(kGreen));  
	bkgsigpdf->plotOn(frame1final,Normalization(fbkgsig->getVal(),RooAbsPdf::Relative),Name("plot_bkgsig_axis1_unbinned"),LineStyle(kDashed),LineColor(kCyan));
	bkgbkgpdf->plotOn(frame1final,Normalization(fbkgbkg->getVal(),RooAbsPdf::Relative),Name("plot_bkgbkg_axis1_unbinned"),LineStyle(kDashed),LineColor(kBlack));  
	frame1final->Draw();
	//    c2->GetPad(1)->SetLogy(1);
	TLegend *leg = new TLegend(0.18,0.74,0.38,0.94);
	leg->AddEntry("data","data","lp");
	leg->AddEntry("fit","fit","l");
	leg->AddEntry("plot_sigsig_axis1_unbinned","prompt-prompt","l");
	leg->AddEntry("plot_sigbkg_axis1_unbinned","prompt-fake","l");
	leg->AddEntry("plot_bkgsig_axis1_unbinned","fake-prompt","l");
	leg->AddEntry("plot_bkgbkg_axis1_unbinned","fake-fake","l");
	leg->SetFillColor(kWhite);
	leg->Draw();
  
	c2->cd(2);
	RooPlot *frame2final = binning_roovar2->frame(Title("Fit axis 2 - binned"));
	dataset->plotOn(frame2final);
	model_2D_uncorrelated->plotOn(frame2final);
	sigsigpdf->plotOn(frame2final,Normalization(fsigsig->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kRed));	  
	sigbkgpdf->plotOn(frame2final,Normalization(fsigbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kGreen));  
	bkgsigpdf->plotOn(frame2final,Normalization(fbkgsig->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kCyan));
	bkgbkgpdf->plotOn(frame2final,Normalization(fbkgbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kBlack));  
	frame2final->Draw();
	//    c2->GetPad(2)->SetLogy(1);
	leg->Draw();
      }

      /*
	c2->cd(1);
	RooPlot *ppnllplot = pp->frame();
	model_2D_uncorrelated_noextended_nll->plotOn(ppnllplot,ShiftToZero());
	ppnllplot->Draw();
      */

      /*
      c2->cd(3);
      RooPlot *j1nllplot = j1->frame();
      model_2D_uncorrelated_noextended_nll->plotOn(j1nllplot,ShiftToZero());
      model_2axes_nll->plotOn(j1nllplot,ShiftToZero(),LineColor(kRed));    
      j1nllplot->Draw();
      */

      /*
      c2->cd(3);
      RooPlot *contourplot = minuit_secondpass->contour(*pp,*j1,1,0,0,0,0,0);
      contourplot->SetAxisRange(out->pp-0.05,out->pp+0.05,"X");
      contourplot->SetAxisRange(j1->getVal()-0.05,j1->getVal()+0.05,"Y");
      contourplot->Draw();
      */

      TH2F *h2d;
      TH2F *h2f;
      TH2F *h2l;
      TH2F *h2h;
      TH1F *h3d;
      TH1F *h3f;
      TH1F *h3l;
      TH1F *h3h;
      h2d = new TH2F("h2d","h2d",n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
      h2f = new TH2F("h2f","h2f",n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
      h2l = new TH2F("h2l","h2l",n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
      h2h = new TH2F("h2h","h2h",n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
      h3d = new TH1F("h3d","Diagonal projection",n_templatebins,0.5*sqrt(2),(0.5+n_templatebins)*sqrt(2));
      h3f = new TH1F("h3f","h3f",n_templatebins,0.5*sqrt(2),(0.5+n_templatebins)*sqrt(2));
      h3l = new TH1F("h3l","h3l",n_templatebins,0.5*sqrt(2),(0.5+n_templatebins)*sqrt(2));
      h3h = new TH1F("h3h","h3h",n_templatebins,0.5*sqrt(2),(0.5+n_templatebins)*sqrt(2));
      h2d->Sumw2();
      h2f->Sumw2();
      h2l->Sumw2();
      h2h->Sumw2();
      h3d->Sumw2();
      h3f->Sumw2();
      h3l->Sumw2();
      h3h->Sumw2();


      for (int i=0; i<dataset->numEntries(); i++){
	float r1 = dataset->get(i)->getRealValue("binning_roovar1");
	float r2 = dataset->get(i)->getRealValue("binning_roovar2");
	float w = dataset->store()->weight(i);
	h2d->Fill(r1,r2,w);
	h3d->Fill((r1+r2)/sqrt(2),w);
      }

      RooDataSet *rand_uncorrelated = model_2D_uncorrelated->generate(RooArgSet(*binning_roovar1,*binning_roovar2),1e6);
      for (int i=0; i<rand_uncorrelated->numEntries(); i++){
	float r1 = rand_uncorrelated->get(i)->getRealValue("binning_roovar1");
	float r2 = rand_uncorrelated->get(i)->getRealValue("binning_roovar2");
	float w = rand_uncorrelated->store()->weight(i);
	h2f->Fill(r1,r2,w);
	h3f->Fill((r1+r2)/sqrt(2),w);
      }
      delete rand_uncorrelated;

      pp->setVal(minpp);
      RooDataSet *rand_uncorrelated_low = model_2D_uncorrelated->generate(RooArgSet(*binning_roovar1,*binning_roovar2),1e6);
      for (int i=0; i<rand_uncorrelated_low->numEntries(); i++){
	float r1 = rand_uncorrelated_low->get(i)->getRealValue("binning_roovar1");
	float r2 = rand_uncorrelated_low->get(i)->getRealValue("binning_roovar2");
	float w = rand_uncorrelated_low->store()->weight(i);
	h2l->Fill(r1,r2,w);
	h3l->Fill((r1+r2)/sqrt(2),w);
      }
      delete rand_uncorrelated_low;

      pp->setVal(maxpp);
      RooDataSet *rand_uncorrelated_high = model_2D_uncorrelated->generate(RooArgSet(*binning_roovar1,*binning_roovar2),1e6);
      for (int i=0; i<rand_uncorrelated_high->numEntries(); i++){
	float r1 = rand_uncorrelated_high->get(i)->getRealValue("binning_roovar1");
	float r2 = rand_uncorrelated_high->get(i)->getRealValue("binning_roovar2");
	float w = rand_uncorrelated_high->store()->weight(i);
	h2h->Fill(r1,r2,w);
	h3h->Fill((r1+r2)/sqrt(2),w);
      }
      delete rand_uncorrelated_high;

      pp->setVal(out->pp); // revert pp to fitted number

      h2d->Print();
      h2f->Print();
      h2l->Print();
      h2h->Print();
      h3d->Print();
      h3f->Print();
      h3l->Print();
      h3h->Print();


//      TCanvas *c3 = new TCanvas("c3","c3",1200,800);
//      c3->Divide(3);
//      c3->cd();

      for (int k=0; k<h2f->GetNbinsX(); k++) for (int l=0; l<h2f->GetNbinsY(); l++) h2f->SetBinError(k+1,l+1,0);
      for (int k=0; k<h2l->GetNbinsX(); k++) for (int l=0; l<h2l->GetNbinsY(); l++) h2l->SetBinError(k+1,l+1,0);
      for (int k=0; k<h2h->GetNbinsX(); k++) for (int l=0; l<h2h->GetNbinsY(); l++) h2h->SetBinError(k+1,l+1,0);
      for (int k=0; k<h3f->GetNbinsX(); k++) h3f->SetBinError(k+1,0);
      for (int k=0; k<h3l->GetNbinsX(); k++) h3l->SetBinError(k+1,0);
      for (int k=0; k<h3h->GetNbinsX(); k++) h3h->SetBinError(k+1,0);

      h2d->Scale(1.0/h2d->Integral());
      h2f->Scale(1.0/h2f->Integral());
      h2l->Scale(1.0/h2l->Integral());
      h2h->Scale(1.0/h2h->Integral());
      h3d->Scale(1.0/h3d->Integral());
      h3f->Scale(1.0/h3f->Integral());
      h3l->Scale(1.0/h3l->Integral());
      h3h->Scale(1.0/h3h->Integral());


      h2d->SetMarkerStyle(20);
      h2f->SetLineWidth(2);
      h2l->SetLineWidth(2);
      h2h->SetLineWidth(2);
      h3d->SetMarkerStyle(20);
      h3f->SetLineWidth(2);
      h3l->SetLineWidth(2);
      h3h->SetLineWidth(2);

      h2d->SetMarkerColor(kBlack);
      h2f->SetLineColor(kBlue);
      h2l->SetLineColor(kRed);
      h2h->SetLineColor(kGreen);
      h3d->SetMarkerColor(kBlack);
      h3f->SetLineColor(kBlue);
      h3l->SetLineColor(kRed);
      h3h->SetLineColor(kGreen);
      
//      c3->cd(1);
//      h2d->ProjectionX()->Draw();
//      h2f->ProjectionX()->Draw("same");
//      c3->cd(2);
//      h2d->ProjectionY()->Draw();
//      h2f->ProjectionY()->Draw("same");
//      c3->cd(3);
//      h3d->Draw();
//      h3f->Draw("same");
//      h3l->Draw("same");
//      h3h->Draw("same");

      {
	c2->cd(4);
	h3d->SetStats(0);
	h3d->GetXaxis()->SetTitle("Diag. projection binned (Iso_{1},Iso_{2})");
	h3d->Draw();
	h3f->Draw("same");
	h3l->Draw("same");
	h3h->Draw("same");
	
	TLegend *leg = new TLegend(0.75,0.7,0.95,0.9);
	leg->AddEntry(h3d,"data","lp");
	leg->AddEntry(h3f,"fit","l");
	leg->AddEntry(h3h,"max pp","l");
	leg->AddEntry(h3l,"min pp","l");
	leg->SetFillColor(kWhite);
	leg->Draw();
      }

      c2->Update();
      //      c3->Update();
      c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));   
      c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));   
      c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));   
      c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));   
      //      c3->SaveAs(Form("plots/fittingplot3_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
   

    } // c2 binned

    delete minuit_secondpass;
    delete minuit_secondpass_constraint;
    delete model_2D_uncorrelated_noextended_nll_constraint;
    delete model_2D_uncorrelated_noextended_nll;
    delete model_2D_uncorrelated_unbinned;
    delete model_2D_uncorrelated;
    delete constraint_gaussian_nll;
    delete constraint_gaussian_j2;
    delete constraint_parameters_set;
    delete constraint_pdf_set;
    delete constraint_gaussian_j1;
    delete fbkgbkg;
    delete fbkgsig;
    delete fsigbkg;
    delete fsigsig;
    delete minuit_firstpass;
    delete model_2axes_nll;
    delete model_axis2_noextended_nll;
    delete model_axis1_noextended_nll;
    delete model_axis2_unbinned;
    delete model_axis1_unbinned;
    delete model_axis2;
    delete model_axis1;
    delete fsig2;
    delete fsig1;
    delete bkgpdf_axis2_unbinned;
    delete sigpdf_axis2_unbinned;
    delete bkgdhist_axis2_unbinned;
    delete sigdhist_axis2_unbinned;
    delete bkgpdf_axis1_unbinned;
    delete sigpdf_axis1_unbinned;
    delete bkgdhist_axis1_unbinned;
    delete sigdhist_axis1_unbinned;
    delete bkgbkgpdf_unbinned;
    delete bkgsigpdf_unbinned;
    delete sigbkgpdf_unbinned;
    delete sigsigpdf_unbinned;
    delete bkgbkgdhist_unbinned;
    delete bkgsigdhist_unbinned;
    delete sigbkgdhist_unbinned;
    delete sigsigdhist_unbinned;

    if (do_syst_string==TString("savepdfMCtrue2D") || do_syst_string==TString("savepdfMCtrue1D")){
      out->pdf_forgen[0] = sigsigpdf;
      out->pdf_forgen[1] = sigbkgpdf;
      out->pdf_forgen[2] = bkgsigpdf;
      out->pdf_forgen[3] = bkgbkgpdf;
      out->pdf_forgen[4] = sigpdf_axis1;
      out->pdf_forgen[5] = bkgpdf_axis1;
      out->pdf_forgen[6] = sigpdf_axis2;
      out->pdf_forgen[7] = bkgpdf_axis2;
      std::cout << "RUN FOR SAVE PDF FINISHED" << std::endl;
      return out;
    }

    if (do_syst_string=="purefitbias" && runcount==0){
      datafr=out;
      continue;
    }

    delete bkgpdf_axis2;
    delete sigpdf_axis2;
    delete bkgdhist_axis2;
    delete sigdhist_axis2;
    delete bkgpdf_axis1;
    delete sigpdf_axis1;
    delete bkgdhist_axis1;
    delete sigdhist_axis1;

    delete bkgbkgpdf;
    delete bkgsigpdf;
    delete sigbkgpdf;
    delete sigsigpdf;
    delete bkgbkgdhist;
    delete bkgsigdhist;
    delete sigbkgdhist;
    delete sigsigdhist;
    delete j2;
    delete j1;
    delete pp;

    if (do_syst_string==TString("templatestatistics")) do_syst_templatestatistics_outputvector.push_back(out);
    if (do_syst_string==TString("purefitbias") && runcount>0) do_syst_purefitbias_outputvector.push_back(out);
    if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE")) do_syst_MCpromptdriven_outputvector.push_back(out);
    if (do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEE")) do_syst_MCfakedriven_outputvector.push_back(out);
    if (do_syst_string==TString("templateshape2frag")) do_syst_2frag_outputvector.push_back(out);

    delete dataset_sigsig   ;
    delete dataset_sigbkg   ;
    delete dataset_bkgsig   ;
    delete dataset_bkgbkg   ;
    delete dataset_sig_axis1;
    delete dataset_bkg_axis1;
    delete dataset_sig_axis2;
    delete dataset_bkg_axis2;
    delete dataset_axis1    ;
    delete dataset_axis2    ;
    delete dataset          ;
      
    dataset_sigsig   =original_dataset_sigsig   ;
    dataset_sigbkg   =original_dataset_sigbkg   ;
    dataset_bkgsig   =original_dataset_bkgsig   ;
    dataset_bkgbkg   =original_dataset_bkgbkg   ;
    dataset_sig_axis1=original_dataset_sig_axis1;
    dataset_bkg_axis1=original_dataset_bkg_axis1;
    dataset_sig_axis2=original_dataset_sig_axis2;
    dataset_bkg_axis2=original_dataset_bkg_axis2;
    dataset_axis1    =original_dataset_axis1    ;
    dataset_axis2    =original_dataset_axis2    ;
    dataset          =original_dataset          ;

  print_mem();
    
  }
  


//  delete dataset_axis2;
//  delete dataset_axis1;
//  delete dataset_bkg_axis2;
//  delete dataset_sig_axis2;
//  delete dataset_bkg_axis1;
//  delete dataset_sig_axis1;
//  delete binning_roovar2;
//  delete binning_roovar1;
//  delete binning_roovar2_threshold;
//  delete binning_roovar1_threshold;
//  delete rooweight;
//  delete dataset;
//  delete dataset_bkgbkg;
//  delete dataset_bkgsig;
//  delete dataset_sigbkg;
//  delete dataset_sigsig;
//  delete dataset_sigsig_orig;
//  delete dataset_sigbkg_orig;
//  delete dataset_bkgsig_orig;
//  delete dataset_bkgbkg_orig;
//  delete dataset_orig;
//  delete roovar1;
//  delete roopt1;
//  delete rooeta1;
//  delete roosieie1;
//  delete roovar2;
//  delete roopt2;
//  delete rooeta2;
//  delete roosieie2;
//  delete roorho;                                                                                                       




  if (do_syst_string==TString("templatestatistics") || do_syst_string==TString("purefitbias") || do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")){

    std::vector<fit_output*> *do_syst_vector = NULL;
    if (do_syst_string==TString("templatestatistics")){
      do_syst_vector = &do_syst_templatestatistics_outputvector;
    }
    else if (do_syst_string==TString("purefitbias")){
      do_syst_vector = &do_syst_purefitbias_outputvector;
    }
    else if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE")){
      do_syst_vector = &do_syst_MCpromptdriven_outputvector;
    }
    else if (do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEE")){
      do_syst_vector = &do_syst_MCfakedriven_outputvector;
    }
    else if (do_syst_string==TString("templateshape2frag")){
      do_syst_vector = &do_syst_2frag_outputvector;
    }
    assert (do_syst_vector);

    float mean=0;
    for (unsigned int i=0; i<do_syst_vector->size(); i++) mean+=do_syst_vector->at(i)->pp;
    mean/=do_syst_vector->size();
    std::cout << "Mean pp " << mean << std::endl;

    float centerval = datafr->pp;
    std::cout << "Center val: " << centerval << std::endl;


    RooRealVar fittedx("fittedx","fittedx",0,1);
    RooRealVar fittedpull("fittedpull","fittedpull",-20,20);    
    RooDataSet dsetx("dsetx","dsetx",fittedx);
    RooDataSet dsetpull("dsetpull","dsetpull",fittedpull);

    RooRealVar meangaus("meangaus","meangaus",0,1);
    RooRealVar sigmagaus("sigmagaus","sigmagaus",0,0.3);
    RooGaussian gaus("gaus","gaus",fittedx,meangaus,sigmagaus);

    meangaus.setVal(0.2);
    sigmagaus.setVal(0.03);

    RooRealVar meangauspull("meangauspull","meangauspull",-20,20);
    RooRealVar sigmagauspull("sigmagauspull","sigmagauspull",0,3);
    RooGaussian gauspull("gauspull","gauspull",fittedpull,meangauspull,sigmagauspull);

    for (unsigned int i=0; i<do_syst_vector->size(); i++){
      fittedx.setVal(do_syst_vector->at(i)->pp);
      fittedpull.setVal((do_syst_vector->at(i)->pp-centerval)/do_syst_vector->at(i)->pp_err);
      std::cout << "Toy nr. " << i << " " << fittedx.getVal() << " " << do_syst_vector->at(i)->pp_err << " " << fittedpull.getVal() << std::endl;
      dsetx.add(fittedx);
      dsetpull.add(fittedpull);
    }

    dsetx.Print();
    dsetpull.Print();
    
    RooFitResult *biasfitresult = gaus.fitTo(dsetx,Save(),Range(0.01,1));
    RooFitResult *biasfitresultpull = gauspull.fitTo(dsetpull,Save(),Range(-10,10));

    fittedx.setRange(mean-0.20,mean+0.20);
    RooPlot *framefittedx = fittedx.frame();
    dsetx.plotOn(framefittedx);
    gaus.plotOn(framefittedx);
    gaus.paramOn(framefittedx);
    fittedpull.setRange(-20,20);
    RooPlot *framefittedpull = fittedpull.frame();
    dsetpull.plotOn(framefittedpull);
    gauspull.plotOn(framefittedpull);
    gauspull.paramOn(framefittedpull);

    TFile *file_bias = new TFile(Form("plots/bias_%s_%s_%s_b%d.root",do_syst_string.Data(),diffvariable.Data(),splitting.Data(),bin),"recreate");
    file_bias->cd();
    framefittedx->Write();
    framefittedpull->Write();
    fittedx.Write();
    fittedpull.Write();
    dsetx.Write();
    dsetpull.Write();
    file_bias->Close();

    TH1F *histo_bias = new TH1F(Form("histo_bias_%s",do_syst_string.Data()),Form("histo_bias_%s",do_syst_string.Data()),bins_to_run,binsdef);
    {
      TString unit = diffvariables_units_list(diffvariable);
      histo_bias->GetXaxis()->SetTitle(Form("%s %s",diffvariables_names_list(diffvariable).Data(),unit!=TString("") ? (TString("(").Append(unit.Append(")"))).Data() : TString("").Data()));
    }
    histo_bias->SetTitle("");
    histo_bias->SetStats(0);
    histo_bias->SetMarkerStyle(20);
    if (do_syst_string==TString("templatestatistics")) {    
      histo_bias->GetYaxis()->SetTitle("Width of fitted purity / data purity");
      histo_bias->SetBinContent(bin+1,sigmagaus.getVal()/centerval);
      histo_bias->SetBinError(bin+1,sigmagaus.getPropagatedError(*biasfitresult)/centerval);
    }
    if (do_syst_string==TString("purefitbias")) {
      histo_bias->GetYaxis()->SetTitle(Form("Fitted mean of the pull (%f pp purity toys)",centerval));
      histo_bias->SetBinContent(bin+1,meangauspull.getVal());
      histo_bias->SetBinError(bin+1,meangauspull.getPropagatedError(*biasfitresultpull));
    }
    if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshape2frag")){
      histo_bias->GetYaxis()->SetTitle(Form("Fitted purity / gen. purity (%f purity toys)",centerval));
      histo_bias->SetBinContent(bin+1,meangaus.getVal()/centerval);
      histo_bias->SetBinError(bin+1,meangaus.getPropagatedError(*biasfitresult)/centerval);
    }
    histo_bias->SaveAs(Form("plots/histo_bias_%s_%s_%s_b%d.root",do_syst_string.Data(),diffvariable.Data(),splitting.Data(),bin));    


  }

  bool writeoutpurity = (do_syst_string==TString("") || do_syst_string==TString("doMCfulldriven") || do_syst_string==TString("newtemplates"));

  if (writeoutpurity){

    TH1F *purity[4];
    TH1F *eventshisto;
    TH1F *redchi2;
    TH1F *probchi2;
    TH1F *fitpull;
    TH1F *fitpull_histo;
    if (bins_to_run>0) {eventshisto = new TH1F("eventshisto","eventshisto",bins_to_run,binsdef); redchi2 = new TH1F("redchi2","redchi2",bins_to_run,binsdef); probchi2 = new TH1F("probchi2","probchi2",bins_to_run,binsdef);}
    else {eventshisto = new TH1F("eventshisto","eventshisto",n_bins,0,n_bins); redchi2 = new TH1F("redchi2","redchi2",n_bins,0,n_bins); probchi2 = new TH1F("probchi2","probchi2",n_bins,0,n_bins);}
    TH1F *overflowremovaleffhisto;
    if (bins_to_run>0) overflowremovaleffhisto = new TH1F("overflowremovaleffhisto","overflowremovaleffhisto",bins_to_run,binsdef);
    else overflowremovaleffhisto = new TH1F("overflowremovaleffhisto","overflowremovaleffhisto",n_bins,0,n_bins);
    fitpull = new TH1F("fitpull","fitpull",100,0,100);
    fitpull_histo = new TH1F("fitpull_histo","fitpull_histo",25,-5,5);

    std::cout << "Output histos created" << std::endl;

    int colors[4] = {kRed, kGreen, kCyan, kBlack};
    
    for (int i=0; i<4; i++){
      TString name = "purity_";
      TString title = Form("Purity - %s category",splitting.Data());
      if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshapeMCfulldriven") || do_syst_string==TString("subtractionZee") || do_syst_string==TString("templateshape2frag")) {name+=do_syst_string; name.Append("_");}
      if (i==0) name.Append("sigsig"); else if (i==1) name.Append("sigbkg"); else if (i==2) name.Append("bkgsig"); else if (i==3) name.Append("bkgbkg");
      if (bins_to_run>0) purity[i] = new TH1F(name.Data(),title.Data(),bins_to_run,binsdef);
      else purity[i] = new TH1F(name.Data(),title.Data(),n_bins,0,n_bins);
      purity[i]->SetMarkerStyle(20);
      purity[i]->SetMarkerColor(colors[i]);
      purity[i]->SetLineColor(colors[i]);
      purity[i]->SetLineWidth(2);
      purity[i]->GetYaxis()->SetRangeUser(0,1);
      purity[i]->GetYaxis()->SetTitle("");
      TString unit = diffvariables_units_list(diffvariable);
      purity[i]->GetXaxis()->SetTitle(Form("%s %s",diffvariables_names_list(diffvariable).Data(),unit!=TString("") ? (TString("(").Append(unit.Append(")"))).Data() : TString("").Data()));
    } 

    std::cout << "Output histos formatted" << std::endl;

    purity[0]->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->pp);
    purity[0]->SetBinError((bin!=n_bins) ? bin+1 : 1,out->pp_err);
    purity[1]->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->pf);
    purity[1]->SetBinError((bin!=n_bins) ? bin+1 : 1,out->pf_err);
    purity[2]->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->fp);
    purity[2]->SetBinError((bin!=n_bins) ? bin+1 : 1,out->fp_err);
    purity[3]->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->ff);
    purity[3]->SetBinError((bin!=n_bins) ? bin+1 : 1,out->ff_err);
    eventshisto->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->tot_events);
    redchi2->SetBinContent((bin!=n_bins) ? bin+1 : 1, out->chi2/out->ndof);
    probchi2->SetBinContent((bin!=n_bins) ? bin+1 : 1, out->probchi2);
    for (int i=0; i<n_templatebins*n_templatebins; i++) fitpull->SetBinContent(i+1,out->fitpulls[i]);
    for (int i=0; i<n_templatebins*n_templatebins; i++) if (!(out->lowstatbin[i])) fitpull_histo->Fill(out->fitpulls[i]);
    overflowremovaleffhisto->SetBinContent((bin!=n_bins) ? bin+1 : 1,out->eff_overflow_removal_pp);

    std::cout << "Output histos filled" << std::endl;

    TString helper("");
    if (do_syst_string==TString("templateshapeMCpromptdrivenEB") || do_syst_string==TString("templateshapeMCfakedrivenEB") || do_syst_string==TString("templateshapeMCpromptdrivenEE") || do_syst_string==TString("templateshapeMCfakedrivenEE") || do_syst_string==TString("templateshapeMCfulldriven") || do_syst_string==TString("subtractionZee") || do_syst_string==TString("templateshape2frag")) {helper = do_syst_string; helper+=TString("_");}
    
    TFile *purityfile = new TFile(Form("plots/histo_purity_%s%s_%s_b%d.root",helper.Data(),diffvariable.Data(),splitting.Data(),bin),"recreate");
    purityfile->cd();
    for (int i=0; i<4; i++) purity[i]->Write();
    eventshisto->Write();
    overflowremovaleffhisto->Write();
    redchi2->Write();
    probchi2->Write();
    fitpull->Write();
    fitpull_histo->Write();
    purityfile->Close();

    std::cout << "Output histos written" << std::endl;

  }


  return out;
    

};

void fit_dataset_allbins(TString diffvariable="", TString splitting="", TString do_syst_string=TString("")){

  int bins_to_run=-1; 
  float *binsdef=NULL;

  if (diffvariable=="invmass"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_invmass_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_invmass_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_invmass_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_invmass_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_invmass_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_invmass_EEEE;
  }
  if (diffvariable=="diphotonpt"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_diphotonpt_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_diphotonpt_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_diphotonpt_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_diphotonpt_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_diphotonpt_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_diphotonpt_EEEE;
  }
  if (diffvariable=="costhetastar"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_costhetastar_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_costhetastar_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_costhetastar_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_costhetastar_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_costhetastar_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_costhetastar_EEEE;
  }
  if (diffvariable=="dphi"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_dphi_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_dphi_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_dphi_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_dphi_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_dphi_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_dphi_EEEE;
  }
  if (diffvariable=="dR"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_dR_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_dR_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_dR_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_dR_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_dR_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_dR_EEEE;
  }
  
  fit_output *fr[n_bins];

  for (int bin=0; bin<bins_to_run; bin++) {
    fr[bin]=fit_dataset(diffvariable,splitting,bin,do_syst_string);
  }

};


void post_process(TString diffvariable="", TString splitting="", bool skipsystematics=false){

  const float intlumi=5.044;

  TH1F *xsec;
  TH1F *xsec_withsyst;
  TH1F *xsec_ngammagammayield = NULL;

  TH1F *fitpull_histo = NULL;

  TH1F *systplot_purefitbias=NULL;
  TH1F *systplot_templatestatistics=NULL;
  TH1F *systplot_templateshapeMCpromptdrivenEB=NULL;
  TH1F *systplot_templateshapeMCfakedrivenEB=NULL;
  TH1F *systplot_templateshapeMCpromptdrivenEE=NULL;
  TH1F *systplot_templateshapeMCfakedrivenEE=NULL;
  TH1F *systplot_templateshape2frag=NULL;
  TH1F *systplot_zee=NULL;
  TH1F *systplot_tot=NULL;
  TH1F *systplot_efficiency=NULL;
  TH1F *systplot_unfolding=NULL;
  TH1F *systplot_totfinal=NULL;

  TH1F *systplot_uncorrelated=NULL;
  TH1F *systplot_1catcorrelated=NULL;
  TH1F *systplot_allcatcorrelated=NULL;

  TH1F *systplot_statistic=NULL;

  int bins_to_run=-1;
  float *binsdef=NULL;

  if (diffvariable=="invmass"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_invmass_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_invmass_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_invmass_EEEE; 
    else bins_to_run+=n_templates_invmass_EEEE;
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_invmass_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_invmass_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_invmass_EEEE;
  }
  if (diffvariable=="diphotonpt"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_diphotonpt_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_diphotonpt_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_diphotonpt_EEEE; 
    else bins_to_run+=n_templates_diphotonpt_EEEE;
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_diphotonpt_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_diphotonpt_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_diphotonpt_EEEE;
  }
  if (diffvariable=="costhetastar"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_costhetastar_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_costhetastar_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_costhetastar_EEEE; 
    else bins_to_run+=n_templates_costhetastar_EEEE;
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_costhetastar_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_costhetastar_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_costhetastar_EEEE;
  }
  if (diffvariable=="dphi"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_dphi_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_dphi_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_dphi_EEEE; 
    else bins_to_run+=n_templates_dphi_EEEE;
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_dphi_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_dphi_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_dphi_EEEE;
  }
  if (diffvariable=="dR"){
    if (splitting=="EBEB")      bins_to_run+=n_templates_dR_EBEB;
    else if (splitting=="EBEE") bins_to_run+=n_templates_dR_EBEE;
    else if (splitting=="EEEE") bins_to_run+=n_templates_dR_EEEE; 
    else bins_to_run+=n_templates_dR_EEEE;
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_dR_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_dR_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_dR_EEEE;
  }



  if (splitting=="inclusive"){
    TString sp[3]={"EBEB","EBEE","EEEE"};
    TFile *axsec_file[3];
    TH1F *axsec[3];
    TH1F *axsec_withsyst[3];
    TH1F *axsec_ngammagammayield[3];
    TH1F *axsec_fitpull_histo[3];
    for (int i=0; i<3; i++) {
      axsec_file[i] = new TFile(Form("plots/histo_xsec_%s_%s.root", diffvariable.Data(),sp[i].Data()));
      axsec_file[i]->GetObject("xsec",axsec[i]);
      axsec_file[i]->GetObject("xsec_withsyst",axsec_withsyst[i]);
      axsec_file[i]->GetObject("xsec_ngammagammayield",axsec_ngammagammayield[i]);
      axsec_file[i]->GetObject("fitpull_histo",axsec_fitpull_histo[i]);
    }
    for (int i=1; i<3; i++) {
      axsec[0]->Add(axsec[i]);
      axsec_withsyst[0]->Add(axsec_withsyst[i]);
      axsec_ngammagammayield[0]->Add(axsec_ngammagammayield[i]);
      axsec_fitpull_histo[0]->Add(axsec_fitpull_histo[i]);
    }
    xsec=axsec[0];
    xsec_withsyst=axsec_withsyst[0];
    xsec_ngammagammayield=axsec_ngammagammayield[0];
    fitpull_histo=axsec_fitpull_histo[0];
  }

  else {

  bool sym = false;
  
  TString s1; TString s2;
  if (splitting=="EBEB") {s1="EB"; s2="EB";}
  else if (splitting=="EEEE") {s1="EE"; s2="EE";}
  else if (splitting=="EBEE") {s1="EB"; s2="EE";}
  sym  = (s1==s2);
  

  TH1F::SetDefaultSumw2(kTRUE);
  
  //  fit_output *fr[n_bins];

  TH1F *purity[4];
  TH1F *eventshisto;
  TH1F *overflowremovaleffhisto;

  TH1F *eff=NULL;
  {
    TFile *eff_file = new TFile("plots/Efficiency_WithSysErr.root");
    std::map<TString,TString> translation;
    translation.insert(std::pair<TString,TString>(TString("invmass"),TString("mgg")));
    translation.insert(std::pair<TString,TString>(TString("diphotonpt"),TString("qtgg")));
    translation.insert(std::pair<TString,TString>(TString("costhetastar"),TString("costhetastar")));
    translation.insert(std::pair<TString,TString>(TString("dphi"),TString("deltaphi")));
    translation.insert(std::pair<TString,TString>(TString("dR"),TString("dR")));
    eff_file->GetObject(Form("h_%s_%s_WithTotErr",translation[diffvariable].Data(),splitting.Data()),eff);
  }
  assert (eff!=NULL);

//  TH1F *eff=NULL;
//  TFile *eff_file = new TFile("efficiencies.root");
//  eff_file->GetObject(Form("w_eff_gg_%s_%s",splitting.Data(),diffvariable.Data()),eff);
//  assert (eff!=NULL);
//  eff->GetYaxis()->SetTitle("selection/ID efficiency");
//  eff->GetXaxis()->SetTitle(diffvariable.Data());


  TH1F *unfoldunc = NULL;

  if (!skipsystematics){
    std::map<TString,TString> translation2;
    translation2.insert(std::pair<TString,TString>(TString("invmass"),TString("mgg")));
    translation2.insert(std::pair<TString,TString>(TString("diphotonpt"),TString("pt")));
    translation2.insert(std::pair<TString,TString>(TString("costhetastar"),TString("costt")));
    translation2.insert(std::pair<TString,TString>(TString("dphi"),TString("phi")));
    translation2.insert(std::pair<TString,TString>(TString("dR"),TString("dR")));
    TFile *unfoldunc_file = new TFile("plots/Unfolding_SysErr.root");
    TString unfoldunc_name = Form("Unfolding_RelativeSysErr_%s_%s",translation2[diffvariable].Data(),splitting.Data());
    unfoldunc_file->GetObject(unfoldunc_name.Data(),unfoldunc);
    assert (unfoldunc!=NULL);
  }


  TFile *purity_file = new TFile(Form("plots/histo_purity_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
  purity_file->GetObject("purity_sigsig",purity[0]);
  purity_file->GetObject("purity_sigbkg",purity[1]);
  purity_file->GetObject("purity_bkgsig",purity[2]);
  purity_file->GetObject("purity_bkgbkg",purity[3]);
  purity_file->GetObject("eventshisto",eventshisto);
  purity_file->GetObject("overflowremovaleffhisto",overflowremovaleffhisto);
  purity_file->GetObject("fitpull_histo",fitpull_histo);

  std::cout << "Purity histos imported" << std::endl;
  for (int i=0; i<4; i++) purity[i]->Print(); eventshisto->Print(); overflowremovaleffhisto->Print(); 

  TH1F *histo_bias_purefitbias = NULL;
  TH1F *histo_bias_templatestatistics = NULL;
  TFile *file_standardsel_dy = NULL;
  TFile *file_pixelrevsel_dy_data = NULL;
  TFile *file_pixelrevsel_dy_mc = NULL;

  TH1F *histo_bias_templateshapeMCpromptdrivenEB = NULL;
  TH1F *histo_bias_templateshapeMCfakedrivenEB = NULL;
  TH1F *histo_bias_templateshapeMCpromptdrivenEE = NULL;
  TH1F *histo_bias_templateshapeMCfakedrivenEE = NULL;
  TH1F *histo_bias_templateshape2frag = NULL;

  bool skipZsubtraction = skipsystematics;

  skipZsubtraction = false; // DEBUG: force Z subtraciton on

  if (!skipZsubtraction){
    file_standardsel_dy = new TFile("outphoton_dy_standard.root");
    assert(file_standardsel_dy);
    file_pixelrevsel_dy_data = new TFile("outphoton_data_pixelrev.root");
    assert(file_pixelrevsel_dy_data);
    file_pixelrevsel_dy_mc = new TFile("outphoton_dy_pixelrev.root");
    assert(file_pixelrevsel_dy_mc);    
  }

  if (!skipsystematics){
    TFile *file_bias_purefitbias  = new TFile(Form("plots/histo_bias_purefitbias_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
    file_bias_purefitbias->GetObject("histo_bias_purefitbias",histo_bias_purefitbias);
    assert(histo_bias_purefitbias);
    TFile *file_bias_templatestatistics  = new TFile(Form("plots/histo_bias_templatestatistics_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
    file_bias_templatestatistics->GetObject("histo_bias_templatestatistics",histo_bias_templatestatistics);
    assert(histo_bias_templatestatistics);
    if (splitting!="EEEE"){
    TFile *file_bias_templateshapepromptEB = new TFile(Form("plots/histo_bias_templateshapeMCpromptdrivenEB_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
    file_bias_templateshapepromptEB->GetObject("histo_bias_templateshapeMCpromptdrivenEB",histo_bias_templateshapeMCpromptdrivenEB);
    TFile *file_bias_templateshapefakeEB = new TFile(Form("plots/histo_bias_templateshapeMCfakedrivenEB_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
    file_bias_templateshapefakeEB->GetObject("histo_bias_templateshapeMCfakedrivenEB",histo_bias_templateshapeMCfakedrivenEB);
    }
    if (splitting!="EBEB"){
    TFile *file_bias_templateshapepromptEE = new TFile(Form("plots/histo_bias_templateshapeMCpromptdrivenEE_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
    file_bias_templateshapepromptEE->GetObject("histo_bias_templateshapeMCpromptdrivenEE",histo_bias_templateshapeMCpromptdrivenEE);
    TFile *file_bias_templateshapefakeEE = new TFile(Form("plots/histo_bias_templateshapeMCfakedrivenEE_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
    file_bias_templateshapefakeEE->GetObject("histo_bias_templateshapeMCfakedrivenEE",histo_bias_templateshapeMCfakedrivenEE);
    }
    TFile *file_bias_templateshape2frag = new TFile(Form("plots/histo_bias_templateshape2frag_%s_%s_allbins.root",diffvariable.Data(),splitting.Data()));
    file_bias_templateshape2frag->GetObject("histo_bias_templateshape2frag",histo_bias_templateshape2frag);
  }
  
    xsec = new TH1F("xsec","xsec",bins_to_run,binsdef);
    xsec->SetMarkerStyle(20);
    xsec->SetMarkerColor(kBlack);
    xsec->SetLineColor(kBlack);
    xsec->SetLineWidth(2);

  
    //    xsec->GetYaxis()->SetTitle("");
    //    xsec->GetXaxis()->SetTitle(diffvariables_names_list(diffvariable).Data());
    {
      TString unit = diffvariables_units_list(diffvariable);
      xsec->GetXaxis()->SetTitle(Form("%s %s",diffvariables_names_list(diffvariable).Data(),unit!=TString("") ? (TString("(").Append(unit.Append(")"))).Data() : TString("").Data()));
    }

    xsec_withsyst = (TH1F*)(xsec->Clone("xsec_withsyst"));
    xsec_withsyst->SetLineColor(kRed);

    xsec_ngammagammayield = (TH1F*)(xsec->Clone("xsec_ngammagammayield"));
    xsec_ngammagammayield->SetMarkerColor(kGreen+2);
    xsec_ngammagammayield->SetLineColor(kGreen+2);

    float uncertainty_scalefactor_Zee_PIXEL=2*0.025; // molt per 2 perche' due leg
    map<TString,pair<float,float> > syst_purity_dy;
    syst_purity_dy.insert(pair<TString,pair<float,float> >("EBEB",pair<float,float>(8.6542e-01,sqrt(pow(4.51e-02,2)+pow(uncertainty_scalefactor_Zee_PIXEL,2)))));
    syst_purity_dy.insert(pair<TString,pair<float,float> >("EBEE",pair<float,float>(7.9537e-01,sqrt(pow(8.22e-02,2)+pow(uncertainty_scalefactor_Zee_PIXEL,2)))));
    syst_purity_dy.insert(pair<TString,pair<float,float> >("EEEE",pair<float,float>(8.5493e-01,sqrt(pow(7.29e-02,2)+pow(uncertainty_scalefactor_Zee_PIXEL,2)))));


    if (!skipsystematics){
      TH1F *systplot = NULL;
      systplot = new TH1F("systplot","Systematic uncertainties",bins_to_run,binsdef);
      systplot->SetLineStyle(kDashed);
      systplot->SetLineWidth(2);
      {
	TString unit = diffvariables_units_list(diffvariable);
	systplot->GetXaxis()->SetTitle(Form("%s %s",diffvariables_names_list(diffvariable).Data(),unit!=TString("") ? (TString("(").Append(unit.Append(")"))).Data() : TString("").Data()));
      }
      systplot_purefitbias=(TH1F*)(systplot->Clone("systplot_purefitbias"));
      systplot_templatestatistics=(TH1F*)(systplot->Clone("systplot_templatestatistics"));
      systplot_templateshapeMCpromptdrivenEB=(TH1F*)(systplot->Clone("systplot_templateshapeMCpromptdrivenEB"));
      systplot_templateshapeMCfakedrivenEB=(TH1F*)(systplot->Clone("systplot_templateshapeMCfakedrivenEB"));
      systplot_templateshapeMCpromptdrivenEE=(TH1F*)(systplot->Clone("systplot_templateshapeMCpromptdrivenEE"));
      systplot_templateshapeMCfakedrivenEE=(TH1F*)(systplot->Clone("systplot_templateshapeMCfakedrivenEE"));
      systplot_templateshape2frag=(TH1F*)(systplot->Clone("systplot_templateshape2frag"));
      systplot_zee=(TH1F*)(systplot->Clone("systplot_zee"));
      systplot_tot=(TH1F*)(systplot->Clone("systplot_tot"));
      systplot_efficiency=(TH1F*)(systplot->Clone("systplot_efficiency"));
      systplot_unfolding=(TH1F*)(systplot->Clone("systplot_unfolding"));
      systplot_totfinal=(TH1F*)(systplot->Clone("systplot_totfinal"));
      systplot_uncorrelated=(TH1F*)(systplot->Clone("systplot_uncorrelated"));
      systplot_1catcorrelated=(TH1F*)(systplot->Clone("systplot_1catcorrelated"));
      systplot_allcatcorrelated=(TH1F*)(systplot->Clone("systplot_allcatcorrelated"));
      systplot_statistic=(TH1F*)(systplot->Clone("systplot_statistic"));
      
    }

    //    std::cout << "start loop: " << std::endl;

    float unfoldingdy_mc_all = (!skipZsubtraction) ? ((RooDataSet*)(file_pixelrevsel_dy_mc->Get(Form("roofit/obs_roodset_%s_%s_b%d",splitting.Data(),diffvariable.Data(),n_bins))))->sumEntries() : 0;

  for (int bin=0; bin<bins_to_run; bin++) {

    //    if (!fr[bin]->fr) continue;

    float pp = 	       purity[0]->GetBinContent(bin+1);
    float pp_err =     purity[0]->GetBinError(bin+1);
//    float pf = 	       purity[1]->GetBinContent(bin+1);
//    float pf_err =     purity[1]->GetBinError(bin+1);
//    float fp = 	       purity[2]->GetBinContent(bin+1);
//    float fp_err =     purity[2]->GetBinError(bin+1);
//    float ff = 	       purity[3]->GetBinContent(bin+1);
//    float ff_err =     purity[3]->GetBinError(bin+1);
    float tot_events = eventshisto->GetBinContent(bin+1);
    float eff_overflow = overflowremovaleffhisto->GetBinContent(bin+1);

    if (!skipZsubtraction) assert((RooDataSet*)(file_standardsel_dy->Get(Form("roofit/obs_roodset_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin))));
    float events_dy = (!skipZsubtraction) ? ((RooDataSet*)(file_standardsel_dy->Get(Form("roofit/obs_roodset_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin))))->sumEntries() : 0; // normalized to 1/fb, xsec normalized to 2475 

    float unfoldingdy_data = (!skipZsubtraction) ? ((RooDataSet*)(file_pixelrevsel_dy_data->Get(Form("roofit/obs_roodset_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin))))->sumEntries()/intlumi : 0;
    float unfoldingdy_mc = (!skipZsubtraction) ? ((RooDataSet*)(file_pixelrevsel_dy_mc->Get(Form("roofit/obs_roodset_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin))))->sumEntries() : 0;

    float purity_dy = syst_purity_dy[splitting].first;
    float purity_dy_err = syst_purity_dy[splitting].second;
    float scale_dy = (!skipZsubtraction) ? unfoldingdy_data/unfoldingdy_mc : 0;

    if (!skipZsubtraction){
      if (unfoldingdy_mc/unfoldingdy_mc_all<0.01) scale_dy=3048./2475.; // fix for low stat bins
    }

    if (splitting=="EBEB") scale_dy*=pow(0.979,2);
    else if (splitting=="EBEE") scale_dy*=0.979*0.985;
    else if (splitting=="EEEE") scale_dy*=pow(0.985,2);
    float rel_error_on_purity_pp = events_dy*purity_dy_err*scale_dy/(pp*tot_events/eff_overflow/intlumi-events_dy*purity_dy*scale_dy);

    std::cout << "bin " << bin << std::endl;
    std::cout << "Data/MC DY " << unfoldingdy_data/unfoldingdy_mc << std::endl;
    std::cout << "UNFOLDING DY FACTOR " << scale_dy << std::endl;
    std::cout << unfoldingdy_data << " " << unfoldingdy_mc << std::endl;
    std::cout << "SUBTRACTION " << events_dy*purity_dy*scale_dy/(pp*tot_events/eff_overflow/intlumi) << std::endl;

    xsec->SetBinContent(bin+1,(pp*tot_events/eff_overflow/intlumi-events_dy*purity_dy*scale_dy)/xsec->GetBinWidth(bin+1));
    xsec_withsyst->SetBinContent(bin+1,xsec->GetBinContent(bin+1));
    xsec_ngammagammayield->SetBinContent(bin+1,pp*tot_events/eff_overflow-events_dy*purity_dy*scale_dy*intlumi);
    std::cout << xsec_ngammagammayield->GetBinContent(bin+1) << std::endl;

    float shapesyst1 = (!skipsystematics && splitting!="EEEE") ? pp*fabs(histo_bias_templateshapeMCpromptdrivenEB->GetBinContent(bin+1)-1) : 0;
    float shapesyst2 = (!skipsystematics && splitting!="EEEE") ? pp*fabs(histo_bias_templateshapeMCfakedrivenEB->GetBinContent(bin+1)-1) : 0;
    float shapesyst3 = (!skipsystematics && splitting!="EBEB") ? pp*fabs(histo_bias_templateshapeMCpromptdrivenEE->GetBinContent(bin+1)-1) : 0;
    float shapesyst4 = (!skipsystematics && splitting!="EBEB") ? pp*fabs(histo_bias_templateshapeMCfakedrivenEE->GetBinContent(bin+1)-1) : 0;
    float shapesyst5 = (!skipsystematics) ? pp*fabs(histo_bias_templateshape2frag->GetBinContent(bin+1)-1) : 0;

    float purity_error_withsyst = pp_err;
    if (!skipsystematics) purity_error_withsyst = sqrt(pow(pp_err,2) + pow(pp*histo_bias_templatestatistics->GetBinContent(bin+1),2) + pow(pp_err*histo_bias_purefitbias->GetBinContent(bin+1),2) + pow(shapesyst1,2) + pow(shapesyst2,2) + pow(shapesyst3,2) + pow(shapesyst4,2) + pow(shapesyst5,2) + pow(pp*rel_error_on_purity_pp,2));


    float errpoiss=1.0/sqrt(tot_events);
    float err=sqrt(pow(pp_err/pp,2)+pow(errpoiss,2));
    float err_withsyst=sqrt(pow(purity_error_withsyst/pp,2)+pow(errpoiss,2));


    xsec->SetBinError(bin+1,err*xsec->GetBinContent(bin+1));
    xsec_withsyst->SetBinError(bin+1,err_withsyst*xsec->GetBinContent(bin+1));
    xsec_ngammagammayield->SetBinError(bin+1,err_withsyst*xsec_ngammagammayield->GetBinContent(bin+1));
    
    if (!skipsystematics){
      systplot_purefitbias->SetBinContent(bin+1,pp_err*fabs(histo_bias_purefitbias->GetBinContent(bin+1))/pp);
      systplot_templatestatistics->SetBinContent(bin+1,histo_bias_templatestatistics->GetBinContent(bin+1));
      systplot_templateshapeMCpromptdrivenEB->SetBinContent(bin+1,shapesyst1/pp);
      systplot_templateshapeMCfakedrivenEB->SetBinContent(bin+1,shapesyst2/pp);
      systplot_templateshapeMCpromptdrivenEE->SetBinContent(bin+1,shapesyst3/pp);
      systplot_templateshapeMCfakedrivenEE->SetBinContent(bin+1,shapesyst4/pp);
      systplot_templateshape2frag->SetBinContent(bin+1,shapesyst5/pp);
      systplot_zee->SetBinContent(bin+1,rel_error_on_purity_pp);
      systplot_tot->SetBinContent(bin+1,sqrt(pow(pp*histo_bias_templatestatistics->GetBinContent(bin+1),2) + pow(pp_err*histo_bias_purefitbias->GetBinContent(bin+1),2) + pow(shapesyst1,2) + pow(shapesyst2,2) + pow(shapesyst3,2) + pow(shapesyst4,2) + pow(shapesyst5,2) + pow(pp*rel_error_on_purity_pp,2))/pp);
      systplot_efficiency->SetBinContent(bin+1,eff->GetBinError(bin+1)/eff->GetBinContent(bin+1));
      systplot_unfolding->SetBinContent(bin+1,unfoldunc->GetBinContent(bin+1));
      systplot_totfinal->SetBinContent(bin+1,sqrt(pow(systplot_tot->GetBinContent(bin+1),2)+pow(systplot_efficiency->GetBinContent(bin+1),2)+pow(systplot_unfolding->GetBinContent(bin+1),2)));
      systplot_totfinal->SetBinError(bin+1,0);
      systplot_statistic->SetBinContent(bin+1,err);
      //      std::cout << systplot_tot->GetBinContent(bin+1) << " " << systplot_totfinal->GetBinContent(bin+1) << std::endl;
    }
    
    //    std::cout << err << " " << err_withsyst << std::endl;

  }


  //  std::cout << "NIENTE DIVIDE(EFF): EFFICIENZA NON CORRETTA!!!" << std::endl;  
  /*
  xsec->Divide(eff);  
  xsec_withsyst->Divide(eff);
  */

  TCanvas *output_canv = new TCanvas("output_canv","output_canv");
  output_canv->cd();

  purity[0]->SetStats(0);
  purity[0]->SetTitle("");
  purity[0]->Draw("e1");

  if (!sym){
    purity[1]->Draw("e1same");
    purity[2]->Draw("e1same");
  }
  else {
    purity[1]->Add(purity[2]);
    purity[1]->Draw("e1same");
  }
  purity[3]->Draw("e1same");
  
  //  purity[0]->GetXaxis()->SetTitle(diffvariables_names_list(diffvariable).Data());
  {
    TString unit = diffvariables_units_list(diffvariable);
    purity[0]->GetXaxis()->SetTitle(Form("%s %s",diffvariables_names_list(diffvariable).Data(),unit!=TString("") ? (TString("(").Append(unit.Append(")"))).Data() : TString("").Data()));
  }

  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  if (sym){
    leg->AddEntry(purity[0],"prompt - prompt","lp");
    leg->AddEntry(purity[1],"prompt - fake","lp");
    leg->AddEntry(purity[3],"fake - fake","lp");
  }
  else {
    leg->AddEntry(purity[0],"prompt - prompt","lp");
    leg->AddEntry(purity[1],Form("prompt %s - fake %s",s1.Data(),s2.Data()),"lp");
    leg->AddEntry(purity[2],Form("fake %s - prompt %s",s1.Data(),s2.Data()),"lp");
    leg->AddEntry(purity[3],"fake - fake","lp");
  }
  leg->SetFillColor(kWhite);
  leg->Draw();

  output_canv->Update();

  output_canv->SaveAs(Form("plots/plot_purity_%s_%s.png", diffvariable.Data(),splitting.Data()));
  output_canv->SaveAs(Form("plots/plot_purity_%s_%s.jpg", diffvariable.Data(),splitting.Data()));
  output_canv->SaveAs(Form("plots/plot_purity_%s_%s.root",diffvariable.Data(),splitting.Data()));
  output_canv->SaveAs(Form("plots/plot_purity_%s_%s.pdf", diffvariable.Data(),splitting.Data()));

  TFile *out1 = new TFile(Form("plots/purityhistos_%s_%s.root",splitting.Data(),diffvariable.Data()),"recreate");
  out1->cd();
  purity[0]->Write();
  purity[1]->Write();
  purity[2]->Write();
  purity[3]->Write();
  out1->Close();

  }



  xsec_withsyst->Print();
  xsec->Print();

  TCanvas *xsec_canv = new TCanvas("xsec_canv","xsec_canv");
  xsec_canv->cd();
  xsec_withsyst->SetStats(0);
  xsec_withsyst->SetTitle(Form("Differential cross section - %s category",splitting.Data()));
  {
    TString unit = diffvariables_units_list(diffvariable);
    xsec_withsyst->GetYaxis()->SetTitle(Form("d#sigma/d%s (fb%s)",diffvariables_names_list(diffvariable).Data(),unit!=TString("") ? (TString("/").Append(unit.Data())).Data() : TString("").Data()));
  }
  xsec_withsyst->SetMinimum(0);
  xsec_withsyst->Draw("e1");
  xsec->Draw("e1same");


  TLegend *legxsec = (diffvariable!="dphi") ? new TLegend(0.7,0.7,0.9,0.9) : new TLegend(0.1,0.7,0.3,0.9);
  legxsec->AddEntry(xsec,"Stat. unc. only","l");
  legxsec->AddEntry(xsec_withsyst,"Stat. + syst. unc.","l");
  legxsec->SetFillColor(kWhite);
  legxsec->Draw();

  xsec_canv->Update();
  
  xsec_canv->SaveAs(Form("plots/plot_xsec_%s_%s.png", diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_%s_%s.jpg", diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_%s_%s.root",diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_%s_%s.pdf", diffvariable.Data(),splitting.Data()));
  
  xsec_canv->SetLogy();
  xsec->GetYaxis()->UnZoom();
  xsec_canv->SaveAs(Form("plots/plot_xsec_log_%s_%s.png", diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_log_%s_%s.jpg", diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_log_%s_%s.root",diffvariable.Data(),splitting.Data()));
  xsec_canv->SaveAs(Form("plots/plot_xsec_log_%s_%s.pdf", diffvariable.Data(),splitting.Data()));


  TFile *xsec_file = new TFile(Form("plots/histo_xsec_%s_%s.root", diffvariable.Data(),splitting.Data()),"recreate");
  xsec_file->cd();
  xsec->Write();
  xsec_withsyst->Write();
  if (xsec_ngammagammayield) xsec_ngammagammayield->Write();
  fitpull_histo->Write();
  xsec_file->Close();  


  TCanvas *fitpull_canv = new TCanvas("fitpull_canv","fitpull_canv");
  fitpull_canv->cd();

  gStyle->SetOptFit(1);

  fitpull_histo->SetStats(1);
  fitpull_histo->SetTitle(Form("Distribution of the fit pull - %s category",splitting.Data()));
  fitpull_histo->Fit("gaus");
  fitpull_histo->Draw();
  fitpull_canv->SaveAs(Form("plots/plot_fitpull_%s_%s.root", diffvariable.Data(),splitting.Data()));
  fitpull_canv->SaveAs(Form("plots/plot_fitpull_%s_%s.pdf", diffvariable.Data(),splitting.Data()));
  fitpull_canv->SaveAs(Form("plots/plot_fitpull_%s_%s.png", diffvariable.Data(),splitting.Data()));

  gStyle->SetOptFit(0);

  if (splitting=="inclusive"){

    TFile *out[3];
    out[0] = new TFile(Form("plots/histo_purity_%s_%s_allbins.root",diffvariable.Data(),"EBEB"));
    out[1] = new TFile(Form("plots/histo_purity_%s_%s_allbins.root",diffvariable.Data(),"EBEE"));
    out[2] = new TFile(Form("plots/histo_purity_%s_%s_allbins.root",diffvariable.Data(),"EEEE"));

    TH1F *h[3][4];
    TH1F *n[3];

    for (int i=0; i<3; i++){
      out[i]->GetObject("purity_sigsig",h[i][0]);
      out[i]->GetObject("purity_sigbkg",h[i][1]);
      out[i]->GetObject("purity_bkgsig",h[i][2]);
      out[i]->GetObject("purity_bkgbkg",h[i][3]);
      out[i]->GetObject("eventshisto",n[i]);
    }

    for (int j=0; j<4; j++) for (int i=0; i<3; i++) h[i][j]->Multiply(n[i]);
    for (int j=0; j<4; j++) for (int i=1; i<3; i++) h[0][j]->Add(h[i][j]);
    for (int i=1; i<3; i++) n[0]->Add(n[i]);  
    for (int j=0; j<4; j++) h[0][j]->Divide(n[0]);
    h[0][1]->Add(h[0][2]);

    TCanvas *canv = new TCanvas("purity_inclusive");
    canv->cd();
    h[0][0]->SetTitle("");
    h[0][0]->SetStats(0);
    h[0][0]->SetMinimum(0);
    h[0][0]->SetMaximum(1);
    {
      TString unit = diffvariables_units_list(diffvariable);
      h[0][0]->GetXaxis()->SetTitle(Form("%s %s",diffvariables_names_list(diffvariable).Data(),unit!=TString("") ? (TString("(").Append(unit.Append(")"))).Data() : TString("").Data()));
      h[0][0]->GetYaxis()->SetTitle("Purity fraction");
    }
    h[0][0]->GetXaxis()->SetLabelSize(0.038);
    h[0][0]->GetXaxis()->SetTitleSize(0.038);
    h[0][0]->GetXaxis()->SetTitleOffset(1.15);
    h[0][0]->Draw();
    h[0][3]->Draw("same");
    h[0][1]->Draw("same");
    h[0][0]->Draw("same");
    TLegend *leg = new TLegend(0.55,0.7,0.9,0.9);
    leg->AddEntry(h[0][0],"prompt - prompt","lp");
    leg->AddEntry(h[0][1],"prompt - fake","lp");
    leg->AddEntry(h[0][3],"fake - fake","lp");
    leg->SetFillColor(kWhite);
    leg->Draw();
    TLatex a;
    a.SetNDC();
    a.SetTextSize(0.03);
    a.DrawLatex(0.13,0.83,"#splitline{CMS Preliminary}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");
    canv->Update();
    canv->SaveAs(Form("plots/plot_purity_%s_%s.root",diffvariable.Data(),"inclusive"));
    canv->SaveAs(Form("plots/plot_purity_%s_%s.pdf",diffvariable.Data(),"inclusive"));
    canv->SaveAs(Form("plots/plot_purity_%s_%s.png",diffvariable.Data(),"inclusive"));


  }

  if (!skipsystematics && splitting!="inclusive"){

    TCanvas *systplot_canv = new TCanvas("systplot_canv","systplot_canv");
    systplot_canv->cd();
    systplot_tot->SetStats(0);
    systplot_tot->SetTitle(Form("Systematic uncertainties on purity - %s category",splitting.Data()));
    systplot_tot->SetMinimum(0);

    systplot_tot->Draw();    
    systplot_templateshapeMCpromptdrivenEB->Draw("same");
    systplot_templateshapeMCfakedrivenEB->Draw("same");
    systplot_templateshapeMCpromptdrivenEE->Draw("same");
    systplot_templateshapeMCfakedrivenEE->Draw("same");
    systplot_templateshape2frag->Draw("same");
    systplot_purefitbias->Draw("same");
    systplot_templatestatistics->Draw("same");
    systplot_zee->Draw("same");

    systplot_templateshapeMCpromptdrivenEB->SetLineColor(kRed);
    systplot_templateshapeMCfakedrivenEB->SetLineColor(kBlue);
    systplot_templateshapeMCpromptdrivenEE->SetLineColor(kRed);
    systplot_templateshapeMCfakedrivenEE->SetLineColor(kBlue);
    systplot_templateshapeMCpromptdrivenEE->SetLineStyle(kDotted);
    systplot_templateshapeMCfakedrivenEE->SetLineStyle(kDotted);
    systplot_templateshape2frag->SetLineColor(kOrange);
    systplot_purefitbias->SetLineColor(kGreen);
    systplot_templatestatistics->SetLineColor(kGray);
    systplot_zee->SetLineColor(kMagenta);
    systplot_tot->SetLineColor(kBlack);

    TLegend *legsystplot = new TLegend(0.6,0.7,0.9,0.9);
    if (splitting!="EEEE"){
      legsystplot->AddEntry(systplot_templateshapeMCpromptdrivenEB,"Prompt template shape EB","l");
      legsystplot->AddEntry(systplot_templateshapeMCfakedrivenEB,"Fakes template shape EB","l");
    }
    if (splitting!="EBEB"){
      legsystplot->AddEntry(systplot_templateshapeMCpromptdrivenEE,"Prompt template shape EE","l");
      legsystplot->AddEntry(systplot_templateshapeMCfakedrivenEE,"Fakes template shape EE","l");
    }
    legsystplot->AddEntry(systplot_templateshape2frag,"Fragmentation description","l");
    legsystplot->AddEntry(systplot_templatestatistics,"Template stat. fluctuation","l");
    legsystplot->AddEntry(systplot_purefitbias,"Fit bias","l");
    legsystplot->AddEntry(systplot_zee,"Zee subtraction ","l");
    legsystplot->AddEntry(systplot_tot,"Total syst. uncertainty","l");
    legsystplot->SetFillColor(kWhite);
    legsystplot->Draw();

    systplot_canv->Update();
    
    systplot_canv->SaveAs(Form("plots/plot_systsummary_%s_%s.png", diffvariable.Data(),splitting.Data()));
    systplot_canv->SaveAs(Form("plots/plot_systsummary_%s_%s.jpg", diffvariable.Data(),splitting.Data()));
    systplot_canv->SaveAs(Form("plots/plot_systsummary_%s_%s.root",diffvariable.Data(),splitting.Data()));
    systplot_canv->SaveAs(Form("plots/plot_systsummary_%s_%s.pdf", diffvariable.Data(),splitting.Data()));

    systplot_tot->SaveAs(Form("plots/histo_systsummary_%s_%s.root", diffvariable.Data(),splitting.Data()));

    TCanvas *systplot_canv2 = new TCanvas("systplot_canv2","systplot_canv2");
    systplot_canv2->cd();
    systplot_totfinal->SetStats(0);
    systplot_totfinal->SetTitle(Form("Systematic uncertainties on cross-section - %s category",splitting.Data()));
    systplot_totfinal->SetMinimum(0);

    systplot_totfinal->Draw();    
    systplot_tot->Draw("same");
    systplot_efficiency->Draw("same");
    systplot_unfolding->Draw("same");

    systplot_tot->SetLineColor(kRed);
    systplot_efficiency->SetLineColor(kRed-7);
    systplot_unfolding->SetLineColor(kCyan);
    systplot_totfinal->SetLineColor(kBlack);

    TLegend *legsystplot2 = new TLegend(0.6,0.7,0.9,0.9);
    legsystplot2->AddEntry(systplot_tot,"Purity uncertainty","l");
    legsystplot2->AddEntry(systplot_efficiency,"Efficiency uncertainty","l");
    legsystplot2->AddEntry(systplot_unfolding,"Unfolding uncertainty","l");
    legsystplot2->AddEntry(systplot_totfinal,"Total syst. uncertainty","l");
    legsystplot2->SetFillColor(kWhite);
    legsystplot2->Draw();

    systplot_canv2->Update();
    
    systplot_canv2->SaveAs(Form("plots/plot_systsummaryfinal_%s_%s.png", diffvariable.Data(),splitting.Data()));
    systplot_canv2->SaveAs(Form("plots/plot_systsummaryfinal_%s_%s.jpg", diffvariable.Data(),splitting.Data()));
    systplot_canv2->SaveAs(Form("plots/plot_systsummaryfinal_%s_%s.root",diffvariable.Data(),splitting.Data()));
    systplot_canv2->SaveAs(Form("plots/plot_systsummaryfinal_%s_%s.pdf", diffvariable.Data(),splitting.Data()));

    TFile *f2 = new TFile(Form("plots/histo_systsummaryfinal_%s_%s.root", diffvariable.Data(),splitting.Data()),"recreate");
    f2->cd();

    systplot_templateshapeMCfakedrivenEB->Write();
    systplot_templateshapeMCpromptdrivenEB->Write();
    systplot_templateshapeMCfakedrivenEE->Write();
    systplot_templateshapeMCpromptdrivenEE->Write();
    systplot_templateshape2frag->Write();
    systplot_purefitbias->Write();
    systplot_templatestatistics->Write();
    systplot_zee->Write();
    systplot_tot->Write();
    systplot_efficiency->Write();
    systplot_unfolding->Write();
    systplot_totfinal->Write();
    systplot_statistic->Write();

    std::vector<TH1F*> toadd_1catcorrelated;
    toadd_1catcorrelated.push_back(systplot_zee);
    toadd_1catcorrelated.push_back(systplot_templatestatistics);
    toadd_1catcorrelated.push_back(systplot_efficiency);
    toadd_1catcorrelated.push_back(systplot_unfolding);

    systplot_1catcorrelated=AddTHInQuadrature(toadd_1catcorrelated,"systplot_1catcorrelated");

    std::vector<TH1F*> toadd_allcatcorrelated;
    toadd_allcatcorrelated.push_back(systplot_templateshapeMCfakedrivenEB);
    toadd_allcatcorrelated.push_back(systplot_templateshapeMCpromptdrivenEB);
    toadd_allcatcorrelated.push_back(systplot_templateshapeMCfakedrivenEE);
    toadd_allcatcorrelated.push_back(systplot_templateshapeMCpromptdrivenEE);
    toadd_allcatcorrelated.push_back(systplot_templateshape2frag);

    systplot_allcatcorrelated=AddTHInQuadrature(toadd_allcatcorrelated,"systplot_allcatcorrelated");

    std::vector<TH1F*> toadd_uncorrelated;
    toadd_uncorrelated.push_back(systplot_purefitbias);

    systplot_uncorrelated=AddTHInQuadrature(toadd_uncorrelated,"systplot_uncorrelated");

    systplot_1catcorrelated->Write();
    systplot_allcatcorrelated->Write();
    systplot_uncorrelated->Write();


    f2->Close();

  }

  if (!skipsystematics && splitting=="inclusive"){

    TH1F *eff[3];
    {
      TFile *eff_file = new TFile("plots/Efficiency_WithSysErr.root");
      std::map<TString,TString> translation;
      translation.insert(std::pair<TString,TString>(TString("invmass"),TString("mgg")));
      translation.insert(std::pair<TString,TString>(TString("diphotonpt"),TString("qtgg")));
      translation.insert(std::pair<TString,TString>(TString("costhetastar"),TString("costhetastar")));
      translation.insert(std::pair<TString,TString>(TString("dphi"),TString("deltaphi")));
      translation.insert(std::pair<TString,TString>(TString("dR"),TString("dR")));
      eff_file->GetObject(Form("h_%s_%s_WithTotErr",translation[diffvariable].Data(),"EBEB"),eff[0]);
      eff_file->GetObject(Form("h_%s_%s_WithTotErr",translation[diffvariable].Data(),"EBEE"),eff[1]);
      eff_file->GetObject(Form("h_%s_%s_WithTotErr",translation[diffvariable].Data(),"EEEE"),eff[2]);
    }


    TFile *fsysts[3];
    fsysts[0] = new TFile(Form("plots/histo_systsummaryfinal_%s_EBEB.root", diffvariable.Data()));;
    fsysts[1] = new TFile(Form("plots/histo_systsummaryfinal_%s_EBEE.root", diffvariable.Data()));;
    fsysts[2] = new TFile(Form("plots/histo_systsummaryfinal_%s_EEEE.root", diffvariable.Data()));;

    const int n_cats = 3;
    const int n_syst_1catcorr = 4;
    const int n_syst_allcatcorr = 5;

    TH1F *hsysts[3];
    TH1F *hsysts_1catcorrelated[3][n_syst_1catcorr];
    TH1F *hsysts_allcatcorrelated[3][n_syst_allcatcorr];
    TH1F *hsysts_uncorrelated[3];
    TH1F *hsysts_statistic[3];

    for (int i=0; i<3; i++) fsysts[i]->GetObject("systplot_totfinal",hsysts[i]);
    for (int i=0; i<3; i++) {
      fsysts[i]->GetObject("systplot_zee",hsysts_1catcorrelated[i][0]);
      fsysts[i]->GetObject("systplot_templatestatistics",hsysts_1catcorrelated[i][1]);
      fsysts[i]->GetObject("systplot_efficiency",hsysts_1catcorrelated[i][2]);
      fsysts[i]->GetObject("systplot_unfolding",hsysts_1catcorrelated[i][3]);
    }
    for (int i=0; i<3; i++) {
      fsysts[i]->GetObject("systplot_templateshapeMCpromptdrivenEB",hsysts_allcatcorrelated[i][0]);
      fsysts[i]->GetObject("systplot_templateshapeMCfakedrivenEB",hsysts_allcatcorrelated[i][1]);
      fsysts[i]->GetObject("systplot_templateshapeMCpromptdrivenEE",hsysts_allcatcorrelated[i][2]);
      fsysts[i]->GetObject("systplot_templateshapeMCfakedrivenEE",hsysts_allcatcorrelated[i][3]);
      fsysts[i]->GetObject("systplot_templateshape2frag",hsysts_allcatcorrelated[i][4]);
    }
    for (int i=0; i<3; i++) fsysts[i]->GetObject("systplot_uncorrelated",hsysts_uncorrelated[i]);
    for (int i=0; i<3; i++) fsysts[i]->GetObject("systplot_statistic",hsysts_statistic[i]);
    for (int i=0; i<3; i++) hsysts[i]->Sumw2();

    TH1F *unfoldnevt[3];
    std::map<TString,TString> translation2;
    translation2.insert(std::pair<TString,TString>(TString("invmass"),TString("mgg")));
    translation2.insert(std::pair<TString,TString>(TString("diphotonpt"),TString("pt")));
    translation2.insert(std::pair<TString,TString>(TString("costhetastar"),TString("costt")));
    translation2.insert(std::pair<TString,TString>(TString("dphi"),TString("phi")));
    translation2.insert(std::pair<TString,TString>(TString("dR"),TString("dR")));
    TFile *unfoldnevt_file = new TFile("plots/Unfolding_SysErr.root");
    unfoldnevt_file->GetObject(Form("Unfolding_Nevt_%s_EBEB",translation2[diffvariable].Data()),unfoldnevt[0]);
    unfoldnevt_file->GetObject(Form("Unfolding_Nevt_%s_EBEE",translation2[diffvariable].Data()),unfoldnevt[1]);
    unfoldnevt_file->GetObject(Form("Unfolding_Nevt_%s_EEEE",translation2[diffvariable].Data()),unfoldnevt[2]);

    TH1F *unfoldnevt_new[3];
    for (int i=0; i<3; i++) {
      unfoldnevt_new[i]=(TH1F*)(hsysts[i]->Clone(Form("unfoldnevt_new_%d",i)));
      unfoldnevt_new[i]->Reset();
      unfoldnevt_new[i]->Sumw2();
      for (int bin=0; bin<bins_to_run; bin++) {
	unfoldnevt_new[i]->SetBinContent(bin+1,unfoldnevt[i]->GetBinContent(bin+1)/eff[i]->GetBinContent(bin+1));
	unfoldnevt_new[i]->SetBinError(bin+1,unfoldnevt[i]->GetBinError(bin+1)/eff[i]->GetBinContent(bin+1));
      }
    }

    TH1F *unfoldnevt_tot =(TH1F*)(unfoldnevt_new[0]->Clone("unfoldnevt_tot"));
    unfoldnevt_tot->Reset();
    unfoldnevt_tot->Sumw2();
//    TH1F *hsyst_tot = (TH1F*)(hsysts[0]->Clone("hsyst_tot"));
//    hsyst_tot->Reset();
//    hsyst_tot->Sumw2();
    for (int i=0; i<3; i++) {
//      hsysts[i]->Multiply(unfoldnevt_new[i]);
//      hsyst_tot->Add(hsysts[i]);
      unfoldnevt_tot->Add(unfoldnevt_new[i]);
    }
    //    hsyst_tot->Divide(unfoldnevt_tot);

    TH1F* systplot_totfinal_inclusive = (TH1F*)(hsysts[0]->Clone("systplot_totfinal_inclusive"));
    systplot_totfinal_inclusive->Reset();
    systplot_totfinal_inclusive->Sumw2();
    systplot_totfinal_inclusive->SetLineStyle(kDashed);
    systplot_totfinal_inclusive->SetMinimum(0);
    systplot_totfinal_inclusive->SetTitle("Total systematic uncertainty on cross-section");

    for (int i=0; i<3; i++) {
      hsysts_statistic[i]->Multiply(unfoldnevt_new[i]);
      hsysts_uncorrelated[i]->Multiply(unfoldnevt_new[i]);
      for (int k=0; k<n_syst_1catcorr; k++) hsysts_1catcorrelated[i][k]->Multiply(unfoldnevt_new[i]);
      for (int k=0; k<n_syst_allcatcorr; k++) hsysts_allcatcorrelated[i][k]->Multiply(unfoldnevt_new[i]);
    }

    
    float a[n_cats][n_bins];
    float b[n_cats][n_bins];
    float c[n_cats][n_syst_1catcorr][n_bins];
    float d[n_cats][n_syst_allcatcorr][n_bins];


    for (int i=0; i<n_cats; i++){
      for (int bin = 0; bin<bins_to_run; bin++){
	a[i][bin]=hsysts_statistic[i]->GetBinContent(bin+1);
	b[i][bin]=hsysts_uncorrelated[i]->GetBinContent(bin+1);
	for (int k=0; k<n_syst_1catcorr; k++) c[i][k][bin]=hsysts_1catcorrelated[i][k]->GetBinContent(bin+1);
	for (int k=0; k<n_syst_allcatcorr; k++) d[i][k][bin]=hsysts_allcatcorrelated[i][k]->GetBinContent(bin+1);
      }
    }



    float stat1bin[n_cats][bins_to_run];
    for (int i=0; i<n_cats; i++){
      for (int bin = 0; bin<bins_to_run; bin++){
	stat1bin[i][bin]=a[i][bin];
      }
    }

    float syst1bin[n_cats][bins_to_run];
    for (int i=0; i<n_cats; i++){
      for (int bin = 0; bin<bins_to_run; bin++){
	float res = 0;
	res+=pow(b[i][bin],2);
	for (int k=0; k<n_syst_1catcorr; k++) res+=pow(c[i][k][bin],2);
	for (int k=0; k<n_syst_allcatcorr; k++) res+=pow(d[i][k][bin],2);
	syst1bin[i][bin]=sqrt(res);
      }
    }




    float statcategory[n_cats];
    for (int i=0; i<n_cats; i++){
      float res = 0;
      for (int bin = 0; bin<bins_to_run; bin++){
	res+=pow(a[i][bin],2);
      }
      statcategory[i]=sqrt(res);
    }

    float systcategory[n_cats];
    for (int i=0; i<n_cats; i++){
      float res = 0;
      float sumc[n_syst_1catcorr];
      float sumd[n_syst_allcatcorr];
      for (int bin = 0; bin<bins_to_run; bin++){
	res+=pow(b[i][bin],2);
      }
      for (int k=0; k<n_syst_1catcorr; k++){
	sumc[k]=0;
	for (int bin = 0; bin<bins_to_run; bin++) sumc[k]+=c[i][k][bin];
      }
      for (int k=0; k<n_syst_allcatcorr; k++){
	sumd[k]=0;
	for (int bin = 0; bin<bins_to_run; bin++) sumd[k]+=d[i][k][bin];
      }
      for (int k=0; k<n_syst_1catcorr; k++) res+=pow(sumc[k],2);
      for (int k=0; k<n_syst_allcatcorr; k++) res+=pow(sumd[k],2);
      systcategory[i]=sqrt(res);
    }


    float statcolumn[bins_to_run];
    for (int bin = 0; bin<bins_to_run; bin++){
      float res = 0;
      for (int i=0; i<n_cats; i++){
	res+=pow(a[i][bin],2);
      }
      statcolumn[bin]=sqrt(res);
    }

    float systcolumn[bins_to_run];
    for (int bin = 0; bin<bins_to_run; bin++){
      float res = 0;
      for (int i=0; i<n_cats; i++){
        res+=pow(b[i][bin],2);
      }
      for (int i=0; i<n_cats; i++){
	for (int k=0; k<n_syst_1catcorr; k++) res+=pow(c[i][k][bin],2);
      }
      for (int k=0; k<n_syst_allcatcorr; k++) {
	float sumd_column=0;
	for (int i=0; i<n_cats; i++) sumd_column+=d[i][k][bin];
	res+=pow(sumd_column,2);
      }
      systcolumn[bin]=sqrt(res);
    }

						   
    float statall = 0;
    {
      float res = 0;
      for (int i=0; i<n_cats; i++){
	for (int bin = 0; bin<bins_to_run; bin++){
	  res+=pow(a[i][bin],2);
	}
      }
      statall=sqrt(res);
    }

    float systall = 0;
    {
      float res = 0;
      float sumc[n_cats][n_syst_1catcorr];
      float sumd[n_cats][n_syst_allcatcorr];
      for (int i=0; i<n_cats; i++){
	for (int bin = 0; bin<bins_to_run; bin++){
	  res+=pow(b[i][bin],2);
	}
	for (int k=0; k<n_syst_1catcorr; k++){
	  sumc[i][k]=0;
	  for (int bin = 0; bin<bins_to_run; bin++) sumc[i][k]+=c[i][k][bin];
	}
	for (int k=0; k<n_syst_allcatcorr; k++){
	  sumd[i][k]=0;
	  for (int bin = 0; bin<bins_to_run; bin++) sumd[i][k]+=d[i][k][bin];
	}
      }
      for (int i=0; i<n_cats; i++) for (int k=0; k<n_syst_1catcorr; k++) res+=pow(sumc[i][k],2);
      for (int k=0; k<n_syst_allcatcorr; k++){
	float sumd_column=0;
	for (int i=0; i<n_cats; i++) sumd_column+=sumd[i][k];
	res+=pow(sumd_column,2);
      }
      systall=sqrt(res);
    }

    for(int bin = 0; bin<bins_to_run; bin++) {
      //      float val = 100*sqrt(pow(systcolumn[bin],2)+pow(statcolumn[bin],2))/unfoldnevt_tot->GetBinContent(bin+1);
      //      (val>=10) ? printf("%.1f\\%%\n\n\n",val) : printf("%.2f\\%%\n\n\n",val);
      systplot_totfinal_inclusive->SetBinContent(bin+1,systcolumn[bin]/unfoldnevt_tot->GetBinContent(bin+1));
      systplot_totfinal_inclusive->SetBinError(bin+1,0);
    }

//    for(int bin = 0; bin<bins_to_run; bin++) {
//      float val = 100*sqrt(pow(systcolumn[bin],2)+pow(statcolumn[bin],2))/unfoldnevt_tot->GetBinContent(bin+1);
//      (val>=10) ? printf("%.1f\\%%\n\n\n",val) : printf("%.2f\\%%\n\n\n",val);
//    }


    for (int i=0; i<n_cats; i++){
      std::cout << std::endl;
      std::cout << "---" << std::endl;
      std::cout << std::endl;

      std::cout << "Category " << i << ":" << std::endl;
      std::cout << std::endl;
      std::cout << "CENTRAL VALUE" << std::endl;
      std::cout << unfoldnevt_new[i]->Integral()/intlumi/1e3 << " pb" << std::endl;
      std::cout << std::endl;
      
      std::cout << "STAT UNCERTAINTY" << std::endl;
      std::cout << statcategory[i]/unfoldnevt_new[i]->Integral() << " relative" << std::endl;
      std::cout << statcategory[i]/intlumi/1e3 << " pb" << std::endl;
      std::cout << std::endl;


      std::cout << "SYST UNCERTAINTY" << std::endl;
      std::cout << systcategory[i]/unfoldnevt_new[i]->Integral() << " relative" << std::endl;
      std::cout << systcategory[i]/intlumi/1e3 << " pb" << std::endl;
      std::cout << std::endl;

      std::cout << "LUMI UNCERTAINTY" << std::endl;
      float lumi_rel = 2.2e-2;
      std::cout << lumi_rel << " relative" << std::endl;
      std::cout << lumi_rel*unfoldnevt_new[i]->Integral()/intlumi/1e3 << " pb" << std::endl;
      std::cout << std::endl;

      std::cout << "TOTAL UNCERTAINTY" << std::endl;
      float e_stat = statcategory[i]/unfoldnevt_new[i]->Integral();
      float e_syst = systcategory[i]/unfoldnevt_new[i]->Integral();
      float e_lumi = lumi_rel;
      std::cout << sqrt(pow(e_stat,2)+pow(e_syst,2)+pow(e_lumi,2)) << " relative" << std::endl;
      std::cout << sqrt(pow(e_stat,2)+pow(e_syst,2)+pow(e_lumi,2))*unfoldnevt_new[i]->Integral()/intlumi/1e3 << " pb" << std::endl;
      std::cout << std::endl;
      
    }


    std::cout << std::endl;
    std::cout << "---" << std::endl;
    std::cout << std::endl;

    std::cout << "Full acceptance:" << std::endl;
    std::cout << std::endl;
    std::cout << "CENTRAL VALUE" << std::endl;
    std::cout << unfoldnevt_tot->Integral()/intlumi/1e3 << " pb" << std::endl;
    std::cout << std::endl;

    std::cout << "STAT UNCERTAINTY" << std::endl;
    std::cout << statall/unfoldnevt_tot->Integral() << " relative" << std::endl;
    std::cout << statall/intlumi/1e3 << " pb" << std::endl;
    std::cout << std::endl;


    std::cout << "SYST UNCERTAINTY" << std::endl;
    std::cout << systall/unfoldnevt_tot->Integral() << " relative" << std::endl;
    std::cout << systall/intlumi/1e3 << " pb" << std::endl;
    std::cout << std::endl;

    std::cout << "LUMI UNCERTAINTY" << std::endl;
    float lumi_rel = 2.2e-2;
    std::cout << lumi_rel << " relative" << std::endl;
    std::cout << lumi_rel*unfoldnevt_tot->Integral()/intlumi/1e3 << " pb" << std::endl;
    std::cout << std::endl;

    std::cout << "TOTAL UNCERTAINTY" << std::endl;
    float e_stat = statall/unfoldnevt_tot->Integral();
    float e_syst = systall/unfoldnevt_tot->Integral();
    float e_lumi = lumi_rel;
    std::cout << sqrt(pow(e_stat,2)+pow(e_syst,2)+pow(e_lumi,2)) << " relative" << std::endl;
    std::cout << sqrt(pow(e_stat,2)+pow(e_syst,2)+pow(e_lumi,2))*unfoldnevt_tot->Integral()/intlumi/1e3 << " pb" << std::endl;
    std::cout << std::endl;

    TCanvas *canv3 = new TCanvas();
    canv3->cd();
    systplot_totfinal_inclusive->Draw();
    canv3->SaveAs(Form("plots/histo_systsummaryfinal_%s_inclusive.jpg", diffvariable.Data()));
    canv3->SaveAs(Form("plots/histo_systsummaryfinal_%s_inclusive.png", diffvariable.Data()));
    canv3->SaveAs(Form("plots/histo_systsummaryfinal_%s_inclusive.root", diffvariable.Data()));
    canv3->SaveAs(Form("plots/histo_systsummaryfinal_%s_inclusive.pdf", diffvariable.Data()));

    TH1F *histo_uncorrelated_allcat;
    {
      std::vector<TH1F*> toadd_uncorrelated_allcat;
      for (int i=0; i<n_cats; i++) toadd_uncorrelated_allcat.push_back(hsysts_uncorrelated[i]);
      histo_uncorrelated_allcat = AddTHInQuadrature(toadd_uncorrelated_allcat,Form("%s_allcat",hsysts_uncorrelated[0]->GetName()));  
      histo_uncorrelated_allcat->Divide(unfoldnevt_tot);
      for (int bin=0; bin<bins_to_run; bin++) histo_uncorrelated_allcat->SetBinError(bin+1,0);
    }
    TH1F *histo_1catcorrelated_allcat[n_syst_1catcorr];
    for (int k=0; k<n_syst_1catcorr; k++){
      std::vector<TH1F*> toadd_1catcorrelated_allcat;
      for (int i=0; i<n_cats; i++) toadd_1catcorrelated_allcat.push_back(hsysts_1catcorrelated[i][k]);
      histo_1catcorrelated_allcat[k] = AddTHInQuadrature(toadd_1catcorrelated_allcat,Form("%s_allcat",hsysts_1catcorrelated[0][k]->GetName()));
      histo_1catcorrelated_allcat[k]->Divide(unfoldnevt_tot);
      for (int bin=0; bin<bins_to_run; bin++) histo_1catcorrelated_allcat[k]->SetBinError(bin+1,0);
    }
    TH1F *histo_allcatcorrelated_allcat[n_syst_allcatcorr];
    for (int k=0; k<n_syst_allcatcorr; k++){
      histo_allcatcorrelated_allcat[k] = (TH1F*)(hsysts_allcatcorrelated[0][k]->Clone(Form("%s_allcat",hsysts_allcatcorrelated[0][k]->GetName())));
      histo_allcatcorrelated_allcat[k]->Reset();
      for (int i=0; i<n_cats; i++) histo_allcatcorrelated_allcat[k]->Add(hsysts_allcatcorrelated[i][k]);
      histo_allcatcorrelated_allcat[k]->Divide(unfoldnevt_tot);
      for (int bin=0; bin<bins_to_run; bin++) histo_allcatcorrelated_allcat[k]->SetBinError(bin+1,0);
    }

    TCanvas *canv3b = new TCanvas();
    canv3b->cd();
    systplot_totfinal_inclusive->SetMinimum(0);
    systplot_totfinal_inclusive->Draw();
    histo_uncorrelated_allcat->Draw("same");
    for (int k=0; k<n_syst_1catcorr; k++) histo_1catcorrelated_allcat[k]->Draw("same");
    for (int k=0; k<n_syst_allcatcorr; k++) histo_allcatcorrelated_allcat[k]->Draw("same");
    TLegend *leg_canv3b = new TLegend(0.6,0.7,0.9,0.9);
    leg_canv3b->AddEntry(histo_uncorrelated_allcat,"Fit bias","l");
    leg_canv3b->AddEntry(histo_1catcorrelated_allcat[0],"Zee subtraction","l");
    leg_canv3b->AddEntry(histo_1catcorrelated_allcat[1],"Template stat. fluctuation","l");
    leg_canv3b->AddEntry(histo_1catcorrelated_allcat[2],"Efficiency uncertainty","l");
    leg_canv3b->AddEntry(histo_1catcorrelated_allcat[3],"Unfolding uncertainty","l");
    leg_canv3b->AddEntry(histo_allcatcorrelated_allcat[0],"Prompt template shape EB","l");
    leg_canv3b->AddEntry(histo_allcatcorrelated_allcat[1],"Fakes template shape EB","l");
    leg_canv3b->AddEntry(histo_allcatcorrelated_allcat[2],"Prompt template shape EE","l");
    leg_canv3b->AddEntry(histo_allcatcorrelated_allcat[3],"Fakes template shape EE","l");
    leg_canv3b->AddEntry(histo_allcatcorrelated_allcat[4],"Fragmentation description","l");
    leg_canv3b->AddEntry(systplot_totfinal_inclusive,"Total syst. uncertainty","l");
    leg_canv3b->Draw();
    systplot_totfinal_inclusive->GetYaxis()->SetRangeUser(0,systplot_totfinal_inclusive->GetBinContent(systplot_totfinal_inclusive->GetMaximumBin())*1.05);
    canv3b->Update();

    canv3b->SaveAs(Form("plots/histo_systsummaryfinal_splitted_%s_inclusive.jpg", diffvariable.Data()));
    canv3b->SaveAs(Form("plots/histo_systsummaryfinal_splitted_%s_inclusive.png", diffvariable.Data()));
    canv3b->SaveAs(Form("plots/histo_systsummaryfinal_splitted_%s_inclusive.root", diffvariable.Data()));
    canv3b->SaveAs(Form("plots/histo_systsummaryfinal_splitted_%s_inclusive.pdf", diffvariable.Data()));



      

    TH1F *histo_finalxs_fortheorycomp = (TH1F*)(systplot_totfinal_inclusive->Clone(Form("histo_finalxs_fortheorycomp_%s",diffvariable.Data())));
    histo_finalxs_fortheorycomp->SetTitle("Cross section (unfolding+efficiency) (stat.+syst.+lumi.)");
    histo_finalxs_fortheorycomp->Reset();
    histo_finalxs_fortheorycomp->GetYaxis()->UnZoom();
    histo_finalxs_fortheorycomp->SetLineStyle(1);
    for (int bin=0; bin<bins_to_run; bin++){
      float xs = unfoldnevt_tot->GetBinContent(bin+1)/intlumi/1e3/unfoldnevt_tot->GetBinWidth(bin+1);
      float relerr = sqrt(pow(systcolumn[bin]/unfoldnevt_tot->GetBinContent(bin+1),2)+pow(statcolumn[bin]/unfoldnevt_tot->GetBinContent(bin+1),2)+pow(lumi_rel,2));
      histo_finalxs_fortheorycomp->SetBinContent(bin+1,xs);
      histo_finalxs_fortheorycomp->SetBinError(bin+1,relerr*xs);
      std::cout << "Bin " << bin << " " << xs << " +/- " << 100*relerr << " % (stat.+syst.+lumi.)" << std::endl;
    }
    TCanvas *canv4 = new TCanvas();
    canv4->cd();
    histo_finalxs_fortheorycomp->Draw("e1");

    histo_finalxs_fortheorycomp->SaveAs(Form("plots/%s.root",histo_finalxs_fortheorycomp->GetName()));

    TH1F *histo_finalxsnolumi_fortheorycomp = (TH1F*)(systplot_totfinal_inclusive->Clone(Form("histo_finalxsnolumi_fortheorycomp_%s",diffvariable.Data())));
    histo_finalxsnolumi_fortheorycomp->SetTitle("Cross section (unfolding+efficiency) (stat.+syst.)");
    histo_finalxsnolumi_fortheorycomp->Reset();
    histo_finalxsnolumi_fortheorycomp->GetYaxis()->UnZoom();
    histo_finalxsnolumi_fortheorycomp->SetLineStyle(1);
    for (int bin=0; bin<bins_to_run; bin++){
      float xs = unfoldnevt_tot->GetBinContent(bin+1)/intlumi/1e3/unfoldnevt_tot->GetBinWidth(bin+1);
      float relerr = sqrt(pow(systcolumn[bin]/unfoldnevt_tot->GetBinContent(bin+1),2)+pow(statcolumn[bin]/unfoldnevt_tot->GetBinContent(bin+1),2));
      histo_finalxsnolumi_fortheorycomp->SetBinContent(bin+1,xs);
      histo_finalxsnolumi_fortheorycomp->SetBinError(bin+1,relerr*xs);
      std::cout << "Bin " << bin << " " << xs << " +/- " << 100*relerr << " % (stat.+syst.)" << std::endl;
    }
    TCanvas *canv5 = new TCanvas();
    canv5->cd();
    histo_finalxsnolumi_fortheorycomp->Draw("e1");

    histo_finalxsnolumi_fortheorycomp->SaveAs(Form("plots/%s.root",histo_finalxsnolumi_fortheorycomp->GetName()));

  }

};

void post_process_all(bool skipsystematics = false, TString var=""){
  for (std::vector<TString>::const_iterator it = diffvariables_list.begin(); it!=diffvariables_list.end(); it++){
    if (var!="" && var!=*it) continue;
    post_process(it->Data(),"EBEB",skipsystematics);
    post_process(it->Data(),"EBEE",skipsystematics);
    post_process(it->Data(),"EEEE",skipsystematics);
    post_process(it->Data(),"inclusive",skipsystematics);
  }
}

void reweight_rhosigma(RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold){

  TH2F *hnum = new TH2F("hnum","hnum",30,0,30,20,0,10);
  TH2F *hden = new TH2F("hden","hden",30,0,30,20,0,10);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue("roorho")),fabs((*dset)->get(i)->getRealValue("roosigma")),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue("roorho")),fabs(dsetdestination->get(i)->getRealValue("roosigma")),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH2F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_rhosigmarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float rho = args.getRealValue("roorho");
    float sigma = args.getRealValue("roosigma");
    float neww = oldw*h->GetBinContent(h->FindBin(rho,sigma));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "RhoSigma2D rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  if (deleteold) delete old_dset;

};

void reweight_pt_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar){



  TH1F *hnum = new TH1F("hnum","hnum",n_ptbins_forreweighting,ptbins_forreweighting);
  TH1F *hden = new TH1F("hden","hden",n_ptbins_forreweighting,ptbins_forreweighting);
  hnum->Sumw2();
  hden->Sumw2();

  const char* ptname=Form("roopt%d",numvar);

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue(ptname)),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue(ptname)),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_ptrew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float pt = args.getRealValue(ptname);
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(pt)));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Pt 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};

void reweight_pt_2d(RooDataSet **dset, RooDataSet *dsetdestination){

  TH2F *hnum = new TH2F("hnum","hnum",n_ptbins_forreweighting,ptbins_forreweighting,n_ptbins_forreweighting,ptbins_forreweighting);
  TH2F *hden = new TH2F("hden","hden",n_ptbins_forreweighting,ptbins_forreweighting,n_ptbins_forreweighting,ptbins_forreweighting);
//  TH2F *hnum = new TH2F("hnum","hnum",30,0,300,30,0,300);
//  TH2F *hden = new TH2F("hden","hden",30,0,300,30,0,300);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue("roopt1")),fabs((*dset)->get(i)->getRealValue("roopt2")),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue("roopt1")),fabs(dsetdestination->get(i)->getRealValue("roopt2")),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH2F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_ptrew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float pt1 = args.getRealValue("roopt1");
    float pt2 = args.getRealValue("roopt2");
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(pt1),fabs(pt2)));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Pt 2d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};

void reweight_eta_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  TH1F *hnum = new TH1F("hnum","hnum",25,0,2.5);
  TH1F *hden = new TH1F("hden","hden",25,0,2.5);
  hnum->Sumw2();
  hden->Sumw2();

  const char* etaname=Form("rooeta%d",numvar);

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue(etaname)),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue(etaname)),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_etarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float eta = args.getRealValue(etaname);
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(eta)));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Eta 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};

void reweight_eta_2d(RooDataSet **dset, RooDataSet *dsetdestination){

  TH2F *hnum = new TH2F("hnum","hnum",25,0,2.5,25,0,2.5);
  TH2F *hden = new TH2F("hden","hden",25,0,2.5,25,0,2.5);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue("rooeta1")),fabs((*dset)->get(i)->getRealValue("rooeta2")),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue("rooeta1")),fabs(dsetdestination->get(i)->getRealValue("rooeta2")),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH2F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_etarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float eta1 = args.getRealValue("rooeta1");
    float eta2 = args.getRealValue("rooeta2");
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(eta1),fabs(eta2)));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Eta 2d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};


void reweight_sigma(RooDataSet **dset, RooDataSet *dsetdestination){

  TH1F *hnum = new TH1F("hnum","hnum",20,0,10);
  TH1F *hden = new TH1F("hden","hden",20,0,10);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill((*dset)->get(i)->getRealValue("roosigma"),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(dsetdestination->get(i)->getRealValue("roosigma"),dsetdestination->store()->weight(i));
  }

  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  //  h->SaveAs("plots/ratio.root");

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_sigmarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float sigma = args.getRealValue("roosigma");
    float neww = oldw*h->GetBinContent(h->FindBin(sigma));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Sigma reweighted: Norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;
  std::cout << "Sigma moving " << old_dset->mean(*roosigma) << " " << newdset->mean(*roosigma)  << std::endl;
  delete old_dset;

};

void reweight_rho(RooDataSet **dset, RooDataSet *dsetdestination, RooPlot *plot){

  TH1F *hnum = new TH1F("hnum","hnum",30,0,30);
  TH1F *hden = new TH1F("hden","hden",30,0,30);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill((*dset)->get(i)->getRealValue("roorho"),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(dsetdestination->get(i)->getRealValue("roorho"),dsetdestination->store()->weight(i));
  }

  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;



  RooPlot *p = plot;
  if (plot){
    (*dset)->plotOn(p);
    dsetdestination->plotOn(p,MarkerColor(kBlue));
  }

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_rhorew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float rho = args.getRealValue("roorho");
    float neww = oldw*h->GetBinContent(h->FindBin(rho));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }



  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Rho reweighted: Norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;
  std::cout << "Rho moving " << old_dset->mean(*roorho) << " " << newdset->mean(*roorho)  << std::endl;
  delete old_dset;

  if (plot)  (*dset)->plotOn(p,MarkerColor(kRed));

};

void plot_template_dependency_axis1(RooDataSet *dset, TString variable, float min, float max, int bins, bool dobinned){

  dset->Print();

  TH1F *histo[bins];
  TString title = Form("templ_dependency_%s_%s",dset->GetName(),variable.Data());
  if (!dobinned) title+=TString("_unbinned"); 
  for (int i=0; i<bins; i++) {
    TString loctitle = Form("%s_bin%d",title.Data(),i);
    if (!dobinned) histo[i] = new TH1F(loctitle.Data(),loctitle.Data(),n_histobins,leftrange,rightrange);
    else histo[i] = new TH1F(loctitle.Data(),loctitle.Data(),n_templatebins,0.5,n_templatebins+0.5);
    histo[i]->Sumw2();
  }

  for (int i=0; i<dset->numEntries(); i++){
    RooArgSet args = *(dset->get(i));
    float w = dset->store()->weight(i);
    float var = 0;
    if (variable=="pt") var = args.getRealValue("roopt1");
    if (variable=="sieie") var = args.getRealValue("roosieie1");
    if (var>=max) var=max-1e-5;
    if (var<min) var=min;
    int bin = (var-min)/(max-min)*bins;
    //    std::cout << var << " " << bin << " " << args.getRealValue("roovar1") << std::endl;
    if (!dobinned) histo[bin]->Fill(args.getRealValue("roovar1"),w);
    else histo[bin]->Fill(args.getRealValue("binning_roovar1"),w);
  }

  Int_t colors[7] = {kBlack,kRed,kBlue,kGreen,kMagenta,kCyan,kOrange};

  std::cout << "Bins:" << std::endl;
  for (int i=0; i<bins; i++) {
    std::cout << "bin" << i << ": " << min+(max-min)/bins*i << " - " << min+(max-min)/bins*(i+1) << std::endl;
    histo[i]->SetMarkerColor(bins<7 ? colors[i] : 40+i);
    histo[i]->SetLineColor(bins<7 ? colors[i] : 40+i);
    histo[i]->SetMarkerStyle(20);
    histo[i]->SetStats(0);
    if (histo[i]->Integral()>0) histo[i]->Scale(1.0/histo[i]->Integral());
  }

  TCanvas *canv_templ = new TCanvas(Form("canv_templ_%s",title.Data()),Form("canv_templ_%s",title.Data()));
  canv_templ->cd();
  histo[0]->Draw("C");
  for (int i=1; i<bins; i++) histo[i]->Draw("Csame");

};

void validate_reweighting(RooDataSet *dset, RooDataSet *dsetdestination, int numvar){

  std::cout << "Validating " << dset->GetName() << " " << dset->sumEntries() << " vs. " << dsetdestination->GetName() << " " << dsetdestination->sumEntries() << std::endl;

  TH1F *test[4];
  TH1F *target[4];

  test[0] = new TH1F("test_rho","test_rho",30,0,30);
  test[1] = new TH1F("test_sigma","test_sigma",20,0,10);
  test[2] = new TH1F("test_pt","test_pt",n_ptbins_forreweighting,ptbins_forreweighting);
  test[3] = new TH1F("test_eta","test_eta",25,0,2.5);
  target[0] = new TH1F("target_rho","target_rho",30,0,30);
  target[1] = new TH1F("target_sigma","target_sigma",20,0,10);
  target[2] = new TH1F("target_pt","target_pt",n_ptbins_forreweighting,ptbins_forreweighting);
  target[3] = new TH1F("target_eta","target_eta",25,0,2.5);

  const char* ptname=Form("roopt%d",numvar);
  const char* etaname=Form("rooeta%d",numvar);

  for (int i=0; i<dset->numEntries(); i++){
    RooArgSet args = *(dset->get(i));
    float oldw = dset->store()->weight(i);
    test[0]->Fill(args.getRealValue("roorho"),oldw);
    test[1]->Fill(args.getRealValue("roosigma"),oldw);
    test[2]->Fill(args.getRealValue(ptname),oldw);
    test[3]->Fill(args.getRealValue(etaname),oldw);
  }    

  for (int i=0; i<dsetdestination->numEntries(); i++){
    target[0]->Fill(dsetdestination->get(i)->getRealValue("roorho"),dsetdestination->store()->weight(i));
    target[1]->Fill(dsetdestination->get(i)->getRealValue("roosigma"),dsetdestination->store()->weight(i));
    target[2]->Fill(dsetdestination->get(i)->getRealValue(ptname),dsetdestination->store()->weight(i));
    target[3]->Fill(dsetdestination->get(i)->getRealValue(etaname),dsetdestination->store()->weight(i));

  }

  for (int i=0; i<4; i++){
    test[i]->Scale(1.0/test[i]->Integral());
    target[i]->Scale(1.0/target[i]->Integral());
    test[i]->SetLineColor(kRed);
    test[i]->SetMarkerColor(kRed);
  }

  TString name(dset->GetName());
  name.Append(Form("_roovar%d",numvar));
  TCanvas *c = new TCanvas(name.Data(),name.Data());
  c->Divide(2,2);

  for (int i=0;i<4; i++){
    c->cd(i+1);
    test[i]->Draw();
    target[i]->Draw("same");
    //    test[i]->SaveAs(Form("plots/test%d.root",i));
  }



  //  c->SaveAs(Form("%s_rew.png",(*dset)->GetName()));

//  for (int i=0;i<4; i++){
//    delete test[i]; delete target[i];
//  }

};

void plot_datasets_axis1(std::vector<plot_dataset_struct> dsets, TString outname, TString legtitle, bool legendup, bool dolin){

  const char* varname = "roovar1";
  const int ndsets = (int)(dsets.size());

  TH1F *h[100];

  for (int j=0; j<ndsets; j++){
    h[j] = new TH1F(Form("histo_%d_rv%d",j,1),Form("histo_%d_rv%d",j,1),n_histobins,leftrange,rightrange);
    for (int i=0; i<(dsets[j].dset)->numEntries(); i++) h[j]->Fill((dsets[j].dset)->get(i)->getRealValue(varname),(dsets[j].dset)->store()->weight(i));
    h[j]->Scale(1.0/h[j]->Integral());
    h[j]->SetLineWidth(2);
    h[j]->SetLineColor(dsets[j].color);
    h[j]->SetMarkerColor(dsets[j].color);
    h[j]->SetTitle("");
    //    if (legdata.leg[j].Index("data",4)!=kNPOS) h[j]->SetMarkerStyle(20);
    h[j]->SetMarkerStyle(kFullCircle);
    if (dsets[j].legend.Index("left",4)!=kNPOS) h[j]->SetMarkerStyle(kOpenTriangleUp);
    if (dsets[j].legend.Index("right",5)!=kNPOS) h[j]->SetMarkerStyle(kOpenSquare);
    h[j]->SetStats(0);
    h[j]->GetXaxis()->SetTitle("Photon PFIso (GeV)");
    //    h[j]->GetXaxis()->SetRangeUser(0,6);
  }

  TCanvas *comp = new TCanvas("shape_comparison");

  float max=0;
  for (int j=0; j<ndsets; j++){
    float thismax = h[j]->GetBinContent(h[j]->GetMaximumBin());
    max = (thismax>max) ? thismax : max;
  }
  h[0]->GetYaxis()->SetRangeUser(TMath::Max(h[0]->GetMinimum(),1e-4),max*1.05);

  if (!dolin) comp->SetLogy(1);

  h[0]->GetYaxis()->SetTitle("a.u.");
  h[0]->Draw();
  for (int j=1; j<ndsets; j++) h[j]->Draw("same");
  h[0]->Draw("same");

  TLegend *leg = (legendup) ? new TLegend(0.6,0.7,0.9,0.9,legtitle.Data()) : new TLegend(0.6,0.15,0.9,0.35,legtitle.Data());
  leg->SetFillColor(kWhite);

  for (int j=0; j<ndsets; j++) {
//    if (legdata.leg[j].Index("data",4)!=kNPOS) leg->AddEntry(h[j],legdata.leg[j].Data(),"p");
//    else leg->AddEntry(h[j],legdata.leg[j].Data(),"lp");
      leg->AddEntry(h[j],dsets[j].legend.Data(),"lp");
  }
  leg->Draw();

  TLatex a;
  a.SetNDC();
  a.SetTextSize(0.03);
  if (legendup)  a.DrawLatex(0.63,0.6,"#splitline{CMS Preliminary}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");
  else a.DrawLatex(0.63,0.85,"#splitline{CMS Preliminary}{#sqrt{s} = 7 TeV L = 5.0 fb^{-1}}");

  comp->SaveAs(Form("%s.%s",outname.Data(),"root"));
  comp->SaveAs(Form("%s.%s",outname.Data(),"pdf"));
  comp->SaveAs(Form("%s.%s",outname.Data(),"png"));

};

void produce_category_binning(RooDataSet **dset, bool deleteold){

  assert ((*dset)->numEntries()>0);
  RooArgSet newargs;
  {
    RooArgSet initialvars = *((*dset)->get(0));
    if (initialvars.find("roovar1")) {(*dset)->addColumn(*binning_roovar1_threshold); newargs.add(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*binning_roovar1));}
    if (initialvars.find("roovar2")) {(*dset)->addColumn(*binning_roovar2_threshold); newargs.add(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*binning_roovar2));}
    newargs.add(RooArgList(*roorho,*roosigma,*rooweight));
  }

  RooDataSet *old_dset = *dset;
  RooDataSet *newdset = new RooDataSet(Form("%s_binned",(*dset)->GetName()),Form("%s_binned",(*dset)->GetName()),newargs,WeightVar(*rooweight));

  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float w = (*dset)->store()->weight(i);
    
    if (args.find("roovar1")){
      roovar1->setVal(args.getRealValue("roovar1"));
      roopt1->setVal(args.getRealValue("roopt1"));
      roosieie1->setVal(args.getRealValue("roosieie1"));
      rooeta1->setVal(args.getRealValue("rooeta1"));
      //      binning_roovar1->setIndex(args.getCatIndex("binning_roovar1_threshold"));
      binning_roovar1->setVal(args.getCatIndex("binning_roovar1_threshold"));
    }
    if (args.find("roovar2")){
      roovar2->setVal(args.getRealValue("roovar2"));
      roopt2->setVal(args.getRealValue("roopt2"));
      roosieie2->setVal(args.getRealValue("roosieie2"));
      rooeta2->setVal(args.getRealValue("rooeta2"));
      //      binning_roovar2->setIndex(args.getCatIndex("binning_roovar2_threshold"));
      binning_roovar2->setVal(args.getCatIndex("binning_roovar2_threshold"));
    }
    roorho->setVal(args.getRealValue("roorho"));
    roosigma->setVal(args.getRealValue("roosigma"));

    newdset->add(newargs,w);
  }


    *dset=newdset;
    TString nametitle = old_dset->GetName();
    old_dset->SetName(Form("%s_OLD",nametitle.Data()));
    old_dset->SetTitle(Form("%s_OLD",nametitle.Data()));
    newdset->SetName(nametitle.Data());
    newdset->SetTitle(nametitle.Data());

    std::cout << "Dataset rebinned from "; old_dset->Print();  std::cout << " to "; newdset->Print();

    if (deleteold) delete old_dset;

};

void randomize_dataset_statistically_binned(RooDataSet **dset){

  bool plot = false; 



  assert ((*dset)->numEntries()>0);  
  RooArgSet initialvars = *((*dset)->get(0));
  assert (initialvars.find("binning_roovar1") || initialvars.find("binning_roovar2"));

  int code = 0;
  if (initialvars.find("binning_roovar1") && initialvars.find("binning_roovar2")) code=3;
  else if (initialvars.find("binning_roovar1")) code=1;
  else code=2;
  assert (code>0);


  TH1F *hnum1d = new TH1F("hnum1d","hnum1d",n_templatebins,0.5,0.5+n_templatebins);
  TH1F *hden1d = NULL;
  hnum1d->Sumw2();
  TH2F *hnum2d = new TH2F("hnum2d","hnum2d",n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
  TH2F *hden2d = NULL;
  hnum2d->Sumw2();

  if (code==3) create_histo_from_dataset_binned(*dset,NULL,&hden2d); else create_histo_from_dataset_binned(*dset,&hden1d,NULL);


  if (code==3){  
    for (int i=0; i<hden2d->GetNbinsX()+1; i++)
      for (int j=0; j<hden2d->GetNbinsY()+1; j++){
	hnum2d->SetBinContent(i,j,hden2d->GetBinContent(i,j)+hden2d->GetBinError(i,j)*_random_generator->Gaus());
	if (hnum2d->GetBinContent(i,j)<0) hnum2d->SetBinContent(i,j,0);
	hnum2d->SetBinError(i,j,0);
      }
    if (plot) {
      TH2F *h2old = (TH2F*)(hden2d->Clone("h2old")); h2old->SetLineColor(kRed); h2old->SetMarkerColor(kRed);
      TH2F *h2new = (TH2F*)(hnum2d->Clone("h2new")); h2new->SetMarkerColor(kBlack); h2new->SetMarkerStyle(20);
      h2old->ProjectionX()->Draw("E1L");
      h2new->ProjectionX()->Draw("PSAME");
    }    
    hnum2d->Divide(hden2d);
  }
  else {
    for (int i=0; i<hden1d->GetNbinsX()+1; i++){
      hnum1d->SetBinContent(i,hden1d->GetBinContent(i)+hden1d->GetBinError(i)*_random_generator->Gaus());
      if (hnum1d->GetBinContent(i)<0) hnum1d->SetBinContent(i,0);
      hnum1d->SetBinError(i,0);
    }
    if (plot) {
      TH1F *h1old = (TH1F*)(hden1d->Clone("h1old")); h1old->SetLineColor(kRed); h1old->SetMarkerColor(kRed);
      TH1F *h1new = (TH1F*)(hnum1d->Clone("h1new")); h1new->SetMarkerColor(kBlack); h1new->SetMarkerStyle(20);
      h1old->Draw("E1L");
      h1new->Draw("PSAME");
    }
    hnum1d->Divide(hden1d);
  }



  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_statfluct",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float neww = 0;
    if (code==3) neww = oldw*hnum2d->GetBinContent(hnum2d->FindBin(args.getRealValue("binning_roovar1"),args.getRealValue("binning_roovar2")));
    else if (code==1) neww = oldw*hnum1d->GetBinContent(hnum1d->FindBin(args.getRealValue("binning_roovar1")));
    else if (code==2) neww = oldw*hnum1d->GetBinContent(hnum1d->FindBin(args.getRealValue("binning_roovar2")));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }

  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum1d; if (hden1d) delete hden1d;
  delete hnum2d; if (hden2d) delete hden2d;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  TString nametitle = old_dset->GetName();
  old_dset->SetName(Form("%s_OLD",nametitle.Data()));
  old_dset->SetTitle(Form("%s_OLD",nametitle.Data()));
  newdset->SetName(nametitle.Data());
  newdset->SetTitle(nametitle.Data());
  
  std::cout << "Dataset randomized from "; old_dset->Print();  std::cout << " to "; newdset->Print();

  delete old_dset;

};


void create_histo_from_dataset_binned(RooDataSet *dset, TH1F** h1out, TH2F** h2out){

  assert (dset->numEntries()>0);  
  RooArgSet initialvars = *(dset->get(0));
  assert (initialvars.find("binning_roovar1") || initialvars.find("binning_roovar2"));

  int code = 0;
  if (initialvars.find("binning_roovar1") && initialvars.find("binning_roovar2")) code=3;
  else if (initialvars.find("binning_roovar1")) code=1;
  else code=2;
  assert (code>0);

  TString nametitle = dset->GetName();
  
  TH1F *h1d = new TH1F(nametitle+TString("_histo1d"),nametitle+TString("_histo1d"),n_templatebins,0.5,0.5+n_templatebins);
  h1d->Sumw2();
  TH2F *h2d = new TH2F(nametitle+TString("_histo2d"),nametitle+TString("_histo2d"),n_templatebins,0.5,0.5+n_templatebins,n_templatebins,0.5,0.5+n_templatebins);
  h2d->Sumw2();

  for (int i=0; i<dset->numEntries(); i++){
    if (code==3) h2d->Fill(dset->get(i)->getRealValue("binning_roovar1"),dset->get(i)->getRealValue("binning_roovar2"),dset->store()->weight(i));
    else if (code==1) h1d->Fill(dset->get(i)->getRealValue("binning_roovar1"),dset->store()->weight(i));
    else if (code==2) h1d->Fill(dset->get(i)->getRealValue("binning_roovar2"),dset->store()->weight(i));
  }

  std::cout << "Produced histo from dset " << nametitle.Data() << std::endl;
  
  //  if (code==3) h2d->Print("v"); else h1d->Print("v");

  if (code==3) {delete h1d; *h2out=h2d;}
  else {delete h2d; *h1out=h1d;}


};

void create_histo_from_dataset_variablebins(RooDataSet *dset, TH1F** h1out, TH2F** h2out){

  assert (dset->numEntries()>0);  
  RooArgSet initialvars = *(dset->get(0));
  assert (initialvars.find("binning_roovar1") || initialvars.find("binning_roovar2"));

  int code = 0;
  if (initialvars.find("binning_roovar1") && initialvars.find("binning_roovar2")) code=3;
  else if (initialvars.find("binning_roovar1")) code=1;
  else code=2;
  assert (code>0);

  TString nametitle = dset->GetName();
  TH1F *h1d = NULL;
  TH2F *h2d = NULL;

  TH1F *hnew1d = (code<3) ? new TH1F(nametitle+TString("_histo1dvb"),nametitle+TString("_histo1dvb"),n_templatebins,templatebinsboundaries) : NULL;
  TH2F *hnew2d = (code==3) ? new TH2F(nametitle+TString("_histo2dvb"),nametitle+TString("_histo2dvb"),n_templatebins,templatebinsboundaries,n_templatebins,templatebinsboundaries) : NULL;


  if (code==3) create_histo_from_dataset_binned(dset,NULL,&h2d); else create_histo_from_dataset_binned(dset,&h1d,NULL);

  if (code==3){  
    for (int i=0; i<h2d->GetNbinsX()+1; i++)
      for (int j=0; j<h2d->GetNbinsY()+1; j++){
	hnew2d->SetBinContent(i,j,h2d->GetBinContent(i,j)/hnew2d->GetXaxis()->GetBinWidth(i)/hnew2d->GetYaxis()->GetBinWidth(j));
	hnew2d->SetBinError(i,j,h2d->GetBinError(i,j)/hnew2d->GetXaxis()->GetBinWidth(i)/hnew2d->GetYaxis()->GetBinWidth(j));
      }
    hnew2d->GetZaxis()->SetTitle("ev. / GeV^{2}");
  }
  else {
    for (int i=0; i<h1d->GetNbinsX()+1; i++){
      hnew1d->SetBinContent(i,h1d->GetBinContent(i)/hnew1d->GetBinWidth(i));
      hnew1d->SetBinError(i,h1d->GetBinError(i)/hnew1d->GetBinWidth(i));
    }
    hnew1d->GetYaxis()->SetTitle("ev. / GeV");
  }

  if (code==3) {delete h2d; *h2out=hnew2d;}
  else {delete h1d; *h1out=hnew1d;}


};


//void generate_toy_dataset_1d(RooDataSet **target, RooAbsPdf *sigpdf, RooAbsPdf *bkgpdf, float fsig1toy){
//  //  std::cout << "TOY GENERATION DEBUG:" << std::endl;
//
//  assert ((*target)->numEntries()>0);  
//  RooArgSet initialvars = *((*target)->get(0));
//  assert (initialvars.find("binning_roovar1") || initialvars.find("binning_roovar2"));
//  int code = 0;
//  if (initialvars.find("binning_roovar1")) code=1;
//  else if (initialvars.find("binning_roovar2")) code=2;
//  assert (code>0);
//
//  RooAddPdf *addpdf = new RooAddPdf("addpdf","addpdf",*sigpdf,*bkgpdf,RooRealConstant::value(fsig1toy));
//
//  RooDataSet *generated = addpdf->generate((code==1) ? RooArgSet(*binning_roovar1) : RooArgSet(*binning_roovar2),_random_generator->Poisson((*target)->sumEntries()),kFALSE,kFALSE,"",kFALSE,kTRUE);
//
////  (*target)->Print();
////  (*target)->Print("v");
//  std::cout << "Generated dataset:" << std::endl;
//  generated->Print();
//  //  generated->Print("v");
//
//  delete *target;
//  *target = generated;
//
//};

void generate_toy_dataset_2d(RooDataSet **target, RooAbsPdf *sigsigpdf, RooAbsPdf *sigbkgpdf, RooAbsPdf *bkgsigpdf, RooAbsPdf *bkgbkgpdf, float pptoy, float pftoy, float fptoy){
  //  std::cout << "TOY GENERATION DEBUG:" << std::endl;

  assert (pptoy+pftoy+fptoy<=1);

  assert ((*target)->numEntries()>0);  
  RooArgSet initialvars = *((*target)->get(0));
  assert (initialvars.find("binning_roovar1") && initialvars.find("binning_roovar2"));

  RooAddPdf *addpdf = new RooAddPdf("addpdf","addpdf",RooArgList(*sigsigpdf,*sigbkgpdf,*bkgsigpdf,*bkgbkgpdf),RooArgList(RooRealConstant::value(pptoy),RooRealConstant::value(pftoy),RooRealConstant::value(fptoy)),kFALSE);

  addpdf->Print();
  (*target)->Print();

  RooDataSet *generated = addpdf->generate(RooArgSet(*binning_roovar1,*binning_roovar2),_random_generator->Poisson((*target)->sumEntries()),kFALSE,kFALSE,"",kFALSE,kTRUE);

  //  (*target)->Print();
  //  (*target)->Print("v");
  std::cout << "Generated dataset:" << std::endl;
  generated->Print();
  //  generated->Print("v");

  delete *target;
  *target = generated;

};

void print_mem(){
  gSystem->GetProcInfo(&procinfo); 
  std::cout << "Resident mem (kB): " << procinfo.fMemResident << std::endl; 
  std::cout << "Virtual mem (kB):  " << procinfo.fMemVirtual << std::endl; 
  gSystem->Sleep(1e3);
};

bool myfunc_sortonfirst(pair<float,float> a, pair<float,float> b) { return (a.first<b.first); };

void find_adaptive_binning(RooDataSet *dset, int *n_found_bins, Double_t *array_bounds, int axis, float threshold){

//  // DEBUG
//  *n_found_bins=3;
//  Double_t templatebinsboundaries_reduced[4]={-3,-1,3,9};
//  for (int i=0; i<4; i++) array_bounds[i] = templatebinsboundaries_reduced[i];
//  return;
//    *n_found_bins=5;
//    Double_t templatebinsboundaries_reduced[6] = {-3,-0.5,0,1,4,9};
//    for (int i=0; i<6; i++) array_bounds[i] = templatebinsboundaries_reduced[i];
//    return;

//  *n_found_bins=48;
//  cout << "DEBUG: FINE BINNING" << endl;
//  for (int i=0; i<49; i++) array_bounds[i]=leftrange+(rightrange-leftrange)/48.*i;
//  return;

  if (threshold<0 && threshold>=-100){
    std::cout << "APPLYING FIXED BIN NUMBER BOUND (10)" << std::endl;
    *n_found_bins=10;
    Double_t templatebinsboundaries_reduced[11] = {-3,-0.75,-0.5,-0.25,0,0.5,1,2,4,6,9};
    for (int i=0; i<11; i++) array_bounds[i] = templatebinsboundaries_reduced[i];
    return;
  }
  if (threshold<-100){
    std::cout << "APPLYING FIXED BIN NUMBER BOUND (5)" << std::endl;
    *n_found_bins=5;
    Double_t templatebinsboundaries_reduced[6] = {-3,-0.5,0,1,4,9};
    for (int i=0; i<6; i++) array_bounds[i] = templatebinsboundaries_reduced[i];
    return;
  }

  assert ((axis==1) || (axis==2));
  TString string = (axis==1) ? TString("roovar1") : TString("roovar2");

  vector<pair<float,float> > values;
  for (int i=0; i<dset->numEntries(); i++){
    values.push_back(pair<float,float>(dset->get(i)->getRealValue(string.Data()),dset->store()->weight(i)));
  }
  sort(values.begin(),values.end(),myfunc_sortonfirst);

  vector<float> boundaries;
  boundaries.push_back(leftrange);
  float totw = 0;
  float totw2 = 0;
  for (int i=0; i<int(values.size()); i++){
    totw+=values.at(i).second;
    totw2+=pow(values.at(i).second,2);
    if (sqrt(totw2)/totw<threshold && totw>0) {
      boundaries.push_back(values.at(i).first);
      totw=0;
      totw2=0;
    }
  }
  if (boundaries.size()>1) boundaries.at(boundaries.size()-1) = rightrange; else boundaries.push_back(rightrange);

  std::cout << "debug adaptive binning:" << std::endl;
  std::cout << boundaries.size()-1 << " bins" << std::endl;
  for (int i=0; i<int(boundaries.size()); i++) std::cout << boundaries.at(i) << " ";
  std::cout << std::endl;

  *n_found_bins = boundaries.size()-1;
  assert (*n_found_bins<=n_templatebins_max);
  for (int i=0; i<int(boundaries.size()); i++) array_bounds[i] = boundaries.at(i);

  if (*n_found_bins<5) {
    std::cout << "APPLYING LOWER BIN NUMBER BOUND AT 5" << std::endl;
    *n_found_bins=5;
    Double_t templatebinsboundaries_reduced[6] = {-3,-0.5,0,1,4,9};
    for (int i=0; i<6; i++) array_bounds[i] = templatebinsboundaries_reduced[i];
  }

};

//RooDataSet** split_in_eta_cats(RooDataSet *dset, int numvar);
//RooDataSet** split_in_eta1eta2_cats(RooDataSet *dset);
//void reweight_pteta(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);

/*
RooDataSet** split_in_eta_cats(RooDataSet *dset, int numvar){

  RooDataSet **output = new RooDataSet*[n_eta_cats];

  std::cout << "Splitting " << dset->GetName() << std::endl;

  for (int k=0; k<n_eta_cats; k++){
    output[k]=(RooDataSet*) (dset->reduce(Cut(Form("TMath::Abs(rooeta%d)>%f && TMath::Abs(rooeta%d)<%f",numvar,etabins[k],numvar,etabins[k+1])),Name(Form("%s_eta%d",dset->GetName(),k)),Title(Form("%s_eta%d",dset->GetName(),k))));
    std::cout << "eta bin " << k << " nentries " << output[k]->sumEntries() << std::endl;
  }

  return output;

};

RooDataSet** split_in_eta1eta2_cats(RooDataSet *dset){

  RooDataSet **output = new RooDataSet*[n_eta1eta2_cats];

  std::cout << "Splitting " << dset->GetName() << std::endl;

  for (int k=0; k<n_eta1eta2_cats; k++){
    int rbin = ((int)k)/((int)n_eta_cats);
    int sbin = ((int)k)%((int)n_eta_cats);
    output[k]=(RooDataSet*) (dset->reduce(Cut(Form("TMath::Abs(rooeta1)>%f && TMath::Abs(rooeta1)<%f && TMath::Abs(rooeta2)>%f && TMath::Abs(rooeta2)<%f",etabins[rbin],etabins[rbin+1],etabins[sbin],etabins[sbin+1])),Name(Form("%s_eta1eta2%d",dset->GetName(),k)),Title(Form("%s_eta1eta2%d",dset->GetName(),k))));
    std::cout << "eta1 bin " << rbin << " eta2 bin " << sbin << " nentries " << output[k]->sumEntries() << std::endl;
  }

  return output;

};
*/

/*
void reweight_pteta(RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  TH2F *hnum = new TH2F("hnum","hnum",55,25,300,25,0,2.5);
  TH2F *hden = new TH2F("hden","hden",55,25,300,25,0,2.5);
  hnum->Sumw2();
  hden->Sumw2();

  const char* ptname=Form("roopt%d",numvar);
  const char* etaname=Form("rooeta%d",numvar);

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill((*dset)->get(i)->getRealValue(ptname),fabs((*dset)->get(i)->getRealValue(etaname)),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(dsetdestination->get(i)->getRealValue(ptname),fabs(dsetdestination->get(i)->getRealValue(etaname)),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH2F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_ptetarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float pt = args.getRealValue(ptname);
    float eta = args.getRealValue(etaname);
    float neww = oldw*h->GetBinContent(h->FindBin(pt,fabs(eta)));
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Pt+Eta rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  //  delete old_dset;

};
*/

fit_output* template_studies_2d_variablebinning(TString diffvariable, TString splitting, int bin, const TString do_syst_string=TString("")){
  setTDRStyle();
  return fit_dataset(diffvariable,splitting,bin,do_syst_string);
};

TH1F* AddTHInQuadrature(std::vector<TH1F*> vector, TString name){

  if (vector.size()==0) return NULL;

  TH1F *tot = (TH1F*)(vector[0]->Clone(name.Data()));
  tot->Reset();
  tot->Sumw2();

  for (unsigned int i=0; i<vector.size(); i++){
    TH1F *h = (TH1F*)(vector[i]->Clone());
    h->Multiply(h);
    tot->Add(h);
    delete h;
  }

  for (int i=0; i<tot->GetNbinsX(); i++) {
    tot->SetBinContent(i+1,sqrt(tot->GetBinContent(i+1)));
    tot->SetBinError(i+1,0);
  }

  return tot;

};
