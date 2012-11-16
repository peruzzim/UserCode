#include <assert.h>

#include "binsdef.h"
#include "RooFitResult.h"
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
//#include "TThread.h"
#include "firstbinpdf.cxx"

using namespace std;
using namespace RooFit;

typedef struct {
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
} fit_output; 

bool study_templates=0;
bool study_templates_plotting=0;

const int numcpu=4;

RooDataSet** split_in_eta_cats(RooDataSet *dset, int numvar);
RooDataSet** split_in_eta1eta2_cats(RooDataSet *dset);
void reweight_pteta(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void reweight_eta(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void reweight_rho(RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_sigma(RooDataSet **dset, RooDataSet *dsetdestination);
void validate_reweighting(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);

RooRealVar *roovar1=NULL;
RooRealVar *roovar2=NULL;
RooRealVar *roopt1=NULL;
RooRealVar *roopt2=NULL;
RooRealVar *rooeta1=NULL;
RooRealVar *rooeta2=NULL;
RooRealVar *roorho=NULL;
RooRealVar *roosigma=NULL;
RooRealVar *rooweight=NULL;

bool ok_for_recycle=0;
float lm_sig1_recycle=0;
float lm_sig2_recycle=0;
float lm_bkg1_recycle=0;
float lm_bkg2_recycle=0;
RooHistPdf *sig1pdf_recycle=NULL;
RooHistPdf *sig2pdf_recycle=NULL;
RooHistPdf *bkg1pdf_recycle=NULL;
RooHistPdf *bkg2pdf_recycle=NULL;


fit_output* fit_dataset(const char* inputfilename_ts, const char* inputfilename_tb, const char* inputfilename_d, TString diffvariable, TString splitting, int bin, int template_bin=-1, int ismc=0){

  bool iststudy=(template_bin>-1);
  if (iststudy) n_rhosigma_cats=0;

  for (int k=0; k<10; k++) std::cout << std::endl;
  std::cout << "Process " << diffvariable.Data() << " bin " << bin << std::endl;
  for (int k=0; k<2; k++) std::cout << std::endl;

//  freopen("/dev/null","w",stdout);
//  freopen("/dev/null","w",stderr);

  fit_output *out = new fit_output();
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

  TFile *inputfile_ts = TFile::Open(inputfilename_ts);
  TFile *inputfile_tb = TFile::Open(inputfilename_tb);
  TFile *inputfile_d = TFile::Open(inputfilename_d);


  //TString data_dir("mc_Tree_standard_sel/");
  //TString data_dir("data_Tree_standard_sel/"); 
  TString data_dir("data_Tree_doublerandomcone_sel/"); 
  
  TString sig_dir("data_Tree_randomcone_signal_template/");
  //  TString sig_dir("mc_Tree_signal_template/");
  //  TString sig_dir("mc_Tree_randomcone_signal_template/");
  if (ismc==1) sig_dir="mc_Tree_signal_template/";
  if (ismc==2) sig_dir="mc_Tree_randomcone_signal_template/";

  TString bkg_dir("data_Tree_sieiesideband_sel/");
  //TString bkg_dir("mc_Tree_sieiesideband_sel/");
  if (ismc==1) bkg_dir="mc_Tree_background_template/"; 
  if (ismc==2) bkg_dir="mc_Tree_sieiesideband_sel/"; 

  if (splitting=="EEEB") splitting="EBEE";

  TH1::SetDefaultSumw2(kTRUE);
  
    TString s1; TString s2;
    if (splitting=="EBEB") {s1="EB"; s2="EB";}
    else if (splitting=="EEEE") {s1="EE"; s2="EE";}
    else if (splitting=="EBEE") {s1="EB"; s2="EE";}
    bool sym  = (s1==s2);
     
    RooWorkspace *wspace_ts=NULL;
    RooWorkspace *wspace_tb=NULL;
    RooWorkspace *wspace_d=NULL;
    inputfile_ts->GetObject(TString(sig_dir).Append("rooworkspace"),wspace_ts);
    inputfile_tb->GetObject(TString(bkg_dir).Append("rooworkspace"),wspace_tb);
    inputfile_d->GetObject(TString(data_dir).Append("rooworkspace"),wspace_d);
    assert(wspace_ts); //wspace_ts->Print();
    assert(wspace_tb); //wspace_tb->Print();
    assert(wspace_d); //wspace_d->Print();  

    if (roovar1==NULL) { roovar1 = wspace_d->var("roovar1"); roovar1->setBins(n_histobins);}
    if (roovar2==NULL) { roovar2 = wspace_d->var("roovar2"); roovar2->setBins(n_histobins);}
    if (roopt1==NULL) { roopt1 = wspace_d->var("roopt1"); }
    if (rooeta1==NULL) { rooeta1 = wspace_d->var("rooeta1"); }
    if (roopt2==NULL) { roopt2 = wspace_d->var("roopt2"); }
    if (rooeta2==NULL) { rooeta2 = wspace_d->var("rooeta2"); }
    if (roorho==NULL) { roorho = wspace_d->var("roorho"); }
    if (roosigma==NULL) { roosigma = wspace_d->var("roosigma"); }
    //if (rooweight==NULL) { rooweight = wspace_d->var("rooweight"); }
    rooweight = new RooRealVar("rooweight","rooweight",0,5);
    assert (roovar1);
    assert (roovar2);
    assert (roopt1);
    assert (rooeta1);
    assert (roopt2);
    assert (rooeta2);
    assert (roorho);
    assert (roosigma);
    assert (rooweight);

    const float pp_init = 0.4;
    const float pf_init = 0.2;
    const float fp_init = 0.1;
    
      // inversione:
      //  pp = pp;
      //  pf = j1 - pp;
      //  fp = j2 - pp;
      //  ff = (1-j1-j2)+pp
      
    RooRealVar *pp = new RooRealVar("pp","pp",pp_init,0,1);
    RooRealVar *j1 = new RooRealVar("j1","j1",pp_init+pf_init,0,1);
    RooRealVar *j2 = new RooRealVar("j2","j2",pp_init+fp_init,0,1);
    

    RooDataSet *sig1dset = (RooDataSet*)(wspace_ts->data(Form("roodset_signal_%s_b%d_rv%d",s1.Data(),(template_bin>-1) ? template_bin : n_bins,1)));
    RooDataSet *bkg1dset = (RooDataSet*)(wspace_tb->data(Form("roodset_background_%s_b%d_rv%d",s1.Data(),(template_bin>-1) ? template_bin : n_bins,1)));
    RooDataSet *sig2dset = (RooDataSet*)(wspace_ts->data(Form("roodset_signal_%s_b%d_rv%d",s2.Data(),(template_bin>-1) ? template_bin : n_bins,2)));
    RooDataSet *bkg2dset = (RooDataSet*)(wspace_tb->data(Form("roodset_background_%s_b%d_rv%d",s2.Data(),(template_bin>-1) ? template_bin : n_bins,2)));
    RooDataSet *dataset = (RooDataSet*)(wspace_d->data(Form("obs_roodset_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin)));    

    
    {
      float etacut = 0.4;
      sig1dset = (RooDataSet*)(sig1dset->reduce( Cut(Form("TMath::Abs(rooeta1)<%f",etacut)) ));
      bkg1dset = (RooDataSet*)(bkg1dset->reduce( Cut(Form("TMath::Abs(rooeta1)<%f",etacut)) ));
      sig2dset = (RooDataSet*)(sig2dset->reduce( Cut(Form("TMath::Abs(rooeta2)<%f",etacut)) ));
      bkg2dset = (RooDataSet*)(bkg2dset->reduce( Cut(Form("TMath::Abs(rooeta2)<%f",etacut)) ));
      dataset =  (RooDataSet*)(dataset->reduce( Cut(Form("TMath::Abs(rooeta1)<%f && TMath::Abs(rooeta2)<%f",etacut,etacut)) ));
    }
//    {
//      sig1dset = (RooDataSet*)(sig1dset->reduce( Cut("TMath::Abs(roosigma)>-1.0 && TMath::Abs(roosigma)<7.0") ));
//      bkg1dset = (RooDataSet*)(bkg1dset->reduce( Cut("TMath::Abs(roosigma)>-1.0 && TMath::Abs(roosigma)<7.0") ));
//      sig2dset = (RooDataSet*)(sig2dset->reduce( Cut("TMath::Abs(roosigma)>-1.0 && TMath::Abs(roosigma)<7.0") ));
//      bkg2dset = (RooDataSet*)(bkg2dset->reduce( Cut("TMath::Abs(roosigma)>-1.0 && TMath::Abs(roosigma)<7.0") ));
//      dataset =  (RooDataSet*)(dataset->reduce( Cut( "TMath::Abs(roosigma)>-1.0 && TMath::Abs(roosigma)<7.0") ));
//    }
    

    RooDataSet *dataset_axis1 = (RooDataSet*)(dataset->reduce(SelectVars(RooArgList(*roovar1,*roopt1,*rooeta1,*roorho,*roosigma))));
    RooDataSet *dataset_axis2 = (RooDataSet*)(dataset->reduce(SelectVars(RooArgList(*roovar2,*roopt2,*rooeta2,*roorho,*roosigma))));



    RooDataSet **sig1dset_eta1binned = split_in_eta_cats(sig1dset,1);
    RooDataSet **bkg1dset_eta1binned = split_in_eta_cats(bkg1dset,1);
    RooDataSet **sig2dset_eta2binned = split_in_eta_cats(sig2dset,2);
    RooDataSet **bkg2dset_eta2binned = split_in_eta_cats(bkg2dset,2);
    RooDataSet **dataset_eta1eta2binned = split_in_eta1eta2_cats(dataset);
    RooDataSet **dataset_eta1binned_axis1 = split_in_eta_cats(dataset_axis1,1);
    RooDataSet **dataset_eta2binned_axis2 = split_in_eta_cats(dataset_axis2,2);


    RooCategory *eta1eta2cat = new RooCategory("eta1eta2cat","eta1eta2cat");
    for (int k=0; k<n_eta1eta2_cats; k++) eta1eta2cat->defineType(Form("cat%d",k));


    map<string,RooDataSet*> mapping_eta1eta2cat;
    for (int k=0; k<n_eta1eta2_cats; k++) mapping_eta1eta2cat.insert(pair<string,RooDataSet*>(Form("cat%d",k),dataset_eta1eta2binned[k]));
    RooDataSet *dataset_combinedeta1eta2binned;
    if (!iststudy) dataset_combinedeta1eta2binned = new RooDataSet("dataset_combinedeta1eta2binned","dataset_combinedeta1eta2binned",RooArgSet(*roovar1,*roovar2),Index(*eta1eta2cat),Import(mapping_eta1eta2cat));


    std::cout << "Before rew" << std::endl;    
    sig1dset->Print();
    bkg1dset->Print();
    sig2dset->Print();
    bkg2dset->Print();
    dataset->Print(); 
    dataset_axis1->Print();
    dataset_axis2->Print();
    for (int k=0; k<n_eta_cats; k++) {
      sig1dset_eta1binned[k]->Print();
      bkg1dset_eta1binned[k]->Print();
      sig2dset_eta2binned[k]->Print();
      bkg2dset_eta2binned[k]->Print();
      std::cout << sig1dset_eta1binned[k]->sumEntries() << std::endl;
      std::cout << bkg1dset_eta1binned[k]->sumEntries() << std::endl;
      std::cout << sig2dset_eta2binned[k]->sumEntries() << std::endl;
      std::cout << bkg2dset_eta2binned[k]->sumEntries() << std::endl;
    }
    for (int k=0; k<n_eta1eta2_cats; k++) {
      dataset_eta1eta2binned[k]->Print();
      std::cout << dataset_eta1eta2binned[k]->sumEntries()  << std::endl;
    }


    
    { // eta reweighting
      reweight_eta(&sig1dset,dataset_axis1,1);
      reweight_eta(&bkg1dset,dataset_axis1,1);
      reweight_eta(&sig2dset,dataset_axis2,2);
      reweight_eta(&bkg2dset,dataset_axis2,2);
      for (int k=0; k<n_eta_cats; k++) {
	reweight_eta(&(sig1dset_eta1binned[k]),dataset_eta1binned_axis1[k],1);
	reweight_eta(&(bkg1dset_eta1binned[k]),dataset_eta1binned_axis1[k],1);
	reweight_eta(&(sig2dset_eta2binned[k]),dataset_eta2binned_axis2[k],2);
	reweight_eta(&(bkg2dset_eta2binned[k]),dataset_eta2binned_axis2[k],2);
      }
    }
//    { // rho reweighting
//      reweight_rho(&sig1dset,dataset_axis1);
//      reweight_rho(&bkg1dset,dataset_axis1);
//      reweight_rho(&sig2dset,dataset_axis2);
//      reweight_rho(&bkg2dset,dataset_axis2);
//      for (int k=0; k<n_eta_cats; k++) {
//	reweight_rho(&(sig1dset_eta1binned[k]),dataset_eta1binned_axis1[k]);
//	reweight_rho(&(bkg1dset_eta1binned[k]),dataset_eta1binned_axis1[k]);
//	reweight_rho(&(sig2dset_eta2binned[k]),dataset_eta2binned_axis2[k]);
//	reweight_rho(&(bkg2dset_eta2binned[k]),dataset_eta2binned_axis2[k]);
//      }
//    }
    { // sigma reweighting
      reweight_sigma(&sig1dset,dataset_axis1);
      reweight_sigma(&bkg1dset,dataset_axis1);
      reweight_sigma(&sig2dset,dataset_axis2);
      reweight_sigma(&bkg2dset,dataset_axis2);
      for (int k=0; k<n_eta_cats; k++) {
	reweight_sigma(&(sig1dset_eta1binned[k]),dataset_eta1binned_axis1[k]);
	reweight_sigma(&(bkg1dset_eta1binned[k]),dataset_eta1binned_axis1[k]);
	reweight_sigma(&(sig2dset_eta2binned[k]),dataset_eta2binned_axis2[k]);
	reweight_sigma(&(bkg2dset_eta2binned[k]),dataset_eta2binned_axis2[k]);
      }
    }



    std::cout << "After rew" << std::endl;    
    sig1dset->Print();
    bkg1dset->Print();
    sig2dset->Print();
    bkg2dset->Print();
    dataset->Print(); 
    dataset_axis1->Print();
    dataset_axis2->Print();
    for (int k=0; k<n_eta_cats; k++) {
      sig1dset_eta1binned[k]->Print();
      bkg1dset_eta1binned[k]->Print();
      sig2dset_eta2binned[k]->Print();
      bkg2dset_eta2binned[k]->Print();
      std::cout << sig1dset_eta1binned[k]->sumEntries() << std::endl;
      std::cout << bkg1dset_eta1binned[k]->sumEntries() << std::endl;
      std::cout << sig2dset_eta2binned[k]->sumEntries() << std::endl;
      std::cout << bkg2dset_eta2binned[k]->sumEntries() << std::endl;
    }
    for (int k=0; k<n_eta1eta2_cats; k++) {
      dataset_eta1eta2binned[k]->Print();
      std::cout << dataset_eta1eta2binned[k]->sumEntries()  << std::endl;
    }



    { // validate reweighting
      validate_reweighting(&sig1dset,dataset_axis1,1);
      validate_reweighting(&bkg1dset,dataset_axis1,1);
      validate_reweighting(&sig2dset,dataset_axis2,2);
      validate_reweighting(&bkg2dset,dataset_axis2,2);
      for (int k=0; k<n_eta_cats; k++) {
	validate_reweighting(&(sig1dset_eta1binned[k]),dataset_eta1binned_axis1[k],1);
	validate_reweighting(&(bkg1dset_eta1binned[k]),dataset_eta1binned_axis1[k],1);
	validate_reweighting(&(sig2dset_eta2binned[k]),dataset_eta2binned_axis2[k],2);
	validate_reweighting(&(bkg2dset_eta2binned[k]),dataset_eta2binned_axis2[k],2);
      }
    }



    if (iststudy){
      TH1* sig1dsethist = sig1dset->createHistogram("sig1dsethist",*roovar1);
      TH1* bkg1dsethist = bkg1dset->createHistogram("bkg1dsethist",*roovar1);
      TH1* sig2dsethist = sig2dset->createHistogram("sig2dsethist",*roovar2);
      TH1* bkg2dsethist = bkg2dset->createHistogram("bkg2dsethist",*roovar2);
      TString n="data"; if (ismc==1) n="mctrue"; if (ismc==2) n="mcdata";
      TFile *outf = new TFile(Form("outtemplates_%s_%s_bin%d.root",splitting.Data(),n.Data(),template_bin),"recreate");
      outf->cd();
      sig1dsethist->Write();
      bkg1dsethist->Write();
      sig2dsethist->Write();
      bkg2dsethist->Write();
      outf->Close();
      return NULL;
    }



    RooHistPdf *sig1pdf=NULL;
    RooHistPdf *sig2pdf=NULL;
    RooHistPdf *bkg1pdf=NULL;
    RooHistPdf *bkg2pdf=NULL;

    RooRealVar *lm_sig1 = new RooRealVar("lm_sig1","lm_sig1",0.5,0,1);
    RooRealVar *lm_sig2 = new RooRealVar("lm_sig2","lm_sig2",0.5,0,1);
    RooRealVar *lm_bkg1 = new RooRealVar("lm_bkg1","lm_bkg1",0.5,0,1);
    RooRealVar *lm_bkg2 = new RooRealVar("lm_bkg2","lm_bkg2",0.5,0,1);

    
    RooHistPdf *sig1pdf_eta1binned[n_eta_cats];
    RooHistPdf *bkg1pdf_eta1binned[n_eta_cats];
    RooHistPdf *sig2pdf_eta2binned[n_eta_cats];
    RooHistPdf *bkg2pdf_eta2binned[n_eta_cats];


    if (!ok_for_recycle){

    RooDataHist *sig1dhist = new RooDataHist("sig1dhist","sig1dhist",*roovar1,*sig1dset);
    RooDataHist *sig2dhist = new RooDataHist("sig2dhist","sig2dhist",*roovar2,*sig2dset);
    RooDataHist *bkg1dhist = new RooDataHist("bkg1dhist","bkg1dhist",*roovar1,*bkg1dset);
    RooDataHist *bkg2dhist = new RooDataHist("bkg2dhist","bkg2dhist",*roovar2,*bkg2dset);
    sig1pdf = (RooHistPdf*)(new RooHistPdf("sig1pdf","sig1pdf",*roovar1,*sig1dhist));
    sig2pdf = (RooHistPdf*)(new RooHistPdf("sig2pdf","sig2pdf",*roovar2,*sig2dhist));
    bkg1pdf = (RooHistPdf*)(new RooHistPdf("bkg1pdf","bkg1pdf",*roovar1,*bkg1dhist));
    bkg2pdf = (RooHistPdf*)(new RooHistPdf("bkg2pdf","bkg2pdf",*roovar2,*bkg2dhist));
    
    RooDataHist *sig1dhist_eta1binned[n_eta_cats];
    RooDataHist *bkg1dhist_eta1binned[n_eta_cats];
    RooDataHist *sig2dhist_eta2binned[n_eta_cats];
    RooDataHist *bkg2dhist_eta2binned[n_eta_cats];

    for (int k=0; k<n_eta_cats; k++){
	sig1dhist_eta1binned[k] = new RooDataHist(Form("sig1dhist_eta1binned_%d",k),Form("sig1dhist_eta1binned_%d",k),*roovar1,*sig1dset_eta1binned[k]);
	bkg1dhist_eta1binned[k] = new RooDataHist(Form("bkg1dhist_eta1binned_%d",k),Form("bkg1dhist_eta1binned_%d",k),*roovar1,*bkg1dset_eta1binned[k]);
	sig2dhist_eta2binned[k] = new RooDataHist(Form("sig2dhist_eta2binned_%d",k),Form("sig2dhist_eta2binned_%d",k),*roovar2,*sig2dset_eta2binned[k]);
	bkg2dhist_eta2binned[k] = new RooDataHist(Form("bkg2dhist_eta2binned_%d",k),Form("bkg2dhist_eta2binned_%d",k),*roovar2,*bkg2dset_eta2binned[k]);
	sig1pdf_eta1binned[k] = (RooHistPdf*)(new RooHistPdf(Form("sig1pdf_eta1binned_%d",k),Form("sig1pdf_eta1binned_%d",k),*roovar1,*sig1dhist_eta1binned[k]));
	bkg1pdf_eta1binned[k] = (RooHistPdf*)(new RooHistPdf(Form("bkg1pdf_eta1binned_%d",k),Form("bkg1pdf_eta1binned_%d",k),*roovar1,*bkg1dhist_eta1binned[k]));
	sig2pdf_eta2binned[k] = (RooHistPdf*)(new RooHistPdf(Form("sig2pdf_eta2binned_%d",k),Form("sig2pdf_eta2binned_%d",k),*roovar2,*sig2dhist_eta2binned[k]));
	bkg2pdf_eta2binned[k] = (RooHistPdf*)(new RooHistPdf(Form("bkg2pdf_eta2binned_%d",k),Form("bkg2pdf_eta2binned_%d",k),*roovar2,*bkg2dhist_eta2binned[k]));
    }


    }
    else {
      /*
      lm_sig1->setVal(lm_sig1_recycle);
      lm_sig2->setVal(lm_sig2_recycle);
      lm_bkg1->setVal(lm_bkg1_recycle);
      lm_bkg2->setVal(lm_bkg2_recycle);
      
      sig1pdf=sig1pdf_recycle;
      sig2pdf=sig2pdf_recycle;
      bkg1pdf=bkg1pdf_recycle;
      bkg2pdf=bkg2pdf_recycle;
      */
    }

    /*
    lm_sig1->setConstant(kTRUE);
    lm_sig2->setConstant(kTRUE);
    lm_bkg1->setConstant(kTRUE);
    lm_bkg2->setConstant(kTRUE);
    */

    {
    TCanvas *c0_incl = new TCanvas(Form("c0_incl"),Form("c0_incl"),1200,800);
    c0_incl->Divide(2,2);
    c0_incl->cd(1);
    RooPlot *frame01sig = roovar1->frame();
    sig1dset->plotOn(frame01sig);
    sig1pdf->plotOn(frame01sig);
    frame01sig->Draw();
    c0_incl->GetPad(1)->SetLogy(1);
    c0_incl->cd(2);
    RooPlot *frame02sig = roovar2->frame();
    sig2dset->plotOn(frame02sig);
    sig2pdf->plotOn(frame02sig);
    frame02sig->Draw();
    //    c0_incl->GetPad(2)->SetLogy(1);
    c0_incl->cd(3);
    RooPlot *frame01bkg = roovar1->frame();
    bkg1dset->plotOn(frame01bkg);
    bkg1pdf->plotOn(frame01bkg);
    frame01bkg->Draw();
    c0_incl->GetPad(3)->SetLogy(1);
    c0_incl->cd(4);
    RooPlot *frame02bkg = roovar2->frame();
    bkg2dset->plotOn(frame02bkg);
    bkg2pdf->plotOn(frame02bkg);
    frame02bkg->Draw();
    //    c0_incl->GetPad(4)->SetLogy(1);
    }

    
    for (int k=0; k<n_eta_cats; k++) {
      TCanvas *c0 = new TCanvas(Form("c0_%d",k),Form("c0_%d",k),1200,800);
    c0->Divide(2,2);
    c0->cd(1);
    RooPlot *frame01sig = roovar1->frame();
    sig1dset_eta1binned[k]->plotOn(frame01sig);
    sig1pdf_eta1binned[k]->plotOn(frame01sig);
    frame01sig->Draw();
    c0->GetPad(1)->SetLogy(1);
    c0->cd(2);
    RooPlot *frame02sig = roovar2->frame();
    sig2dset_eta2binned[k]->plotOn(frame02sig);
    sig2pdf_eta2binned[k]->plotOn(frame02sig);
    frame02sig->Draw();
    //    c0->GetPad(2)->SetLogy(1);
    c0->cd(3);
    RooPlot *frame01bkg = roovar1->frame();
    bkg1dset_eta1binned[k]->plotOn(frame01bkg);
    bkg1pdf_eta1binned[k]->plotOn(frame01bkg);
    frame01bkg->Draw();
    c0->GetPad(3)->SetLogy(1);
    c0->cd(4);
    RooPlot *frame02bkg = roovar2->frame();
    bkg2dset_eta2binned[k]->plotOn(frame02bkg);
    bkg2pdf_eta2binned[k]->plotOn(frame02bkg);
    frame02bkg->Draw();
    //    c0->GetPad(4)->SetLogy(1);

    }
    

//    c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
//    c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));
//    c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));
//    c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));






    RooFormulaVar *fsig1 = new RooFormulaVar("fsig1","fsig1","j1",RooArgList(*j1));
    //    RooFormulaVar *fbkg1 = new RooFormulaVar("fbkg1","fbkg1","1-fsig1",RooArgList(*fsig1));
    RooFormulaVar *fsig2 = new RooFormulaVar("fsig2","fsig2","@0", (sym) ? RooArgList(*j1) : RooArgList(*j2) );
    //    RooFormulaVar *fbkg2 = new RooFormulaVar("fbkg2","fbkg2","1-fsig2",RooArgList(*fsig2));

    RooFormulaVar *fsigsig = new RooFormulaVar("fsigsig","fsigsig","pp",RooArgList(*pp));
    RooFormulaVar *fsigbkg = new RooFormulaVar("fsigbkg","fsigbkg","fsig1-pp",RooArgList(*fsig1,*pp));  
    RooFormulaVar *fbkgsig = new RooFormulaVar("fbkgsig","fbkgsig","fsig2-pp",RooArgList(*fsig2,*pp));
    RooFormulaVar *fbkgbkg = new RooFormulaVar("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg-fbkgsig",RooArgList(*fsigsig,*fsigbkg,*fbkgsig));
    
    RooProdPdf *sigsigpdf = new RooProdPdf("sigsigpdf","sigsigpdf",RooArgList(*sig1pdf,*sig2pdf));
    RooProdPdf *sigbkgpdf = new RooProdPdf("sigbkgpdf","sigbkgpdf",RooArgList(*sig1pdf,*bkg2pdf));
    RooProdPdf *bkgsigpdf = new RooProdPdf("bkgsigpdf","bkgsigpdf",RooArgList(*bkg1pdf,*sig2pdf));
    RooProdPdf *bkgbkgpdf = new RooProdPdf("bkgbkgpdf","bkgbkgpdf",RooArgList(*bkg1pdf,*bkg2pdf));

    RooProdPdf *sigsigpdf_eta1eta2binned[n_eta1eta2_cats];
    RooProdPdf *sigbkgpdf_eta1eta2binned[n_eta1eta2_cats];
    RooProdPdf *bkgsigpdf_eta1eta2binned[n_eta1eta2_cats];
    RooProdPdf *bkgbkgpdf_eta1eta2binned[n_eta1eta2_cats];

    for (int k1=0; k1<n_eta_cats; k1++){
      for (int k2=0; k2<n_eta_cats; k2++){
	int k = k1*n_eta_cats+k2;
	sigsigpdf_eta1eta2binned[k] = new RooProdPdf(Form("sigsigpdf_%d",k),Form("sigsigpdf_%d",k),RooArgList(*sig1pdf_eta1binned[k1],*sig2pdf_eta2binned[k2]));
	sigbkgpdf_eta1eta2binned[k] = new RooProdPdf(Form("sigbkgpdf_%d",k),Form("sigbkgpdf_%d",k),RooArgList(*sig1pdf_eta1binned[k1],*bkg2pdf_eta2binned[k2]));
	bkgsigpdf_eta1eta2binned[k] = new RooProdPdf(Form("bkgsigpdf_%d",k),Form("bkgsigpdf_%d",k),RooArgList(*bkg1pdf_eta1binned[k1],*sig2pdf_eta2binned[k2]));
	bkgbkgpdf_eta1eta2binned[k] = new RooProdPdf(Form("bkgbkgpdf_%d",k),Form("bkgbkgpdf_%d",k),RooArgList(*bkg1pdf_eta1binned[k1],*bkg2pdf_eta2binned[k2]));
      }
    }
    
    //    RooRealVar *nevents = new RooRealVar("nevents","nevents",dataset->sumEntries(),0,1e6);


    RooAddPdf *model_axis1 = new RooAddPdf("model_axis1","model_axis1",RooArgList(*sig1pdf,*bkg1pdf),RooArgList(*fsig1));
    RooAddPdf *model_axis2 = new RooAddPdf("model_axis2","model_axis2",RooArgList(*sig2pdf,*bkg2pdf),RooArgList(*fsig2));

    RooFitResult *firstpass;




      //    RooExtendPdf *model_axis1_extended = new RooExtendPdf("model_axis1_extended","model_axis1_extended",*model_axis1,*nevents);
      //    RooNLLVar *model_axis1_extended_nll = new RooNLLVar("model_axis1_extended_nll","model_axis1_extended_nll",*model_axis1_extended,*dataset_axis1,NumCPU(numcpu/2));
    RooNLLVar *model_axis1_noextended_nll = new RooNLLVar("model_axis1_noextended_nll","model_axis1_noextended_nll",*model_axis1,*dataset_axis1,NumCPU(numcpu/2));
    RooNLLVar *model_axis1_nll = model_axis1_noextended_nll;
    //    RooExtendPdf *model_axis2_extended = new RooExtendPdf("model_axis2_extended","model_axis2_extended",*model_axis2,*nevents);
    //    RooNLLVar *model_axis2_extended_nll = new RooNLLVar("model_axis2_extended_nll","model_axis2_extended_nll",*model_axis2_extended,*dataset_axis2,NumCPU(numcpu/2));
    RooNLLVar *model_axis2_noextended_nll = new RooNLLVar("model_axis2_noextended_nll","model_axis2_noextended_nll",*model_axis2,*dataset_axis2,NumCPU(numcpu/2));
    RooNLLVar *model_axis2_nll = model_axis2_noextended_nll;
    
    RooAddition *model_2axes_nll = new RooAddition("model_2axes_nll","model_2axes_nll",RooArgSet(*model_axis1_nll,*model_axis2_nll));
    //    RooAddition *model_2D_nll_constrained = new RooAddition("model_2D_nll_constrained","model_2D_nll_constrained",RooArgSet(*model_2D_extended_nll,*constrain_asym_nll));
      
      
      RooMinimizer *minuit_firstpass = new RooMinimizer(*model_2axes_nll);
      minuit_firstpass->migrad();
      minuit_firstpass->hesse();
      //      minuit_firstpass->minos();
      firstpass = minuit_firstpass->save("firstpass","firstpass");
      firstpass->Print();
      
      /*
      j1->setVal(6.00750e-01);
      j2->setVal(4.77380e-01);
      */
      delete model_axis1_nll;
      delete model_axis2_nll;



    TCanvas *c1 = new TCanvas("c1","c1",1200,800);
    c1->Divide(2);
    c1->cd(1);
    RooPlot *frame1bla = roovar1->frame();
    dataset_axis1->plotOn(frame1bla);
    model_axis1->plotOn(frame1bla);
    model_axis1->plotOn(frame1bla,Components("sig1pdf"),LineStyle(kDashed),LineColor(kRed));
    model_axis1->plotOn(frame1bla,Components("bkg1pdf"),LineStyle(kDashed),LineColor(kBlack));
    frame1bla->Draw();
    c1->GetPad(1)->SetLogy(1);
    c1->cd(2);
    RooPlot *frame2bla = roovar2->frame();
    dataset_axis2->plotOn(frame2bla);
    model_axis2->plotOn(frame2bla);
    model_axis2->plotOn(frame2bla,Components("sig2pdf"),LineStyle(kDashed),LineColor(kRed));
    model_axis2->plotOn(frame2bla,Components("bkg2pdf"),LineStyle(kDashed),LineColor(kBlack));
    frame2bla->Draw();
    c1->GetPad(2)->SetLogy(1);


    c1->SaveAs(Form("plots/fittingplot1_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
    c1->SaveAs(Form("plots/fittingplot1_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));
    c1->SaveAs(Form("plots/fittingplot1_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));
    c1->SaveAs(Form("plots/fittingplot1_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));




    float lowerbounds[4]={0,fsig1->getVal()-1,fsig2->getVal()-1,fsig1->getVal()+fsig2->getVal()-1};
    float upperbounds[4]={1,fsig1->getVal(),fsig2->getVal(),fsig1->getVal()+fsig2->getVal()};


    float minpp = TMath::MaxElement(4,lowerbounds);
    float maxpp = TMath::MinElement(4,upperbounds);
    pp->setVal((minpp+maxpp)/2);


    // DEBUUUUUUUUUUUUUUUUUUUUG
    //    std::cout << "j1 " << j1->getVal() << " " << j1->getPropagatedError(*firstpass) << std::endl;
    //    std::cout << "j2 " << j2->getVal() << " " << j2->getPropagatedError(*firstpass) << std::endl;

    
    {

      float nsigma_tolerance = 3;

      std::cout << "setting constrain pp val at " << pp->getVal() << " between " << minpp << " and " << maxpp << std::endl;
      pp->setRange(minpp,maxpp);

      std::cout << "setting constrain j1 val at " << j1->getVal() << " between " << j1->getVal()-nsigma_tolerance*j1->getPropagatedError(*firstpass) << " and " << j1->getVal()+nsigma_tolerance*j1->getPropagatedError(*firstpass) << std::endl;
      j1->setRange(j1->getVal()-nsigma_tolerance*j1->getPropagatedError(*firstpass),j1->getVal()+nsigma_tolerance*j1->getPropagatedError(*firstpass));
      
      if (!sym){
	std::cout << "setting constrain j2 val at " << j2->getVal() << " between " << j2->getVal()-nsigma_tolerance*j2->getPropagatedError(*firstpass) << " and " << j2->getVal()+nsigma_tolerance*j2->getPropagatedError(*firstpass) << std::endl;
	j2->setRange(j2->getVal()-nsigma_tolerance*j2->getPropagatedError(*firstpass),j2->getVal()+nsigma_tolerance*j2->getPropagatedError(*firstpass));
      }
     
    }
    

    RooAddPdf *model_2D_uncorrelated = new RooAddPdf("model_2D_uncorrelated","model_2D_uncorrelated",RooArgList(*sigsigpdf,*sigbkgpdf,*bkgsigpdf,*bkgbkgpdf),RooArgList(*fsigsig,*fsigbkg,*fbkgsig,*fbkgbkg),kFALSE);
    //    RooExtendPdf *model_2D_uncorrelated_extended = new RooExtendPdf("model_2D_uncorrelated_extended","model_2D_uncorrelated_extended",*model_2D_uncorrelated,*nevents);
    //    RooNLLVar *model_2D_uncorrelated_extended_nll = new RooNLLVar("model_2D_uncorrelated_extended_nll","model_2D_uncorrelated_extended_nll",*model_2D_uncorrelated_extended,*dataset,NumCPU(numcpu));
    RooNLLVar *model_2D_uncorrelated_noextended_nll = new RooNLLVar("model_2D_uncorrelated_noextended_nll","model_2D_uncorrelated_noextended_nll",*model_2D_uncorrelated,*dataset,NumCPU(numcpu));


    RooSimultaneous *model_2D_combinedeta1eta2binned = new RooSimultaneous("model_2D_eta1eta2binned","model_2D_eta1eta2binned",*eta1eta2cat);
    RooAddPdf *model_2D_eta1eta2binned[n_eta1eta2_cats];
    for (int k=0; k<n_eta1eta2_cats; k++) {
      model_2D_eta1eta2binned[k] = new RooAddPdf(Form("model_2D_eta1eta2binned_%d",k),Form("model_2D_eta1eta2binned_%d",k),RooArgList(*sigsigpdf_eta1eta2binned[k],*sigbkgpdf_eta1eta2binned[k],*bkgsigpdf_eta1eta2binned[k],*bkgbkgpdf_eta1eta2binned[k]),RooArgList(*fsigsig,*fsigbkg,*fbkgsig,*fbkgbkg),kFALSE);
      model_2D_combinedeta1eta2binned->addPdf(*model_2D_eta1eta2binned[k],Form("cat%d",k));
    }
    RooNLLVar *model_2D_combinedeta1eta2binned_noextended_nll = new RooNLLVar("model_2D_combinedeta1eta2binned_noextended_nll","model_2D_combinedeta1eta2binned_noextended_nll",*model_2D_combinedeta1eta2binned,*dataset_combinedeta1eta2binned,NumCPU(numcpu));
    RooNLLVar *model_2D_combinedeta1eta2binned_nll = model_2D_combinedeta1eta2binned_noextended_nll;


    RooArgList sigsig_eta1eta2binned_pdflist;
    RooArgList sigbkg_eta1eta2binned_pdflist;
    RooArgList bkgsig_eta1eta2binned_pdflist;
    RooArgList bkgbkg_eta1eta2binned_pdflist;
    RooArgList eta1eta2binned_coeflist;
    for (int k=0; k<n_eta1eta2_cats; k++) sigsig_eta1eta2binned_pdflist.add(*(sigsigpdf_eta1eta2binned[k]));
    for (int k=0; k<n_eta1eta2_cats; k++) sigbkg_eta1eta2binned_pdflist.add(*(sigbkgpdf_eta1eta2binned[k]));
    for (int k=0; k<n_eta1eta2_cats; k++) bkgsig_eta1eta2binned_pdflist.add(*(bkgsigpdf_eta1eta2binned[k]));
    for (int k=0; k<n_eta1eta2_cats; k++) bkgbkg_eta1eta2binned_pdflist.add(*(bkgbkgpdf_eta1eta2binned[k]));
    for (int k=0; k<n_eta1eta2_cats; k++) eta1eta2binned_coeflist.add(RooRealConstant::value(dataset_eta1eta2binned[k]->sumEntries()));
    RooAddPdf *sigsigpdf_summedeta1eta2binned = new RooAddPdf("sigsigpdf_summedeta1eta2binned","sigsigpdf_summedeta1eta2binned",sigsig_eta1eta2binned_pdflist,eta1eta2binned_coeflist);
    RooAddPdf *sigbkgpdf_summedeta1eta2binned = new RooAddPdf("sigbkgpdf_summedeta1eta2binned","sigbkgpdf_summedeta1eta2binned",sigbkg_eta1eta2binned_pdflist,eta1eta2binned_coeflist);
    RooAddPdf *bkgsigpdf_summedeta1eta2binned = new RooAddPdf("bkgsigpdf_summedeta1eta2binned","bkgsigpdf_summedeta1eta2binned",bkgsig_eta1eta2binned_pdflist,eta1eta2binned_coeflist);
    RooAddPdf *bkgbkgpdf_summedeta1eta2binned = new RooAddPdf("bkgbkgpdf_summedeta1eta2binned","bkgbkgpdf_summedeta1eta2binned",bkgbkg_eta1eta2binned_pdflist,eta1eta2binned_coeflist);
    RooAddPdf *model_2D_summedeta1eta2binned = new RooAddPdf("model_2D_summedeta1eta2binned","model_2D_summedeta1eta2binned",RooArgList(*sigsigpdf_summedeta1eta2binned,*sigbkgpdf_summedeta1eta2binned,*bkgsigpdf_summedeta1eta2binned,*bkgbkgpdf_summedeta1eta2binned),RooArgList(*fsigsig,*fsigbkg,*fbkgsig,*fbkgbkg),kFALSE);

    RooNLLVar *model_2D_summedeta1eta2binned_noextended_nll = new RooNLLVar("model_2D_summedeta1eta2binned_noextended_nll","model_2D_summedeta1eta2binned_noextended_nll",*model_2D_summedeta1eta2binned,*dataset,NumCPU(numcpu));
    RooNLLVar *model_2D_summedeta1eta2binned_nll = model_2D_summedeta1eta2binned_noextended_nll;


    //    RooNLLVar *model_2D_nll = model_2D_uncorrelated_noextended_nll;  // THIS IS THE MODEL THAT IS GOING TO BE USED IN PHASE 2
    RooNLLVar *model_2D_nll = model_2D_summedeta1eta2binned_nll;  // THIS IS THE MODEL THAT IS GOING TO BE USED IN PHASE 2

    
    RooMinimizer *minuit_secondpass = new RooMinimizer(*model_2D_nll);
    //    minuit_secondpass->migrad();
    //    minuit_secondpass->hesse();
    RooFitResult *secondpass;
    //    secondpass = minuit_secondpass->save("secondpass","secondpass");
    //    secondpass->Print();


    
//    TCanvas *c2fits = new TCanvas("c2fits","c2fits",1200,800);
//    c2fits->Divide(abs(sqrt(n_eta1eta2_cats))+1,abs(sqrt(n_eta1eta2_cats))+1);
//    //    c2fits->Divide(9);
//    for (int k=0; k<n_eta1eta2_cats; k++){
//      c2fits->cd(k+1);
//      RooPlot *framek = roovar1->frame();
//      eta1eta2cat->setLabel(Form("cat%d",k));
//      dataset_eta1eta2binned[k]->plotOn(framek); 
//      model_2D_eta1eta2binned[k]->plotOn(framek);
//      model_2D_eta1eta2binned[k]->plotOn(framek,Components(Form("sigsigpdf_%d",k)),LineStyle(kDashed),LineColor(kRed));
//      model_2D_eta1eta2binned[k]->plotOn(framek,Components(Form("sigbkgpdf_%d",k)),LineStyle(kDashed),LineColor(kGreen));
//      model_2D_eta1eta2binned[k]->plotOn(framek,Components(Form("bkgsigpdf_%d",k)),LineStyle(kDashed),LineColor(kGreen+2));
//      model_2D_eta1eta2binned[k]->plotOn(framek,Components(Form("bkgbkgpdf_%d",k)),LineStyle(kDashed),LineColor(kBlack));
//      framek->Draw();
//    }
//    c2fits->Update();



    TCanvas *c2 = new TCanvas("c2","c2",1200,800);
    c2->Divide(2,2);
    
    c2->cd(1);
    RooPlot *frame1final = roovar1->frame();
    dataset->plotOn(frame1final);
    model_2D_summedeta1eta2binned->plotOn(frame1final);
    model_2D_summedeta1eta2binned->plotOn(frame1final,Components("sigsigpdf_summedeta1eta2binned"),LineStyle(kDashed),LineColor(kRed));
    model_2D_summedeta1eta2binned->plotOn(frame1final,Components("sigbkgpdf_summedeta1eta2binned"),LineStyle(kDashed),LineColor(kGreen));
    model_2D_summedeta1eta2binned->plotOn(frame1final,Components("bkgsigpdf_summedeta1eta2binned"),LineStyle(kDashed),LineColor(kGreen+2));
    model_2D_summedeta1eta2binned->plotOn(frame1final,Components("bkgbkgpdf_summedeta1eta2binned"),LineStyle(kDashed),LineColor(kBlack));
    frame1final->Draw();
    //    c2->GetPad(1)->SetLogy(1);

    c2->cd(2);
    RooPlot *frame2final = roovar2->frame();
    dataset->plotOn(frame2final);
    model_2D_summedeta1eta2binned->plotOn(frame2final);
    model_2D_summedeta1eta2binned->plotOn(frame2final,Components("sigsigpdf_summedeta1eta2binned"),LineStyle(kDashed),LineColor(kRed));
    model_2D_summedeta1eta2binned->plotOn(frame2final,Components("sigbkgpdf_summedeta1eta2binned"),LineStyle(kDashed),LineColor(kGreen));
    model_2D_summedeta1eta2binned->plotOn(frame2final,Components("bkgsigpdf_summedeta1eta2binned"),LineStyle(kDashed),LineColor(kGreen+2));
    model_2D_summedeta1eta2binned->plotOn(frame2final,Components("bkgbkgpdf_summedeta1eta2binned"),LineStyle(kDashed),LineColor(kBlack));
    frame2final->Draw();
    //    c2->GetPad(2)->SetLogy(1);

    /*
    c2->cd(3);
    RooPlot *ppnllplot = pp->frame();
    model_2D_nll->plotOn(ppnllplot,ShiftToZero());
    ppnllplot->Draw();
    */

    c2->cd(4);

    int colors[3] = {kBlack, kRed, kGreen};
    TH2F *h2[3];
    TH1F *h3[3];
    for (int i=0; i<3; i++) h2[i] = new TH2F(Form("h2%d",i),Form("h2%d",i),60,leftrange,rightrange,60,leftrange,rightrange);
    for (int i=0; i<3; i++) h3[i] = new TH1F(Form("h3%d",i),Form("h3%d",i),60,leftrange*sqrt(2),rightrange*sqrt(2));


    TRandom3 *randomgen = new TRandom3(0);
    RooRealVar *diagvar = new RooRealVar("diagvar","diagvar",0,leftrange*sqrt(2),rightrange*sqrt(2));
    //    diagvar->setBins(n_histobins);

    RooDataSet *diag_dataset = new RooDataSet("diag_dataset","diag_dataset",*diagvar);
    for (int i=0; i<dataset->numEntries(); i++){
      float r1 = dataset->get(i)->getRealValue("roovar1");
      float r2 = dataset->get(i)->getRealValue("roovar2");
      float w = dataset->store()->weight(i);
      diagvar->setVal((r1+r2)/sqrt(2));
      diag_dataset->add(*diagvar,w);
      h2[0]->Fill(r1,r2,w);
    }

    RooDataSet *rand_uncorrelated = model_2D_uncorrelated->generate(RooArgSet(*roovar1,*roovar2,*rooweight),1e6);
    RooDataSet *rand_diag_uncorrelated = new RooDataSet("rand_diag_uncorrelated","rand_diag_uncorrelated",RooArgList(*diagvar,*rooweight),WeightVar(*rooweight));
    std::cout << "rand_diag_uncorrelated" << std::endl;
    ofstream file_uncorr("file_uncorr.txt");
    for (int i=0; i<rand_uncorrelated->numEntries(); i++){
      float r1 = rand_uncorrelated->get(i)->getRealValue("roovar1");
      float r2 = rand_uncorrelated->get(i)->getRealValue("roovar2");
      float w = rand_uncorrelated->store()->weight(i);
      rooweight->setVal(w);
      diagvar->setVal((r1+r2)/sqrt(2));
      rand_diag_uncorrelated->add(*diagvar,w);
      h2[1]->Fill(r1,r2,w);
      file_uncorr << r1 << " " << r2 << " " << (r1+r2)/sqrt(2) << " " << w << std::endl;
    }
    RooDataHist *randhist_diag_uncorrelated = new RooDataHist("randhist_diag_uncorrelated","randhist_diag_uncorrelated",*diagvar,*rand_diag_uncorrelated);
    RooHistPdf *randpdf_diag_uncorrelated = new RooHistPdf("randpdf_diag_uncorrelated","randpdf_diag_uncorrelated",*diagvar,*randhist_diag_uncorrelated);
    file_uncorr.close();    

    RooDataSet *rand_correlated = model_2D_summedeta1eta2binned->generate(RooArgSet(*roovar1,*roovar2),1e6);
    RooDataSet *rand_diag_correlated = new RooDataSet("rand_diag_correlated","rand_diag_correlated",RooArgList(*diagvar,*rooweight),WeightVar(*rooweight));
    std::cout << "rand_diag_correlated" << std::endl;
    ofstream file_corr("file_corr.txt");
    for (int i=0; i<rand_correlated->numEntries(); i++){
      float r1 = rand_correlated->get(i)->getRealValue("roovar1");
      float r2 = rand_correlated->get(i)->getRealValue("roovar2");
      float w = rand_correlated->store()->weight(i);
      rooweight->setVal(w);
      diagvar->setVal((r1+r2)/sqrt(2));
      rand_diag_correlated->add(*diagvar,w);
      h2[2]->Fill(r1,r2,w);
      file_corr << r1 << " " << r2 << " " << (r1+r2)/sqrt(2) << " " << w << std::endl;
    }
    RooDataHist *randhist_diag_correlated = new RooDataHist("randhist_diag_correlated","randhist_diag_correlated",*diagvar,*rand_diag_correlated);
    RooHistPdf *randpdf_diag_correlated = new RooHistPdf("randpdf_diag_correlated","randpdf_diag_correlated",*diagvar,*randhist_diag_correlated);
    file_corr.close();

    /*
    RooPlot *diagplot = diagvar->frame();
    diag_dataset->plotOn(diagplot);
    randpdf_diag_uncorrelated->plotOn(diagplot,LineColor(kRed));
    randpdf_diag_correlated->plotOn(diagplot,LineColor(kGreen));
    diagplot->Draw();
    */

    for (int i=0; i<3; i++) {

      for (int k=0; k<h2[i]->GetNbinsX(); k++) for (int l=0; l<h2[i]->GetNbinsY(); l++) h2[i]->SetBinError(k+1,l+1,0);
      h2[i]->Scale(1.0/h2[i]->Integral());
      h2[i]->SetLineColor(colors[i]);
      h2[i]->SetMarkerColor(colors[i]);
      h2[i]->SetMarkerStyle(20);

      for (int k=0; k<h2[i]->GetNbinsX(); k++) for (int l=0; l<h2[i]->GetNbinsY(); l++) {
	  float xcenter = h2[i]->GetBinCenter(k);
	  float ycenter = h2[i]->GetBinCenter(l);
	  h3[i]->Fill((xcenter+ycenter)/sqrt(2),h2[i]->GetBinContent(k,l));
	    }
      for (int k=0; k<h3[i]->GetNbinsX(); k++) h3[i]->SetBinError(k+1,0);
      h3[i]->Scale(1.0/h3[i]->Integral());
      h3[i]->SetLineColor(colors[i]);
      h3[i]->SetMarkerColor(colors[i]);
      h3[i]->SetMarkerStyle(20);

    }



//    randpdf_diag_sigsig->plotOn(diagplot,Normalization(1,RooAbsPdf::NumEvent),LineColor(kRed));
//    randpdf_diag_sigbkg->plotOn(diagplot,Normalization(1,RooAbsPdf::NumEvent),LineColor(kGreen));
//    randpdf_diag_bkgsig->plotOn(diagplot,Normalization(1,RooAbsPdf::NumEvent),LineColor(kGreen+2));
//    randpdf_diag_bkgbkg->plotOn(diagplot,Normalization(1,RooAbsPdf::NumEvent),LineColor(kBlack));
//    randpdf_diag_sigsig->plotOn(diagplot,Normalization(fsigsig->getVal(),RooAbsPdf::Relative),LineColor(kRed));
//    randpdf_diag_sigbkg->plotOn(diagplot,Normalization(fsigbkg->getVal(),RooAbsPdf::Relative),LineColor(kGreen));
//    randpdf_diag_bkgsig->plotOn(diagplot,Normalization(fbkgsig->getVal(),RooAbsPdf::Relative),LineColor(kGreen+2));
//    randpdf_diag_bkgbkg->plotOn(diagplot,Normalization(fbkgbkg->getVal(),RooAbsPdf::Relative),LineColor(kBlack));

    c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
    c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.pdf",splitting.Data(),diffvariable.Data(),bin));
    c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.jpg",splitting.Data(),diffvariable.Data(),bin));
    c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));

    TCanvas *c3 = new TCanvas("c3","c3",1200,800);
    c3->Divide(3);
    c3->cd(1);
    h2[0]->ProjectionX()->Draw();
    h2[1]->ProjectionX()->Draw("same");
    h2[2]->ProjectionX()->Draw("same");
    c3->cd(2);
    h2[0]->ProjectionY()->Draw();
    h2[1]->ProjectionY()->Draw("same");
    h2[2]->ProjectionY()->Draw("same");
    c3->cd(3);
    h3[0]->Draw();
    h3[1]->Draw("same");
    h3[2]->Draw("same");
    c3->Update();
   

    return;

//    std::cout << "expecting purities = " << pp_init << " " << pf_init << " " << fp_init << " " << 1-pp_init-pf_init-fp_init << std::endl;
    std::cout << "pp " << fsigsig->getVal() << " " << fsigsig->getPropagatedError(*secondpass) << std::endl; 
    std::cout << "pf " << fsigbkg->getVal() << " " << fsigbkg->getPropagatedError(*secondpass) << std::endl; 
    std::cout << "fp " << fbkgsig->getVal() << " " << fbkgsig->getPropagatedError(*secondpass) << std::endl; 
    std::cout << "ff " << fbkgbkg->getVal() << " " << fbkgbkg->getPropagatedError(*secondpass) << std::endl; 


//    out->fr=fr;
//    out->tot_events=nevents->getVal();




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


    RooWorkspace *wspace = new RooWorkspace("fittingwspace","fittingwspace");
    wspace->import(*firstpass);
    wspace->import(*secondpass);
    wspace->writeToFile(Form("plots/fittingwspace_%s_%s_b%d.root",splitting.Data(),diffvariable.Data(),bin));

    delete model_2D_nll;

    return out;

    

//  RooMCStudy *mcstudy = NULL;
//
//  if (study_templates){
//
//    rf1.setVal(0.6);
//    rf2.setVal(1./2);
//    int howmanyevents = 2e+3;
//    int howmanytoys = 500;
//
//    if (!study_templates_plotting){
//    if (dosingle) mcstudy = new RooMCStudy(*model_extended,RooArgSet(roovar1),Silence(),Extended(),FitOptions(Save(kTRUE),PrintEvalErrors(0)));
//    else mcstudy = new RooMCStudy(*model_extended,RooArgSet(roovar1,roovar2),Silence(),Extended(),FitOptions(Save(kTRUE),PrintEvalErrors(0)));
//    mcstudy->generateAndFit(howmanytoys, howmanyevents);
//
//    TCanvas *c1 = new TCanvas();
//    c1->SetWindowSize(800,600);
//    c1->Divide(3);
//
//    c1->cd(1);    
//    RooPlot* frame3 = mcstudy->plotPull(rf1,FitGauss(kTRUE));
//    frame3->Draw();
//
//    TVirtualPad* subpad = c1->cd(3);
//    subpad->Divide(2,2);
//
//    subpad->cd(1);    
//    RooPlot* frame1 = mcstudy->plotParam(rf1);
//    frame1->Draw();
//    subpad->cd(2);    
//    RooPlot* frame2 = mcstudy->plotError(rf1);
//    frame2->Draw();
//
//    if (!dosingle){
//      c1->cd(2);    
//      RooPlot* frame6 = mcstudy->plotPull(rf2,FitGauss(kTRUE));
//      frame6->Draw();
//      subpad->cd(3);    
//      RooPlot* frame4 = mcstudy->plotParam(rf2);
//      frame4->Draw();
//      subpad->cd(4);    
//      RooPlot* frame5 = mcstudy->plotError(rf2);
//      frame5->Draw();
//    }
//
//    c1->cd();
//    c1->SaveAs(Form("plots/biasstudy_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
//
//    return out;
//  }
//
//    if (study_templates_plotting){
//    RooDataHist *datahist_toy;
//    datahist_toy = model_extended->generateBinned(RooArgSet(roovar1,roovar2),howmanyevents,Name("Toy_dataset"));
//    datahist=datahist_toy;
//    rf1.setVal(0);
//    rf2.setVal(0);
//    }
//
//    
//
//  }
//
//

/*
  RooFitResult *fitres = model_extended->fitTo(*datahist,Save());

  model_extended->Print();

  std::cout << "---------------------------" << std::endl;
  std::cout << rf1.getVal() << " " << rf2.getVal() << std::endl;

  out->fr=fitres;

  //  if (dosingle) out->tot_events=h_datahist_1D->Integral();
  //  else out->tot_events=h_datahist_2D->Integral();
  out->tot_events=nevents->getVal();

  TCanvas *canv = new TCanvas();
  canv->SetName(Form("fittingplot_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin));
  canv->SetTitle(Form("fittingplot_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin));

  RooPlot *varframe[2];
  varframe[0]=roovar1.frame();
  varframe[1]=roovar2.frame();
  canv->Divide(2);
  for (int i=0; i<2; i++){
    if (dosingle && i==1) continue;
    canv->cd(i+1);
    datahist->plotOn(varframe[i]);
    model_extended->plotOn(varframe[i]);
    if (dosingle){
      model_extended->plotOn(varframe[i],Components("sig1pdf"),LineStyle(kDashed),LineColor(kRed));
      model_extended->plotOn(varframe[i],Components("bkg1pdf"),LineStyle(kDashed),LineColor(kBlack));
    }
    else {
      model_extended->plotOn(varframe[i],Components("sigsigpdf"),LineStyle(kDashed),LineColor(kRed));
      model_extended->plotOn(varframe[i],Components("sigbkgpdf"),LineStyle(kDashed),LineColor(kGreen));
      model_extended->plotOn(varframe[i],Components("bkgbkgpdf"),LineStyle(kDashed),LineColor(kBlack));
    }
    varframe[i]->Draw();
    //    canv->GetPad(i+1)->SetLogy(1);
  }

  canv->SaveAs(Form("plots/fittingplot_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));

  return out;
*/
};




void run_fits(TString inputfilename_ts="templates_sig.root", TString inputfilename_tb="templates_bkg.root",TString inputfilename_d="tobefitted.root", TString diffvariable="", TString splitting=""){

  TH1F::SetDefaultSumw2(kTRUE);

  int bins_to_run=0; 
  float *binsdef=NULL;

  if (diffvariable=="invmass"){
    if (splitting=="EBEB")      bins_to_run=n_templates_invmass_EBEB;
    else if (splitting=="EBEE") bins_to_run=n_templates_invmass_EBEE;
    else if (splitting=="EEEE") bins_to_run=n_templates_invmass_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_invmass_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_invmass_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_invmass_EEEE;
  }
  if (diffvariable=="diphotonpt"){
    if (splitting=="EBEB")      bins_to_run=n_templates_diphotonpt_EBEB;
    else if (splitting=="EBEE") bins_to_run=n_templates_diphotonpt_EBEE;
    else if (splitting=="EEEE") bins_to_run=n_templates_diphotonpt_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_diphotonpt_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_diphotonpt_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_diphotonpt_EEEE;
  }
  if (diffvariable=="costhetastar"){
    if (splitting=="EBEB")      bins_to_run=n_templates_costhetastar_EBEB;
    else if (splitting=="EBEE") bins_to_run=n_templates_costhetastar_EBEE;
    else if (splitting=="EEEE") bins_to_run=n_templates_costhetastar_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_costhetastar_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_costhetastar_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_costhetastar_EEEE;
  }
  if (diffvariable=="dphi"){
    if (splitting=="EBEB")      bins_to_run=n_templates_dphi_EBEB;
    else if (splitting=="EBEE") bins_to_run=n_templates_dphi_EBEE;
    else if (splitting=="EEEE") bins_to_run=n_templates_dphi_EEEE; 
    if (splitting=="EBEB")      binsdef=binsdef_diphoton_dphi_EBEB;
    else if (splitting=="EBEE") binsdef=binsdef_diphoton_dphi_EBEE;
    else if (splitting=="EEEE") binsdef=binsdef_diphoton_dphi_EEEE;
  }
//  if (diffvariable=="dosingle"){
//    dosingle=1;
//    if (splitting=="EB")      bins_to_run=n_templates_EB;
//    else if (splitting=="EE") bins_to_run=n_templates_EE;
//    if (splitting=="EB")      binsdef=binsdef_single_gamma_EB_eta;
//    else if (splitting=="EE") binsdef=binsdef_single_gamma_EE_eta;
//  }

  
  fit_output *fr[n_bins];

  TH1F *purity[4];
  TH1F *eff=NULL;
  TH1F *xsec;

  TFile *eff_file = new TFile("efficiencies.root");
  eff_file->GetObject(Form("w_eff_gg_%s_%s",splitting.Data(),diffvariable.Data()),eff);
  assert (eff!=NULL);
  eff->GetYaxis()->SetTitle("selection/ID efficiency");
  eff->GetXaxis()->SetTitle(diffvariable.Data());

  int colors[4] = {kRed, kGreen, kGreen+2, kBlack};

  for (int i=0; i<4; i++){
    TString name = "purity_";
    if (i==0) name.Append("sigsig"); else if (i==1) name.Append("sigbkg"); else if (i==2) name.Append("bkgsig"); else if (i==3) name.Append("bkgbkg");
    purity[i] = new TH1F(name.Data(),name.Data(),bins_to_run,binsdef);
    purity[i]->SetMarkerStyle(20);
    purity[i]->SetMarkerColor(colors[i]);
    purity[i]->SetLineColor(colors[i]);
    purity[i]->SetLineWidth(2);
    purity[i]->GetYaxis()->SetRangeUser(0,1);
    purity[i]->GetYaxis()->SetTitle("purity");
    purity[i]->GetXaxis()->SetTitle(diffvariable.Data());
  }

    xsec = new TH1F("xsec","xsec",bins_to_run,binsdef);
    xsec->SetMarkerStyle(20);
    xsec->SetMarkerColor(kGreen);
    xsec->SetLineColor(kGreen);
    xsec->SetLineWidth(2);


    //    xsec->GetYaxis()->SetTitle("");
    xsec->GetXaxis()->SetTitle(diffvariable.Data());



  for (int bin=0; bin<bins_to_run; bin++) {

    fr[bin]=fit_dataset(inputfilename_ts.Data(),inputfilename_tb.Data(),inputfilename_d.Data(),diffvariable,splitting,bin);
    
    if (!fr[bin]->fr) continue;

    float intlumi=4519.0;

    float pp = fr[bin]->pp;
    float pp_err = fr[bin]->pp_err;
    float pf = fr[bin]->pf;
    float pf_err = fr[bin]->pf_err;
    float fp = fr[bin]->fp;
    float fp_err = fr[bin]->fp_err;
    float ff = fr[bin]->ff;
    float ff_err = fr[bin]->ff_err;

    purity[0]->SetBinContent(bin+1,pp);
    purity[0]->SetBinError(bin+1,pp_err);
    purity[1]->SetBinContent(bin+1,pf);
    purity[1]->SetBinError(bin+1,pf_err);
    purity[2]->SetBinContent(bin+1,fp);
    purity[2]->SetBinError(bin+1,fp_err);
    purity[3]->SetBinContent(bin+1,ff);
    purity[3]->SetBinError(bin+1,ff_err);
      
    //      xsec->SetBinContent(bin+1,pp->getVal()*fr[bin]->tot_events/xsec->GetBinWidth(bin+1)/intlumi);
    //      xsec->SetBinContent(bin+1,pp->getVal()*fr[bin]->tot_events);
    xsec->SetBinContent(bin+1,pp*fr[bin]->tot_events/xsec->GetBinWidth(bin+1)/intlumi);

    float err1=purity[0]->GetBinError(bin+1)/purity[0]->GetBinContent(bin+1);
    float err2=1.0/sqrt(fr[bin]->tot_events);
    float err=sqrt(err1*err1+err2*err2);
    xsec->SetBinError(bin+1,err*xsec->GetBinContent(bin+1));
    
    //for (int k=0; k<100; k++) std::cout << "WARNING: NO EFFICIENCY CORRECTION!!!" << std::endl;
    xsec->Divide(eff);

  }


  

  TCanvas *output_canv = new TCanvas("output_canv","output_canv");
  output_canv->cd();

  purity[0]->Draw("e1");
  purity[1]->Draw("e1same");
  purity[2]->Draw("e1same");
  purity[3]->Draw("e1same");
  
  output_canv->Update();


  TCanvas *xsec_canv = new TCanvas("xsec_canv","xsec_canv");
  xsec_canv->cd();

  xsec->Draw("e1");

  xsec_canv->Update();


//  TCanvas *eff_canv = new TCanvas("eff_canv","eff_canv");
//  eff_canv->cd();
//  eff->Draw("e1");
//  eff_canv->Update();

  TFile *out1 = new TFile(Form("plots/purity_%s_%s.root",splitting.Data(),diffvariable.Data()),"recreate");
  out1->cd();
  purity[0]->Write();
  purity[1]->Write();
  purity[2]->Write();
  purity[3]->Write();

  TFile *outfile = new TFile(Form("plots/xsec_%s_%s.root",splitting.Data(),diffvariable.Data()),"recreate");
  outfile->cd();
  xsec->Write();
  

};



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

void reweight_pteta(RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  TH2F *hnum = new TH2F("hnum","hnum",55,25,300,25,0,2.5);
  TH2F *hden = new TH2F("hden","hden",55,25,300,25,0,2.5);
  hnum->Sumw2();
  hden->Sumw2();

  const char* ptname=Form("roopt%d",numvar);
  const char* etaname=Form("rooeta%d",numvar);

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill((*dset)->get(i)->getRealValue(ptname),fabs((*dset)->get(i)->getRealValue(etaname)),(*dset)->weight());
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(dsetdestination->get(i)->getRealValue(ptname),fabs(dsetdestination->get(i)->getRealValue(etaname)),dsetdestination->weight());
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH2F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_ptetarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->weight();
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

void reweight_eta(RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  TH1F *hnum = new TH1F("hnum","hnum",25,0,2.5);
  TH1F *hden = new TH1F("hden","hden",25,0,2.5);
  hnum->Sumw2();
  hden->Sumw2();

  const char* etaname=Form("rooeta%d",numvar);

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue(etaname)),(*dset)->weight());
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue(etaname)),dsetdestination->weight());
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_etarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->weight();
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
  std::cout << "Eta rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  //  delete old_dset;

};

void reweight_rho(RooDataSet **dset, RooDataSet *dsetdestination){

  TH1F *hnum = new TH1F("hnum","hnum",30,0,30);
  TH1F *hden = new TH1F("hden","hden",30,0,30);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill((*dset)->get(i)->getRealValue("roorho"),(*dset)->weight());
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(dsetdestination->get(i)->getRealValue("roorho"),dsetdestination->weight());
  }

  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_rhorew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->weight();
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
  //  delete old_dset;

};

void reweight_sigma(RooDataSet **dset, RooDataSet *dsetdestination){

  TH1F *hnum = new TH1F("hnum","hnum",20,0,10);
  TH1F *hden = new TH1F("hden","hden",20,0,10);
  hnum->Sumw2();
  hden->Sumw2();

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill((*dset)->get(i)->getRealValue("roosigma"),(*dset)->weight());
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(dsetdestination->get(i)->getRealValue("roosigma"),dsetdestination->weight());
  }

  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_sigmarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->weight();
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
  //  delete old_dset;

};

void validate_reweighting(RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  TH1F *test[4];
  TH1F *target[4];

  test[0] = new TH1F("test_rho","test_rho",30,0,30);
  test[1] = new TH1F("test_sigma","test_sigma",20,0,10);
  test[2] = new TH1F("test_pt","test_pt",55,25,300);
  test[3] = new TH1F("test_eta","test_eta",25,0,2.5);
  target[0] = new TH1F("target_rho","target_rho",30,0,30);
  target[1] = new TH1F("target_sigma","target_sigma",20,0,10);
  target[2] = new TH1F("target_pt","target_pt",55,25,300);
  target[3] = new TH1F("target_eta","target_eta",25,0,2.5);

  const char* ptname=Form("roopt%d",numvar);
  const char* etaname=Form("rooeta%d",numvar);

  for (int i=0; i<(*dset)->numEntries(); i++){
    test[0]->Fill((*dset)->get(i)->getRealValue("roorho"),(*dset)->weight());
    test[1]->Fill((*dset)->get(i)->getRealValue("roosigma"),(*dset)->weight());
    test[2]->Fill((*dset)->get(i)->getRealValue(ptname),(*dset)->weight());
    test[3]->Fill((*dset)->get(i)->getRealValue(etaname),(*dset)->weight());
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    target[0]->Fill(dsetdestination->get(i)->getRealValue("roorho"),dsetdestination->weight());
    target[1]->Fill(dsetdestination->get(i)->getRealValue("roosigma"),dsetdestination->weight());
    target[2]->Fill(dsetdestination->get(i)->getRealValue(ptname),dsetdestination->weight());
    target[3]->Fill(dsetdestination->get(i)->getRealValue(etaname),dsetdestination->weight());

  }

  for (int i=0; i<4; i++){
    test[i]->Scale(1.0/test[i]->Integral());
    target[i]->Scale(1.0/target[i]->Integral());
    test[i]->SetLineColor(kRed);
  }

  TCanvas *c = new TCanvas((*dset)->GetName(),(*dset)->GetName());
  c->Divide(2,2);

  for (int i=0;i<4; i++){
    c->cd(i+1);
    test[i]->Draw();
    target[i]->Draw("same");
  }

  c->SaveAs(Form("%s_rew.png",(*dset)->GetName()));

//  for (int i=0;i<4; i++){
//    delete test[i]; delete target[i];
//  }

};


    /*
    const float rho1 = (s1=="EB") ? 0.9 : 0.5;
    const float rho2 = (s2=="EB") ? 0.9 : 0.5;

    RooDataSet *sig1dset_nozero = (RooDataSet*)(sig1dset->reduce(Cut("roovar1>0.1")));
    firstbinpdf *firstbin_sig1 = new firstbinpdf("firstbin_sig1","firstbin_sig1",*roovar1);
    RooNDKeysPdf *kpdf_sig1 = new RooNDKeysPdf("kpdf_sig1","kpdf_sig1",*roovar1,*sig1dset_nozero,"av",rho1);
    kpdf_sig1->fixShape(1);
    sig1pdf = new RooAddPdf("sig1pdf","sig1pdf",RooArgList(*firstbin_sig1,*kpdf_sig1),*lm_sig1);

    RooDataSet *sig2dset_nozero = (RooDataSet*)(sig2dset->reduce(Cut("roovar2>0.1")));
    firstbinpdf *firstbin_sig2 = new firstbinpdf("firstbin_sig2","firstbin_sig2",*roovar2);
    //RooAbsPdf *firstbin_sig2 = RooClassFactory::makePdfInstance("firstbin_sig2","(roovar2<0.1)",RooArgSet(*roovar2));
    RooNDKeysPdf *kpdf_sig2 = new RooNDKeysPdf("kpdf_sig2","kpdf_sig2",*roovar2,*sig2dset_nozero,"av",rho2);
    kpdf_sig2->fixShape(1);
    sig2pdf = new RooAddPdf("sig2pdf","sig2pdf",RooArgList(*firstbin_sig2,*kpdf_sig2),*lm_sig2);

    RooDataSet *bkg1dset_nozero = (RooDataSet*)(bkg1dset->reduce(Cut("roovar1>0.1")));
    firstbinpdf *firstbin_bkg1 = new firstbinpdf("firstbin_bkg1","firstbin_bkg1",*roovar1);
    //RooAbsPdf *firstbin_bkg1 = RooClassFactory::makePdfInstance("firstbin_bkg1","(roovar1<0.1)",RooArgSet(*roovar1));
    RooNDKeysPdf *kpdf_bkg1 = new RooNDKeysPdf("kpdf_bkg1","kpdf_bkg1",*roovar1,*bkg1dset_nozero,"av",rho1);
    kpdf_bkg1->fixShape(1);
    bkg1pdf = new RooAddPdf("bkg1pdf","bkg1pdf",RooArgList(*firstbin_bkg1,*kpdf_bkg1),*lm_bkg1);

    RooDataSet *bkg2dset_nozero = (RooDataSet*)(bkg2dset->reduce(Cut("roovar2>0.1")));
    firstbinpdf *firstbin_bkg2 = new firstbinpdf("firstbin_bkg2","firstbin_bkg2",*roovar2);
    //RooAbsPdf *firstbin_bkg2 = RooClassFactory::makePdfInstance("firstbin_bkg2","(roovar2<0.1)",RooArgSet(*roovar2));
    RooNDKeysPdf *kpdf_bkg2 = new RooNDKeysPdf("kpdf_bkg2","kpdf_bkg2",*roovar2,*bkg2dset_nozero,"av",rho2);
    kpdf_bkg2->fixShape(1);
    bkg2pdf = new RooAddPdf("bkg2pdf","bkg2pdf",RooArgList(*firstbin_bkg2,*kpdf_bkg2),*lm_bkg2);

    
    
    sig1pdf->fitTo(*sig1dset,NumCPU(numcpu));
    bkg1pdf->fitTo(*bkg1dset,NumCPU(numcpu));
    if (sym) {
      lm_sig2->setVal(lm_sig1->getVal());
      lm_bkg2->setVal(lm_bkg1->getVal());
    }
    else {
      sig2pdf->fitTo(*sig2dset);
      bkg2pdf->fitTo(*bkg2dset);        
    }
    ok_for_recycle=1;
    lm_sig1_recycle=lm_sig1->getVal();
    lm_sig2_recycle=lm_sig2->getVal();
    lm_bkg1_recycle=lm_bkg1->getVal();
    lm_bkg2_recycle=lm_bkg2->getVal();

    sig1pdf_recycle=sig1pdf;
    sig2pdf_recycle=sig2pdf;
    bkg1pdf_recycle=bkg1pdf;
    bkg2pdf_recycle=bkg2pdf;
    */
