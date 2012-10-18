#include <assert.h>

#include "binsdef.h"
#include "RooFitResult.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
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

const int numcpu=12;

//typedef struct {
//  RooAbsPdf *pdf;
//  RooAbsData *data;
//} necessary_for_fit;
//
//void *threaded_fitter(void *arg){
//  necessary_for_fit *input = (necessary_for_fit*)arg;
//  return input->pdf->fitTo(*(input->data),Save());
//};

fit_output* fit_dataset(const char* inputfilename_ts, const char* inputfilename_tb, const char* inputfilename_d, TString diffvariable, TString splitting, int bin){

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
  TString data_dir("data_Tree_standard_sel/"); 
  
  TString sig_dir("data_Tree_randomcone_signal_template/");
  //  TString sig_dir("mc_Tree_signal_template/");
  //  TString sig_dir("mc_Tree_randomcone_signal_template/");

  TString bkg_dir("data_Tree_sieiesideband_sel/");
  //TString bkg_dir("mc_Tree_sieiesideband_sel/");

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
    assert(wspace_ts); wspace_ts->Print();
    assert(wspace_tb); wspace_tb->Print();
    assert(wspace_d); wspace_d->Print();  

    RooRealVar *roovar1 = wspace_d->var("roovar1");
    RooRealVar *roovar2 = wspace_d->var("roovar2");
    //    RooRealVar *roovar_helper = wspace_ts->var("roovar_helper");
    assert (roovar1); roovar1->setBins(200); roovar1->Print();    
    assert (roovar2); roovar2->setBins(200); roovar2->Print();    

    const float pp_init = 0.4;
    const float pf_init = 0.1;
    const float fp_init = 0.1;
    
      // inversione:
      //  pp = pp;
      //  pf = j1 - pp;
      //  fp = j2 - pp;
      //  ff = (1-j1-j2)+pp
      
    RooRealVar *pp = new RooRealVar("pp","pp",pp_init,0,1);
    RooRealVar *j1 = new RooRealVar("j1","j1",pp_init+pf_init,0,1);
    RooRealVar *j2 = new RooRealVar("j2","j2",pp_init+fp_init,0,1);
    

    RooDataSet *sig1dset = (RooDataSet*)(wspace_ts->data(Form("roodset_signal_%s_b%d_rv%d",s1.Data(),n_bins,1)));
    RooDataSet *bkg1dset = (RooDataSet*)(wspace_tb->data(Form("roodset_background_%s_b%d_rv%d",s1.Data(),n_bins,1)));
    RooDataSet *sig2dset = (RooDataSet*)(wspace_ts->data(Form("roodset_signal_%s_b%d_rv%d",s2.Data(),n_bins,2)));
    RooDataSet *bkg2dset = (RooDataSet*)(wspace_tb->data(Form("roodset_background_%s_b%d_rv%d",s2.Data(),n_bins,2)));
    RooDataSet *dataset = (RooDataSet*)(wspace_d->data(Form("obs_roodset_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin)));    

    sig1dset->Print();
    bkg1dset->Print();
    sig2dset->Print();
    bkg2dset->Print();
    dataset->Print(); 


    /*
    const int nev = 1e4; // DEBUUUUUUUUUUUUUUUUUUUUUUUUG!!!!!!!!!!!!!!!!
    RooDataSet *sig1dset_red=(RooDataSet*)(sig1dset->reduce(EventRange(0,nev))); // DEBUUUUUUUUUUUUUUUUUUUUUUUUG!!!!!!!!!!!!!!!!
    RooDataSet *bkg1dset_red=(RooDataSet*)(bkg1dset->reduce(EventRange(0,nev))); // DEBUUUUUUUUUUUUUUUUUUUUUUUUG!!!!!!!!!!!!!!!!
    RooDataSet *sig2dset_red=(RooDataSet*)(sig2dset->reduce(EventRange(0,nev))); // DEBUUUUUUUUUUUUUUUUUUUUUUUUG!!!!!!!!!!!!!!!!
    RooDataSet *bkg2dset_red=(RooDataSet*)(bkg2dset->reduce(EventRange(0,nev))); // DEBUUUUUUUUUUUUUUUUUUUUUUUUG!!!!!!!!!!!!!!!!
    RooDataSet *dataset_red=(RooDataSet*)(dataset->reduce(EventRange(0,nev))); // DEBUUUUUUUUUUUUUUUUUUUUUUUUG!!!!!!!!!!!!!!!!
    sig1dset=sig1dset_red;
    bkg1dset=bkg1dset_red;
    sig2dset=sig2dset_red;
    bkg2dset=bkg2dset_red;
    dataset=dataset_red;
    */
    sig1dset->Print();
    bkg1dset->Print();
    sig2dset->Print();
    bkg2dset->Print();
    dataset->Print(); 
    

    RooDataSet *dataset_axis1 = (RooDataSet*)(dataset->reduce(SelectVars(*roovar1)));
    RooDataSet *dataset_axis2 = (RooDataSet*)(dataset->reduce(SelectVars(*roovar2)));


    RooDataSet *sig1dset_nozero = (RooDataSet*)(sig1dset->reduce(Cut("roovar1>0.1")));
    firstbinpdf *firstbin_sig1 = new firstbinpdf("firstbin_sig1","firstbin_sig1",*roovar1);
    RooNDKeysPdf *kpdf_sig1 = new RooNDKeysPdf("kpdf_sig1","kpdf_sig1",*roovar1,*sig1dset_nozero,"av",1.0);
    kpdf_sig1->fixShape(1);
    RooRealVar *lm_sig1 = new RooRealVar("lm_sig1","lm_sig1",0.5,0,1);
    RooAddPdf *sig1pdf = new RooAddPdf("sig1pdf","sig1pdf",RooArgList(*firstbin_sig1,*kpdf_sig1),*lm_sig1);

    RooDataSet *sig2dset_nozero = (RooDataSet*)(sig2dset->reduce(Cut("roovar2>0.1")));
    firstbinpdf *firstbin_sig2 = new firstbinpdf("firstbin_sig2","firstbin_sig2",*roovar2);
    //RooAbsPdf *firstbin_sig2 = RooClassFactory::makePdfInstance("firstbin_sig2","(roovar2<0.1)",RooArgSet(*roovar2));
    RooNDKeysPdf *kpdf_sig2 = new RooNDKeysPdf("kpdf_sig2","kpdf_sig2",*roovar2,*sig2dset_nozero,"av",1.0);
    kpdf_sig2->fixShape(1);
    RooRealVar *lm_sig2 = new RooRealVar("lm_sig2","lm_sig2",0.5,0,1);
    RooAddPdf *sig2pdf = new RooAddPdf("sig2pdf","sig2pdf",RooArgList(*firstbin_sig2,*kpdf_sig2),*lm_sig2);

    RooDataSet *bkg1dset_nozero = (RooDataSet*)(bkg1dset->reduce(Cut("roovar1>0.1")));
    firstbinpdf *firstbin_bkg1 = new firstbinpdf("firstbin_bkg1","firstbin_bkg1",*roovar1);
    //RooAbsPdf *firstbin_bkg1 = RooClassFactory::makePdfInstance("firstbin_bkg1","(roovar1<0.1)",RooArgSet(*roovar1));
    RooNDKeysPdf *kpdf_bkg1 = new RooNDKeysPdf("kpdf_bkg1","kpdf_bkg1",*roovar1,*bkg1dset_nozero,"av",1.0);
    kpdf_bkg1->fixShape(1);
    RooRealVar *lm_bkg1 = new RooRealVar("lm_bkg1","lm_bkg1",0.5,0,1);
    RooAddPdf *bkg1pdf = new RooAddPdf("bkg1pdf","bkg1pdf",RooArgList(*firstbin_bkg1,*kpdf_bkg1),*lm_bkg1);

    RooDataSet *bkg2dset_nozero = (RooDataSet*)(bkg2dset->reduce(Cut("roovar2>0.1")));
    firstbinpdf *firstbin_bkg2 = new firstbinpdf("firstbin_bkg2","firstbin_bkg2",*roovar2);
    //RooAbsPdf *firstbin_bkg2 = RooClassFactory::makePdfInstance("firstbin_bkg2","(roovar2<0.1)",RooArgSet(*roovar2));
    RooNDKeysPdf *kpdf_bkg2 = new RooNDKeysPdf("kpdf_bkg2","kpdf_bkg2",*roovar2,*bkg2dset_nozero,"av",1.0);
    kpdf_bkg2->fixShape(1);
    RooRealVar *lm_bkg2 = new RooRealVar("lm_bkg2","lm_bkg2",0.5,0,1);
    RooAddPdf *bkg2pdf = new RooAddPdf("bkg2pdf","bkg2pdf",RooArgList(*firstbin_bkg2,*kpdf_bkg2),*lm_bkg2);




    


//    necessary_for_fit *input_thrsig1 = new necessary_for_fit(); input_thrsig1->pdf=sig1pdf; input_thrsig1->data=sig1dset;
//    TThread *thrsig1 = new TThread("thrsig1", threaded_fitter, (void*)input_thrsig1);
//    necessary_for_fit *input_thrsig2 = new necessary_for_fit(); input_thrsig2->pdf=sig2pdf; input_thrsig2->data=sig2dset;
//    TThread *thrsig2 = new TThread("thrsig2", threaded_fitter, (void*)input_thrsig2);
//    necessary_for_fit *input_thrbkg1 = new necessary_for_fit(); input_thrbkg1->pdf=bkg1pdf; input_thrbkg1->data=bkg1dset;
//    TThread *thrbkg1 = new TThread("thrbkg1", threaded_fitter, (void*)input_thrbkg1);
//    necessary_for_fit *input_thrbkg2 = new necessary_for_fit(); input_thrbkg2->pdf=bkg2pdf; input_thrbkg2->data=bkg2dset;
//    TThread *thrbkg2 = new TThread("thrbkg2", threaded_fitter, (void*)input_thrbkg2);
//
//    thrsig1->Run();
//    thrsig2->Run();
//    thrbkg1->Run();
//    thrbkg2->Run();
//
//    TThread::Ps();
//
//    thrsig1->Join();
//    thrsig2->Join();
//
//    TThread::Ps();
//
//    thrbkg1->Join();
//    thrbkg2->Join();
//
//    TThread::Ps();



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

    lm_sig1->setConstant(kTRUE);
    lm_sig2->setConstant(kTRUE);
    lm_bkg1->setConstant(kTRUE);
    lm_bkg2->setConstant(kTRUE);

    
    TCanvas *c0 = new TCanvas("c0","c0",800,600);
    c0->Divide(2,2);
    c0->cd(1);
    RooPlot *frame01sig = roovar1->frame();
    sig1dset->plotOn(frame01sig);
    sig1pdf->plotOn(frame01sig);
    frame01sig->Draw();
    c0->GetPad(1)->SetLogy(1);
    c0->cd(2);
    RooPlot *frame02sig = roovar2->frame();
    sig2dset->plotOn(frame02sig);
    sig2pdf->plotOn(frame02sig);
    frame02sig->Draw();
    c0->GetPad(2)->SetLogy(1);
    c0->cd(3);
    RooPlot *frame01bkg = roovar1->frame();
    bkg1dset->plotOn(frame01bkg);
    bkg1pdf->plotOn(frame01bkg);
    frame01bkg->Draw();
    c0->GetPad(3)->SetLogy(1);
    c0->cd(4);
    RooPlot *frame02bkg = roovar2->frame();
    bkg2dset->plotOn(frame02bkg);
    bkg2pdf->plotOn(frame02bkg);
    frame02bkg->Draw();
    c0->GetPad(4)->SetLogy(1);

    
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

    
    //    RooRealVar *nevents = new RooRealVar("nevents","nevents",dataset->sumEntries(),0,1e6);

    RooAddPdf *model_2D = new RooAddPdf("model_2D","model_2D",RooArgList(*sigsigpdf,*sigbkgpdf,*bkgsigpdf,*bkgbkgpdf),RooArgList(*fsigsig,*fsigbkg,*fbkgsig,*fbkgbkg),kFALSE);
    //    RooExtendPdf *model_2D_extended = new RooExtendPdf("model_2D_extended","model_2D_extended",*model_2D,*nevents);
    //    RooNLLVar *model_2D_extended_nll = new RooNLLVar("model_2D_extended_nll","model_2D_extended_nll",*model_2D_extended,*dataset,NumCPU(numcpu));
    RooNLLVar *model_2D_noextended_nll = new RooNLLVar("model_2D_noextended_nll","model_2D_noextended_nll",*model_2D,*dataset,NumCPU(numcpu));
    RooNLLVar *model_2D_nll = model_2D_noextended_nll;

    RooAddPdf *model_axis1 = new RooAddPdf("model_axis1","model_axis1",RooArgList(*sig1pdf,*bkg1pdf),RooArgList(*fsig1));
    //    RooExtendPdf *model_axis1_extended = new RooExtendPdf("model_axis1_extended","model_axis1_extended",*model_axis1,*nevents);
    //    RooNLLVar *model_axis1_extended_nll = new RooNLLVar("model_axis1_extended_nll","model_axis1_extended_nll",*model_axis1_extended,*dataset_axis1,NumCPU(numcpu));
    RooNLLVar *model_axis1_noextended_nll = new RooNLLVar("model_axis1_noextended_nll","model_axis1_noextended_nll",*model_axis1,*dataset_axis1,NumCPU(numcpu));
    RooNLLVar *model_axis1_nll = model_axis1_noextended_nll;

    RooAddPdf *model_axis2 = new RooAddPdf("model_axis2","model_axis2",RooArgList(*sig2pdf,*bkg2pdf),RooArgList(*fsig2));
    //    RooExtendPdf *model_axis2_extended = new RooExtendPdf("model_axis2_extended","model_axis2_extended",*model_axis2,*nevents);
    //    RooNLLVar *model_axis2_extended_nll = new RooNLLVar("model_axis2_extended_nll","model_axis2_extended_nll",*model_axis2_extended,*dataset_axis2,NumCPU(numcpu));
    RooNLLVar *model_axis2_noextended_nll = new RooNLLVar("model_axis2_noextended_nll","model_axis2_noextended_nll",*model_axis2,*dataset_axis2,NumCPU(numcpu));
    RooNLLVar *model_axis2_nll = model_axis2_noextended_nll;

    //    RooFormulaVar *asym = new RooFormulaVar("asym","asym","TMath::Abs(j1-j2)",RooArgSet(*j1,*j2));
    //    RooGaussian *constrain_asym = new RooGaussian("constrain_asym","constrain_asym",*asym,RooRealConstant::value(0),RooRealConstant::value(0.001));
    //    RooConstraintSum *constrain_asym_nll = new RooConstraintSum("constrain_asym_nll","constrain_asym_nll",*constrain_asym,RooArgList(*asym));

    RooAddition *model_2axes_nll = new RooAddition("model_2axes_nll","model_2axes_nll",RooArgSet(*model_axis1_nll,*model_axis2_nll));
    //    RooAddition *model_2D_nll_constrained = new RooAddition("model_2D_nll_constrained","model_2D_nll_constrained",RooArgSet(*model_2D_extended_nll,*constrain_asym_nll));


    

    RooFitResult *firstpass;


    RooMinimizer *minuit_firstpass = new RooMinimizer(*model_2axes_nll);
    minuit_firstpass->migrad();
    minuit_firstpass->hesse();
    firstpass = minuit_firstpass->save("firstpass","firstpass");
    firstpass->Print();


    /*
    j1->setVal(7.17932e-01);
    j2->setVal(7.17913e-01);
    */

    TCanvas *c1 = new TCanvas("c1","c1",800,600);
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






    RooMinimizer *minuit_secondpass = new RooMinimizer(*model_2D_nll);
    minuit_secondpass->migrad();
    minuit_secondpass->hesse();
    RooFitResult *secondpass = minuit_secondpass->save("secondpass","secondpass");
    secondpass->Print();

    TCanvas *c2 = new TCanvas("c2","c2",800,600);
    c2->Divide(2,2);
    c2->cd(1);
    RooPlot *frame1final = roovar1->frame();
    dataset->plotOn(frame1final);
    model_2D->plotOn(frame1final);
    model_2D->plotOn(frame1final,Components("sigsigpdf"),LineStyle(kDashed),LineColor(kRed));
    model_2D->plotOn(frame1final,Components("sigbkgpdf"),LineStyle(kDashed),LineColor(kGreen));
    model_2D->plotOn(frame1final,Components("bkgsigpdf"),LineStyle(kDashed),LineColor(kGreen+2));
    model_2D->plotOn(frame1final,Components("bkgbkgpdf"),LineStyle(kDashed),LineColor(kBlack));
    frame1final->Draw();
    c2->GetPad(1)->SetLogy(1);
    c2->cd(2);
    RooPlot *frame2final = roovar2->frame();
    dataset->plotOn(frame2final);
    model_2D->plotOn(frame2final);
    model_2D->plotOn(frame2final,Components("sigsigpdf"),LineStyle(kDashed),LineColor(kRed));
    model_2D->plotOn(frame2final,Components("sigbkgpdf"),LineStyle(kDashed),LineColor(kGreen));
    model_2D->plotOn(frame2final,Components("bkgsigpdf"),LineStyle(kDashed),LineColor(kGreen+2));
    model_2D->plotOn(frame2final,Components("bkgbkgpdf"),LineStyle(kDashed),LineColor(kBlack));
    frame2final->Draw();
    c2->GetPad(2)->SetLogy(1);
//    c2->cd(3);
//    RooPlot *ppnllplot = pp->frame();
//    model_2D_extended_nll_constrained->plotOn(ppnllplot);
//    ppnllplot->Draw();



//    std::cout << "expecting purities = " << pp_init << " " << pf_init << " " << fp_init << " " << 1-pp_init-pf_init-fp_init << std::endl;
//    std::cout << "pp " << fsigsig->getVal() << " " << fsigsig->getPropagatedError(*fr) << " " << fabs(pp_init-fsigsig->getVal())/fsigsig->getPropagatedError(*fr) << std::endl;
//    std::cout << "pf " << fsigbkg->getVal() << " " << fsigbkg->getPropagatedError(*fr) << " " << fabs(pf_init-fsigbkg->getVal())/fsigbkg->getPropagatedError(*fr) << std::endl;
//    std::cout << "fp " << fbkgsig->getVal() << " " << fbkgsig->getPropagatedError(*fr) << " " << fabs(fp_init-fbkgsig->getVal())/fbkgsig->getPropagatedError(*fr) << std::endl;
//    std::cout << "ff " << fbkgbkg->getVal() << " " << fbkgbkg->getPropagatedError(*fr) << " " << fabs(1-pp_init-pf_init-fp_init-fbkgbkg->getVal())/fbkgbkg->getPropagatedError(*fr) << std::endl;


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




    c0->SaveAs(Form("plots/fittingplot0_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
    c1->SaveAs(Form("plots/fittingplot1_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));
    c2->SaveAs(Form("plots/fittingplot2_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));


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

//  TFile *eff_file = new TFile("efficiencies.root");
//  eff_file->GetObject(Form("w_eff_gg_%s_%s",splitting.Data(),diffvariable.Data()),eff);
//  assert (eff!=NULL);
//  eff->GetYaxis()->SetTitle("selection/ID efficiency");
//  eff->GetXaxis()->SetTitle(diffvariable.Data());

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
    
    for (int k=0; k<100; k++) std::cout << "WARNING: NO EFFICIENCY CORRECTION!!!" << std::endl;
    //    xsec->Divide(eff);

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



