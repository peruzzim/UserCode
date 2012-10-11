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

using namespace std;
using namespace RooFit;

typedef struct {
  RooFitResult *fr;
  float tot_events;
} fit_output; 

bool study_templates=0;

fit_output* fit_dataset(const char* inputfilename_t, const char* inputfilename_d, TString diffvariable, TString splitting, int bin){

  bool dosingle=0;

  if (diffvariable==TString("dosingle")) {dosingle=1; diffvariable=TString("singlegamma_eta");}

  fit_output *out = new fit_output();
  out->fr=NULL;
  out->tot_events=0;

  TFile *inputfile_t = TFile::Open(inputfilename_t);
  TFile *inputfile_d = TFile::Open(inputfilename_d);


  //TString data_dir("mc_Tree_standard_sel/");
  TString data_dir("data_Tree_standard_sel/"); 
  
  TString sig_dir("data_Tree_randomcone_signal_template/");
  //  TString sig_dir("mc_Tree_signal_template/");
  //  TString sig_dir("mc_Tree_randomcone_signal_template/");

  TString bkg_dir("data_Tree_sieiesideband_sel/");
  //TString bkg_dir("mc_Tree_sieiesideband_sel/");

  if (splitting=="EEEB") splitting="EBEE";

  TH1F *h_datahist_1D;
  TH2F *h_datahist_2D;

  TH1F *h_sig1hist;
  TH1F *h_sig2hist;
  TH1F *h_bkg1hist;
  TH1F *h_bkg2hist;

  TH1::SetDefaultSumw2(kTRUE);
  
    TString s1; TString s2;
    if (splitting=="EBEB") {s1="EB"; s2="EB";}
    else if (splitting=="EEEE") {s1="EE"; s2="EE";}
    else if (splitting=="EBEE") {s1="EB"; s2="EE";}

    if (dosingle){
      inputfile_t->GetObject(TString(sig_dir).Append(Form("template_signal_%s_b%d",splitting.Data(),n_bins)),h_sig1hist);
      inputfile_t->GetObject(TString(bkg_dir).Append(Form("template_background_%s_b%d",splitting.Data(),n_bins)),h_bkg1hist);
      inputfile_d->GetObject(TString(data_dir).Append(Form("obs_hist_single_%s_b%d",splitting.Data(),bin)),h_datahist_1D);
      std::cout << inputfilename_d << std::endl;
      std::cout << TString(data_dir).Append(Form("obs_hist_single_%s_b%d",splitting.Data(),bin)) << std::endl;
      assert (h_sig1hist!=NULL);
      assert (h_bkg1hist!=NULL);
      assert (h_datahist_1D!=NULL);
    }
    else {
      inputfile_t->GetObject(TString(sig_dir).Append(Form("template_signal_%s_b%d",s1.Data(),n_bins)),h_sig1hist);
      inputfile_t->GetObject(TString(bkg_dir).Append(Form("template_background_%s_b%d",s1.Data(),n_bins)),h_bkg1hist);
      inputfile_t->GetObject(TString(sig_dir).Append(Form("template_signal_%s_b%d",s2.Data(),n_bins)),h_sig2hist);
      inputfile_t->GetObject(TString(bkg_dir).Append(Form("template_background_%s_b%d",s2.Data(),n_bins)),h_bkg2hist);
      inputfile_d->GetObject(TString(data_dir).Append(Form("obs_hist_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin)),h_datahist_2D);
      std::cout << inputfilename_d << std::endl;
      std::cout << TString(data_dir).Append(Form("obs_hist_%s_%s_b%d",splitting.Data(),diffvariable.Data(),bin)) << std::endl;
      assert (h_sig1hist!=NULL);
      assert (h_bkg1hist!=NULL);
      assert (h_sig2hist!=NULL);
      assert (h_bkg2hist!=NULL);
      assert (h_datahist_2D!=NULL);
    }


    
  
  RooRealVar roovar1("roovar1","roovar1",0,5);
  roovar1.setBins(n_histobins);
  RooRealVar roovar2("roovar2","roovar2",0,5);
  roovar2.setBins(n_histobins);

  RooRealVar rf1("rf1","rf1",1./3,0,1);
  RooRealVar rf2("rf2","rf2",1./2,0,1);

  RooDataHist *sig1hist;
  RooDataHist *sig2hist;
  RooDataHist *bkg1hist;
  RooDataHist *bkg2hist;
  RooDataHist *datahist;

  RooHistPdf *sig1pdf;
  RooHistPdf *sig2pdf;
  RooHistPdf *bkg1pdf;
  RooHistPdf *bkg2pdf;

  RooProdPdf *sigsigpdf;
  RooProdPdf *sigbkgpdf_order;
  RooProdPdf *bkgsigpdf_order;
  RooAddPdf *sigbkgpdf_noorder;
  RooProdPdf *bkgbkgpdf;

  RooFormulaVar *fsigsig;
  RooFormulaVar *fsigbkg;
  RooFormulaVar *fbkgbkg;
  RooFormulaVar *fsig;
  RooFormulaVar *fbkg;
  
  RooAddPdf *model;

  if (dosingle){
    sig1hist = new RooDataHist("sighist","sighist",RooArgList(roovar1),h_sig1hist,1.0/h_sig1hist->Integral());
    bkg1hist = new RooDataHist("bkghist","bkghist",RooArgList(roovar1),h_bkg1hist,1.0/h_bkg1hist->Integral());
    datahist = new RooDataHist("datahist","datahist",roovar1,h_datahist_1D);

    sig1pdf = new RooHistPdf("sig1pdf","sig1pdf",RooArgList(roovar1),*sig1hist);
    bkg1pdf = new RooHistPdf("bkg1pdf","bkg1pdf",RooArgList(roovar1),*bkg1hist);
    fsig = new RooFormulaVar("fsig","fsig","rf1",RooArgList(rf1));
    fbkg = new RooFormulaVar("fbkg","fbkg","1-fsigsig",RooArgList(*fsig));
    model = new RooAddPdf("model","model",RooArgList(*sig1pdf,*bkg1pdf),RooArgList(rf1),kTRUE);
  }
  else {
    sig1hist = new RooDataHist("sighist","sighist",RooArgList(roovar1),h_sig1hist,1.0/h_sig1hist->Integral());
    sig2hist = new RooDataHist("sighist","sighist",RooArgList(roovar2),h_sig2hist,1.0/h_sig2hist->Integral());
    bkg1hist = new RooDataHist("bkghist","bkghist",RooArgList(roovar1),h_bkg1hist,1.0/h_bkg1hist->Integral());
    bkg2hist = new RooDataHist("bkghist","bkghist",RooArgList(roovar2),h_bkg2hist,1.0/h_bkg2hist->Integral());
    datahist = new RooDataHist("datahist","datahist",RooArgList(roovar1,roovar2),h_datahist_2D);

    sig1pdf = new RooHistPdf("sig1pdf","sig1pdf",RooArgList(roovar1),*sig1hist);
    sig2pdf = new RooHistPdf("sig2pdf","sig2pdf",RooArgList(roovar2),*sig2hist);
    bkg1pdf = new RooHistPdf("bkg1pdf","bkg1pdf",RooArgList(roovar1),*bkg1hist);
    bkg2pdf = new RooHistPdf("bkg2pdf","bkg2pdf",RooArgList(roovar2),*bkg2hist);

    sigsigpdf = new RooProdPdf("sigsigpdf","sigsigpdf",RooArgList(*sig1pdf,*sig2pdf));
    sigbkgpdf_order = new RooProdPdf("sigbkgpdf_order","sigbkgpdf_order",RooArgList(*sig1pdf,*bkg2pdf));
    bkgsigpdf_order = new RooProdPdf("bkgsigpdf_order","bkgsigpdf_order",RooArgList(*bkg1pdf,*sig2pdf));
    sigbkgpdf_noorder = new RooAddPdf("sigbkgpdf","sigbkgpdf",RooArgList(*sigbkgpdf_order,*bkgsigpdf_order),RooArgList(RooRealConstant::value(0.5),RooRealConstant::value(0.5)));
    bkgbkgpdf = new RooProdPdf("bkgbkgpdf","bkgbkgpdf",RooArgList(*bkg1pdf,*bkg2pdf));

    fsigsig = new RooFormulaVar("fsigsig","fsigsig","rf1",RooArgList(rf1));
    fsigbkg = new RooFormulaVar("fsigbkg","fsigbkg","(1-rf1)*rf2",RooArgList(rf1,rf2));
    fbkgbkg = new RooFormulaVar("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg",RooArgList(*fsigsig,*fsigbkg));

    model = new RooAddPdf("model","model",RooArgList(*sigsigpdf,*sigbkgpdf_noorder,*bkgbkgpdf),RooArgList(rf1,rf2),kTRUE);
  }

    RooMCStudy *mcstudy = NULL;

  if (study_templates){

    rf1.setVal(0.2);
    rf2.setVal(1./2);

    if (dosingle) mcstudy = new RooMCStudy(*model,RooArgSet(roovar1),Silence(),FitOptions(Save(kTRUE),PrintEvalErrors(0)));
    else mcstudy = new RooMCStudy(*model,RooArgSet(roovar1,roovar2),Silence(),FitOptions(Save(kTRUE),PrintEvalErrors(0)));
    mcstudy->generateAndFit(500, 2e+3);

    TCanvas *c1 = new TCanvas();
    c1->SetWindowSize(800,600);
    c1->Divide(3);

    c1->cd(1);    
    RooPlot* frame3 = mcstudy->plotPull(rf1,FitGauss(kTRUE));
    frame3->Draw();

    TVirtualPad* subpad = c1->cd(3);
    subpad->Divide(2,2);

    subpad->cd(1);    
    RooPlot* frame1 = mcstudy->plotParam(rf1);
    frame1->Draw();
    subpad->cd(2);    
    RooPlot* frame2 = mcstudy->plotError(rf1);
    frame2->Draw();

    if (!dosingle){
      c1->cd(2);    
      RooPlot* frame6 = mcstudy->plotPull(rf2,FitGauss(kTRUE));
      frame6->Draw();
      subpad->cd(3);    
      RooPlot* frame4 = mcstudy->plotParam(rf2);
      frame4->Draw();
      subpad->cd(4);    
      RooPlot* frame5 = mcstudy->plotError(rf2);
      frame5->Draw();
    }

    c1->cd();
    c1->SaveAs(Form("plots/biasstudy_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));

    /*
    RooDataHist *datahist_toy;
    datahist_toy = model->generateBinned(RooArgSet(roovar1,roovar2),1e4,Name("Toy_dataset"));
    datahist=datahist_toy;
    rf1.setVal(0);
    rf2.setVal(0);
    */

    return out;

  }


//  roovar1.setRange("fullrange",0,5);
//  roovar2.setRange("fullrange",0,5);
//  roovar1.setRange("fitrange",2,4);
//  roovar2.setRange("fitrange",2,4);

//  RooAbsReal* nll = model->createNLL(*datahist,SumCoefRange("fitrange"),Range("fitrange"));
//  RooMinuit m(*nll) ;
//  m.setVerbose(kTRUE) ;
//  m.migrad() ;
//  m.setVerbose(kFALSE) ;
//  m.hesse() ;
//  RooFitResult *fitres = m.save();

  RooFitResult *fitres = model->fitTo(*datahist,Save());

  model->Print();

  std::cout << "---------------------------" << std::endl;
  std::cout << rf1.getVal() << " " << rf2.getVal() << std::endl;

  out->fr=fitres;

  if (dosingle) out->tot_events=h_datahist_1D->Integral();
  else out->tot_events=h_datahist_2D->Integral();

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
    model->plotOn(varframe[i]);
    if (dosingle){
      model->plotOn(varframe[i],Components("sigpdf"),LineStyle(kDashed),LineColor(kRed));
      model->plotOn(varframe[i],Components("bkgpdf"),LineStyle(kDashed),LineColor(kBlack));
    }
    else {
      model->plotOn(varframe[i],Components("sigsigpdf"),LineStyle(kDashed),LineColor(kRed));
      model->plotOn(varframe[i],Components("sigbkgpdf"),LineStyle(kDashed),LineColor(kGreen));
      model->plotOn(varframe[i],Components("bkgbkgpdf"),LineStyle(kDashed),LineColor(kBlack));
    }
    varframe[i]->Draw();
    canv->GetPad(i+1)->SetLogy(1);
  }

  canv->SaveAs(Form("plots/fittingplot_%s_%s_b%d.png",splitting.Data(),diffvariable.Data(),bin));

  return out;

};




void run_fits(TString inputfilename_t="templates.root", TString inputfilename_d="tobefitted.root", TString diffvariable="", TString splitting=""){

  TH1F::SetDefaultSumw2(kTRUE);

  bool dosingle=0;

  int bins_to_run=0; 
  float *binsdef;

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
  if (diffvariable=="dosingle"){
    dosingle=1;
    if (splitting=="EB")      bins_to_run=n_templates_EB;
    else if (splitting=="EE") bins_to_run=n_templates_EE;
    if (splitting=="EB")      binsdef=binsdef_single_gamma_EB_eta;
    else if (splitting=="EE") binsdef=binsdef_single_gamma_EE_eta;
    diffvariable="singlegamma_eta";
  }

  
  fit_output *fr[n_bins];

  TH1F *purity[3];
  TH1F *eff;
  TH1F *xsec;

  TFile *eff_file = new TFile("efficiencies.root");
  if (dosingle) eff_file->GetObject(Form("w_eff_sg_%s",splitting.Data()),eff);
  eff_file->GetObject(Form("w_eff_gg_%s_%s",splitting.Data(),diffvariable.Data()),eff);

  assert (eff!=NULL);

  eff->GetYaxis()->SetTitle("selection/ID efficiency");
  eff->GetXaxis()->SetTitle(diffvariable.Data());

  int colors[3] = {kRed, kGreen, kBlack};

  for (int i=0; i<3; i++){
    TString name = "purity_";
    if (i==0) name.Append("sigsig"); else if (i==1) name.Append("sigbkg"); else if (i==2) name.Append("bkgbkg");
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




  RooRealVar *rf1;
  RooRealVar *rf2;

  if (study_templates) bins_to_run=1;

  for (int bin=0; bin<bins_to_run; bin++) {

    fr[bin]=fit_dataset(inputfilename_t.Data(),inputfilename_d.Data(),diffvariable,splitting,bin);
    
    if (!fr[bin]->fr) continue;

    float intlumi=4.519;

    if (dosingle) {

	rf1=(RooRealVar*)(fr[bin]->fr->floatParsFinal().find("rf1"));
	RooFormulaVar fsig("fsig","fsig","rf1",RooArgList(*rf1));
	RooFormulaVar fbkg("fbkg","fbkg","1-fsig",RooArgList(fsig));
	std::cout << "------------------------------" << std::endl;
	std::cout << bin << " " << fsig.getVal() << " " << fbkg.getVal() << std::endl;

	purity[0]->SetBinContent(bin+1,fsig.getVal());
	purity[0]->SetBinError(bin+1,fsig.getPropagatedError(*(fr[bin]->fr)));
	purity[2]->SetBinContent(bin+1,fbkg.getVal());
	purity[2]->SetBinError(bin+1,fbkg.getPropagatedError(*(fr[bin]->fr)));

	xsec->SetBinContent(bin+1,fsig.getVal()*fr[bin]->tot_events/xsec->GetBinWidth(bin+1)/intlumi);
	float err1=purity[0]->GetBinError(bin+1)/purity[0]->GetBinContent(bin+1);
	float err2=1.0/sqrt(fr[bin]->tot_events);
	float err=sqrt(err1*err1+err2*err2);
	xsec->SetBinError(bin+1,err*xsec->GetBinContent(bin+1));
	
	xsec->Divide(eff);

    }

    else {
	rf1=(RooRealVar*)(fr[bin]->fr->floatParsFinal().find("rf1"));
	rf2=(RooRealVar*)(fr[bin]->fr->floatParsFinal().find("rf2"));
	RooFormulaVar fsigsig("fsigsig","fsigsig","rf1",RooArgList(*rf1));
	RooFormulaVar fsigbkg("fsigbkg","fsigbkg","(1-rf1)*rf2",RooArgList(*rf1,*rf2));
	RooFormulaVar fbkgbkg("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg",RooArgList(fsigsig,fsigbkg));
	std::cout << "------------------------------" << std::endl;
	std::cout << bin << " " << fsigsig.getVal() << " " << fsigbkg.getVal() << " " << fbkgbkg.getVal() << std::endl;

	purity[0]->SetBinContent(bin+1,fsigsig.getVal());
	purity[0]->SetBinError(bin+1,fsigsig.getPropagatedError(*(fr[bin]->fr)));
	purity[1]->SetBinContent(bin+1,fsigbkg.getVal());
	purity[1]->SetBinError(bin+1,fsigbkg.getPropagatedError(*(fr[bin]->fr)));
	purity[2]->SetBinContent(bin+1,fbkgbkg.getVal());
	purity[2]->SetBinError(bin+1,fbkgbkg.getPropagatedError(*(fr[bin]->fr)));

	xsec->SetBinContent(bin+1,fsigsig.getVal()*fr[bin]->tot_events/xsec->GetBinWidth(bin+1)/intlumi);
	float err1=purity[0]->GetBinError(bin+1)/purity[0]->GetBinContent(bin+1);
	float err2=1.0/sqrt(fr[bin]->tot_events);
	float err=sqrt(err1*err1+err2*err2);
	xsec->SetBinError(bin+1,err*xsec->GetBinContent(bin+1));
	
	xsec->Divide(eff);
    }


  }

  if (study_templates) return;


  TCanvas *output_canv = new TCanvas("output_canv","output_canv");
  output_canv->cd();

  purity[0]->Draw("e1");
  if (!dosingle) purity[1]->Draw("e1same");
  purity[2]->Draw("e1same");

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
  if (!dosingle) purity[1]->Write();
  purity[2]->Write();

  TFile *outfile = new TFile(Form("plots/xsec_%s_%s.root",splitting.Data(),diffvariable.Data()),"recreate");
  outfile->cd();
  xsec->Write();
  

};


