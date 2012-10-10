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

using namespace std;
using namespace RooFit;

typedef struct {
  RooFitResult *fr;
  float tot_events;
} fit_output; 

bool study_templates=0;

fit_output* fit_dataset(const char* inputfilename_t, const char* inputfilename_d, TString diffvariable, TString splitting, float leftrange, float rightrange, int bin, TCanvas *canv=NULL, int rbin=1){

fit_output *out = new fit_output();

  TFile *inputfile_t = TFile::Open(inputfilename_t);
  TFile *inputfile_d = TFile::Open(inputfilename_d);

  //  TString data_dir("data_Tree_standard_sel/");
  //  if (!isdata) data_dir="mc_Tree_standard_sel/";
  TString data_dir("data_Tree_standard_sel/"); //testing without mc
  
  TString sig_dir("data_Tree_randomcone_signal_template/");
  //  TString bkg_dir("data_Tree_impinging_track_template/");
  //  TString sig_dir("mc_Tree_randomcone_signal_template/");
  //  TString bkg_dir("mc_Tree_impinging_track_template/");
  //  TString sig_dir("mc_Tree_signal_template/");
  TString bkg_dir("data_Tree_sieiesideband_sel/");

  if (splitting=="EEEB") splitting="EBEE";

  
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

    h_sig1hist->Rebin(rbin);
    h_bkg1hist->Rebin(rbin);
    h_datahist_2D->Rebin2D(rbin,rbin);

  
  RooRealVar roovar("roovar","roovar",-5,5);
  RooRealVar roovar1("roovar1","roovar1",-5,5);
  RooRealVar roovar2("roovar2","roovar2",-5,5);

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
  
  RooAddPdf *model;


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


  if (study_templates){
    RooDataHist *datahist_toy;
    datahist_toy = model->generateBinned(RooArgSet(roovar1,roovar2),100000,Name("Toy dataset"));
    datahist=datahist_toy;
  }

  roovar.setRange("fullrange",-5,5);
  roovar1.setRange("fullrange",-5,5);
  roovar2.setRange("fullrange",-5,5);

  roovar.setRange("fitrange",leftrange,rightrange);
  roovar1.setRange("fitrange",leftrange,rightrange);
  roovar2.setRange("fitrange",leftrange,rightrange);

  RooFitResult *fitres = model->fitTo(*datahist,Save(),Range("fitrange",kFALSE));
  model->Print();

  out->fr=fitres;

  out->tot_events=h_datahist_2D->Integral();


  //  canv=NULL;
  if (canv!=NULL){

    RooPlot *varframe[2];
    varframe[0]=roovar1.frame();
    varframe[1]=roovar2.frame();
    int i=0;
    canv->cd(bin+i+1);
    datahist->plotOn(varframe[i]);
    model->plotOn(varframe[i],Range("fitrange"),NormRange("fullrange"));
    model->plotOn(varframe[i],Components("sigsigpdf"),LineStyle(kDashed),LineColor(kRed),Range("fitrange"),NormRange("fullrange"));
    model->plotOn(varframe[i],Components("sigbkgpdf"),LineStyle(kDashed),LineColor(kGreen),Range("fitrange"),NormRange("fullrange"));
    model->plotOn(varframe[i],Components("bkgbkgpdf"),LineStyle(kDashed),LineColor(kBlack),Range("fitrange"),NormRange("fullrange"));
    varframe[i]->Draw();

  }



  return out;

};




void run_fits(TString inputfilename_t="templates.root", TString inputfilename_d="tobefitted.root", TString diffvariable="", TString splitting="", float leftrange=-5, float rightrange=5, int rebin=1){

  TH1F::SetDefaultSumw2(kTRUE);

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

  
  fit_output *fr[n_bins];

  TCanvas *fits_canv;

  TH1F *purity[3];
  TH1F *eff;
  TH1F *xsec;

  TFile *eff_file = new TFile("efficiencies.root");
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
    xsec->GetXaxis()->SetTitle("diffvariable.Data()");


    fits_canv = new TCanvas("fits","fits");
    if (bins_to_run==4) fits_canv->Divide(2,3);
    if (bins_to_run==6) fits_canv->Divide(2,3);
    if (bins_to_run==8) fits_canv->Divide(2,4);
    if (bins_to_run==9) fits_canv->Divide(2,5);
    if (bins_to_run==12) fits_canv->Divide(2,6);
    if (bins_to_run==15) fits_canv->Divide(2,8);




  RooRealVar *rf1;
  RooRealVar *rf2;

  for (int bin=0; bin<bins_to_run; bin++) {

    fr[bin]=fit_dataset(inputfilename_t.Data(),inputfilename_d.Data(),diffvariable,splitting,leftrange,rightrange,bin,fits_canv,rebin);
    
    float intlumi=4.519;

	rf1=(RooRealVar*)(fr[bin]->fr->floatParsFinal().find("rf1"));
	rf2=(RooRealVar*)(fr[bin]->fr->floatParsFinal().find("rf2"));
	RooFormulaVar fsigsig("fsigsig","fsigsig","rf1",RooArgList(*rf1));
	RooFormulaVar fsigbkg("fsigbkg","fsigbkg","(1-rf1)*rf2",RooArgList(*rf1,*rf2));
	RooFormulaVar fbkgbkg("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg",RooArgList(fsigsig,fsigbkg));
	//	std::cout << bin << " " << fsigsig.getVal() << " " << fsigbkg.getVal() << " " << fbkgbkg.getVal() << std::endl;

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

  fits_canv->Update();

  TCanvas *output_canv = new TCanvas("output_canv","output_canv");
  output_canv->cd();

  purity[0]->Draw("e1");
  purity[1]->Draw("e1same");
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

  TFile *out1 = new TFile(Form("purity_%s_data.root",splitting.Data()),"recreate");
  out1->cd();
  purity[0]->Write();
  purity[1]->Write();
  purity[2]->Write();

  TFile *outfile = new TFile(Form("xsec_%s_data.root",splitting.Data()),"recreate");
  outfile->cd();
  xsec->Write();
  

};


