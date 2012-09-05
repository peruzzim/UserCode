#include <iostream>
#include <vector>
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "TStyle.h"


using namespace RooFit;

std::vector<TString> gethistoname(TString pref1, TString dset1, TString temp1, TString reg="EB"){

  std::vector<TString> out(2,"");

  TString file1=pref1;
  file1.Append("_");
  file1.Append(dset1);
  file1.Append("_");
  file1.Append(temp1);
  file1.Append(".root");

  out[0]=file1;

  TString name1="";
  if (dset1=="data") name1.Append("data_Tree_"); else name1.Append("mc_Tree_");
  if (temp1=="bkg") name1.Append("background_template/template_background_");
  if (temp1=="sig") name1.Append("signal_template/template_signal_");
  if (temp1=="rcone") name1.Append("randomcone_signal_template/template_signal_");
  if (temp1=="impinging") name1.Append("impinging_track_template/template_background_");
  if (temp1=="sieiesideband") name1.Append("sieiesideband_sel/template_background_");
  if (temp1=="combisosideband") name1.Append("combisosideband_sel/template_background_");
  name1.Append(reg);
  name1.Append("_b9");

  out[1]=name1;

  return out;

};

void doBiasStudy(\
		 TString pref1, TString dset1, TString temp1,	\
		 TString pref2, TString dset2, TString temp2,	\
		 TString pref3, TString dset3, TString temp3,	\
		 TString pref4, TString dset4, TString temp4,	\
		 TString reg="EB", int nsamples = 10, int nevents_persample = 100){
  
  std::vector<TString> names_gen[2];
  std::vector<TString> names_fit[2];

  names_gen[0]=gethistoname(pref1,dset1,temp1,reg);
  names_gen[1]=gethistoname(pref2,dset2,temp2,reg);
  names_fit[0]=gethistoname(pref3,dset3,temp3,reg);
  names_fit[1]=gethistoname(pref4,dset4,temp4,reg);

  TFile *f1_gen = new TFile(names_gen[0][0].Data(),"read");
  TFile *f2_gen = new TFile(names_gen[1][0].Data(),"read");
  TH1F *hs_gen;
  TH1F *hb_gen;
  f1_gen->GetObject(names_gen[0][1].Data(),hs_gen);
  assert (hs_gen!=NULL);
  f2_gen->GetObject(names_gen[1][1].Data(),hb_gen);
  assert (hb_gen!=NULL);

  TFile *f1_fit = new TFile(names_fit[0][0].Data(),"read");
  TFile *f2_fit = new TFile(names_fit[1][0].Data(),"read");
  TH1F *hs_fit;
  TH1F *hb_fit;
  f1_fit->GetObject(names_fit[0][1].Data(),hs_fit);
  assert (hs_fit!=NULL);
  f2_fit->GetObject(names_fit[1][1].Data(),hb_fit);
  assert (hb_fit!=NULL);


  RooRealVar *photoniso = new RooRealVar("photoniso","PFphotonIso (GeV)",0,5);
  //  photoniso->setBins(40);
  RooRealVar sigfrac("sigfrac","sigfrac",0.5,0,1);
  RooRealVar nevents("nevents","nevents",1,0,1e+6);
  RooFormulaVar nsig("nsig","nevents*sigfrac",RooArgList(sigfrac,nevents));
  RooFormulaVar nbkg("nbkg","nevents*(1-sigfrac)",RooArgList(sigfrac,nevents));


  RooDataHist datahs_gen("datahs_gen","datahs_gen",*photoniso,Import(*hs_gen));
  RooDataHist datahb_gen("datahb_gen","datahb_gen",*photoniso,Import(*hb_gen));
  RooHistPdf pdfhs_gen("pdfhs_gen","pdfhs_gen",*photoniso,datahs_gen);
  RooHistPdf pdfhb_gen("pdfhb_gen","pdfhb_gen",*photoniso,datahb_gen);
  RooAddPdf genmodel("genmodel","genmodel",RooArgList(pdfhs_gen,pdfhb_gen),RooArgList(nsig,nbkg));


  RooDataHist datahs_fit("datahs_fit","datahs_fit",*photoniso,Import(*hs_fit));
  RooDataHist datahb_fit("datahb_fit","datahb_fit",*photoniso,Import(*hb_fit));
  RooHistPdf pdfhs_fit("pdfhs_fit","pdfhs_fit",*photoniso,datahs_fit);
  RooHistPdf pdfhb_fit("pdfhb_fit","pdfhb_fit",*photoniso,datahb_fit);
  RooAddPdf fitmodel("fitmodel","fitmodel",RooArgList(pdfhs_fit,pdfhb_fit),RooArgList(nsig,nbkg));


  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);


//  TCanvas *c0 = new TCanvas();
//  c0->cd();
//  RooPlot *f0 = photoniso->frame(Title("photonisoframe"));


  RooMCStudy *mcstudy = new RooMCStudy(genmodel,*photoniso,FitModel(fitmodel),Silence(),Extended(),FitOptions(Save(kTRUE),PrintEvalErrors(0)));
  mcstudy->generateAndFit(nsamples, nevents_persample);

  RooPlot* frame1 = mcstudy->plotParam(sigfrac);
  RooPlot* frame2 = mcstudy->plotError(sigfrac);
  RooPlot* frame3 = mcstudy->plotPull(sigfrac,Bins(40),FitGauss(kTRUE));


  TCanvas* c = new TCanvas();
  c->Divide(3);
  c->cd(1); frame1->Draw();
  c->cd(2); frame2->Draw();
  c->cd(3); frame3->Draw();

  //  Bool_t generateAndFit(Int_t nSamples, Int_t nEvtPerSample = 0, Bool_t keepGenData = kFALSE, const char* asciiFilePat = 0)
  //  RooAbsData* data = mcstudy->genData(0);


//  mcstudy->genData(0)->plotOn(f0);
//  genmodel.plotOn(f0,LineColor(kGreen));//,DataError(RooAbsData::SumW2));
//  genmodel.plotOn(f0,Components(pdfhb),LineColor(kRed));//,DataError(RooAbsData::SumW2));
//  genmodel.plotOn(f0,Components(pdfhs),LineColor(kBlue));//,DataError(RooAbsData::SumW2));
//  mcstudy->genData(0)->plotOn(f0);
//
//
//  f0->Draw();
//



  gDirectory->Add(mcstudy);
  gDirectory->Add(photoniso);

};
