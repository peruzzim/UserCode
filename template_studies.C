#include <assert.h>

using namespace std;
using namespace RooFit;

typedef struct {
  RooFitResult *fr;
  float purity;
  float purity5;
  float efficiency;
  float histo_population;
  float histo_population5;
  float sig_events;
} fit_output; 

  bool study_templates=false;

void run_fits(TString inputfilename="out_NEW.root", TString splitting, float leftrange=-5, float rightrange=35){

  bool single_gamma;
  if (splitting=="EB" || splitting=="EE") single_gamma=true; else single_gamma=false;

  const int n_bins=5;
  int bins_to_run=3; if (single_gamma) bins_to_run=1;

  fit_output *fr[n_bins];

  TGraphErrors *out = new TGraphErrors(bins_to_run);
  out->SetMarkerStyle(20);
  out->SetMarkerColor(kRed);
  out->SetLineColor(kRed);
  out->SetLineWidth(2);

  TGraphErrors *out2 = new TGraphErrors(bins_to_run);
  out2->SetMarkerStyle(20);
  out2->SetMarkerColor(kBlue);
  out2->SetLineColor(kBlue);
  out2->SetLineWidth(2);


  TCanvas *fits_canv = new TCanvas("fits","fits");
  fits_canv->Divide(2,bins_to_run);

  RooRealVar *rf1;
  RooRealVar *rf2;

  for (int bin=0; bin<bins_to_run; bin++) {
    fr[bin]=fit_dataset(inputfilename.Data(),splitting.Data(),leftrange,rightrange,bin,fits_canv);

    if (!single_gamma){  
      rf1=(RooRealVar*)(fr[bin]->fr->floatParsFinal().find("rf1"));
      rf2=(RooRealVar*)(fr[bin]->fr->floatParsFinal().find("rf2"));
      RooFormulaVar fsigsig("fsigsig","fsigsig","rf1",RooArgList(*rf1));
      RooFormulaVar fsigbkg("fsigbkg","fsigbkg","(1-rf1)*rf2",RooArgList(*rf1,*rf2));
      RooFormulaVar fbkgbkg("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg",RooArgList(fsigsig,fsigbkg));
      std::cout << bin << " " << fsigsig.getVal() << " " << fsigbkg.getVal() << " " << fbkgbkg.getVal() << std::endl;
      out->SetPoint(bin,bin,fsigsig.getVal());
      out->SetPointError(bin,0,fsigsig.getPropagatedError(*(fr[bin]->fr)));
    }

    if (single_gamma){
      rf1=(RooRealVar*)(fr[bin]->fr->floatParsFinal().find("rf1"));
      RooFormulaVar fsig("fsig","fsig","rf1",RooArgList(*rf1));
      RooFormulaVar fbkg("fbkg","fbkg","1-rf1",RooArgList(*rf1));
      out->SetPoint(bin,bin,fsig.getVal());
      //      out->SetPoint(bin,bin,fr[bin]->purity5);
      out->SetPointError(bin,0,fsig.getPropagatedError(*(fr[bin]->fr)));
      out2->SetPoint(bin,bin,fr[bin]->purity5);
    }
  }

  fits_canv->Update();

  if (!study_templates){
  TCanvas *output_canv = new TCanvas("output_canv","output_canv");
  output_canv->cd();
  out->GetYaxis()->SetRangeUser(0,1);
  out->Draw("AP");
  out2->Draw("P");
  output_canv->Update();
  }
		

};


fit_output* fit_dataset(const char* inputfilename, TString splitting, float leftrange, float rightrange, int bin, TCanvas *canv=NULL){
 
  bool single_gamma;

  TFile *inputfile = TFile::Open(inputfilename);
  TFile *inputfile_noselection = TFile::Open("test_signal_noselection.root");


  TString data_dir("data_Tree_standard_sel/");
  TString mc_dir("mc_Tree_standard_sel/");

  fit_output *out = new fit_output();


  if (splitting=="EB" || splitting=="EE") single_gamma=true; else single_gamma=false;
  if (splitting=="EEEB") splitting="EBEE";

  TH1F *h_sighist;
  TH1F *h_sighist_noselection;
  TH1F *h_bkghist;
  TH1F *h_datahist;
  TH2F *h_datahist_2D;

  TH1F *h_sig1hist;
  TH1F *h_sig2hist;
  TH1F *h_bkg1hist;
  TH1F *h_bkg2hist;
  //  TH1F *h_datahist;

  
  if (single_gamma) {
    inputfile->GetObject(TString(mc_dir).Append(Form("template_signal_%s_b%d",splitting.Data(),bin)),h_sighist);
    inputfile_noselection->GetObject(TString(mc_dir).Append(Form("template_signal_%s_b%d",splitting.Data(),bin)),h_sighist_noselection);
    inputfile->GetObject(TString(mc_dir).Append(Form("template_background_%s_b%d",splitting.Data(),bin)),h_bkghist);
    inputfile->GetObject(TString(data_dir).Append(Form("obs_hist_single_%s_b%d",splitting.Data(),bin)),h_datahist);
    assert (h_sighist!=NULL);
    assert (h_bkghist!=NULL);
    assert (h_datahist!=NULL);
  }
  if (!single_gamma) {
    TString s1; TString s2;
    if (splitting=="EBEB") {s1="EB"; s2="EB";}
    else if (splitting=="EEEE") {s1="EE"; s2="EE";}
    else if (splitting=="EBEE") {s1="EB"; s2="EE";}
    inputfile->GetObject(TString(mc_dir).Append(Form("template_signal_%s_b%d",s1.Data(),bin)),h_sig1hist);
    inputfile->GetObject(TString(mc_dir).Append(Form("template_background_%s_b%d",s1.Data(),bin)),h_bkg1hist);
    inputfile->GetObject(TString(mc_dir).Append(Form("template_signal_%s_b%d",s2.Data(),bin)),h_sig2hist);
    inputfile->GetObject(TString(mc_dir).Append(Form("template_background_%s_b%d",s2.Data(),bin)),h_bkg2hist);
    inputfile->GetObject(TString(data_dir).Append(Form("obs_hist_%s_b%d",splitting.Data(),bin)),h_datahist_2D);
    assert (h_sig1hist!=NULL);
    assert (h_bkg1hist!=NULL);
    assert (h_sig2hist!=NULL);
    assert (h_bkg2hist!=NULL);
    assert (h_datahist_2D!=NULL);
  }


  RooRealVar roovar("roovar","roovar",-5,35);
  RooRealVar roovar1("roovar1","roovar1",-5,35);
  RooRealVar roovar2("roovar2","roovar2",-5,35);

  RooRealVar rf1("rf1","rf1",1,0,1);
  RooRealVar rf2("rf2","rf2",1,0,1);
  RooRealVar rf3("rf3","rf3",1,0,1);

  if (!single_gamma){
      rf1.setVal(1./3);
      rf2.setVal(1./2);
  }

  RooDataHist *sighist;
  RooDataHist *bkghist;
  RooDataHist *datahist;


  RooHistPdf *sigpdf;
  RooHistPdf *bkgpdf;

  RooFormulaVar *fsig;
  RooFormulaVar *fbkg;

  RooAddPdf *model;

if (single_gamma){
  
  sighist   = new RooDataHist("sighist","sighist",RooArgList(roovar),h_sighist,1.0/h_sighist->Integral());
  bkghist   = new RooDataHist("bkghist","bkghist",RooArgList(roovar),h_bkghist,1.0/h_bkghist->Integral());
  datahist   = new RooDataHist("datahist","datahist",RooArgList(roovar),h_datahist);
  
  sigpdf = new RooHistPdf("sigpdf","sigpdf",roovar,*sighist);
  bkgpdf = new RooHistPdf("bkgpdf","bkgpdf",roovar,*bkghist);
  
  fsig = new RooFormulaVar("fsig","fsig","rf1",RooArgList(rf1));
  fbkg = new RooFormulaVar("fbkg","fbkg","1-rf1",RooArgList(rf1));

  model = new RooAddPdf("model","model",RooArgList(*sigpdf,*bkgpdf),RooArgList(rf1),kTRUE);

 }


 RooDataHist *sig1hist;
 RooDataHist *sig2hist;
 RooDataHist *bkg1hist;
 RooDataHist *bkg2hist;
 // RooDataHist *datahist;

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

    //  RooAddPdf *model;

  if (!single_gamma){

    sig1hist = new RooDataHist("sighist","sighist",RooArgList(roovar1),h_sig1hist);
    sig2hist = new RooDataHist("sighist","sighist",RooArgList(roovar2),h_sig2hist);
    bkg1hist = new RooDataHist("bkghist","bkghist",RooArgList(roovar1),h_bkg1hist);
    bkg2hist = new RooDataHist("bkghist","bkghist",RooArgList(roovar2),h_bkg2hist);
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


//  if (study_templates){
//    RooDataSet *roodset_toy = model->generate(RooArgList(*(roovar[0]),*(roovar[1])),Name("toy_model"),NumEvents(1000000),AutoBinned(kTRUE));
//    roodset=roodset_toy;
//  }

  RooFitResult *fitres = model->fitTo(*datahist,Save(),Range(leftrange,rightrange,kFALSE));

  out->fr=fitres;
  out->purity=rf1.getVal();

  float integral5_tot;
  float integral5_sig;

  //  if (single_gamma){

    //    roovar.setRange("range5",-5,5);
    //    integral5_tot=model->createIntegral(RooArgSet(roovar),Range("range5"))->getVal();
    //    integral5_sig=sigpdf->createIntegral(RooArgSet(roovar),Range("range5"))->getVal()*(rf1.getVal());

    //    integral5_tot=model->createIntegral(RooArgSet(roovar),Range("range5"))->getVal();
    //    integral5_sig=sigpdf->createIntegral(RooArgSet(roovar),Range("range5"))->getVal();
//    out->purity5=integral5_sig/integral5_tot;
//    
//    out->histo_population=h_datahist->Integral();
//    out->histo_population5=h_datahist->Integral(h_datahist->FindBin(-5),h_datahist->FindBin(5));
//
//    out->efficiency=h_sighist->Integral()/h_sighist_noselection->Integral();

    //  }


  //  model->Print();

    roovar.setRange("range5",-5,35);
  cout << model->createIntegral(RooArgSet(roovar),Range("range5"))->getVal() << std::endl;
  cout << sigpdf->createIntegral(RooArgSet(roovar),Range("range5"))->getVal() << std::endl;
  cout << bkgpdf->createIntegral(RooArgSet(roovar),Range("range5"))->getVal() << std::endl;

  cout << h_datahist->Integral() << endl;
  cout << h_bkghist->Integral() << endl;
  cout << h_sighist->Integral() << endl;

  cout << "-" << endl;

  std::cout << out->purity << std::endl;
  std::cout << out->histo_population << std::endl;
  std::cout << out->histo_population5 << std::endl;
  std::cout << integral5_sig << std::endl;
  std::cout << integral5_tot << std::endl;
  std::cout << out->purity5 << std::endl;
  std::cout << out->efficiency << std::endl;

  if (canv!=NULL){

    RooPlot *varframe[2];

    if (single_gamma){
      canv->cd(2*bin+1);
      varframe[0]=roovar.frame("","");
      datahist->plotOn(varframe[0]);
      model->plotOn(varframe[0]);
      model->plotOn(varframe[0],Components("sigpdf"),LineStyle(kDashed),LineColor(kRed));
      model->plotOn(varframe[0],Components("bkgpdf"),LineStyle(kDashed),LineColor(kBlack));
      varframe[0]->Draw();
    }

    if (!single_gamma){
      varframe[0]=roovar1.frame("","");
      varframe[1]=roovar2.frame("","");
      for (int i=0; i<2; i++){
      canv->cd(2*bin+i+1);
      datahist->plotOn(varframe[i]);
      model->plotOn(varframe[i]);
      model->plotOn(varframe[i],Components("sigsigpdf"),LineStyle(kDashed),LineColor(kRed));
      model->plotOn(varframe[i],Components("sigbkgpdf"),LineStyle(kDashed),LineColor(kGreen));
      model->plotOn(varframe[i],Components("bkgbkgpdf"),LineStyle(kDashed),LineColor(kBlack));
      varframe[i]->Draw();
      }
    }

  }



  return out;
  
};


