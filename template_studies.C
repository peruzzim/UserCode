#include <assert.h>

using namespace std;
using namespace RooFit;

typedef struct {
  RooFitResult *fr;
  float purity;
  float efficiency;
  float sig_events;
  float tot_events;
  float purity5;
  float efficiency5;
  float sig_events5;
  float tot_events5;
} fit_output; 

bool study_templates=false;

static const int n_bins=15;

float binsdef_single_gamma[n_bins+1]={30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180};
float binsdef_diphoton[n_bins+1]={80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230};

void run_fits(TString inputfilename_t="out_NEW.root", TString inputfilename_d_mc="out_NEW_mc.root", TString inputfilename_d_data="out_NEW_data.root", TString splitting, float leftrange=-5, float rightrange=35, int rebin=1){

  TH1F::SetDefaultSumw2(kTRUE);

  bool single_gamma;
  if (splitting=="EB" || splitting=="EE") single_gamma=true; else single_gamma=false;


  int bins_to_run=15; if (single_gamma) bins_to_run=15;

  fit_output *fr[n_bins][2];

  TCanvas *fits_canv[2];

  TH1F *purity[2];
  TH1F *purity5[2];
  TH1F *eff[2];
  TH1F *xsec[2];


  for (int i=0; i<2; i++){
    TString name;

    name="purity";
    purity[i] = new TH1F(Form("%s_%d",name.Data(),i),Form("%s_%d",name.Data(),i),n_bins,single_gamma ? binsdef_single_gamma : binsdef_diphoton);
    purity[i]->SetMarkerStyle(20);
    purity[i]->SetMarkerColor(kRed);
    purity[i]->SetLineColor(kRed);
    purity[i]->SetLineWidth(2);
    purity[i]->GetYaxis()->SetRangeUser(0,1);
    if (i==0) purity[i]->SetLineStyle(kDashed);

    purity[i]->GetYaxis()->SetTitle("purity  (red) / purity < 5 GeV (blue)");
    if (single_gamma)  purity[i]->GetXaxis()->SetTitle("p_{T}");
    if (!single_gamma) purity[i]->GetXaxis()->SetTitle("m_{#gamma #gamma}");

    name="purity5";
    purity5[i] = new TH1F(Form("%s_%d",name.Data(),i),Form("%s_%d",name.Data(),i),n_bins,single_gamma ? binsdef_single_gamma : binsdef_diphoton);
    purity5[i]->SetMarkerStyle(20);
    purity5[i]->SetMarkerColor(kBlue);
    purity5[i]->SetLineColor(kBlue);
    purity5[i]->SetLineWidth(2);
    if (i==0) purity5[i]->SetLineStyle(kDashed);

    if (single_gamma)  purity5[i]->GetXaxis()->SetTitle("p_{T}");
    if (!single_gamma) purity5[i]->GetXaxis()->SetTitle("m_{#gamma #gamma}");

    name="efficiency";
    eff[i] = new TH1F(Form("%s_%d",name.Data(),i),Form("%s_%d",name.Data(),i),n_bins,single_gamma ? binsdef_single_gamma : binsdef_diphoton);
    eff[i]->SetMarkerStyle(20);
    eff[i]->SetLineWidth(2);
    eff[i]->GetYaxis()->SetRangeUser(0,1);
    if (i==0) eff[i]->SetLineStyle(kDashed);

    eff[i]->GetYaxis()->SetTitle("selection/ID efficiency (H/e, sieie)");
    if (single_gamma)  eff[i]->GetXaxis()->SetTitle("p_{T}");
    if (!single_gamma) eff[i]->GetXaxis()->SetTitle("m_{#gamma #gamma}");

    name="events";
    xsec[i] = new TH1F(Form("%s_%d",name.Data(),i),Form("%s_%d",name.Data(),i),n_bins,single_gamma ? binsdef_single_gamma : binsdef_diphoton);
    xsec[i]->SetMarkerStyle(20);
    xsec[i]->SetMarkerColor(kGreen);
    xsec[i]->SetLineColor(kGreen);
    xsec[i]->SetLineWidth(2);
    if (i==0) xsec[i]->SetLineStyle(kDashed);

    xsec[i]->GetYaxis()->SetTitle("events after efficiency corr.");
    if (single_gamma)  xsec[i]->GetXaxis()->SetTitle("p_{T}");
    if (!single_gamma) xsec[i]->GetXaxis()->SetTitle("m_{#gamma #gamma}");

    name="fits";
    fits_canv[i] = new TCanvas(Form("%s_%d",name.Data(),i),Form("%s_%d",name.Data(),i));
    if (bins_to_run==6) fits_canv[i]->Divide(2,3);
    if (bins_to_run==9) fits_canv[i]->Divide(2,5);
    if (bins_to_run==12) fits_canv[i]->Divide(2,6);
    if (bins_to_run==15) fits_canv[i]->Divide(2,8);

  }

  

  RooRealVar *rf1;
  RooRealVar *rf2;

  for (int bin=0; bin<bins_to_run; bin++) {

    fr[bin][0]=fit_dataset(inputfilename_t.Data(),inputfilename_d_mc.Data(),splitting.Data(),leftrange,rightrange,bin,fits_canv[0],rebin,n_bins);
    fr[bin][1]=fit_dataset(inputfilename_t.Data(),inputfilename_d_data.Data(),splitting.Data(),leftrange,rightrange,bin,fits_canv[1],rebin,n_bins);
    

    for (int i=0; i<2; i++){
      if (!single_gamma){  
	rf1=(RooRealVar*)(fr[bin][i]->fr->floatParsFinal().find("rf1"));
	rf2=(RooRealVar*)(fr[bin][i]->fr->floatParsFinal().find("rf2"));
	RooFormulaVar fsigsig("fsigsig","fsigsig","rf1",RooArgList(*rf1));
	RooFormulaVar fsigbkg("fsigbkg","fsigbkg","(1-rf1)*rf2",RooArgList(*rf1,*rf2));
	RooFormulaVar fbkgbkg("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg",RooArgList(fsigsig,fsigbkg));
	std::cout << bin << " " << fsigsig.getVal() << " " << fsigbkg.getVal() << " " << fbkgbkg.getVal() << std::endl;
	purity[i]->SetBinContent(bin+1,fsigsig.getVal());
	purity[i]->SetBinError(bin+1,fsigsig.getPropagatedError(*(fr[bin][i]->fr)));
	purity5[i]->SetBinContent(bin+1,fr[bin][i]->purity5);
	xsec[i]->SetBinContent(bin+1,fr[bin][i]->sig_events/fr[bin][i]->efficiency);
	eff[i]->SetBinContent(bin+1,fr[bin][i]->efficiency);
      }

      if (single_gamma){
	rf1=(RooRealVar*)(fr[bin][i]->fr->floatParsFinal().find("rf1"));
	RooFormulaVar fsig("fsig","fsig","rf1",RooArgList(*rf1));
	RooFormulaVar fbkg("fbkg","fbkg","1-rf1",RooArgList(*rf1));
	purity[i]->SetBinContent(bin+1,fsig.getVal());
	purity[i]->SetBinError(bin+1,fsig.getPropagatedError(*(fr[bin][i]->fr)));
	purity5[i]->SetBinContent(bin+1,fr[bin][i]->purity5);
	xsec[i]->SetBinContent(bin+1,fr[bin][i]->sig_events/fr[bin][i]->efficiency);
	eff[i]->SetBinContent(bin+1,fr[bin][i]->efficiency);
      }
    }
  }

  fits_canv[0]->Update();
  fits_canv[1]->Update();


  TCanvas *output_canv = new TCanvas("output_canv","output_canv");
  output_canv->cd();

  purity[0]->Draw("e1");
  purity5[0]->Draw("e1same");
  purity[1]->Draw("e1same");
  purity5[1]->Draw("e1same");

  output_canv->Update();

		
   
  TCanvas *xsec_canv = new TCanvas("xsec_canv","xsec_canv");
  xsec_canv->cd();
  xsec[0]->Draw("e1");
  xsec[1]->Draw("e1same");
  xsec_canv->Update();

  TCanvas *eff_canv = new TCanvas("eff_canv","eff_canv");
  eff_canv->cd();
  eff[0]->Draw("e1");
  eff[1]->Draw("e1same");
  eff_canv->Update();

};

fit_output* bla(){

  return NULL;

};

fit_output* fit_dataset(const char* inputfilename_t, const char* inputfilename_d, TString splitting, float leftrange, float rightrange, int bin, TCanvas *canv=NULL, int rbin, int n_bins){

fit_output *out = new fit_output();

  bool single_gamma;

  TFile *inputfile_t = TFile::Open(inputfilename_t);
  TFile *inputfile_d = TFile::Open(inputfilename_d);
  TFile *inputfile_noselection = TFile::Open("test_signal_noselection.root");


  TString data_dir("data_Tree_standard_sel/");
  TString mc_dir("mc_Tree_standard_sel/");

  //  TString data_dir("data_Tree_DY_sel/");
  //  TString mc_dir("mc_Tree_DY_sel/");




  if (splitting=="EB" || splitting=="EE") single_gamma=true; else single_gamma=false;
  if (splitting=="EEEB") splitting="EBEE";

  TH1F *h_sighist;
  TH1F *h_sighist_selection;
  TH1F *h_sighist_noselection;
  TH1F *h_bkghist;
  TH1F *h_datahist;
  TH2F *h_datahist_2D;

  TH1F *h_sig1hist;
  TH1F *h_sig2hist;
  TH1F *h_bkg1hist;
  TH1F *h_bkg2hist;
  //  TH1F *h_datahist;

  TH2F *h_sigsighist_selection;
  TH2F *h_sigsighist_noselection;

  
  if (single_gamma) {
    inputfile_t->GetObject(TString(mc_dir).Append(Form("template_signal_%s_b%d",splitting.Data(),bin)),h_sighist_selection);
    inputfile_noselection->GetObject(TString(mc_dir).Append(Form("template_signal_%s_b%d",splitting.Data(),bin)),h_sighist_noselection);

    inputfile_t->GetObject(TString(mc_dir).Append(Form("template_signal_%s_b%d",splitting.Data(),n_bins)),h_sighist);
    inputfile_t->GetObject(TString(mc_dir).Append(Form("template_background_%s_b%d",splitting.Data(),n_bins)),h_bkghist);

    inputfile_d->GetObject(TString(data_dir).Append(Form("obs_hist_single_%s_b%d",splitting.Data(),bin)),h_datahist);
    assert (h_sighist!=NULL);
    assert (h_bkghist!=NULL);
    assert (h_datahist!=NULL);
    h_sighist->Rebin(rbin);
    h_sighist_noselection->Rebin(rbin);
    h_bkghist->Rebin(rbin);
    h_datahist->Rebin(rbin);
  }
  if (!single_gamma) {
    TString s1; TString s2;
    if (splitting=="EBEB") {s1="EB"; s2="EB";}
    else if (splitting=="EEEE") {s1="EE"; s2="EE";}
    else if (splitting=="EBEE") {s1="EB"; s2="EE";}
    inputfile_t->GetObject(TString(mc_dir).Append(Form("template_signal_%s_b%d",s1.Data(),n_bins)),h_sig1hist);
    inputfile_t->GetObject(TString(mc_dir).Append(Form("template_background_%s_b%d",s1.Data(),n_bins)),h_bkg1hist);
    inputfile_t->GetObject(TString(mc_dir).Append(Form("template_signal_%s_b%d",s2.Data(),n_bins)),h_sig2hist);
    inputfile_t->GetObject(TString(mc_dir).Append(Form("template_background_%s_b%d",s2.Data(),n_bins)),h_bkg2hist);
    inputfile_t->GetObject(TString(mc_dir).Append(Form("template_sigsig_%s_b%d",splitting.Data(),bin)),h_sigsighist_selection);
    inputfile_noselection->GetObject(TString(mc_dir).Append(Form("template_sigsig_%s_b%d",splitting.Data(),bin)),h_sigsighist_noselection);
    inputfile_d->GetObject(TString(data_dir).Append(Form("obs_hist_%s_b%d",splitting.Data(),bin)),h_datahist_2D);
    assert (h_sig1hist!=NULL);
    assert (h_bkg1hist!=NULL);
    assert (h_sig2hist!=NULL);
    assert (h_bkg2hist!=NULL);
    assert (h_datahist_2D!=NULL);
    h_sig1hist->Rebin(rbin);
    h_bkg1hist->Rebin(rbin);
    //    h_sig2hist->Rebin(rbin);
    //    h_bkg2hist->Rebin(rbin);
    cout << "a" << endl;
    h_datahist_2D->Rebin2D(rbin,rbin);
  }


  RooRealVar roovar("roovar","roovar",-5,35);
  RooRealVar roovar1("roovar1","roovar1",-5,35);
  RooRealVar roovar2("roovar2","roovar2",-5,35);

  RooRealVar rf1("rf1","rf1",0.5,0,1);
  RooRealVar rf2("rf2","rf2",0.5,0,1);
  RooRealVar rf3("rf3","rf3",0.5,0,1);

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


  if (study_templates){
    RooDataHist *datahist_toy;
    if (!single_gamma) datahist_toy = model->generateBinned(RooArgSet(roovar1,roovar2),100000,Name("Toy dataset"));
    if (single_gamma) datahist_toy = model->generateBinned(RooArgSet(roovar),100000,Name("Toy dataset"));
    datahist=datahist_toy;
  }

  roovar.setRange("fullrange",-5,35);
  roovar1.setRange("fullrange",-5,35);
  roovar2.setRange("fullrange",-5,35);


  roovar.setRange("fitrange",leftrange,rightrange);
  roovar1.setRange("fitrange",leftrange,rightrange);
  roovar2.setRange("fitrange",leftrange,rightrange);

  RooFitResult *fitres = model->fitTo(*datahist,Save(),Range("fitrange",kFALSE));
  model->Print();

  out->fr=fitres;

  if (!single_gamma){

    out->purity=fsigsig->getVal();

    out->tot_events=h_datahist_2D->Integral();
    out->sig_events=h_datahist_2D->Integral()*out->purity;

    roovar1.setRange("range5",-5,5);
    roovar2.setRange("range5",-5,5);

    out->sig_events5=out->sig_events*sigsigpdf->createIntegral(RooArgSet(roovar1,roovar2),Range("range5"))->getVal();

    cout << sigsigpdf->createIntegral(RooArgSet(roovar1,roovar2),Range("range5"))->getVal() << endl;
    cout << sigsigpdf->createIntegral(RooArgSet(roovar1,roovar2))->getVal() << endl;

    float sigbkg_events5=h_datahist_2D->Integral()*fsigbkg->getVal()*sigbkgpdf_noorder->createIntegral(RooArgSet(roovar1,roovar2),Range("range5"))->getVal();
    float bkgbkg_events5=h_datahist_2D->Integral()*fbkgbkg->getVal()*bkgbkgpdf->createIntegral(RooArgSet(roovar1,roovar2),Range("range5"))->getVal();

    out->tot_events5=out->sig_events5+sigbkg_events5+bkgbkg_events5;

    out->purity5=out->sig_events5/out->tot_events5;
  

    out->efficiency=h_sigsighist_selection->Integral()/h_sigsighist_noselection->Integral();

    int xbl=h_sigsighist_selection->GetXaxis()->FindBin(-5);
    int xbh=h_sigsighist_selection->GetXaxis()->FindBin(5);
    int ybl=h_sigsighist_selection->GetYaxis()->FindBin(-5);
    int ybh=h_sigsighist_selection->GetYaxis()->FindBin(5);

    out->efficiency5=h_sigsighist_selection->Integral(xbl,xbh,ybl,ybh)/h_sigsighist_noselection->Integral(xbl,xbh,ybl,ybh);


  }


  if (single_gamma){

    out->purity=fsig->getVal();

    out->tot_events=h_datahist->Integral();
    out->sig_events=h_datahist->Integral()*out->purity;
  
    out->sig_events5=out->sig_events*h_sighist->Integral(h_datahist->FindBin(-5),h_datahist->FindBin(5))/h_sighist->Integral();
    out->tot_events5=out->sig_events5+out->tot_events*(1-out->purity)*h_bkghist->Integral(h_datahist->FindBin(-5),h_datahist->FindBin(5))/h_bkghist->Integral();

    out->purity5=out->sig_events5/out->tot_events5;
  
    out->efficiency=h_sighist_selection->Integral()/h_sighist_noselection->Integral();
    out->efficiency5=h_sighist_selection->Integral(h_datahist->FindBin(-5),h_datahist->FindBin(5))/h_sighist_noselection->Integral(h_datahist->FindBin(-5),h_datahist->FindBin(5));


  }


  cout << "-" << endl;
  std::cout << out->purity << std::endl;
  std::cout << out->tot_events << std::endl;
  std::cout << out->sig_events << std::endl;
  std::cout << out->efficiency << std::endl;
  std::cout << out->purity5 << std::endl;
  std::cout << out->tot_events5 << std::endl;
  std::cout << out->sig_events5 << std::endl;
  std::cout << out->efficiency5 << std::endl;
  cout << "-" << endl;


  //  canv=NULL;

  if (canv!=NULL){

    RooPlot *varframe[2];

    if (single_gamma){
      canv->cd(bin+1);
      varframe[0]=roovar.frame();
      datahist->plotOn(varframe[0]);
      model->plotOn(varframe[0]);
      model->plotOn(varframe[0],Components("sigpdf"),LineStyle(kDashed),LineColor(kRed));
      model->plotOn(varframe[0],Components("bkgpdf"),LineStyle(kDashed),LineColor(kBlack));
      varframe[0]->Draw();
    }

    if (!single_gamma){

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

  }



  return out;

};


