
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

void doBiasStudy(TString pref1, TString pref2, TString dset1, TString dset2, TString temp1, TString temp2, TString reg="EB"){

  std::vector<TString> names[2];

  names[0]=gethistoname(pref1,dset1,temp1,reg);
  names[1]=gethistoname(pref2,dset2,temp2,reg);

  TFile *f1 = new TFile(names[0][0].Data(),"read");
  TFile *f2 = new TFile(names[1][0].Data(),"read");

  TH1F *hs;
  TH1F *hb;

  f1->GetObject(names[0][1].Data(),hs);
  assert (hs!=NULL);

  f2->GetObject(names[1][1].Data(),hb);
  assert (hb!=NULL);

  RooRealVar *photoniso = new RooRealVar("photoniso","photoniso",0,5);
  photoniso->setBins(40);

  RooDataHist datahs("datahs","datahs",*photoniso,Import(*hs));
  RooDataHist datahb("datahb","datahb",*photoniso,Import(*hb));
  RooHistPdf pdfhs("pdfhs","pdfhs",*photoniso,datahs);
  RooHistPdf pdfhb("pdfhb","pdfhb",*photoniso,datahb);


  RooRealVar sigfrac("sigfrac","sigfrac",0.5,0,1);
  RooRealVar nevents("nevents","nevents",1,0,1e+6);
  RooFormulaVar psig("psig","nevents*sigfrac",RooArgList(sigfrac,nevents));
  RooFormulaVar pbkg("pbkg","nevents*(1-sigfrac)",RooArgList(sigfrac,nevents));

  RooAddPdf genmodel("genmodel","genmodel",RooArgList(pdfhs,pdfhb),RooArgList(psig,pbkg));
  RooAddPdf fitmodel(genmodel,"fitmodel");

  TCanvas *c0 = new TCanvas();
  c0->cd();
  RooPlot *f0 = photoniso->frame(Title("photonisoframe"));
  genmodel.plotOn(f0,DataError(RooAbsData::SumW2));
  genmodel.plotOn(f0,Components(pdfhb),LineColor(kRed),DataError(RooAbsData::SumW2));
  genmodel.plotOn(f0,Components(pdfhs),LineColor(kGreen),DataError(RooAbsData::SumW2));
  f0->Draw();

  RooMCStudy *mcstudy = new RooMCStudy(genmodel,*photoniso,FitModel(fitmodel),Silence(),Extended(),FitOptions(Save(kTRUE),PrintEvalErrors(0)));

  mcstudy->generateAndFit(1000,1000);

  RooPlot* frame1 = mcstudy->plotParam(sigfrac,Bins(50));
  RooPlot* frame2 = mcstudy->plotError(sigfrac,Bins(50));
  RooPlot* frame3 = mcstudy->plotPull(sigfrac,Bins(50),FitGauss(1));

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","c",900,900);
  c->Divide(3);

  c->cd(1);
  frame1->Draw();
  c->cd(2);
  frame2->Draw();
  c->cd(3);
  frame3->Draw();

  gDirectory->Add(mcstudy);
  gDirectory->Add(photoniso);

};
