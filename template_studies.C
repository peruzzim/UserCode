using namespace std;
using namespace RooFit;

void compare_sidebands(TString varname){

   TFile *datafile_side = TFile::Open("dataEB_side.root");
   TFile *mcfile_side = TFile::Open("mcEB_side.root");
   TFile *mcfile = TFile::Open("mcEB.root");

  // TFile *datafile_side = TFile::Open("dataEE.root");
  // TFile *mcfile_side = TFile::Open("dataEE_electron.root");
  // TFile *mcfile = TFile::Open("DYEE_electron.root");

  RooDataHist *datahist_side, *mchist_bkg, *mchist_side_all;

  datafile_side->GetObject("sieie_all",datahist_side);
  mcfile_side->GetObject("sieie_all",mchist_side_all);
  mcfile->GetObject("sieie_all",mchist_bkg);

  RooRealVar sieievar("sieievar","sieievar",0,0.05);
  RooPlot *sieievarframe = sieievar.frame(Name("sietaieta"),Title("sietaieta"));
  
  mchist_side_all->plotOn(sieievarframe,  LineColor(kRed),  MarkerColor(kRed),    DrawOption("P"),Rescale(1.0/mchist_side_all->sum(kTRUE)));
  mchist_bkg->plotOn(sieievarframe,LineColor(kGreen),MarkerColor(kGreen), DrawOption("P"),Rescale(1.0/mchist_bkg->sum(kTRUE)));
  datahist_side->plotOn(sieievarframe,    LineColor(kBlue), MarkerColor(kBlue),     DrawOption("P"),Rescale(1.0/datahist_side->sum(kTRUE)));
  
  std::cout << mchist_side_all->sum(kFALSE) << " " << datahist_side->sum(kFALSE) << " " << mchist_bkg->sum(kFALSE) << std::endl;
  
  TCanvas *c1 = new TCanvas();
  sieievarframe->Draw();
  c1->Update();
  c1->SaveAs("compare.png");
  
};

void fit_dataset(const char* datafilename, const char* sigfilename, const char* bkgfilename, TString varname, TString splitting, float leftrange, float rightrange){

  RooRealVar *roovar[2];
  RooDataHist *roohist[2][2]; // roohist[sig,bkg][x,y]
  RooDataSet *roodset;

  // hsig, hbkg 0:barrel 1:endcap

  TFile *datafile = TFile::Open(datafilename);
  TFile *sigfile = TFile::Open(sigfilename);
  TFile *bkgfile = TFile::Open(bkgfilename);

  {
    template_production helper;

    if (splitting=="EBEB") {
      datafile->GetObject(helper.get_roovar_name(varname,0,0).Data(),roovar[0]);
      datafile->GetObject(helper.get_roovar_name(varname,0,1).Data(),roovar[1]);
      sigfile->GetObject(helper.get_roohist_name(varname,TString("sig"),TString("EB"),TString("1")).Data(),roohist[0][0]);
      sigfile->GetObject(helper.get_roohist_name(varname,TString("sig"),TString("EB"),TString("2")).Data(),roohist[0][1]);
      bkgfile->GetObject(helper.get_roohist_name(varname,TString("bkg"),TString("EB"),TString("1")).Data(),roohist[1][0]);
      bkgfile->GetObject(helper.get_roohist_name(varname,TString("bkg"),TString("EB"),TString("2")).Data(),roohist[1][1]);
    }
    else if (splitting=="EBEE"){
      datafile->GetObject(helper.get_roovar_name(varname,0,0),roovar[0]);
      datafile->GetObject(helper.get_roovar_name(varname,1,1),roovar[1]);
      sigfile->GetObject(helper.get_roohist_name(varname,TString("sig"),TString("EB"),TString("1")).Data(),roohist[0][0]);
      sigfile->GetObject(helper.get_roohist_name(varname,TString("sig"),TString("EE"),TString("2")).Data(),roohist[0][1]);
      bkgfile->GetObject(helper.get_roohist_name(varname,TString("bkg"),TString("EB"),TString("1")).Data(),roohist[1][0]);
      bkgfile->GetObject(helper.get_roohist_name(varname,TString("bkg"),TString("EE"),TString("2")).Data(),roohist[1][1]);
    }
    else if (splitting=="EEEE"){
      datafile->GetObject(helper.get_roovar_name(varname,1,0),roovar[0]);
      datafile->GetObject(helper.get_roovar_name(varname,1,1),roovar[1]);
      sigfile->GetObject(helper.get_roohist_name(varname,TString("sig"),TString("EE"),TString("1")).Data(),roohist[0][0]);
      sigfile->GetObject(helper.get_roohist_name(varname,TString("sig"),TString("EE"),TString("2")).Data(),roohist[0][1]);
      bkgfile->GetObject(helper.get_roohist_name(varname,TString("bkg"),TString("EE"),TString("1")).Data(),roohist[1][0]);
      bkgfile->GetObject(helper.get_roohist_name(varname,TString("bkg"),TString("EE"),TString("2")).Data(),roohist[1][1]);
    }

    datafile->GetObject(helper.get_roodset_name(varname,splitting).Data(),roodset);
    cout << "using dset " << helper.get_roodset_name(varname,splitting).Data() << endl;

  }
 

  RooHistPdf sig1pdf("sig1pdf","sig1pdf",*(roovar[0]),*(roohist[0][0]));
  RooHistPdf sig2pdf("sig2pdf","sig2pdf",*(roovar[1]),*(roohist[0][1]));
  RooHistPdf bkg1pdf("bkg1pdf","bkg1pdf",*(roovar[0]),*(roohist[1][0]));
  RooHistPdf bkg2pdf("bkg2pdf","bkg2pdf",*(roovar[1]),*(roohist[1][1]));

  RooProdPdf sigsigpdf("sigsigpdf","sigsigpdf",RooArgList(sig1pdf,sig2pdf));
  RooProdPdf sigbkgpdf("sigbkgpdf","sigbkgpdf",RooArgList(sig1pdf,bkg2pdf));
  RooProdPdf bkgsigpdf("bkgsigpdf","bkgsigpdf",RooArgList(bkg1pdf,sig2pdf));
  RooProdPdf bkgbkgpdf("bkgbkgpdf","bkgbkgpdf",RooArgList(bkg1pdf,bkg2pdf));
  
  RooRealVar nsigsig("nsigsig","nsigsig",0,100000);
  RooRealVar nsigbkg("nsigbkg","nsigbkg",0,100000);
  RooRealVar nbkgsig("nbkgsig","nbkgsig",0,100000);
  RooRealVar nbkgbkg("nbkgbkg","nbkgbkg",0,100000);

  RooAddPdf model("model","model",RooArgList(sigsigpdf,sigbkgpdf,bkgsigpdf,bkgbkgpdf),RooArgList(nsigsig,nsigbkg,nbkgsig,nbkgbkg));

  RooFitResult *fitres = model.fitTo(*roodset);

  model.Print();

  RooPlot *varframe[2];

  for (int i=0; i<2; i++){
    varframe[i] = roovar[i]->frame("","");
    roodset->plotOn(varframe[i]);
    model.plotOn(varframe[i],LineColor(kBlue));
    if (i==0){
      sig1pdf.plotOn(varframe[i],LineColor(kRed),Normalization(nsigsig.getVal()+nsigbkg.getVal(),RooAbsReal::NumEvent));
      bkg1pdf.plotOn(varframe[i],LineColor(kGreen),Normalization(nbkgbkg.getVal()+nbkgsig.getVal(),RooAbsReal::NumEvent));
    }
    else if (i==1){
    sig2pdf.plotOn(varframe[i],LineColor(kRed),Normalization(nsigsig.getVal()+nbkgsig.getVal(),RooAbsReal::NumEvent));
    bkg2pdf.plotOn(varframe[i],LineColor(kGreen),Normalization(nbkgbkg.getVal()+nsigbkg.getVal(),RooAbsReal::NumEvent));
    }
    
  }


  Float_t ntot=nsigsig.getVal()+nsigbkg.getVal()+nbkgsig.getVal()+nbkgbkg.getVal();
  cout << "sigsig: " << nsigsig.getVal()/ntot << endl;
  cout << "bkgsig: " << nbkgsig.getVal()/ntot << endl;
  cout << "sigbkg: " << nsigbkg.getVal()/ntot << endl;
  cout << "bkgbkg: " << nbkgbkg.getVal()/ntot << endl;


  TCanvas *c1[2];
  for (int i=0; i<2; i++){
    c1[i] = new TCanvas();
    c1[i]->cd();
    varframe[i]->Draw();
    c1[i]->SaveAs(Form("fit_%d.png",i));
    c1[i]->Close();
  }



};


