using namespace RooFit;

void fit_templates(){

gROOT->ProcessLine(".L get_sieie_template.C+O");

  RooRealVar sieievar("sieievar","sieievar",0,0.05);
  RooPlot *sieievarframe = sieievar.frame(Name("sietaieta"),Title("sietaieta"));

  TFile *datafile = TFile::Open("data_inclusive.root");
  TFile *mcfile = TFile::Open("mc_inclusive.root");

  TTree *treedata;
  RooDataHist *hdata;
  datafile->GetObject("Tree",treedata);
  {
    get_sieie_template helperdata(treedata,true);
    helperdata.Loop();
    hdata = new RooDataHist(*(helperdata.hsieie[0]));
  }

  hdata->plotOn(sieievarframe);
  sieievarframe->Draw();

  TTree *treemc;
  RooDataHist *hmcbkg, *hmcsig;
  mcfile->GetObject("Tree",treemc);
  {
    get_sieie_template helpermc(treemc,false);
    helpermc.Loop();
    hmcbkg = new RooDataHist(*(helpermc.hsieie[0]));
    hmcsig = new RooDataHist(*(helpermc.hsieie[1]));
  }



  RooRealVar nsig("nsig","nsig",100,0,1e+6);
  RooRealVar nbkg("nbkg","nbkg",100,0,1e+6);
  
  RooHistPdf pdfmcsig("pdfmcsig","pdfmcsig",sieievar,*hmcsig);
  RooHistPdf pdfmcbkg("pdfmcbkg","pdfmcbkg",sieievar,*hmcbkg);

  RooAddPdf pdfmc("pdfmc","pdfmc",RooArgList(pdfmcsig,pdfmcbkg),RooArgList(nsig,nbkg));

  RooFitResult *fitres = pdfmc.fitTo(*hdata,Save());

  fitres->Print();



  pdfmc.plotOn(sieievarframe,LineColor(kBlue));
  pdfmc.plotOn(sieievarframe,Components(pdfmcsig),LineColor(kRed),LineStyle(kDashed));
  pdfmc.plotOn(sieievarframe,Components(pdfmcbkg),LineColor(kGreen),LineStyle(kDashed));
  sieievarframe->Draw();

}
