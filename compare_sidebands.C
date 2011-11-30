using namespace RooFit;
using namespace std;



void compare_sidebands(){

  TFile *datafile_side = TFile::Open("dataEE_side.root");
  TFile *mcfile_side = TFile::Open("mcEE_side.root");
  TFile *mcfile_noside = TFile::Open("mcEE.root");


  RooDataHist *datahist_side, *mchist_noside_bkg, *mchist_side_all;

  datafile_side->GetObject("sieie_all",datahist_side);
  mcfile_noside->GetObject("sieie_bkg",mchist_noside_bkg);
  mcfile_side->GetObject("sieie_all",mchist_side_all);


  RooRealVar sieievar("sieievar","sieievar",0,0.05);
  RooPlot *sieievarframe = sieievar.frame(Name("sietaieta"),Title("sietaieta"));
  
  mchist_side_all->plotOn(sieievarframe,LineColor(kRed),DrawOption("L"),Rescale(1.0/mchist_side_all->sum(kTRUE)));
  datahist_side->plotOn(sieievarframe,LineColor(kBlue),DrawOption("L"),Rescale(1.0/datahist_side->sum(kTRUE)));
  mchist_noside_bkg->plotOn(sieievarframe,LineColor(kGreen),DrawOption("L"),Rescale(1.0/mchist_noside_bkg->sum(kTRUE)));
 
  sieievarframe->Draw();

}
