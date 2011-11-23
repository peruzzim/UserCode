#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;
using namespace RooFit;
using namespace TMath;

RooFitResult* fitter_Zee(TH1D *hist){

  RooRealVar  Zmassvar("Zmassvar","Zmassvar", 82, 100);
  RooDataHist *datahist = new RooDataHist("data","Z Mass",Zmassvar,hist);
  RooPlot *Zmassvarframe = Zmassvar.frame(Name("Zmassvarframe"),Title(hist->GetTitle())) ;
  datahist->plotOn(Zmassvarframe);
  
  RooRealVar alpha  ("alpha"  ,        "alpha" , 0.005,0.001,0.1); 
  RooRealVar n      ("n"      ,            "n" , 1,0.001,10); 
  RooRealVar cbmean ("cbmean" ,       "cbmean" , 1, 0.8, 1.2);
  RooRealVar cbsigma("cbsigma",      "cbsigma" , 0.01, 0.001, 0.2);

  RooRealVar bwmean("bwmean","bwmean",91,85,95);
  RooRealVar bwsigma("bwsigma","bwsigma",3,2,4);
  RooRealVar expoconst("expoconst","expoconst",-0.1,-0.5,0);

  RooCBShape cball  ("cball"  , "crystal ball" , Zmassvar, cbmean, cbsigma, alpha, n);
  RooBreitWigner bw("bw","breit wigner",Zmassvar,bwmean,bwsigma);

  RooFFTConvPdf cballXbw("cballXbw","cball (X) bw",Zmassvar,bw,cball);
  RooExponential expo("expo", "exponential", Zmassvar, expoconst);

  RooRealVar frac("frac","frac",0.1,0.001,0.2);

  RooAddPdf Zshapemodel("Zshapemodel","expo + cball (X) bw",RooArgList(expo,cballXbw),frac);


  RooFitResult *fitres =Zshapemodel.fitTo(*datahist,Range(82,100),Save());   
  Zshapemodel.plotOn(Zmassvarframe,LineColor(kBlue));

  Zmassvarframe->Draw();

  fitres->Print();

  return fitres;

};

void fitmacro_Zee(){

  const int NCorrStyles = 2;
  int corrstyles[NCorrStyles]={0,15};

  TH1D *hist[NCorrStyles];
  TCanvas *canv[NCorrStyles];

  RooFitResult *fitres[NCorrStyles];

  for (int i=0; i<NCorrStyles;i++){
    TString histname("eeInvMass");
    histname+=corrstyles[i];
    std::cout << histname.Data(); std::endl;
    gDirectory->GetObject(histname,hist[i]);
  }
  
  for (int i=0; i<NCorrStyles;i++){
    canv[i] = new TCanvas();
    fitres[i]=fitter_Zee(hist[i]); 
  }

 for (int i=0; i<NCorrStyles;i++){
   fitres[i]->Print();
 }


};
