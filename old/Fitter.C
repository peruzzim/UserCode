#include <iostream>
#include <stdio.h>
#include <stdlib.h>


using namespace RooFit;
using namespace std;

struct fitresult {
    RooPlot* plot;
    RooFitResult* res;
};


void ApplyTDRStyle(){
	
	
	// For the canvas:
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetCanvasDefH(600); //Height of canvas
	gStyle->SetCanvasDefW(800); //Width of canvas
	gStyle->SetCanvasDefX(0);   //POsition on screen
	gStyle->SetCanvasDefY(0);
	
	// For the Pad:
	gStyle->SetPadBorderMode(0);
	// gStyle->SetPadBorderSize(Width_t size = 1);
	gStyle->SetPadColor(kWhite);
	gStyle->SetPadGridX(false);
	gStyle->SetPadGridY(false);
	gStyle->SetGridColor(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	
	// For the frame:
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(1);
	gStyle->SetFrameFillColor(0);
	gStyle->SetFrameFillStyle(0);
	gStyle->SetFrameLineColor(1);
	gStyle->SetFrameLineStyle(1);
	gStyle->SetFrameLineWidth(1);
	
	gStyle->SetHistLineColor(1);
	gStyle->SetHistLineStyle(0);
	gStyle->SetHistLineWidth(1);
	
	gStyle->SetMarkerStyle(20);
	
	//For the fit/function:
	gStyle->SetOptFit(1);
	gStyle->SetFitFormat("5.4g");
	gStyle->SetFuncColor(2);
	gStyle->SetFuncStyle(1);
	gStyle->SetFuncWidth(1);
	
	//gStyle->SetOptFile(0);
	gStyle->SetOptStat(1); // To display the mean and RMS:   SetOptStat("mr");
	gStyle->SetStatColor(kWhite);
	gStyle->SetStatFont(42);
	gStyle->SetStatFontSize(0.02);
	gStyle->SetStatTextColor(1);
	gStyle->SetStatFormat("6.4g");
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.15);
	
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetPadLeftMargin(0.13);
	gStyle->SetPadRightMargin(0.05);
	
	gStyle->SetOptTitle(1);
	gStyle->SetTitleFont(42);
	gStyle->SetTitleColor(1);
	gStyle->SetTitleTextColor(1);
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleFontSize(0.05);
	
	
	
	gStyle->SetTitleColor(1, "XYZ");
	gStyle->SetTitleFont(42, "XYZ");
	gStyle->SetTitleSize(0.06, "XYZ");
	// gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
	// gStyle->SetTitleYSize(Float_t size = 0.02);
	gStyle->SetTitleXOffset(0.9);
	gStyle->SetTitleYOffset(1.05);
	// gStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
	
	// For the axis labels:
	
	gStyle->SetLabelColor(1, "XYZ");
	gStyle->SetLabelFont(42, "XYZ");
	gStyle->SetLabelOffset(0.007, "XYZ");
	gStyle->SetLabelSize(0.05, "XYZ");
	
	// For the axis:
	
	gStyle->SetAxisColor(1, "XYZ");
	gStyle->SetStripDecimals(kTRUE);
	//gStyle->SetTickLength(0.03, "XYZ");
	//gStyle->SetNdivisions(510, "XYZ");
	gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
	gStyle->SetPadTickY(1);
	
	
	return;	
}


void ZToEEFit(const char* datamc, int code, int cat, float fixedA=-1 , float fixedN=-1, fitresult* fit, RooPlot *outplot=NULL, const char* plotOpt = "NEU", const char* TitleMass="title_foo") {
	
    bool isData;
    
    if (datamc=="data") isData=true;
    else if (datamc=="mc") isData=false;
    
  Float_t massMin(70), massMax(110);
  //    Float_t massMin(80), massMax(100);
	Int_t nBins(50);
	
	
	RooRealVar  mass("mass","M(e^{+}e^{-})", massMin, massMax,"GeV/c^{2}");
	mass.setBins(10000) ;
	
	//  Parameters for Bifurcated Gaussian
	RooRealVar  mu("m_{0}", "Bias", 0,-5.0,5.0,"GeV/c^{2}"); 
	RooRealVar  sigL("#sigma_ {L}","Left Width", 4.0,0.1,10.0,"GeV/c^{2}"); 
	RooRealVar  sigR("#sigma_{R}","Right Width", 2.0,0.1,10.0,"GeV/c^{2}"); 
	
	// Gaussian
	RooRealVar  sig("#sigma_{G}","Core Width", 1.0,0.5,2.0,"GeV/c^{2}"); 
	
	//  Parameters for Crystal Ball Lineshape 
	RooRealVar  m0("#Delta m_{0}", "Bias", 0.0, -5.0, 8.0);//,"GeV/c^{2}"); 
	RooRealVar  sigma("#sigma_{CB}","Width", 1.8,0.0,5.0);//,"GeV/c^{2}"); 
	//	RooRealVar  sigma("#sigma_{CB}","Width", 1.8,0.5,5.0);//,"GeV/c^{2}"); 
	//  RooRealVar  cut("#alpha","Cut", 0.959,0.8,4.0); 
	//	RooRealVar  cut("#alpha","Cut", 1.0,0.8,4.0); // should be 0.8, 4.0
	RooRealVar  cut("#alpha","Cut", 1.298,0.2,4.0); // should be 0.8, 4.0
	RooRealVar  power("n","Power", 3.92, 0.5, 25.0); 
	//	cut.setConstant();
	if (fixedA!=-1) {
        cut.setVal(fixedA);
        cut.setConstant();
    }
    if (fixedN!=-1) {
        power.setVal(fixedN);
        power.setConstant();
    }
    
	//  Parameters for Breit-Wigner Distribution
	RooRealVar  mRes("M_{Z^{0}}", "Z0 Resonance  Mass", 91.188, 85.0, 95.0);//,"GeV/c^{2}"); 
	RooRealVar  Gamma("#Gamma", "#Gamma", 2.4952, 2.0,3.0);//,"GeV/c^{2}"); 
	mRes.setConstant();
	Gamma.setConstant();
	
	RooRealVar  bgtau("a_{BG}", "Backgroung Shape", -0.15, -1.0, 0.0, "1/GeV/c^{2}");
//	RooRealVar  frac("frac", "Signal Fraction", 1.0,0.0,1.0);
	
//	if (){
		RooRealVar  nsig("N_{S}", "#signal events", 5000,0,1000000.);
		RooRealVar  nbkg("N_{B}", "#background events", 420,0,50000.);
//	}
//	else{
//		RooRealVar  nsig("N_{S}", "#signal events", 50000,0.,1000000.);
//		RooRealVar  nbkg("N_{B}", "#background events", 1000,0,500000);
//	}
//	frac.setConstant();
	//  bgtau.setConstant();
//	RooRealVar  a0("a_{0}", "Backgroung coeff #0", 100, 1, 1000);
//	RooRealVar  a1("a_{1}", "Backgroung coeff #1", -0.1, -1, -0.05);
	
	//  Introduce a resolution model
	//  RooBifurGauss  res("res", "A Bifurcated Gaussian Distribution", deltam, mu,sigL,sigR);
	//RooGaussian    resG("resG",   "A  Gaussian Lineshape",     mass, m0,sig);
	RooCBShape     resCB("resCB", "A  Crystal Ball Lineshape", mass, m0,sigma, cut, power);
	//RooRealVar     fracG("f_{G}",  "Gaussian Fraction",        0.0,0.0,1.0);
	//RooAddPdf      res("res",     "Resolution Model",          resG, resCB, fracG); 
	//fracG.setConstant();
	
	//  Breit-Wigner Lineshape 
	RooBreitWigner bw("bw","A Breit-Wigner Distribution",mass,mRes,Gamma);
	
	
	//  Convolution p.d.f. using numeric convolution operator
	//  RooNumConvPdf bw_res("bw_res","Convolution", mass, bw, res);
	//  bw_res.setConvolutionWindow(m0,sigma,10);
	//  Convolution p.d.f. using numeric convolution operator based on Fourier Transforms
	//RooFFTConvPdf bw_res("bw_res","Convolution", mass, bw, res);
	RooFFTConvPdf bw_res("bw_res","Convolution", mass, bw, resCB);
	
	// Background  p.d.f.
	RooExponential bg("bg", "Backgroung Distribution", mass, bgtau);
	// RooPolynomial bg("bg", "Backgroung Distribution", mass, RooArgList(a0,a1),0);
	// RooAbsReal *powlaw = bindFunction(TMath::Power, mass, a0);
	

	// Fit Model
	//  RooAddPdf      model("model", "Signal + Background", bw_res, bg, frac);
	RooAddPdf      model("model", "Di-photon mass model", RooArgList(bw_res, bg), RooArgList(nsig, nbkg));
	
	// Read data set 


	//	RooDataSet *data = RooDataSet::read(filename, mass);

	RooDataSet *data = new RooDataSet("dset","dset",mass);
	{

	  ifstream in;

        TString filename("output_");
        filename.Append(datamc);
        filename.Append("_no2D2D_nophietacracksexcl/masspoints_ee_");
        filename+=code;
        filename.Append("_cat");
        filename+=cat;
        filename.Append(".txt");
        
	  in.open(filename.Data());
	  int n;
	  n=0;
	  //        while(n<30000){
	  //        while(n<98069){
	  //	  while(n<53929){
	    	  while(n<9755){
	    float x;
	    float w;
	    in >> x;
	    in >> w;
	    //	    std::cout << x << " " << w << std::endl;
	    if (!in.good()) break;
	    if (massMin<x && x<massMax){
	      mass.setVal(x);
	      data->add(mass,w);
	    }
	    n++;
	  }
	  in.close();
	  
	}
	
	TStopwatch t ;
	t.Start() ;

	//RooFitResult *fitres=bw_res.fitTo(*data,Save(),Optimize(0),Timer(1));
	(*fit).res=model.fitTo(*data,Save(),Optimize(0),Timer(1));

    (*fit).plot=mass.frame(Title(TitleMass),Range(massMin,massMax),Bins(nBins));
    
    RooPlot* plot[2];
    plot[0]=(*fit).plot;
    plot[1]=outplot;

    for (int i=0;i<2;i++){
        if (plot[i]==NULL) continue;
        Int_t color = (isData) ? kBlue : kRed;
	data->plotOn(plot[i],MarkerColor(color));
	model.plotOn(plot[i],LineColor(color));
	if (isData) model.paramOn(plot[i], Format(plotOpt,AutoPrecision(2)), Parameters(RooArgSet(m0, sigma, cut, power)), Layout(.61, 0.92, 0.92), ShowConstants(kTRUE) );
	else model.paramOn(plot[i], Format(plotOpt,AutoPrecision(2)), Parameters(RooArgSet(m0, sigma, cut, power)), Layout(.15, 0.46, 0.92), ShowConstants(kTRUE) );
	}
    
    
	model.Print();
	

}


void Fitter(){
  
  ApplyTDRStyle();

  const int ncats=1;
  const int cats[ncats]={6};
  const int ncodes=2;
  const int codes[ncodes]={15,20};

    
    
    fitresult fitresmc[ncodes][ncats];
    fitresult fitresdata[ncodes][ncats];

    
    for (int i=0; i<ncodes; i++){
        for (int j=0; j<ncats; j++){
            fitresmc[i][j].plot=NULL;
            fitresmc[i][j].res=NULL;
            fitresdata[i][j].plot=NULL;
            fitresdata[i][j].res=NULL;
        }   
    }
    
    

  TCanvas *canv[ncodes][ncats];

  for (int i=0; i<ncodes; i++) for (int j=0; j<ncats; j++){

    TString title("Zee ");
    if (cats[j]==1) title.Append("EB-EB ");
    if (cats[j]==2) title.Append("EB-EE ");
      if (cats[j]==3) title.Append("EE-EE ");
      if (cats[j]==4) title.Append("central-central EB ");
      if (cats[j]==5) title.Append("central-outer EB ");
      if (cats[j]==6) title.Append("outer-outer EB ");
    if (codes[i]==15) title.Append("new scheme");
    if (codes[i]==20) title.Append("old scheme");

    ZToEEFit("mc",codes[i],cats[j],-1,-1,&(fitresmc[i][j]),NULL,"NEU",title);
      RooFitResult *res = fitresmc[i][j].res;
      float alpha=((RooRealVar*)res->floatParsFinal().find("#alpha"))->getVal();
      float cbn=((RooRealVar*)res->floatParsFinal().find("n"))->getVal();
    ZToEEFit("data",codes[i],cats[j],alpha,cbn,&(fitresdata[i][j]),fitresmc[i][j].plot,"NEU",title);

    canv[i][j]=new TCanvas();
    canv[i][j]->cd();
    fitresmc[i][j].plot->Draw();

    TPaveText *pl = new TPaveText(0.15, 0.52 ,0.30, 0.62,"BRNDC");
    pl->SetFillColor(0);
    pl->SetTextColor(2);
    pl->AddText("MC");
    pl->Draw();
    
    TPaveText *pl2 = new TPaveText(0.77, 0.52 ,0.92, 0.62,"BRNDC");
    pl2->SetFillColor(0);
    pl2->SetTextColor(4);
    pl2->AddText("DATA");
    pl2->Draw();

    canv[i][j]->Update();
    //    canv[i][j]->Draw();


    title.Append(".png");
    canv[i][j]->SaveAs(title.Data());

  }

//  TFile *outfile = new TFile("fit_results_mcdata.root","recreate");
//  for (int i=0; i<ncodes; i++) for (int j=0; j<ncats; j++){
//    TString mc("mc_");
//    TString data("data_");
//    mc+=codes[i];
//    data+=codes[i];
//    mc.Append("_cat");
//    data.Append("_cat");
//    mc+=cats[j];
//    data+=cats[j];
//    outfile->WriteObject(fitresmc[i][j],mc.Data());
//    outfile->WriteObject(fitresdata[i][j],data.Data());
//  }
//  outfile->Close();


  /*
  for (int i=0; i<ncodes; i++) for (int j=0; j<ncats; j++){
    std::cout << "FIT MC " << codes[i] << " cat " << cats[j] << std::endl;
    fitresmc[i][j]->Print();
    std::cout << "FIT DATA " << codes[i] << " cat " << cats[j] << std::endl;
    fitresdata[i][j]->Print();
  }
  */


	

}





