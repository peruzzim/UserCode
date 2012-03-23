using namespace std;
using namespace RooFit;


void fit_dataset(const char* _varname, const char* inputfilename, TString splitting, float leftrange, float rightrange){

  TString varname(_varname);

  RooRealVar *roovar[2];
  RooDataHist *roohist[2][2]; // roohist[sig,bkg][x,y]
  RooDataSet *roodset;

  roodset=NULL;
  for (int i=0; i<2; i++){
    roovar[i]=NULL;
    for (int j=0; j<2; j++){
      roohist[i][j]=NULL;
    }
  }

  // roovar[EB,EE][1,2]
  // roohist[sig,bkg,all][EB,EE][1,2]
  // templatehist[sig,bkg,all][EB,EE]
  // roodataset[EBEB,EBEE,EEEB,EEEE]

  TFile *inputfile = TFile::Open(inputfilename);

  TString data_dir("data_Tree_standard_sel/");
  TString mc_dir("mc_Tree_standard_sel/");

  template_production *helper = new template_production(NULL);

  {


    if (splitting=="EBEB") {
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,0).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,1).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("1")).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("2")).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("1")).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("2")).Data()),roohist[1][1]);
    }
    else if (splitting=="EBEE"){
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,0).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,1).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("1")).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("2")).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("1")).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("2")).Data()),roohist[1][1]);
    }
    else if (splitting=="EEEB"){
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,0).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,1).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("1")).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("2")).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("1")).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("2")).Data()),roohist[1][1]);
    }
    else if (splitting=="EEEE"){
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,0).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,1).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("1")).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("2")).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("1")).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("2")).Data()),roohist[1][1]);
    }

    inputfile->GetObject(TString(data_dir).Append(helper->get_roodset_name(varname,splitting).Data()),roodset);
    cout << "using dset " << TString(data_dir).Append(helper->get_roodset_name(varname,splitting).Data()) << endl;

    bool wrong=false;

    if (roodset==NULL) wrong=true;
    for (int i=0; i<2; i++){
      if (roovar[i]==NULL) wrong=true;
      for (int j=0; j<2; j++){
	if (roohist[i][j]==NULL) wrong=true;
      }
    }
    
    if (wrong==true) {
      std::cout << "something wrong with initialization; exiting." << std::endl;
      return;
    }

  }


  RooHistPdf sig1pdf("sig1pdf","sig1pdf",*(roovar[0]),*(roohist[0][0]));
  RooHistPdf sig2pdf("sig2pdf","sig2pdf",*(roovar[1]),*(roohist[0][1]));
  RooHistPdf bkg1pdf("bkg1pdf","bkg1pdf",*(roovar[0]),*(roohist[1][0]));
  RooHistPdf bkg2pdf("bkg2pdf","bkg2pdf",*(roovar[1]),*(roohist[1][1]));

  RooProdPdf sigsigpdf("sigsigpdf","sigsigpdf",RooArgList(sig1pdf,sig2pdf));
  RooProdPdf sigbkgpdf("sigbkgpdf","sigbkgpdf",RooArgList(sig1pdf,bkg2pdf));
  RooProdPdf bkgsigpdf("bkgsigpdf","bkgsigpdf",RooArgList(bkg1pdf,sig2pdf));
  RooProdPdf bkgbkgpdf("bkgbkgpdf","bkgbkgpdf",RooArgList(bkg1pdf,bkg2pdf));
  
  RooRealVar nsigsig("nsigsig","nsigsig",100,0,100000);
  RooRealVar nsigbkg("nsigbkg","nsigbkg",100,0,100000);
  RooRealVar nbkgsig("nbkgsig","nbkgsig",100,0,100000);
  RooRealVar nbkgbkg("nbkgbkg","nbkgbkg",100,0,100000);

  RooAddPdf model("model","model",RooArgList(sigsigpdf,sigbkgpdf,bkgsigpdf,bkgbkgpdf),RooArgList(nsigsig,nsigbkg,nbkgsig,nbkgbkg));

  RooFitResult *fitres = model.fitTo(*roodset,Range(leftrange,rightrange,kFALSE));

  model.Print();

  Float_t ntot=nsigsig.getVal()+nsigbkg.getVal()+nbkgsig.getVal()+nbkgbkg.getVal();
  cout << "sigsig: " << nsigsig.getVal()/ntot << endl;
  cout << "bkgsig: " << nbkgsig.getVal()/ntot << endl;
  cout << "sigbkg: " << nsigbkg.getVal()/ntot << endl;
  cout << "bkgbkg: " << nbkgbkg.getVal()/ntot << endl;

  RooPlot *varframe[2];
  varframe[0] = roovar[0]->frame("","");
  sigsigpdf.plotOn(varframe[0]);
  varframe[0]->Draw();


  return;




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


  TCanvas *c1[2];
  for (int i=0; i<2; i++){
    c1[i] = new TCanvas();
    c1[i]->cd();
    varframe[i]->Draw();
    c1[i]->SaveAs(Form("fit_%d.png",i));
    c1[i]->Close();
  }



};


