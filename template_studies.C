using namespace std;
using namespace RooFit;


void fit_dataset(const char* _varname, const char* inputfilename, TString splitting, float leftrange, float rightrange){

  TString varname(_varname);

  RooRealVar *roovar[2];
  RooDataHist *roohist[2][2]; // roohist[sig,bkg][x,y]
  RooDataSet *roodset;

  RooFormulaVar *cut[2];      

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

  int bin=0;

  {
    if (splitting=="EBEB") {
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,0).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,1).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("1"),bin).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("2"),bin).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("1"),bin).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("2"),bin).Data()),roohist[1][1]);
    }
    else if (splitting=="EBEE"){
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,0).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,1).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("1"),bin).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("2"),bin).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("1"),bin).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("2"),bin).Data()),roohist[1][1]);
    }
    else if (splitting=="EEEB"){
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,0).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,1).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("1"),bin).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("2"),bin).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("1"),bin).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("2"),bin).Data()),roohist[1][1]);
    }
    else if (splitting=="EEEE"){
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,0).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,1).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("1"),bin).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("2"),bin).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("1"),bin).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("2"),bin).Data()),roohist[1][1]);
    }

      cut[0] = new RooFormulaVar("cut_0","cut_0", Form("((%s>%f) && (%s<%f))",roovar[0]->getTitle().Data(),leftrange,roovar[0]->getTitle().Data(),rightrange), RooArgList(*(roovar[0])));
      cut[1] = new RooFormulaVar("cut_1","cut_1", Form("((%s>%f) && (%s<%f))",roovar[1]->getTitle().Data(),leftrange,roovar[1]->getTitle().Data(),rightrange), RooArgList(*(roovar[1])));
      std::cout << "Using cuts: " << std::endl;
      cut[0]->Print();
      cut[1]->Print();

    inputfile->GetObject(TString(data_dir).Append(helper->get_roodset_name(varname,splitting,bin).Data()),roodset);
    cout << "using dset " << TString(data_dir).Append(helper->get_roodset_name(varname,splitting,bin).Data()) << endl;

    RooFormulaVar *bothcut = new RooFormulaVar("bothcut","bothcut","cut_0 && cut_1",RooArgList(*(cut[0]),*(cut[1])));
    roodset=(RooDataSet*)(roodset->reduce(Cut(*bothcut)));
    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
	roohist[i][j]=(RooDataHist*)(roohist[i][j]->reduce(Cut(*bothcut)));

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

  roovar[0]->setBins(50);
  roovar[1]->setBins(50);




  RooHistPdf sig1pdf("sig1pdf","sig1pdf",*(roovar[0]),*(roohist[0][0]));
  RooHistPdf sig2pdf("sig2pdf","sig2pdf",*(roovar[1]),*(roohist[0][1]));
  RooHistPdf bkg1pdf("bkg1pdf","bkg1pdf",*(roovar[0]),*(roohist[1][0]));
  RooHistPdf bkg2pdf("bkg2pdf","bkg2pdf",*(roovar[1]),*(roohist[1][1]));

  RooProdPdf sigsigpdf("sigsigpdf","sigsigpdf",RooArgList(sig1pdf,sig2pdf));
  RooProdPdf sigbkgpdf("sigbkgpdf","sigbkgpdf",RooArgList(sig1pdf,bkg2pdf));
  RooProdPdf bkgsigpdf("bkgsigpdf","bkgsigpdf",RooArgList(bkg1pdf,sig2pdf));
  RooProdPdf bkgbkgpdf("bkgbkgpdf","bkgbkgpdf",RooArgList(bkg1pdf,bkg2pdf));


  RooRealVar rf1("rf1","rf1",1e-1,0,1);
  RooRealVar rf2("rf2","rf2",1e-1,0,1);
  RooRealVar rf3("rf3","rf3",1e-1,0,1);

  
  RooFormulaVar fsigsig("fsigsig","fsigsig","rf1",RooArgList(rf1));
  RooFormulaVar fsigbkg("fsigbkg","fsigbkg","(1-rf1)*rf2",RooArgList(rf1,rf2));
  RooFormulaVar fbkgsig("fbkgsig","fbkgsig","(1-rf1)*(1-rf2)*rf3",RooArgList(rf1,rf2,rf3));

  RooFormulaVar fbkgbkg("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg-fbkgsig",RooArgList(fsigsig,fsigbkg,fbkgsig));
  RooFormulaVar fsigbkg_noorder("fsigbkg_noorder","fsigbkg_noorder","fsigbkg+fbkgsig",RooArgList(fsigbkg,fbkgsig));

  RooAddPdf model("model","model",RooArgList(sigsigpdf,sigbkgpdf,bkgsigpdf,bkgbkgpdf),RooArgList(rf1,rf2,rf3),kTRUE);

  RooFitResult *fitres = model.fitTo(*roodset,Save());

  model.Print();

  cout << "sigsig: " << fsigsig.getVal() << " +/- " << fsigsig.getPropagatedError(*fitres) << endl;
  cout << "sigbkg_noorder: " << fsigbkg_noorder.getVal() << " +/- " << fsigbkg_noorder.getPropagatedError(*fitres) << endl;
  cout << "bkgbkg: " << fbkgbkg.getVal() << " +/- " << fbkgbkg.getPropagatedError(*fitres) << endl;
  cout << endl;
  cout << "bkgsig: " << fbkgsig.getVal() << " +/- " << fbkgsig.getPropagatedError(*fitres) << endl;
  cout << "sigbkg: " << fsigbkg.getVal() << " +/- " << fsigbkg.getPropagatedError(*fitres) << endl;




  RooPlot *varframe[2];
  varframe[0] = roovar[0]->frame("","");
  roodset->plotOn(varframe[0]);
  model.plotOn(varframe[0]);
  model.plotOn(varframe[0],Components("sigsigpdf"),LineStyle(kDashed),LineColor(kRed));
  model.plotOn(varframe[0],Components("sigbkgpdf"),LineStyle(kDashed),LineColor(kGreen));
  model.plotOn(varframe[0],Components("bkgsigpdf"),LineStyle(kDashed),LineColor(kBlue));
  model.plotOn(varframe[0],Components("bkgbkgpdf"),LineStyle(kDashed),LineColor(kBlack));
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


