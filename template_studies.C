#include <assert.h>

using namespace std;
using namespace RooFit;

void run_fits(){
  TString varname("PhoIso04");
  TString inputfilename("out_NEW.root");

  const int n_bins=3;

  RooFitResult *fr[n_bins];

  TGraphErrors *out = new TGraphErrors(n_bins);
  out->SetMarkerStyle(20);
  out->SetMarkerColor(kRed);
  out->SetLineColor(kRed);
  out->SetLineWidth(2);

  TCanvas *fits_canv = new TCanvas("fits","fits");
  fits_canv->Divide(2,3);

  RooRealVar *rf1;
  RooRealVar *rf2;
  RooRealVar *rf3;

  for (int bin=0; bin<n_bins; bin++) {
    fr[bin]=fit_dataset(varname.Data(),inputfilename.Data(),"EBEB",-2,5,bin,fits_canv);
    rf1=(RooRealVar*)(fr[bin]->floatParsFinal().find("rf1"));
    rf2=(RooRealVar*)(fr[bin]->floatParsFinal().find("rf2"));
    rf3=(RooRealVar*)(fr[bin]->floatParsFinal().find("rf3"));
    RooFormulaVar fsigsig("fsigsig","fsigsig","rf1",RooArgList(*rf1));
    RooFormulaVar fsigbkg("fsigbkg","fsigbkg","(1-rf1)*rf2",RooArgList(*rf1,*rf2));
    RooFormulaVar fbkgsig("fbkgsig","fbkgsig","(1-rf1)*(1-rf2)*rf3",RooArgList(*rf1,*rf2,*rf3));
    RooFormulaVar fbkgbkg("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg-fbkgsig",RooArgList(fsigsig,fsigbkg,fbkgsig));
    RooFormulaVar fsigbkg_noorder("fsigbkg_noorder","fsigbkg_noorder","fsigbkg+fbkgsig",RooArgList(fsigbkg,fbkgsig));
    out->SetPoint(bin,bin,fsigsig.getVal());
    out->SetPointError(bin,0,fsigsig.getPropagatedError(*(fr[bin])));
  }

  fits_canv->Update();

  TCanvas *output_canv = new TCanvas("output_canv","output_canv");
  output_canv->cd();
  out->Draw("AP");
  output_canv->Update();

		

};


RooFitResult* fit_dataset(const char* _varname, const char* inputfilename, TString splitting, float leftrange, float rightrange, int bin, TCanvas *canv=NULL){

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

  // roovar[EB,EE][1,2][no,ord]
  // roohist[sig,bkg,all][EB,EE][1,2][temp][no,ord]
  // templatehist[sig,bkg,all][EB,EE][temp]
  // roodataset[EBEB,EBEE,EEEB,EEEE][temp][no,ord]

  TFile *inputfile = TFile::Open(inputfilename);

  TString data_dir("data_Tree_standard_sel/");
  TString mc_dir("mc_Tree_standard_sel/");

  template_production *helper = new template_production(NULL);

  TString ord="noorder";

  {
    if (splitting=="EBEB") {
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,0,ord).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,1,ord).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("1"),bin,ord).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("2"),bin,ord).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("1"),bin,ord).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("2"),bin,ord).Data()),roohist[1][1]);
    }
    else if (splitting=="EBEE"){
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,0,ord).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,1,ord).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("1"),bin,ord).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("2"),bin,ord).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("1"),bin,ord).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("2"),bin,ord).Data()),roohist[1][1]);
    }
    else if (splitting=="EEEB"){
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,0,ord).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,1,ord).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("1"),bin,ord).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("2"),bin,ord).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("1"),bin,ord).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("2"),bin,ord).Data()),roohist[1][1]);
    }
    else if (splitting=="EEEE"){
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,0,ord).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,1,ord).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("1"),bin,ord).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("2"),bin,ord).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("1"),bin,ord).Data()),roohist[1][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("2"),bin,ord).Data()),roohist[1][1]);
    }

    for (int i=0; i<2; i++){
      assert (roovar[i]!=NULL);
          for (int j=0; j<2; j++){
	    assert (roohist[i][j]!=NULL);
	  }
    }

      cut[0] = new RooFormulaVar("cut_0","cut_0", Form("((%s>%f) && (%s<%f))",roovar[0]->getTitle().Data(),leftrange,roovar[0]->getTitle().Data(),rightrange), RooArgList(*(roovar[0])));
      cut[1] = new RooFormulaVar("cut_1","cut_1", Form("((%s>%f) && (%s<%f))",roovar[1]->getTitle().Data(),leftrange,roovar[1]->getTitle().Data(),rightrange), RooArgList(*(roovar[1])));
      std::cout << "Using cuts: " << std::endl;
      cut[0]->Print();
      cut[1]->Print();

      inputfile->GetObject(TString(data_dir).Append(helper->get_roodset_name(varname,splitting,bin,ord).Data()),roodset);
      cout << "using dset " << TString(data_dir).Append(helper->get_roodset_name(varname,splitting,bin,ord).Data()) << endl;

    RooFormulaVar *bothcut = new RooFormulaVar("bothcut","bothcut","cut_0 && cut_1",RooArgList(*(cut[0]),*(cut[1])));
    roodset=(RooDataSet*)(roodset->reduce(Cut(*bothcut)));
    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
	roohist[i][j]=(RooDataHist*)(roohist[i][j]->reduce(Cut(*(cut[j]))));

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


    for (int i=0; i<2; i++){
      roovar[i]->Print();
      for (int j=0; j<2; j++){
	roohist[i][j]->Print();
      }
    }
    roodset->Print();



  }




  roovar[0]->setRange(leftrange,rightrange);
  roovar[1]->setRange(leftrange,rightrange);

//  const int n_variable_bins=25;
//  roovar[0]->setBins(n_variable_bins);
//  roovar[1]->setBins(n_variable_bins);

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


  if (canv!=NULL){
    RooPlot *varframe[2];
    for (int i=0; i<2; i++){
      varframe[i] = roovar[i]->frame("","");
      roodset->plotOn(varframe[i]);
      model.plotOn(varframe[i]);
      model.plotOn(varframe[i],Components("sigsigpdf"),LineStyle(kDashed),LineColor(kRed));
      model.plotOn(varframe[i],Components("sigbkgpdf"),LineStyle(kDashed),LineColor(kGreen));
      model.plotOn(varframe[i],Components("bkgsigpdf"),LineStyle(kDashed),LineColor(kBlue));
      model.plotOn(varframe[i],Components("bkgbkgpdf"),LineStyle(kDashed),LineColor(kBlack));
      canv->cd(2*bin+i+1);
      varframe[i]->Draw();
    }
  }

  return fitres;
  
};


