#include <assert.h>

using namespace std;
using namespace RooFit;

typedef struct {
  RooFitResult *fr;
  float purity;
  float purity5;
} fit_output; 

  bool study_templates=false;

void run_fits(TString varname="PhoIso04", TString inputfilename="out_NEW.root", TString splitting, bool single_gamma=false, bool do_order=false, float leftrange=-5, float rightrange=25){

  const int n_bins=3;


  fit_output *fr[n_bins];

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
    fr[bin]=fit_dataset(varname.Data(),inputfilename.Data(),splitting.Data(),leftrange,rightrange,bin,fits_canv,single_gamma,do_order);

    if (!single_gamma){  
      rf1=(RooRealVar*)(fr[bin]->fr->floatParsFinal().find("rf1"));
      rf2=(RooRealVar*)(fr[bin]->fr->floatParsFinal().find("rf2"));
      RooFormulaVar fsigsig("fsigsig","fsigsig","rf1",RooArgList(*rf1));
      RooFormulaVar fsigbkg("fsigbkg","fsigbkg","(1-rf1)*rf2",RooArgList(*rf1,*rf2));
      RooFormulaVar fbkgbkg("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg",RooArgList(fsigsig,fsigbkg));
      std::cout << bin << " " << fsigsig.getVal() << " " << fsigbkg.getVal() << " " << fbkgbkg.getVal() << std::endl;
      out->SetPoint(bin,bin,fsigsig.getVal());
      out->SetPointError(bin,0,fsigsig.getPropagatedError(*(fr[bin]->fr)));
    }

    if (single_gamma){
      rf1=(RooRealVar*)(fr[bin]->fr->floatParsFinal().find("rf1"));
      RooFormulaVar fsig("fsig","fsig","rf1",RooArgList(*rf1));
      RooFormulaVar fbkg("fbkg","fbkg","1-rf1",RooArgList(*rf1));
      out->SetPoint(bin,bin,fsig.getVal());
      //      out->SetPoint(bin,bin,fr[bin]->purity5);
      out->SetPointError(bin,0,fsig.getPropagatedError(*(fr[bin]->fr)));
    }
  }

  fits_canv->Update();

  if (!study_templates){
  TCanvas *output_canv = new TCanvas("output_canv","output_canv");
  output_canv->cd();
  out->Draw("AP");
  output_canv->Update();
  }
		

};


fit_output* fit_dataset(const char* _varname, const char* inputfilename, TString splitting, float leftrange, float rightrange, int bin, TCanvas *canv=NULL, bool single_gamma=false, bool do_order=false){
 
  fit_output *out = new fit_output();

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

  //TString data_dir("mc_Tree_standard_sel/");
  TString data_dir("data_Tree_standard_sel/");
  TString mc_dir("mc_Tree_standard_sel/");

  template_production *helper = new template_production(NULL);

  TString ord= (do_order) ? "order" : "noorder";
  if (single_gamma) ord="noorder";

  if (single_gamma){

  if (splitting=="EB"){
    inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,2,ord).Data()),roovar[0]);
    inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("both"),bin,ord).Data()),roohist[0][0]);
    inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("both"),bin,ord).Data()),roohist[1][0]);
  }
  else if (splitting=="EE"){
    inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,1,2,TString("noorder")).Data()),roovar[0]);
    inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EE"),TString("both"),bin,ord).Data()),roohist[0][0]);
    inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EE"),TString("both"),bin,ord).Data()),roohist[1][0]);
  }

  assert (roovar[0]!=NULL);
  assert (roohist[0][0]!=NULL);
  assert (roohist[1][0]!=NULL);

  }

  if (!single_gamma && !do_order && splitting=="EEEB") splitting="EBEE";

  if (!single_gamma){
  
    if (splitting=="EBEB") {
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,0,ord).Data()),roovar[0]);
      inputfile->GetObject(TString(data_dir).Append(helper->get_roovar_name(varname,0,1,ord).Data()),roovar[1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("1"),bin,ord).Data()),roohist[0][0]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("sig"),TString("EB"),TString("2"),bin,ord).Data()),roohist[0][1]);
      inputfile->GetObject(TString(mc_dir).Append(helper->get_roohist_name(varname,TString("bkg"),TString("EB"),TString("1"),bin,ord).Data()),roohist[1][0]);
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


    for (int i=0; i<2; i++){
      assert (roovar[i]!=NULL);
      std::cout << "using roovar named " << roovar[i]->getTitle().Data() << std::endl; 
          for (int j=0; j<2; j++){
	    assert (roohist[i][j]!=NULL);
	    std::cout << "using roohist named " << roohist[i][j]->GetTitle() << std::endl; 
	  }
    }

  }
  
  RooFormulaVar *bothcut;

      cut[0] = new RooFormulaVar("cut_0","cut_0", Form("((%s>%f) && (%s<%f))",roovar[0]->getTitle().Data(),leftrange,roovar[0]->getTitle().Data(),rightrange), RooArgList(*(roovar[0])));
      if (!single_gamma){
	cut[1] = new RooFormulaVar("cut_1","cut_1", Form("((%s>%f) && (%s<%f))",roovar[1]->getTitle().Data(),leftrange,roovar[1]->getTitle().Data(),rightrange), RooArgList(*(roovar[1])));
	bothcut = new RooFormulaVar("bothcut","bothcut","cut_0 && cut_1",RooArgList(*(cut[0]),*(cut[1])));
      }
//      std::cout << "Using cuts: " << std::endl;
//      cut[0]->Print();
//      cut[1]->Print();

      TString dset_name="";
      if (!single_gamma){
	dset_name=TString(data_dir).Append(helper->get_roodset_name(varname,splitting,bin,ord).Data());
      }
      if (single_gamma){
	dset_name=TString(data_dir).Append(helper->get_roodset_name_single(varname,splitting,bin).Data());
      }
      inputfile->GetObject(dset_name.Data(),roodset);
      cout << "using dset " << dset_name.Data() << std::endl;

      if (!single_gamma){
    roodset=(RooDataSet*)(roodset->reduce(Cut(*bothcut)));
    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
	roohist[i][j]=(RooDataHist*)(roohist[i][j]->reduce(Cut(*(cut[j]))));
      }
      if (single_gamma){
	roodset=(RooDataSet*)(roodset->reduce(Cut(*(cut[0]))));
      for (int i=0; i<2; i++)
	roohist[i][0]=(RooDataHist*)(roohist[i][0]->reduce(Cut(*(cut[0]))));
      }

    bool wrong=false;

    if (roodset==NULL) wrong=true;
    if (!single_gamma){
    for (int i=0; i<2; i++){
      if (roovar[i]==NULL) wrong=true;
      for (int j=0; j<2; j++){
	if (roohist[i][j]==NULL) wrong=true;
      }
    }
    }
    if (single_gamma){
      if (roovar[0]==NULL) wrong=true;
      for (int i=0; i<2; i++){
        if (roohist[i][0]==NULL) wrong=true;
      }
    }

    if (wrong==true) {
      std::cout << "something wrong with initialization; exiting." << std::endl;
      return;
    }


//    for (int i=0; i<2; i++){
//      roovar[i]->Print();
//      for (int j=0; j<2; j++){
//	roohist[i][j]->Print();
//      }
//    }
//    roodset->Print();

  roovar[0]->setRange(leftrange,rightrange);
  if (!single_gamma) roovar[1]->setRange(leftrange,rightrange);

//  const int n_variable_bins=25;
//  roovar[0]->setBins(n_variable_bins);
//  roovar[1]->setBins(n_variable_bins);

  RooRealVar rf1("rf1","rf1",0.5,0,1);
  RooRealVar rf2("rf2","rf2",0.5,0,1);
  RooRealVar rf3("rf3","rf3",0.5,0,1);

  if (!single_gamma){
    if (do_order) {
      rf1.setVal(1./4);
      rf2.setVal(1./3);
      rf3.setVal(1./2);
    }
    if (!do_order) {
      rf1.setVal(1./3);
      rf2.setVal(1./2);
    }
  }

    RooHistPdf *sig1pdf;
    RooHistPdf *sig2pdf;
    RooHistPdf *bkg1pdf;
    RooHistPdf *bkg2pdf;
    RooProdPdf *sigsigpdf;
    RooProdPdf *sigbkgpdf;
    RooProdPdf *bkgsigpdf;
    RooProdPdf *sigbkgpdf_order;
    RooProdPdf *bkgsigpdf_order;
    RooAddPdf *sigbkgpdf_noorder;
    RooProdPdf *bkgbkgpdf;
    RooFormulaVar *fsigsig;
    RooFormulaVar *fsigbkg;
    RooFormulaVar *fbkgsig;
    RooFormulaVar *fbkgbkg;

    RooHistPdf *sigpdf;
    RooHistPdf *bkgpdf;
    RooFormulaVar *fsig;
    RooFormulaVar *fbkg;

  if (!single_gamma){

sig1pdf = new RooHistPdf("sig1pdf","sig1pdf",*(roovar[0]),*(roohist[0][0]));
sig2pdf = new RooHistPdf("sig2pdf","sig2pdf",*(roovar[1]),*(roohist[0][1]));
bkg1pdf = new RooHistPdf("bkg1pdf","bkg1pdf",*(roovar[0]),*(roohist[1][0]));
bkg2pdf = new RooHistPdf("bkg2pdf","bkg2pdf",*(roovar[1]),*(roohist[1][1]));

 sigsigpdf = new RooProdPdf("sigsigpdf","sigsigpdf",RooArgList(*sig1pdf,*sig2pdf));
 bkgbkgpdf = new RooProdPdf("bkgbkgpdf","bkgbkgpdf",RooArgList(*bkg1pdf,*bkg2pdf));

fsigsig = new RooFormulaVar("fsigsig","fsigsig","rf1",RooArgList(rf1));

 if (do_order) {
   sigbkgpdf = new RooProdPdf("sigbkgpdf","sigbkgpdf",RooArgList(*sig1pdf,*bkg2pdf));
   bkgsigpdf = new RooProdPdf("bkgsigpdf","bkgsigpdf",RooArgList(*bkg1pdf,*sig2pdf));
   fsigbkg = new RooFormulaVar("fsigbkg","fsigbkg","(1-rf1)*rf2",RooArgList(rf1,rf2));
   fbkgsig = new RooFormulaVar("fbkgsig","fbkgsig","(1-rf1)*(1-rf2)*rf3",RooArgList(rf1,rf2,rf3));
   fbkgbkg = new RooFormulaVar("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg-fbkgsig",RooArgList(*fsigsig,*fsigbkg,*fbkgsig));
 }
 if (!do_order){
   sigbkgpdf_order = new RooProdPdf("sigbkgpdf_order","sigbkgpdf_order",RooArgList(*sig1pdf,*bkg2pdf));
   bkgsigpdf_order = new RooProdPdf("bkgsigpdf_order","bkgsigpdf_order",RooArgList(*bkg1pdf,*sig2pdf));
   sigbkgpdf_noorder = new RooAddPdf("sigbkgpdf","sigbkgpdf",RooArgList(*sigbkgpdf_order,*bkgsigpdf_order),RooArgList(RooRealConstant::value(0.5),RooRealConstant::value(0.5)));
   fsigbkg = new RooFormulaVar("fsigbkg","fsigbkg","(1-rf1)*rf2",RooArgList(rf1,rf2));
   fbkgbkg = new RooFormulaVar("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg",RooArgList(*fsigsig,*fsigbkg));
 }



  }

  if (single_gamma){
sigpdf = new RooHistPdf("sigpdf","sigpdf",*(roovar[0]),*(roohist[0][0]));
bkgpdf = new RooHistPdf("bkgpdf","bkgpdf",*(roovar[0]),*(roohist[1][0]));
fsig = new RooFormulaVar("fsig","fsig","rf1",RooArgList(rf1));
fbkg = new RooFormulaVar("fbkg","fbkg","1-rf1",RooArgList(rf1));
  }

  RooAddPdf *model;

  if (single_gamma) model = new RooAddPdf("model","model",RooArgList(*sigpdf,*bkgpdf),RooArgList(rf1),kTRUE);
  if (!single_gamma && do_order) model = new RooAddPdf("model","model",RooArgList(*sigsigpdf,*sigbkgpdf,*bkgsigpdf,*bkgbkgpdf),RooArgList(rf1,rf2,rf3),kTRUE);
  if (!single_gamma && !do_order) model = new RooAddPdf("model","model",RooArgList(*sigsigpdf,*sigbkgpdf_noorder,*bkgbkgpdf),RooArgList(rf1,rf2),kTRUE);



  if (study_templates){
    RooDataSet *roodset_toy = model->generate(RooArgList(*(roovar[0]),*(roovar[1])),Name("toy_model"),NumEvents(1000000),AutoBinned(kTRUE));
    roodset=roodset_toy;
  }

  RooFitResult *fitres = model->fitTo(*roodset,Save());

  out->fr=fitres;
  out->purity=rf1.getVal();
  if (single_gamma){
    roovar[0]->setRange("range5",-1,5);
  float integral5_tot=model->createIntegral(RooArgSet(*(roovar[0])),Range("range5"))->getVal();
  float integral5_sig=sigpdf->createIntegral(RooArgSet(*(roovar[0])),Range("range5"))->getVal()*(rf1.getVal());
  out->purity5=integral5_sig/integral5_tot;
  cout << "blablabla" << endl;
  cout << integral5_sig << " " << integral5_tot << " " << out->purity5 << " " << out->purity << endl;
  }

  model->Print();


  if (canv!=NULL){
    RooPlot *varframe[2];
    int a;
    if (single_gamma) a=1; else a=2;
    for (int i=0; i<a; i++){
      varframe[i] = roovar[i]->frame("","");
      if(!study_templates){
      roodset->plotOn(varframe[i]);
      model->plotOn(varframe[i]);
      }
      if (!single_gamma){
      model->plotOn(varframe[i],Components("sigsigpdf"),LineStyle(kDashed),LineColor(kRed));
      model->plotOn(varframe[i],Components("sigbkgpdf"),LineStyle(kDashed),LineColor(kGreen));
      if (do_order) model->plotOn(varframe[i],Components("bkgsigpdf"),LineStyle(kDashed),LineColor(kBlue));
      model->plotOn(varframe[i],Components("bkgbkgpdf"),LineStyle(kDashed),LineColor(kBlack));
      }
      if (single_gamma){
      model->plotOn(varframe[i],Components("sigpdf"),LineStyle(kDashed),LineColor(kRed));
      model->plotOn(varframe[i],Components("bkgpdf"),LineStyle(kDashed),LineColor(kBlack));
      }
      canv->cd(2*bin+i+1);
      varframe[i]->Draw();
    }
  }



  return out;
  
};


