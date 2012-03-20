using namespace RooFit;

void compare_sidebands(TVirtualPad* c1, const char* reg){

  if (c1==NULL) return;
  if (TString(reg)=="") return;

  TFile *file = TFile::Open("out_PhoIso04.root");

  TH1F *datahist_side, *mchist_noside_bkg, *mchist_side_all;


  file->GetObject(TString("data_Tree_sideband_sel/templatehist_PhoIso04_all_").Append(reg).Data(),datahist_side);
  file->GetObject(TString("mc_Tree_sideband_sel/templatehist_PhoIso04_all_").Append(reg).Data(),mchist_side_all);
  file->GetObject(TString("mc_Tree_standard_sel/templatehist_PhoIso04_bkg_").Append(reg).Data(),mchist_noside_bkg);

  datahist_side->SetLineWidth(2);
  mchist_noside_bkg->SetLineWidth(2);
  mchist_side_all->SetLineWidth(2);
  datahist_side->SetLineColor(kBlue);
  mchist_noside_bkg->SetLineColor(kGreen);
  mchist_side_all->SetLineColor(kRed);

  mchist_noside_bkg->SetMarkerColor(kGreen);
  datahist_side->SetMarkerColor(kBlue);
  mchist_side_all->SetMarkerColor(kRed);

  if (TString(reg)=="EE") {
    datahist_side->Rebin(5);
    mchist_noside_bkg->Rebin(5);
    mchist_side_all->Rebin(5);    
  }
 else if (TString(reg)=="EB") {
    datahist_side->Rebin(2);
    mchist_noside_bkg->Rebin(2);
    mchist_side_all->Rebin(2);    
  }

  datahist_side->Scale(1.0/datahist_side->Integral());
  mchist_noside_bkg->Scale(1.0/mchist_noside_bkg->Integral());
  mchist_side_all->Scale(1.0/mchist_side_all->Integral());

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(datahist_side,"DATA ALL SIDE","l");
  leg->AddEntry(mchist_noside_bkg,"MC BKG STD","l");
  leg->AddEntry(mchist_side_all,"MC ALL SIDE","l");

  //  TCanvas *c1 = new TCanvas();

  mchist_noside_bkg->Draw();
  mchist_side_all->Draw("same");
  datahist_side->Draw("same");

  leg->Draw();

  c1->Update();
  //  c1->SaveAs("compare.png");
  
}

void run_compare_sidebands(){

  TCanvas *canv = new TCanvas();
  canv->Divide(2,1);
  TVirtualPad *c;
  c=canv->cd(1);
  compare_sidebands(c,"EB");
  c=canv->cd(2);
  compare_sidebands(c,"EE");

};

//
//void fit_templates(){
//
//  RooRealVar sieievar("sieievar","sieievar",0,0.05);
//  RooPlot *sieievarframe = sieievar.frame(Name("sietaieta"),Title("sietaieta"));
//
//  RooDataHist *hdata, *hbkgtmp, *hsigtmp;
//
//  TFile *datafile = TFile::Open("dataEE.root");
//  datafile->GetObject("sieie_all",hdata);
//
//  TFile *sigtemplatefile = TFile::Open("mcEE.root");
//  sigtemplatefile->GetObject("sieie_sig",hsigtmp);
//
//  //TFile *bkgtemplatefile = TFile::Open("mcEE.root");
//  //bkgtemplatefile->GetObject("sieie_bkg",hbkgtmp);
//  //TFile *bkgtemplatefile = TFile::Open("dataEE_side.root");
//  //bkgtemplatefile->GetObject("sieie_all",hbkgtmp);
//  TFile *bkgtemplatefile = TFile::Open("mcEE_side.root");
//  bkgtemplatefile->GetObject("sieie_all",hbkgtmp);
//
//  RooRealVar fsig("fsig","fsig",0.8,0.5,1.0);
//  
//  RooHistPdf pdfmcsig("pdfmcsig","pdfmcsig",sieievar,*hsigtmp);
//  RooHistPdf pdfmcbkg("pdfmcbkg","pdfmcbkg",sieievar,*hbkgtmp);
//
//  RooAddPdf pdfmc("pdfmc","pdfmc",pdfmcsig,pdfmcbkg,fsig);
//
//  RooFitResult *fitres = pdfmc.fitTo(*hdata,Save());
//
//  fitres->Print();
//
//
//  TCanvas *c1 = new TCanvas();
//  hdata->plotOn(sieievarframe);
//  pdfmc.plotOn(sieievarframe,LineColor(kBlue));
//  pdfmc.plotOn(sieievarframe,Components(pdfmcsig),LineColor(kRed),LineStyle(kDashed));
//  pdfmc.plotOn(sieievarframe,Components(pdfmcbkg),LineColor(kGreen),LineStyle(kDashed));
//  pdfmc.paramOn(sieievarframe);
//  sieievarframe->Draw();
//  c1->SetLogy();
//  c1->Update();
//  c1->SaveAs("fit.png");
//
//  cout << "Purity: " << fsig.getVal() << endl;
//
//}
//
