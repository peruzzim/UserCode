#include "../../binning.h"

int fitCorr_brEta(){ 
  
  gROOT->SetStyle("Plain");
  
  TFile *f = new TFile("./histCorrections_brEta.root","read");

  Double_t xcorr[nBinsEta]; // it contains the first bin, not fitted by the function ftest
  
  Double_t par0[nBinsEta];
  Double_t par1[nBinsEta];
  Double_t par2[nBinsEta];
  Double_t par3[nBinsEta];
  Double_t par4[nBinsEta];

  TF1 * ftest = new TF1("ftest","pol1",1.2,10);
  ftest->SetParameters(1,-1);

  TF1 * ftest2 = new TF1("ftest2","pol2",1.2,10);
  ftest2->SetParameters(1,-1,-1);
  
  //  TF1 * ftest = new TF1("ftest","[0]*(1.-exp(-(x-[4])/[1]))*[2]*x + [3]",1.2,10);
  //  TF1 * ftest2 = new TF1("ftest2","([0]+x*[3])*(1-[2]*exp(x/[1]))",0,10);
  //ftest2->SetParameters(1,-1,1,-0.1);
  
  for (Int_t i = 0 ; i< nBinsEta; ++i){
    //  for (Int_t i = 0 ; i< 7; ++i){
    
    TH1F *h = (TH1F*)f->Get(Form("h_corr_%d",i));    
    
    xcorr[i] = (h->GetBinContent(2*0+1) + h->GetBinContent(2*1+1)) / 2.;
    cout <<     xcorr[i]  << " " << h->GetBinContent(2*0+1)  << " " << h->GetBinContent(2*1+1) << endl;

    TCanvas *c = new TCanvas("c","c",600,600);
    c->cd();
    h->GetYaxis()->SetRangeUser(0.90,1.01);
    if (i>6) h->GetYaxis()->SetRangeUser(0.85,1.01);
    h->GetXaxis()->SetRangeUser(0,5);
    h->Draw();
    TLatex l(3.5,0.95,Form("h_corr_%d",i));
    l.Draw();
    
    if ( i < 7){
      ftest->SetParameters(1,-1.);
      h->Fit("ftest","","same",1.2,5.0);
      par0[i] = ftest->GetParameter(0);
      par1[i] = ftest->GetParameter(1);
      ftest->Draw("same");    
    }
    else{
      ftest2->SetParameters(-1,-1,-1);
      h->Fit("ftest2","","same",0.8,5.0);
      par0[i] = ftest2->GetParameter(0);
      par1[i] = ftest2->GetParameter(1);
      par2[i] = ftest2->GetParameter(2);
      ftest2->Draw("same");    
    }    
    
    c->Update();
    getchar();
  }
  
  // check printouts
  for (Int_t i = 0 ; i<7; ++i){
    //    cout << par0[i] << "\t" <<par1[i] << "\t" <<par2[i] << "\t" <<par3[i] << "\t" <<par4[i] << endl;
    cout << "TF1 ftest_"<< i << "(\"ftest\",\"pol1\",0.8,10); " << endl;
    cout << "ftest_" << i << "->SetParameters(" << par0[i] << ", " << par1[i] << " );"<< endl;
  }
  for (Int_t i = 7 ; i<nBinsEta; ++i){
    cout << "TF1 ftest_"<< i << "(\"ftest\",\"pol2\",0.8,10); " << endl;
    cout << "ftest_" << i << "->SetParameters(" << par0[i] << ", " << par1[i] << ", " << par2[i] << " );"<< endl;
  }

  for (Int_t i = 0 ; i<nBinsEta; ++i){
    cout << i << " " <<   xcorr[i]  << endl;
  }


  // test plot functions
  TCanvas *cf = new TCanvas("cf","cf",600,600);
  cf->cd();
  
  TF1 * ftest_0 = new TF1("ftest","pol1",0.8,10); 
  ftest_0->SetParameters(par0[0],par1[0]);
  TF1 * ftest_1 = new TF1("ftest","pol1",0.8,10); 
  ftest_1->SetParameters(par0[1],par1[1]);
  TF1 * ftest_2 = new TF1("ftest","pol1",0.8,10); 
  ftest_2->SetParameters(par0[2],par1[2]);
  TF1 * ftest_3 = new TF1("ftest","pol1",0.8,10); 
  ftest_3->SetParameters(par0[3],par1[3]);
  TF1 * ftest_4 = new TF1("ftest","pol1",0.8,10); 
  ftest_4->SetParameters(par0[4],par1[4]);
  TF1 * ftest_5 = new TF1("ftest","pol1",0.8,10); 
  ftest_5->SetParameters(par0[5],par1[5]);
  TF1 * ftest_6 = new TF1("ftest","pol1",0.8,10); 
  ftest_6->SetParameters(par0[6],par1[6]);

  TF1 * ftest_7 = new TF1("ftest","pol2",0.8,10); 
  ftest_7->SetParameters(par0[7],par1[7],par2[7]);
  TF1 * ftest_8 = new TF1("ftest","pol2",0.8,10); 
  ftest_8->SetParameters(par0[8],par1[8],par2[8]);
  TF1 * ftest_9 = new TF1("ftest","pol2",0.8,10); 
  ftest_9->SetParameters(par0[9],par1[9],par2[9]);
  TF1 * ftest_10 = new TF1("ftest","pol2",0.8,10); 
  ftest_10->SetParameters(par0[10],par1[10],par2[10]);
  TF1 * ftest_11 = new TF1("ftest","pol2",0.8,10); 
  ftest_11->SetParameters(par0[11],par1[11],par2[11]);
  TF1 * ftest_12 = new TF1("ftest","pol2",0.8,10); 
  ftest_12->SetParameters(par0[12],par1[12],par2[12]);
  TF1 * ftest_13 = new TF1("ftest","pol2",0.8,10); 
  ftest_13->SetParameters(par0[13],par1[13],par2[13]);
  
  TH2F *hh = new TH2F("hh","hh",100,0,5,100,0.7,1.1);
  hh->Draw();
  ftest_0 ->Draw("same");
  ftest_1 ->Draw("same");
  ftest_2 ->Draw("same");
  ftest_3 ->Draw("same");
  ftest_4 ->Draw("same");
  ftest_5 ->Draw("same");
  ftest_6 ->Draw("same");
  ftest_7 ->Draw("same");
  ftest_8 ->Draw("same");
  ftest_9 ->Draw("same");
  ftest_10->Draw("same");
  ftest_11->Draw("same");
  ftest_12->Draw("same");
  ftest_13->Draw("same");
  
  c->Update();

  cout << "================================================================================"<< endl;

  ofstream outfile;  
  TString filename = "./corrections_brEta.C";
  outfile.open(filename);
  //  outfile << "#include \"TF1.h\"" << endl;
  outfile << "#include \"TMath.h\"" << endl;
  outfile << endl;
  outfile << "Double_t applyScCorrectionsBrEta(Double_t eta, Double_t sigmaPhiSigmaEta){ " << endl;
  outfile << endl;
  outfile << "  bool DBG = false;" << endl;
  outfile << endl;
  outfile << "  // eta binning ------------------------------------------------------------------------------- " << endl;
  outfile << "  // " << endl;
  outfile << "  const Double_t etaCrackMin = 1.44; " << endl;
  outfile << "  const Double_t etaCrackMax = 1.56; " << endl;
  outfile << "   " << endl;
  outfile << "  //STD " << endl;
  outfile << "  const Int_t    nBinsEta              = 14; " << endl;
  outfile << "  Double_t       leftEta  [nBinsEta]   = { 0.02, 0.25, 0.46, 0.81, 0.91, 1.01, 1.16,           etaCrackMax,  1.653,  1.8, 2.0, 2.2, 2.3, 2.4 };  " << endl;
  outfile << "  Double_t       rightEta [nBinsEta]   = { 0.25, 0.42, 0.77, 0.91, 1.01, 1.13, etaCrackMin,    1.653,        1.8  ,  2.0, 2.2, 2.3, 2.4, 2.5 };  " << endl;
  outfile << endl;
  outfile << "  Double_t xcorr[nBinsEta];" << endl;
  for (Int_t i = 0; i < nBinsEta; ++i){
    outfile << "  xcorr[" << i << "]=" << xcorr[i] << ";"<< endl;
  }
  outfile << endl;  
  outfile << "  Double_t par0[nBinsEta]; " << endl;
  outfile << "  Double_t par1[nBinsEta]; " << endl;
  outfile << "  Double_t par2[nBinsEta]; " << endl;  
  outfile << endl;
  // check printouts
  for (Int_t i = 0 ; i<7; ++i){
    //     outfile << "  fcorr["<<i<<"] = TF1(\"ftest_"<<i<<"\",\"pol1\",1,10); " <<endl;
    //     outfile << "  fcorr[" << i << "].SetParameters(" << par0[i] << ", " << par1[i] <<  ");"<< endl;
    outfile << "  par0[" << i << "] = " << par0[i] << ";" << endl;
    outfile << "  par1[" << i << "] = " << par1[i] << ";" << endl;
    outfile << "  par2[" << i << "] = 0 ;" << endl;
    outfile << endl;
  }  
  for (Int_t i = 7 ; i<nBinsEta; ++i){
    outfile << "  par0[" << i << "] = " << par0[i] << ";" << endl;
    outfile << "  par1[" << i << "] = " << par1[i] << ";" << endl;
    outfile << "  par2[" << i << "] = " << par2[i] << ";" << endl;
    outfile << endl;
  }
  //  outfile << "  TF1 fcorr[" << nBinsEta << "]; " << endl;
  //   for (Int_t i = 7 ; i<nBinsEta; ++i){
  //     outfile << "  fcorr["<<i<<"] = TF1(\"ftest_"<<i<<"\",\"pol2\",1,10); " <<endl;
  //     outfile << "  fcorr[" << i << "].SetParameters(" << par0[i] << ", " << par1[i] << ", " << par2[i] <<  " );"<< endl;
  //   }  
  outfile << endl;
  outfile << "  // extra protections																					   " << endl;
  outfile << "  // fix sigmaPhiSigmaEta boundaries " <<endl;
  outfile << "  if (sigmaPhiSigmaEta < 0.8)  sigmaPhiSigmaEta = 0.8; "<< endl;
  outfile << "  if (sigmaPhiSigmaEta > 5  )  sigmaPhiSigmaEta = 5; " << endl;
  outfile << endl;
  outfile << "  // eta = 0																						   " << endl;
  outfile << "  if (TMath::Abs(eta)  <  leftEta[0]            ) { eta = " <<  leftEta[0] << " ; }																   " << endl;
  outfile << "  // outside acceptance																					   " << endl;
  outfile << "  if (TMath::Abs(eta)  >=  rightEta[nBinsEta-1] ) { eta = " << rightEta[nBinsEta - 1]-0.01<< "; if (DBG) std::cout << \" WARNING [applyScCorrections]: TMath::Abs(eta)  >=  rightEta[nBinsEta-1] \" << std::endl;}  " << endl;
  outfile << "  																								   " << endl;
  outfile << "  Int_t tmpEta = -1;                                                                                                                                                                         " << endl;
  outfile << "  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								              								      	   " << endl;
  outfile << "    if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] ){				       									      	   " << endl;
  outfile << "      tmpEta = iEta;											       										   " << endl;
  outfile << "    }													       										   " << endl;
  outfile << "  }													       										           " << endl; 
  outfile << endl;
  outfile << "  // Interpolation																					         " << endl;
  outfile << "  Double_t tmpInter = 1;																				         " << endl;
  outfile << "  // In eta cracks/gaps 																				         " << endl;
  outfile << "  if (tmpEta == -1 ) { // need to interpolate    " << endl;
  outfile << "    for (Int_t iEta = 0; iEta < nBinsEta-1; ++iEta){								       								         " << endl;
  outfile << "      if (rightEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[iEta+1] ){													         " << endl;
  outfile << "	if (sigmaPhiSigmaEta >= 1.2)  tmpInter = ( par0[iEta] + sigmaPhiSigmaEta*par1[iEta] + sigmaPhiSigmaEta*sigmaPhiSigmaEta*par2[iEta] +  " << endl;
  outfile << "						   par0[iEta+1] + sigmaPhiSigmaEta*par1[iEta+1] + sigmaPhiSigmaEta*sigmaPhiSigmaEta*par2[iEta+1]) / 2. ; " << endl;
  outfile << "	else tmpInter = (xcorr[iEta] + xcorr[iEta+1])/2.; " << endl;
  outfile << "      }																						         " << endl;
  outfile << "    }																						         " << endl;
  outfile << "    return tmpInter;																					         " << endl;
  outfile << "  }  																							         " << endl;
  outfile << "  if (sigmaPhiSigmaEta >= 1.2) return par0[tmpEta] + sigmaPhiSigmaEta*par1[tmpEta] + sigmaPhiSigmaEta*sigmaPhiSigmaEta*par2[tmpEta]; " << endl;
  outfile << "  else return xcorr[tmpEta]; " << endl;
  outfile << "} " << endl;

//  outfile << "  // Interpolation																					         " << endl;
//  outfile << "  Double_t tmpInter = 1;																				         " << endl;
//  outfile << "  // In eta cracks/gaps 																				         " << endl;
//  outfile << "  if (tmpEta == -1 ) { // need to interpolate    " << endl;
//  outfile << "    for (Int_t iEta = 0; iEta < nBinsEta-1; ++iEta){								       								         " << endl;
//  outfile << "      if (rightEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[iEta+1] ){													         " << endl;
//  outfile << "        if (sigmaPhiSigmaEta >= 1.2)  tmpInter = ( fcorr[iEta]  .Eval(sigmaPhiSigmaEta) + 															          " << endl;
//  outfile << "	         	    		         fcorr[iEta+1].Eval(sigmaPhiSigmaEta) ) / 2. ;														          " << endl;
//  outfile << "        else tmpInter = (xcorr[iEta] + xcorr[iEta+1])/2.; " << endl;
//  outfile << "      }																						         " << endl;
//  outfile << "    }																						         " << endl;
//  outfile << "    return tmpInter;																					         " << endl;
//  outfile << "  }  																							         " << endl;
//  outfile << "  if (sigmaPhiSigmaEta >= 1.2) return fcorr[tmpEta].Eval(sigmaPhiSigmaEta);						       						       " << endl;
//  outfile << "  else return xcorr[tmpEta]; " << endl;
//  outfile << "} " << endl;													       										         

  outfile << endl;
  outfile << endl;
  outfile << "Double_t F_applyScCorrectionsBrEta(Double_t *xx, Double_t *pp){  " << endl;
  outfile << " " << endl;
  outfile << "  Double_t eta  = xx[0] ;                      " << endl;
  outfile << "  Double_t brem = xx[1] ;                      " << endl;
  outfile << "  return applyScCorrectionsBrEta(eta, brem);    		     " << endl;
  outfile << "} " << endl;

  return 0;
}

