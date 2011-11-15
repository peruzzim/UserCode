
int fitCorr_E(){
 
  TFile *f = new TFile("./histCorrections_E.root","read");  

  TH1F *hEB = (TH1F*)f->Get("h_CBE_EB");
  TH1F *hEE = (TH1F*)f->Get("h_CBE_EE");
  
  gROOT->SetStyle("Plain");
  
  TCanvas *c = new TCanvas("c","c",1200,600);
  c->Divide(2,1);
//   c->cd(1);
//   hEB->SetMinimum(0.97);
//   hEB->Draw();
//   TF1 * fEB = new TF1("fEB","([0]+x*[3])*(1-[2]*exp(x/[1]))",0,250);
//   fEB->SetParameters(1,-1,1,1);
//   hEB->Fit("fEB","","same",10,80);  
//   hEB->Fit("fEB","","same",10,200);  

  // Write out the correction for the EB
  ofstream outfile;  
  TString filename = "./corrections_E.C";
  outfile.open(filename);

//   outfile << " Double_t applyScCorrectionsE_EB(Double_t E){      " << endl;    
//   outfile << "   				 	  " << endl;
//   outfile << "   if (E > 200) E =200;   		  " << endl;  
//   outfile << "   if (             E <    5 ) return         1.;  " << endl;
//   outfile << "   if (  5 <= E && E <   10 ) return         " << hEB->GetBinContent(0*1+1) << ";  " << endl;
//   outfile << "   if ( 10 <= E && E <= 200 ) return         (" << fEB->GetParameter(0) << " + E*" << fEB->GetParameter(3) << ")*(1-"<< fEB->GetParameter(2) 
// 	  <<"*exp(E/"<< fEB->GetParameter(1)<<  "));" << endl;
//   outfile << " 						  " << endl;
//   outfile << " }                                           " << endl;
//   outfile << "                					        " << endl;
//   outfile << " Double_t F_applyScCorrectionsE_EB(Double_t *xx, Double_t *pp){  " << endl;        
//   outfile << "   Double_t x = xx[0] ;                     " << endl;
//   outfile << "   return applyScCorrectionsE_EB(x);    		  " << endl;
//   outfile << " }		                          " << endl;
//   outfile << "  		                          " << endl;
  
  c->cd(2);
  hEE->SetMaximum(1.02);
  hEE->SetMinimum(0.92);
  hEE->Draw();
  
  TF1 * fEE = new TF1("fEE","pol1",0,5000);
  fEE->SetParameters(1,1);
  hEE->Fit("fEE","","same",10,80);  
  hEE->Fit("fEE","","same",10,200);  
  hEE->Fit("fEE","","same",10,850);  
  
  Double_t upperLimit =  (1- fEE->GetParameter(0) ) /  fEE->GetParameter(1);
  
  outfile << "Double_t applyScCorrectionsE_EE(Double_t E){      " << endl;    
  outfile << "  				 	  " << endl;
  outfile << " Double_t par0 = 850;               " << endl; 
  outfile << " Double_t par1 = " << fEE->GetParameter(0) << " ;	  " << endl;
  outfile << " Double_t par2 = " << fEE->GetParameter(1) << " ;     " << endl;
  outfile << "  				 	  " << endl;
  outfile << "  if (E  > par0 ) E = par0 ;   		  " << endl;  
  outfile << "  if (            E <   0     ) return      1.;  " << endl;
  outfile << "  if (  0 <= E && E <=  par0  ) return      par1 + E*par2; " << endl;
  outfile << "						   " << endl;
  outfile << "}                                           " << endl;
  outfile << "               					        " << endl;
  outfile << "Double_t F_applyScCorrectionsE_EE(Double_t *xx, Double_t *pp){  " << endl;        
  outfile << "  Double_t x = xx[0] ;                     " << endl;
  outfile << "  return applyScCorrectionsE_EE(x);    		  " << endl;
  outfile << "}		                          " << endl;
  outfile.close();
  
  return 0;
}
