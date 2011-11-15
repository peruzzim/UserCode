
int fitCorr_ET(){
  
  TFile *f = new TFile("./histCorrections_ET.root","read");  

  TH1F *hEB = (TH1F*)f->Get("h_CBET_EB");
  TH1F *hEE = (TH1F*)f->Get("h_CBET_EE");
  
  gROOT->SetStyle("Plain");
  
  TCanvas *c = new TCanvas("c","c",1200,600);
  c->Divide(2,1);
  c->cd(1);
  hEB->SetMinimum(0.97);
  hEB->Draw();
  TF1 * fEB = new TF1("fEB","([0]+x*[3])*(1-[2]*exp(x/[1]))",0,250);
  fEB->SetParameters(1,-1,1,1);
  hEB->Fit("fEB","","same",10,80);  
  hEB->Fit("fEB","","same",10,200);  

  // Write out the correction for the EB
  ofstream outfile;  
  TString filename = "./corrections_ET.C";
  outfile.open(filename);
  outfile << "Double_t applyScCorrectionsET_EB(Double_t ET){      " << endl;    
  outfile << "   				 	  " << endl;
  outfile << "  Double_t par0 = " << hEB->GetBinContent(0*1+1) << "; " << endl; 
  outfile << "  Double_t par1 = " << fEB->GetParameter(0)      << "; " << endl;
  outfile << "  Double_t par2 = " << fEB->GetParameter(3)      << "; " << endl;
  outfile << "  Double_t par3 = " << fEB->GetParameter(2)      << "; " << endl;
  outfile << "  Double_t par4 = " << fEB->GetParameter(1)      << "; " << endl;
  outfile << endl;
  outfile << "  if (ET > 200) ET =200;   		  " << endl;  
  outfile << "  if (             ET <    5 ) return         1.;  " << endl;
  outfile << "  if (  5 <= ET && ET <   10 ) return         par0 ;  " << endl;
  outfile << "  if ( 10 <= ET && ET <= 200 ) return         (par1  + ET*par2)*(1- par3*exp(ET/ par4));" << endl;
  outfile << " 						  " << endl;
  outfile << "}                                           " << endl;
  outfile << "                					        " << endl;
  outfile << "Double_t F_applyScCorrectionsET_EB(Double_t *xx, Double_t *pp){  " << endl;        
  outfile << "  Double_t x = xx[0] ;                     " << endl;
  outfile << "  return applyScCorrectionsET_EB(x);    		  " << endl;
  outfile << "}		                          " << endl;
  outfile << "  		                          " << endl;
  
  c->cd(2);
  hEE->SetMaximum(1.02);
  hEE->SetMinimum(0.92);
  hEE->Draw();
  
  TF1 * fEE = new TF1("fEE","([0]+x*[3])*(1-[2]*exp(x/[1]))",0,250);
  fEE->SetParameters(1,-1,1,1);
  hEE->Fit("fEE","","same",10,80);  
  hEE->Fit("fEE","","same",10,200);  
  
  outfile << "Double_t applyScCorrectionsET_EE(Double_t ET){      " << endl;    
  outfile << "   				 	  " << endl;
  outfile << "  Double_t par0 = " << hEE->GetBinContent(0*1+1) << "; " << endl; 
  outfile << "  Double_t par1 = " << fEE->GetParameter(0)      << "; " << endl;
  outfile << "  Double_t par2 = " << fEE->GetParameter(3)      << "; " << endl;
  outfile << "  Double_t par3 = " << fEE->GetParameter(2)      << "; " << endl;
  outfile << "  Double_t par4 = " << fEE->GetParameter(1)      << "; " << endl;
  outfile << endl;
  outfile << "  if (ET > 200) ET =200;   		  " << endl;  
  outfile << "  if (             ET <    5 ) return         1.;  " << endl;
  outfile << "  if (  5 <= ET && ET <   10 ) return         par0;  " << endl;
  outfile << "  if ( 10 <= ET && ET <= 200 ) return         ( par1  + ET*par2)*(1-par3*exp(ET/par4));" << endl;
  outfile << " 						  " << endl;
  outfile << "}                                           " << endl;
  outfile << "                					        " << endl;
  outfile << "Double_t F_applyScCorrectionsET_EE(Double_t *xx, Double_t *pp){  " << endl;        
  outfile << "  Double_t x = xx[0] ;                     " << endl;
  outfile << "  return applyScCorrectionsET_EE(x);    		  " << endl;
  outfile << "}		                          " << endl;
  outfile.close();
  
  return 0;
}
