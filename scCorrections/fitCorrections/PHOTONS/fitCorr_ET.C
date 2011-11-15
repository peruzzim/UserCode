//#include "corrections_ET.C"


int fitCorr_ET(){
 
  TFile *f = new TFile("./histCorrections_ET.root","read");  



  TH1F *hEB = (TH1F*)f->Get("h_CBET_EB");
  TH1F *hEE = (TH1F*)f->Get("h_CBET_EE");
  
  gROOT->SetStyle("Plain");
  
  TCanvas *c = new TCanvas("c","c",1200,600);
  c->Divide(2,1);
  c->cd(1);
  hEB->Draw();
  TF1 * fEB1 = new TF1("fEB1","pol1[0]",20,140);
  hEB->Fit("fEB1","","same",20,140);  
  TF1 * fEB0 = new TF1("fEB0","pol0[0]",140,250);
  hEB->Fit("fEB0","","same",140,250);  
  
  // Write out the correction for the EB
  ofstream outfile;  
  TString filename = "./corrections_ET.C";
  outfile.open(filename);
  outfile << "Double_t applyScCorrectionsET_EB(Double_t ET){      " << endl;    
  outfile << endl;
  outfile << "  Double_t par0 =  " << hEB->GetBinContent(0*1+1) << "; " << endl;
  outfile << "  Double_t par1 =  " << hEB->GetBinContent(2*1+1) << "; " << endl;
  outfile << "  Double_t par2 =  " << fEB1->GetParameter(0) << ";	   " << endl;
  outfile << "  Double_t par3 =  " << fEB1->GetParameter(1) << ";	   " << endl;
  outfile << "  Double_t par4 =  " << fEB0->GetParameter(0) << ";     " << endl;
  outfile << endl;  
  outfile << "  if (             ET <   5 ) return         1.;  " << endl;
  outfile << "  if (  5 <= ET && ET <  10 ) return         par0 ;  " << endl;
  outfile << "  if ( 10 <= ET && ET <  20 ) return         par1 ;  " << endl;
  outfile << "  if ( 20 <= ET && ET < 140 ) return         par2 + par3*ET ;  " << endl;
  outfile << "  if (140 <= ET             ) return         par4;  " << endl;
  outfile << " 						  " << endl;
  outfile << "}                                           " << endl;
  outfile << "                					        " << endl;
  outfile << "Double_t F_applyScCorrectionsET_EB(Double_t *xx, Double_t *pp){  " << endl;        
  outfile << "  Double_t x = xx[0] ;                     " << endl;
  outfile << "  return applyScCorrectionsET_EB(x);    		  " << endl;
  outfile << "}		                          " << endl;
  
  c->cd(2);
  hEE->Draw();
  
  TF1 * fEE2 = new TF1("fEE2","pol2[0]",30,250);
  hEE->Fit("fEE2","","same",30,250);  
  
  outfile << "Double_t applyScCorrectionsET_EE(Double_t ET){      " << endl;    
  outfile << "   					  " << endl;
  outfile << "  Double_t par0 =  " << hEE->GetBinContent(2*0+1)  << "; " << endl;
  outfile << "  Double_t par1 =  " << hEE->GetBinContent(2*1+1)  << "; " << endl;
  outfile << "  Double_t par2 =  " << hEE->GetBinContent(2*2+1)  << ";	   " << endl;
  outfile << "  Double_t par3 =  " << fEE2->GetParameter(0)      << ";	   " << endl;
  outfile << "  Double_t par4 =  " << fEE2->GetParameter(1)      << ";     " << endl;
  outfile << "  Double_t par5 =  " << fEE2->GetParameter(2)      << ";     " << endl;
  outfile << "  Double_t par6 =  " << hEE->GetBinContent(hEE->GetNbinsX())  << ";     " << endl;
  outfile << endl;
  outfile << "  if (             ET <   5 ) return         1.;  " << endl;
  outfile << "  if (  5 <= ET && ET <  10 ) return          par0 ;  " << endl;
  outfile << "  if ( 10 <= ET && ET <  20 ) return          par1 ;  " << endl;
  outfile << "  if ( 20 <= ET && ET <  30 ) return          par2 ;  " << endl;
  outfile << "  if ( 30 <= ET && ET < 200 ) return          par3 + par4 *ET + par5 *ET*ET ;  " << endl;
  outfile << "  if ( 200 <= ET            ) return          par6 ;"<< endl;
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
