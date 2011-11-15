//================================================================================
// SuperCluster corrections: 
//  Implementation of the parameterization f(sigma_phi/sigma_eta, eta) * F(ET)
//
// Author List:
//      Mauro Donega
//================================================================================
//
// To run ROOFIT add these lines in your rootlogon.C :
//   gSystem->Load("libRooFit.so") ;
//   gSystem->Load("libRooFitCore.so") ;
//     
// To extract f(sigma_phi/sigma_eta, eta)            : root -l run_scCorrections_brEta_PHOTONS.C
// To extract f(sigma_phi/sigma_eta, eta)*F(ET)      : root -l run_scCorrections_brEtaET_PHOTONS.C
// Closure test of f(sigma_phi/sigma_eta, eta)       : root -l run_applyscCorrections_brEta_PHOTONS.C
// Closure test of f(sigma_phi/sigma_eta, eta)*F(ET) : root -l run_applyscCorrections_brEtaET_PHOTONS.C
// (same for ELECTRONS)
//
// The extracted corrections are written out as a simple macro: applyScCorrections.C
//  ( plus the correction histograms histCorrections_ET.root ; histCorrections_brEta.root )
//
// NOTES:
// 0- modify scCorrections.h to represent your NT structure
// 1- to work on a different range/granularity just change the binning
// 2- minimum number of entries  set to minEntries = 100 in Erec/Egen histograms
// 3- CB fit is asymmetric (parameters from the driving gaussian): sigmaL, sigmaR
// 4- there's a min ET cut at minETCut = 5 GeV because of the ET range of the ntuples used
// 5- there's a fiducial range of 2 degrees in phi around the cracks
// 6- eta corrections are symmetric i.e. abs(eta)
//
// By default the corrections are taken from the correction histograms (histograms histCorrections_ET.root ; histCorrections_brEta.root )
// Uncomment the MDDB to read the corrections from file.
//
// The following tests were implemented:
// - R9 instead of sigma_phi/sigma_eta
// - 3D : sigma_phi/sigma_eta, eta, ET in one step
// - f(sigma_phi/sigma_eta, eta) * F(E) : energy instead of ET
//
//================================================================================
#define  scCorrections_cxx
#include "scCorrections.h"
#include "Utils.C"
#include <iostream>
#include <fstream>
#include <sstream>
#include "TROOT.h"
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include "TMath.h"
#include <TStyle.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TBox.h>
#include "TFile.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooPlot.h"

#include "binning.h"

//MDDB #include "applyScCorrections.C" //MDDB

// #include "./200/PHOTONS/corrections_brEta.C"
// #include "./200/PHOTONS/corrections_ET.C"

// #include "ELECTRONS/corrections_brEta.C"
// #include "ELECTRONS/corrections_ET.C"

using std::cout;
using std::endl;

using namespace RooFit;

bool plots         =  false;   //   true;      //   
TString partType   = "photon"; //  "electron"; //   

// !!! protections added only on the standard ET / E methods !!!

//ELECTRONS  !!!  - minR9  - min ET

Float_t minR9_EB       = 0.94; // 9999; // 
Float_t minR9_EE       = 0.95; // 9999; // 
Float_t minETCut       = 10.;  //   5.; // 

Float_t fiducialPhiCut = 2.; // degrees around the crack

// --------------------------------------------------------------------------------

// Binning for theE_rec/Egen histograms
Double_t nBinsErecEgen = 650;
Double_t ErecEgenMin   = 0.0;
Double_t ErecEgenMax   = 1.3;

// matrix (eta,br) of histograms used for the Erec/Egen CB fits
TH1F *h_ErecEgen[nBinsEta][nBinsBr]; 

// 2D histogram containing the results of the CB fits
TH2F *hh_CBetabr    = new TH2F("hh_CBetabr"    ,"hh_CBetabr"    ,nBinsBr*2-1, brbins, nBinsEta*2-1, etabins);

// 1D to fit the corrections
TH1F *h_CBetabr[nBinsEta];

// 2D histogram containing the chi2 of the CB fits
TH2F *hh_Chi2etabr  = new TH2F("hh_Chi2etabr"  ,"hh_Chi2etabr"  ,nBinsBr*2-1, brbins, nBinsEta*2-1, etabins);

// vector ET dependence 
TH1F *h_ETcorrBrEta_EB[nBinsET];
TH1F *h_CBET_EB    = new TH1F("h_CBET_EB"    ,"h_CBET_EB"    ,nBinsET*2-1, ETBins);
TH1F *h_Chi2ET_EB  = new TH1F("h_Chi2ET_EB"  ,"h_Chi2ET_EB"  ,nBinsET*2-1, ETBins);
TH1F *h_ETcorrBrEta_EE[nBinsET];
TH1F *h_CBET_EE    = new TH1F("h_CBET_EE"    ,"h_CBET_EE"    ,nBinsET*2-1, ETBins);
TH1F *h_Chi2ET_EE  = new TH1F("h_Chi2ET_EE"  ,"h_Chi2ET_EE"  ,nBinsET*2-1, ETBins);

// Fraction of candidates with reconstructed energy less than EgenFrac
Double_t EgenFrac = 0.3;
TH2F *hh_EgenFracEtaBr  = new TH2F("hh_EgenFracEtaBr"  ,"hh_EgenFracEtaBr"  ,nBinsBr*2-1, brbins, nBinsEta*2-1, etabins);
TH1F *h_EgenFracET_EB      = new TH1F("h_EgenFracET_EB"      ,"h_EgenFracET_EB"      ,nBinsET*2-1, ETBins);
TH1F *h_EgenFracET_EE      = new TH1F("h_EgenFracET_EE"      ,"h_EgenFracET_EE"      ,nBinsET*2-1, ETBins);




//--------------------------------------------------------------------------------

// Extract f(sigma_phi,sigma_eta, eta)
void scCorrections::run_brEta()
{    
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  initHistograms();
  correctionType Corr  = none;
  fillHistograms(Corr);
  bool apply = false;
  getErecEgenFits_BrEta(Corr, apply);
  deriveCorrections_BrEta();
}

// Apply f(sigma_phi/sigma_eta, eta) (closure test)
void scCorrections::run_apply_brEta()
{  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  initHistograms();
  correctionType Corr  = brEta;
  fillHistograms(Corr);
  bool apply = true;
  getErecEgenFits_BrEta(Corr, apply);
}

// Extract the ET correction after applying f(sigma_phi,sigma_eta, eta)
void scCorrections::run_brEtaET()
{  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  // extract the BrEta corrections
  run_brEta();

  // apply the BrEta correction
  run_apply_brEta();

  // write out the correction macro
  deriveCorrections_BrEtaET();  
}

// Apply f(sigma_phi/sigma_eta, eta)*F(ET) (closure test)
void scCorrections::run_apply_brEtaET()
{  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  //
  initHistograms();
  correctionType Corr  = brEtaET;
  fillHistograms(Corr);
  bool apply = true;
  getErecEgenFits_BrEta(Corr, apply);
}

//
void scCorrections::initHistograms()
{
  // Br
  for (Int_t i = 0; i<nBinsEta; ++i){
    for (Int_t j = 0; j<nBinsBr; ++j){
      h_ErecEgen[i][j] = new TH1F(Form("h_ErecEgen_%d_%d",i,j),Form("h_ErecEgen_%d_%d",i,j),nBinsErecEgen,ErecEgenMin, ErecEgenMax);
    }
  }
  for (Int_t i = 0; i<nBinsEta; ++i){
    h_CBetabr[i] = new TH1F(Form("h_CBetabr_%d",i), Form("h_CBetabr_%d",i),nBinsBr*2-1, brbins);
  }
  // ET - Br
  for (Int_t i = 0; i<nBinsET; ++i){
    h_ETcorrBrEta_EB[i] = new TH1F(Form("h_ETcorrBrEta_EB_%d",i), Form("h_ETcorrBrEta_EB_%d",i), nBinsErecEgen,ErecEgenMin, ErecEgenMax);
    h_ETcorrBrEta_EE[i] = new TH1F(Form("h_ETcorrBrEta_EE_%d",i), Form("h_ETcorrBrEta_EE_%d",i), nBinsErecEgen,ErecEgenMin, ErecEgenMax);
  }
  cout << "initHistograms(): done" << endl;
}

//
void scCorrections::fillHistograms(correctionType applyCorrections )
{
  // read the Correction parameters
  //
  // f(sigma_phi/sigma_eta, eta ) from histogram 
  //
  TFile *histFile_brEta = 0;
  TH1F *h_corr[nBinsEta];
  if (applyCorrections == brEta || applyCorrections == brEtaET){
    histFile_brEta = new TFile("./histCorrections_brEta.root","READ");
    if (histFile_brEta->IsZombie()) {cout << "ERROR: ./histCorrections_brEta.root does not exist " << endl; exit(-1);}
    //
    for (Int_t i = 0; i<nBinsEta; ++i){
      h_corr[i] = (TH1F*)histFile_brEta->Get(Form("h_corr_%d",i));
    }
  }  
  //
  // F(ET) correction from histogram ( *** BARREL / ENDCAP *** )
  //
  TFile *histFile_ET = 0;
  TH1F *h_corrET_EB = 0;
  TH1F *h_corrET_EE = 0;
  if (applyCorrections == brEtaET){
    histFile_ET = new TFile("./histCorrections_ET.root","READ");
    if (histFile_ET->IsZombie()) {cout << "ERROR: ./histCorrections_ET.root does not exist " << endl; exit(-1);}
    //    
    h_corrET_EB = (TH1F*)histFile_ET->Get("h_CBET_EB");
    h_corrET_EE = (TH1F*)histFile_ET->Get("h_CBET_EE");
  }

  // To speed up the loop read only the used variables from the NTuple
  // Add ALL variables used, else they are set to 0 and the compiler won't notice it
  //
  if (fChain == 0) return;
  fChain->SetBranchStatus("*",1);

  // Loop over all entries:
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry % 1000 == 0) cout << jentry << endl;

    // Loop over the reconstructed objects (electrons or photons)
    Int_t nobj = 0;
    if      (partType == "electron") nobj = NEles;
    else if (partType == "photon"  ) nobj = NPhotons;

    // skip pathological events
    if (nobj < 1 || nobj > 2) continue;

    {
      bool goodevt = true;
      if      (partType == "electron") { 
	if (ElSCindex[0]==-1 || ElSCindex[1]==-1) goodevt=false;
	if (ElGenE[0]<0 || ElGenE[1]<0) goodevt=false;
      }
      else if (partType == "photon"  ) { 
	if (PhotSCindex[0]==-1 || PhotSCindex[1]==-1) goodevt=false;
	if (PhoMCmatchexitcode[0]<=0 || PhoMCmatchexitcode[1]<=0) goodevt=false;
      }
      if (!goodevt) continue;
    }
    
    for (Int_t iobj = 0; iobj < nobj; ++iobj) {
      
      Double_t nt_em_eta        = 0;     
      Double_t nt_em_phi        = 0;
      Double_t nt_emCorrEta_e   = 0;
      Double_t nt_emCorrEta_et  = 0;
      Double_t nt_em_br1        = 0;
      Double_t nt_em_r9         = 0;
      Double_t nt_mc_e	        = 0; 
      
      // Use the correct naming of your ntuple
      if (partType == "photon"){
	int scindex = PhotSCindex[iobj];
	nt_em_eta        = SCEta[scindex];        
	nt_em_phi        = SCPhi[scindex];        
	nt_emCorrEta_e   = (fabs(nt_em_eta)<1.4442) ? SCRaw[scindex]*getEtaCorrectionBarrel(SCEta[scindex]) : SCRaw[scindex]+SCPre[scindex]; 
	nt_emCorrEta_et  = nt_emCorrEta_e/cosh(nt_em_eta);
	nt_em_br1        = SCBrem[scindex];     
	nt_em_r9         = PhoR9[iobj]; 
	nt_mc_e          = GenPhotonPt[PhoMCmatchindex[iobj]]*cosh(GenPhotonEta[PhoMCmatchindex[iobj]]);
	//
	// Cheap truth match:
	// the photons/electrons are back to back, just need to compare one reco-photon to the two mc-photons
	// 	Float_t dr0 = TMath::Sqrt(pow(phoSCEta[0]-mcEta[0],2) + pow(phoSCPhi[0]-mcPhi[0],2));
	// 	Float_t dr1 = TMath::Sqrt(pow(phoSCEta[0]-mcEta[1],2) + pow(phoSCPhi[0]-mcPhi[1],2));
	// 	if (dr0 < dr1 ) nt_mc_e = mcE[0];            
	// 	else nt_mc_e = mcE[1];  
	// 	cout << phoSCEta[0] << " " << mcEta[0] <<  " " << phoSCPhi[0] << " " << mcPhi[0] <<  " " << mcE[0]  << endl;
	// 	cout << phoSCEta[0] << " " << mcEta[1] <<  " " << phoSCPhi[0] << " " << mcPhi[1] <<  " " << mcE[1]  << endl;
	// 	cout << nt_mc_e << endl;
      }
      else if (partType == "electron"){
	int scindex = ElSCindex[iobj];
	nt_em_eta        = SCEta[scindex];        
	nt_em_phi        = SCPhi[scindex];        
	nt_emCorrEta_e   = (fabs(nt_em_eta)<1.4442) ? SCRaw[scindex]*getEtaCorrectionBarrel(SCEta[scindex]) : SCRaw[scindex]+SCPre[scindex]; 
	nt_emCorrEta_et  = nt_emCorrEta_e/cosh(nt_em_eta);
	nt_em_br1        = SCBrem[scindex];     
	nt_em_r9         = SCR9[scindex]; 
	nt_mc_e          = ElGenE[iobj];
	//
	// Cheap truth match:
	// the photons/electrons are back to back, just need to compare one reco-photon to the two mc-photons
	// 	Float_t dr0 = TMath::Sqrt(pow(eleSCEta[0]-mcEta[0],2) + pow(eleSCPhi[0]-mcPhi[0],2));
	// 	Float_t dr1 = TMath::Sqrt(pow(eleSCEta[0]-mcEta[1],2) + pow(eleSCPhi[0]-mcPhi[1],2));
	// 	if (dr0 < dr1 ) nt_mc_e = mcE[0];            
	// 	else nt_mc_e = mcE[1];  
	// 	cout << eleSCEta[0] << " " << mcEta[0] <<  " " << eleSCPhi[0] << " " << mcPhi[0] <<  " " << mcE[0]  << endl;
	// 	cout << eleSCEta[0] << " " << mcEta[1] <<  " " << eleSCPhi[0] << " " << mcPhi[1] <<  " " << mcE[1]  << endl;
	// 	cout << nt_mc_e << endl;
      }

      // main event selection
      //---------------------------------
      if ( 
	  // select showering objects: 
	  ( (TMath::Abs(nt_em_eta) < etaCrackMin && (nt_em_r9 < minR9_EB)) || (TMath::Abs(nt_em_eta) > etaCrackMax && (nt_em_r9 < minR9_EE)) ) 
	  
	  &&
	  
	  // remove the lowest ET objects
	  nt_emCorrEta_et > minETCut    
	  
	  && 
	  
	  // remove the phi gaps between the super modules (only for the barrel):: local corrections will take care of them
	  // for the endcap there are no local corrections, so don't exclude any object
	  !isInPhiCracks(nt_em_phi, nt_em_eta,  fiducialPhiCut) 	  
	  
	   ){ 
	
	// Find the correct (br,eta) bin for the h_ErecEgen[nBinsEta][nBinsBr] and h_ETcorrBrEta_EB[nBinsET];
	Int_t iEta    = -1;
	Int_t iBr     = -1;
	Int_t iET     = -1;
	Int_t iETcorr = -1;	
	//
	for (Int_t iiEta = 0; iiEta < nBinsEta; ++iiEta){ 
	  if ( leftEta[iiEta] <= TMath::Abs(nt_em_eta) && TMath::Abs(nt_em_eta) <rightEta[iiEta]) {
	    iEta = iiEta; 
	    // cout << "ETA " << leftEta[iiEta]  << " " << TMath::Abs(nt_em_eta) << " " << rightEta[iiEta] << " " << iEta << endl;
	  }
	}
	for (Int_t iiBr  = 0; iiBr  < nBinsBr;  ++iiBr ){ 
	  if ( leftBr [iiBr]  <= nt_em_br1             &&             nt_em_br1 <rightBr [iiBr] ) {
	    iBr  = iiBr;  
	    // cout << "BR " << leftBr[iiBr]  << " " << nt_em_br1 << " " << rightBr[iiBr] << " " << iBr << endl;
	  }
	}
	// Find the correct (ET) bin for the h_ETcorrBrEta_EX[nBinsET]
	for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 
	  if ( leftET [iiET]  <= nt_emCorrEta_et      &&       nt_emCorrEta_et < rightET[iiET] ) 
	    iET  = iiET;  
	  // cout << "ET " << leftBr[iiET]  << " " << nt_emCorrEta_et << " " << rightBr[iiET] << " " << nt_emCorrEta_et << endl;
	}		
	if (iEta == -1 || iBr == -1 || iET == -1) continue;


	// Initialize the correction factors
	Double_t corrBrEta = 1.; 
	Double_t corrEt    = 1.; 	
	
 	// get the corrections f(br,eta) and F(ET) from histogram                                                                                         //MDDB
 	//																		  //MDDB
 	if (applyCorrections == brEta || applyCorrections == brEtaET) {											  //MDDB
 	  corrBrEta  = h_corr[iEta]->GetBinContent(2*iBr+1);                               								  //MDDB
 	  //																		  //MDDB
 	  // Find the correct (ET corrected with f(br,Eta) ) bin for the correction histogram 								  //MDDB
 	  for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 												  //MDDB
 	    if ( leftET [iiET]  <= (nt_emCorrEta_et / corrBrEta)      &&       (nt_emCorrEta_et / corrBrEta) < rightET[iiET] ) {			  //MDDB
 	      iETcorr  = iiET;  															  //MDDB
 	      // cout << "ET " << leftET [iiET] << " " <<  nt_emCorrEta_et / corrBrEta  << " " <<  rightET[iiET] <<  " " << iETcorr <<endl;		  //MDDB
 	    }																		  //MDDB
 	  }																		  //MDDB
 	  if (iETcorr == -1) continue; // always skip the event, there is something patological	    							  //MDDB
 	  if (applyCorrections == brEtaET) {														  //MDDB
 	    //																		  //MDDB
 	    // histogram correction ET *** BARREL, ENDCAP ***												  //MDDB
 	    if (TMath::Abs(nt_em_eta) < etaCrackMin ){                    										  //MDDB
 	      corrEt = h_corrET_EB    ->GetBinContent(2*iETcorr+1);											  //MDDB
 	      //	    cout << " corr " << iETcorr << " " << h_corrET_EB    ->GetBinContent(2*iETcorr+1) <<endl;					  //MDDB
 	    }																		  //MDDB
 	    if (TMath::Abs(nt_em_eta) > etaCrackMax ) {													  //MDDB
 	      corrEt = h_corrET_EE    ->GetBinContent(2*iETcorr+1);											  //MDDB
 	      //	    cout << "corr " << iETcorr << " " << h_corrET_EE    ->GetBinContent(2*iETcorr+1) <<endl;					  //MDDB
 	    }																		  //MDDB
 	  }																		  //MDDB
 	}																		  //MDDB

//MDDB 	// get the corrections f(br,eta) and F(ET) from File                                                                                                  
//MDDB 	//																		      
//MDDB 	if (applyCorrections == brEta || applyCorrections == brEtaET) {											      
//MDDB 	  corrBrEta = applyScCorrectionsBrEta(nt_em_eta,  nt_em_br1);                           							      
//MDDB 	  																		      
//MDDB 	  // Find the correct (ET corrected with f(br,Eta) ) bin for the correction histogram 								      
//MDDB 	  for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 												      
//MDDB 	    if ( leftET [iiET]  <= (nt_emCorrEta_et / corrBrEta)      &&       (nt_emCorrEta_et / corrBrEta) < rightET[iiET] ) {			      
//MDDB 	      iETcorr  = iiET;  															      
//MDDB 	      // cout << "ET " << leftET [iiET] << " " <<  nt_emCorrEta_et / corrBrEta  << " " <<  rightET[iiET] <<  " " << iETcorr <<endl;		      
//MDDB 	    }																		      
//MDDB 	  }																		      
//MDDB 	  if (iETcorr == -1) continue; // always skip the event, there is something patological	   							      
//MDDB 	  if (applyCorrections == brEtaET) {														      
//MDDB 	    if (TMath::Abs(nt_em_eta) < etaCrackMin ) corrEt = applyScCorrectionsET_EB( nt_emCorrEta_et / corrBrEta );   				      
//MDDB 	    if (TMath::Abs(nt_em_eta) > etaCrackMax ) corrEt = applyScCorrectionsET_EE( nt_emCorrEta_et / corrBrEta );   				      
//MDDB 	  }																		      
//MDDB 	}																		      

	// fill the matrix of histograms h_ErecEgen[nBinsEta][nBinsBr] 
	// ===========================================================
	//
	// *** f(sigma_phi/sigma_eta,eta) ***
	if (applyCorrections == brEta) {
	  // fill
	  if (nt_mc_e!=0 && corrBrEta!=0) h_ErecEgen[iEta][iBr]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta);	  
	}	
	// *** f(sigma_phi/sigma_eta,eta)*F(ET) ***
	else if (applyCorrections == brEtaET) {
	  // fill
	  if (nt_mc_e!=0 && corrEt!=0) h_ErecEgen[iEta][iBr]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta /corrEt);	  
	}	
	// *** no corrections ***
	else { 
	  // fill
	  if (nt_mc_e!=0 )  h_ErecEgen[iEta][iBr]->Fill(nt_emCorrEta_e / nt_mc_e );     
	}
	
	// fill the vector of ET dependence histograms
	// ===========================================
	//	
	// *** f(sigma_phi/sigma_eta,eta) ***
	if (applyCorrections == brEta) {
	  if (TMath::Abs(nt_em_eta) < etaCrackMin ){
	    // fill
	    if (nt_mc_e!=0 && corrBrEta!=0)    h_ETcorrBrEta_EB[iET]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta);
	  }
	  else if (TMath::Abs(nt_em_eta) > etaCrackMax ) {
	    // fill
	    if (nt_mc_e!=0 && corrBrEta!=0)    h_ETcorrBrEta_EE[iET]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta);		
	  }
	}
	// *** f(sigma_phi/sigma_eta,eta)*F(ET) ***
	else if (applyCorrections == brEtaET) {		
	  // bin position for the br,eta corrected ET
	  if (TMath::Abs(nt_em_eta) < etaCrackMin ){
	    
	    // fill
	    if (nt_mc_e!=0 && corrEt!=0) h_ETcorrBrEta_EB[iETcorr]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta / corrEt );	   
	  }
	  if (TMath::Abs(nt_em_eta) > etaCrackMax ) {
	    // fill
	    if (nt_mc_e!=0 && corrEt!=0) h_ETcorrBrEta_EE[iETcorr]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta / corrEt);	   
	  }	      
	}	    
	// *** no corrections ***
	else{
	  if (TMath::Abs(nt_em_eta) < etaCrackMin ){
	    h_ETcorrBrEta_EB[iET]->Fill(nt_emCorrEta_e / nt_mc_e );
	  }
	  if (TMath::Abs(nt_em_eta) > etaCrackMax ) { 
	    h_ETcorrBrEta_EE[iET]->Fill(nt_emCorrEta_e / nt_mc_e );
	  }
	}

      } // if selected event
    } // loop over obj
  } // loop over entries	
  
  if (applyCorrections == brEta || applyCorrections == brEtaET) histFile_brEta->Close();    
  if (applyCorrections == brEtaET) histFile_ET->Close();
  
  return;
  
}


void scCorrections::fitCB(TH1F *h, Double_t &mean, Double_t &emean, Double_t &chi2, Double_t& badReco  )
{
  // This is to quantify the fraction of events where the reconstructed energy is below EgenFrac
  badReco = h->Integral(1,h->FindBin(EgenFrac)) / (Double_t)h->GetEntries();
  
  // driving gaussian. Used to estimate the fit range
  TF1 *gtmp = 0;  
  if (plots){
    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    c1->cd();
    h->Draw();

    // more solid peak position with a temporary histogram filled only with the bins above 50% of the peak max
    TH1F *htmp = new TH1F("htmp","htmp",nBinsErecEgen,ErecEgenMin, ErecEgenMax);
    // protect the fits from small statistics rebinning
    if (h->GetMaximum()< 50.) { h->Rebin(2); htmp->Rebin(2); }
    if (h->GetMaximum()< 25.) { h->Rebin(2); htmp->Rebin(2); }
    for (Int_t i=0; i<h->GetNbinsX(); ++i) {if (h->GetBinContent(i) >h->GetMaximum()/2.) htmp->SetBinContent(i,h->GetBinContent(i));}
    htmp->Draw();
    c1->Update();
    getchar();
    
    // driving gaussian
    gtmp = new TF1("gtmp","gaus(0)",0,2.0);    
    // use the mean and RMS of the Erec/Egen as initial parameters for the driving gaussian
    // gtmp->SetParameters(10,1,0.05);
    cout << h->GetMaximum() << " " <<  h->GetMean() << " " <<  h->GetRMS() << " " <<  htmp->GetMean() << " " <<  htmp->GetRMS() << endl;
    gtmp->SetParameters(10, htmp->GetMean(), htmp->GetRMS());
    h->Fit("gtmp","","",htmp->GetMean()-5*htmp->GetRMS(),htmp->GetMean()+5*htmp->GetRMS()); // 0.8,1.1);
    h->Fit("gtmp","","",
	   gtmp->GetParameter(1)-2*TMath::Abs(gtmp->GetParameter(2)),
	   gtmp->GetParameter(1)+2*TMath::Abs(gtmp->GetParameter(2)));
    cout << "***** GTMP " 
	 << gtmp->GetParameter(0) << " " << gtmp->GetParameter(1) << " " << gtmp->GetParameter(2) 
	 << endl;
    c1->Update();
    getchar();
    delete htmp;
  }
  else{ 
    // more solid peak position with a temporary histogram filled only with the bins above 50% of the peak max
    TH1F *htmp = new TH1F("htmp","htmp",nBinsErecEgen,ErecEgenMin, ErecEgenMax);
    // protect the fits from small statistics rebinning
    if (h->GetMaximum()< 50.) { h->Rebin(2); htmp->Rebin(2); }
    if (h->GetMaximum()< 25.) { h->Rebin(2); htmp->Rebin(2); }
    for (Int_t i=0; i<h->GetNbinsX(); ++i) {if (h->GetBinContent(i) >h->GetMaximum()/2.) htmp->SetBinContent(i,h->GetBinContent(i));}
    htmp->Draw();



    // driving gaussian
    gtmp = new TF1("gtmp","gaus(0)",ErecEgenMin, ErecEgenMax);
    // use the mean and RMS of the Erec/Egen as initial parameters for the driving gaussian
    //    gtmp->SetParameters(10,1,0.01);
    cout << h->GetMaximum() << " " <<  h->GetMean() << " " <<  h->GetRMS() << " " <<  htmp->GetMean() << " " <<  htmp->GetRMS() << endl;
    gtmp->SetParameters(10, htmp->GetMean(), htmp->GetRMS());
    h->Fit("gtmp","","",htmp->GetMean()-5*htmp->GetRMS(),htmp->GetMean()+5*htmp->GetRMS()); // 0.8,1.1);
    h->Fit("gtmp","","",
	   gtmp->GetParameter(1)-2*TMath::Abs(gtmp->GetParameter(2)),
	   gtmp->GetParameter(1)+2*TMath::Abs(gtmp->GetParameter(2)));
    cout << "***** GTMP " 
	 << gtmp->GetParameter(0) << " " << gtmp->GetParameter(1) << " " << gtmp->GetParameter(2) 
	 << endl;
     delete htmp;
 }
     
  // RooDataHist from TH1F
  RooRealVar  x("x","x", ErecEgenMin, ErecEgenMax);  
  RooDataHist data("data","Ereco/Egen",x,h); 
  RooPlot*    xframe = x.frame(Name("xframe"),Title("E^{reco} / E^{gen}")) ;
  data.plotOn(xframe);

  // Initialize CB parameters
  Double_t  alphaStart =   0.5;
  Double_t  alphaMin   =   0.1;
  Double_t  alphaMax   =   2.0;
  Double_t  nStart     =  50.0;
  Double_t  nMin       =   0.0;
  Double_t  nMax       = 200.0;
  if ( gtmp->GetParameter(1) < 0.95 && gtmp->GetParameter(2) > 0.03 ) {
    alphaStart =   0.4;
    alphaMin   =   0.;
    alphaMax   =   1.0;
    nStart     =  20.;
    nMin       =   0.;
    nMax       = 100.;
  }
  
  // fit function
  RooRealVar alpha  ("alpha"  ,        "alpha" , alphaStart, alphaMin, alphaMax); 
  RooRealVar n      ("n"      ,            "n" , nStart, nMin, nMax); 
  RooRealVar cbmean ("cbmean" ,       "cbmean" , gtmp->GetParameter(1) ,0.5, 1.2); // ErecEgenMin, ErecEgenMax);
  RooRealVar cbsigma("cbsigma",      "cbsigma" , 0.01, 0.001,0.1) ;
  RooCBShape cball  ("cball"  , "crystal ball" , x, cbmean, cbsigma, alpha, n);
  
  // other functions
  // CB (x) Gaussian
  //   RooRealVar mg("mg","mg",1) ;
  //   RooRealVar sg("sg","sg",0.1,0.01,0.2) ;
  //   RooGaussian gauss("gauss","gauss",x,mg,sg) ; 
  //   RooFFTConvPdf CBgaus("CBgaus","CB (X) gauss",x,cball,gauss);
  //
  // CB + Gaussian
  //   RooRealVar ng("ng","ng",100,1,10000);
  //   RooRealVar nc("nc","nc",10, 1,10000);
  //   RooAddPdf CBgaus("CBgaus","CB + gauss",RooArgList(cball,gauss),RooArgList(nc,ng));
  
  // Fit Range using as units the sigma of the driving gaussian
  //
  Float_t nsigmaL = 10;
  Float_t nsigmaR = 2;
  //
  // Fit      
  RooFitResult* fitres =cball.fitTo(data,
				    Range(gtmp->GetParameter(1)-nsigmaL*TMath::Abs(gtmp->GetParameter(2)),
					  gtmp->GetParameter(1)+nsigmaR*TMath::Abs(gtmp->GetParameter(2))),
				    Minos(kFALSE));    
  mean  = cbmean.getVal();
  emean = cbmean.getError();
  cball.plotOn (xframe,LineColor(kBlue));
  cball.paramOn(xframe,Layout(0.12,0.64,0.88));
  data. statOn (xframe,Layout(0.12,0.64,0.55));
  chi2 = xframe->chiSquare(4); // dof = 4 = number fof floating parametes
  //  cout << "***** CHI2 " << xframe->chiSquare(4) << endl; getchar();


  // plot
  //  gStyle->SetOptFit(1);
  if (plots){
    TCanvas *cpull = new TCanvas("cpull","cpull",600,700);
    TPad *cpull_1 = new TPad("cpull_1", "cpull_1",0,0.,1,1);
    cpull_1->Draw();
    cpull_1->cd();
    cpull_1->SetFillColor(0);
    cpull_1->SetBorderMode(0);
    cpull_1->SetBorderSize(2);
    cpull_1->SetFrameBorderMode(0);
    cpull_1->SetFrameBorderMode(0);
    cpull_1->SetBottomMargin(0.25);
    TString chi2txt = "#chi^{2}=" + floatToString(xframe->chiSquare(4)) ; // dof = 4 = number dof floating parametes
    TLatex* txt = new TLatex(0.72,10,chi2txt) ;
    txt->SetTextSize(0.04) ;
    txt->SetTextColor(kRed) ;
    xframe->addObject(txt) ;
    xframe->GetXaxis()->SetTitle("E^{reco}/E^{gen}");
    xframe->Draw();
    TPad *cpull_2 = new TPad("cpull_2", "cpull_2",0,0.,1.0,0.17);
    cpull_2->Draw();
    cpull_2->cd();
    cpull_2->SetFillColor(0);
    cpull_2->SetBorderMode(0);
    cpull_2->SetBorderSize(0);
    cpull_2->SetBottomMargin(0.45);
    cpull_2->SetFrameBorderMode(0);
    cpull_2->SetFrameBorderMode(0);
    // pulls    
    //     RooHist* hpull      = xframe->pullHist() ;
    //     RooPlot* xframePull = x.frame(Title("Pull Distribution")) ;
    //     xframePull->addPlotable(hpull,"P") ;
    //     xframePull->GetXaxis()->SetTitle("E^{reco}/E^{gen}");
    //     xframePull->GetXaxis()->SetTitleSize(0.2);
    //     xframePull->GetXaxis()->SetLabelSize(0.2);
    //     xframePull->GetYaxis()->SetRangeUser(-3,3);
    //     xframePull->GetYaxis()->SetTitle("Pull");
    //     xframePull->GetYaxis()->SetTitleSize(0.2);
    //     xframePull->GetYaxis()->SetLabelSize(0.1);
    //     xframePull->Draw();
    cpull->Update();  
    getchar();
  }      
  return;
}

void scCorrections::deriveCorrections_BrEta()
{
  // fit the eta slices
  //
  //   TCanvas *ctmp = new TCanvas("ctmp","ctmp",600,600);
  //   ctmp->cd();
  //   TF1 *f[nBinsEta];
  //   for (Int_t i = 0; i<nBinsEta ; ++i){
  
  //     //polynomials
  //     f[i] = new TF1("f","pol2",0.,10.);
  //     f[i]->SetParameters(1,1,-1);    
  //     h_CBetabr[i]->Fit("f","","same",1.5,6.0);    
  
  //   }
  //   delete ctmp;
  //   // save the fit parameters
  //   ofstream outfile;  
  //   outfile.open ("fitParametes.txt");
  //   for (Int_t i = 0; i<nBinsEta ; ++i){
  //     outfile << f[i]->GetParameter(0) << " " << f[i]->GetParameter(1) << " " << f[i]->GetParameter(2) << endl;
  //   }
  //   outfile.close();
  
  // Write out the correction macro
  ofstream outfile2;  
  TString filename = "./applyScCorrections.C";
  outfile2.open(filename);
  outfile2 << "#ifndef scCorrectionsBrETAET_h                                           " << endl;
  outfile2 << "#define scCorrectionsBrETAET_h                                           " << endl;
  outfile2 << "#include \"TFile.h\"                                                     " << endl;
  outfile2 << "#include \"TH1F.h\"				                        " << endl;
  outfile2 << "#include <iostream>  			                                " << endl;
  outfile2 << "  			                                                " << endl;
  outfile2 << "bool DBG = false;			                                " << endl;
  outfile2 << "  	                                                                " << endl;
  outfile2 << "Double_t applyScCorrectionsBrEta(Double_t eta, Double_t sigmaPhiSigmaEta){  " << endl;
  outfile2 << "  	                                                                " << endl;
  outfile2 << "  const Int_t    nBinsEta = " << nBinsEta << ";                          " << endl;                                  
  outfile2 << "  Double_t       leftEta  ["  << nBinsEta << "];                         " << endl;
  outfile2 << "  Double_t       rightEta ["  << nBinsEta << "];                         " << endl;
  outfile2 << "  const Int_t    nBinsBr  = " << nBinsBr << ";                           " << endl;                                    
  outfile2 << "  Double_t       leftBr  [" << nBinsBr << "];                            " << endl;
  outfile2 << "  Double_t       rightBr [" << nBinsBr << "];                            " << endl;
  outfile2 << "  Double_t brbins  [" << nBinsBr*2 << "];                                " << endl;
  outfile2 << "  TH1F *h_corr["<< nBinsEta << "];                                       " << endl;
  outfile2 << "  	                                                                " << endl;
  for (Int_t i = 0; i< nBinsEta; ++i){ 
    outfile2 << "  leftEta[" << i << "] =  " << leftEta[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsEta; ++i){ 
    outfile2 << "  rightEta[" << i<< "] =  " << rightEta[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsBr; ++i){ 
    outfile2 << "  leftBr[" << i << "] =  " << leftBr[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsBr; ++i){ 
    outfile2 << "  rightBr["<< i << "] =  " << rightBr[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< 2*nBinsBr; ++i){ 
    outfile2 << "  brbins[" << i << "] =  " << brbins[i] << " ; " << endl;
  }
  outfile2 << "  for (Int_t i = 0; i<nBinsEta; ++i){                                                   " << endl;
  outfile2 << "    h_corr[i] = new TH1F(Form(\"h_corr_%d\",i),Form(\"h_corr_%d\",i),nBinsBr*2-1, brbins);  " << endl;
  outfile2 << "  }                                                                                     " << endl;
  for (Int_t i = 0; i<nBinsEta ; ++i){
    for (Int_t j = 0; j< nBinsBr; ++j){
      outfile2 << "  h_corr[" << i << "]->SetBinContent(" <<  2*j+1 << "," << h_CBetabr[i]->GetBinContent(2*j+1) << " );"<< endl;
    }
    outfile2 << endl;
  }


  outfile2 << "     Int_t tmpEta = -1;                                                                                                                                                                        " << endl;
  outfile2 << "     Int_t tmpBr  = -1;																					      " << endl;
  outfile2 << "   																							      " << endl;
  outfile2 << "     // Extra protections										       										      " << endl;
  outfile2 << "      																							      " << endl;
  outfile2 << "     if (TMath::Abs(eta)  <   leftEta[0]           ) { tmpEta = 0;          if (DBG) std::cout << \" WARNING [applyScCorrections]: (TMath::Abs(eta)  <   leftEta[0]          \" << std::endl;}   " << endl;
  outfile2 << "     if (TMath::Abs(eta)  >=  rightEta[nBinsEta-1] ) { tmpEta = nBinsEta-1; if (DBG) std::cout << \" WARNING [applyScCorrections]: TMath::Abs(eta)  >=  rightEta[nBinsEta-1] \" << std::endl;}   " << endl;
  outfile2 << "   																							      " << endl;
  outfile2 << "     if (sigmaPhiSigmaEta <  leftBr [0]            ) {tmpBr = 0;            if (DBG) std::cout << \" WARNING [applyScCorrections]: sigmaPhiSigmaEta <  leftBr [0]            \" << std::endl;}   " << endl;
  outfile2 << "     if (sigmaPhiSigmaEta >= rightBr[nBinsBr]      ) {tmpBr = nBinsBr  -1;  if (DBG) std::cout << \" WARNING [applyScCorrections]: sigmaPhiSigmaEta >= rightBr[nBinsBr]      \" << std::endl;}   " << endl;    
  outfile2 << "                                                                                                                                                                                               " << endl;
  outfile2 << "     for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								              								      " << endl;
  outfile2 << "       if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] ){				       									      " << endl;
  outfile2 << "         tmpEta = iEta;											       										      " << endl;
  outfile2 << "       }													       										      " << endl;
  outfile2 << "     }													         									      " << endl;
  outfile2 << "   																							      " << endl;
  outfile2 << "     for (Int_t iSigmaPhiSigmaEta = 0; iSigmaPhiSigmaEta < nBinsBr; ++iSigmaPhiSigmaEta){			       									      " << endl;
  outfile2 << "       if (leftBr [iSigmaPhiSigmaEta]  <= sigmaPhiSigmaEta && sigmaPhiSigmaEta <rightBr [iSigmaPhiSigmaEta] ){      									      " << endl;
  outfile2 << "         tmpBr = iSigmaPhiSigmaEta;										       									      " << endl;
  outfile2 << "       }													       										      " << endl;
  outfile2 << "     }													       										      " << endl;
  outfile2 << "     																							      " << endl;
  outfile2 << "     // Interpolation																					      " << endl;
  outfile2 << "     Double_t tmpInter = 1;																				      " << endl;
  outfile2 << "     // In eta cracks/gaps 																				      " << endl;
  outfile2 << "     if (tmpEta == -1 && tmpBr != -1 ) { // need to interpolate only eta, if br is out of range skip this										      " << endl;
  outfile2 << "       																							      " << endl;
  outfile2 << "       if (TMath::Abs(eta) >=rightEta[nBinsEta-1] ) { // out of ECAL boundary														      " << endl;
  outfile2 << "         for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];    															      " << endl;
  outfile2 << "         return 1; // don't correct																			      " << endl;
  outfile2 << "       }																							      " << endl;
  outfile2 << "       																							      " << endl;
  outfile2 << "       // central bin at eta = 0																				      " << endl;
  outfile2 << "       if (0 <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[0] ) {															      " << endl;
  outfile2 << "         																						      " << endl;
  outfile2 << "         tmpInter = h_corr[0]->GetBinContent(2*tmpBr+1);																	      " << endl;
  outfile2 << "         																						      " << endl;
  outfile2 << "       }																							      " << endl;
  outfile2 << "       else { // all other crack-bins																			      " << endl;
  outfile2 << "         																						      " << endl;
  outfile2 << "         for (Int_t iEta = 0; iEta < nBinsEta-1; ++iEta){								       								      " << endl;
  outfile2 << "   	if (rightEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[iEta+1] ){													      " << endl;
  outfile2 << "   	  tmpInter = ( h_corr[iEta]  ->GetBinContent(2*tmpBr+1) + 															      " << endl;
  outfile2 << "   		       h_corr[iEta+1]->GetBinContent(2*tmpBr+1) ) / 2. ;														      " << endl;
  outfile2 << "   	}																						      " << endl;
  outfile2 << "         }																						      " << endl;
  outfile2 << "       }																							      " << endl;
  outfile2 << "       for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];    															      " << endl;
  outfile2 << "       return tmpInter;																					      " << endl;
  outfile2 << "     }  																							      " << endl;
  outfile2 << "     // end interpolation                                                                                                                                                                      " << endl;

//   outfile2 << "   Int_t tmpEta = -1;                                                                                             " << endl;         
//   outfile2 << "   for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								       " << endl;
//   outfile2 << "     if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] ){				       " << endl;
//   outfile2 << "       tmpEta = iEta;											       " << endl;
//   outfile2 << "     }													       " << endl;
//   outfile2 << "   }													       " << endl;
//   outfile2 << "   													       " << endl;
//   outfile2 << "   Int_t tmpBr = -1;											       " << endl;
//   outfile2 << "   for (Int_t iSigmaPhiSigmaEta = 0; iSigmaPhiSigmaEta < nBinsBr; ++iSigmaPhiSigmaEta){			       " << endl;
//   outfile2 << "     if (leftBr [iSigmaPhiSigmaEta]  <= sigmaPhiSigmaEta && sigmaPhiSigmaEta <rightBr [iSigmaPhiSigmaEta] ){      " << endl;	   
//   outfile2 << "       tmpBr = iSigmaPhiSigmaEta;										       " << endl;
//   outfile2 << "     }													       " << endl;
//   outfile2 << "   }													       " << endl;

  outfile2 << "   if (tmpEta == -1 || tmpBr == -1){									       " << endl;
  outfile2 << "     for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  		       " << endl;
  outfile2 << "     return  1;	// don't correct												       " << endl;
  outfile2 << "   }													       " << endl;
  outfile2 << "   													       " << endl;
  outfile2 << "   Double_t tmp = h_corr[tmpEta]->GetBinContent(2*tmpBr+1);				   		       " << endl;
  outfile2 << "   for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  		       " << endl;
  outfile2 << "   return  tmp;                                                                                                   " << endl;                         
  outfile2 << "}													   " << endl;  

  //  outfile2 << "  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								   " << endl;
  //  outfile2 << "    for (Int_t iSigmaPhiSigmaEta = 0; iSigmaPhiSigmaEta < nBinsBr; ++iSigmaPhiSigmaEta){			   " << endl;
  //  outfile2 << "      													   " << endl;
  //  outfile2 << "      // select the right (eta,br) bin in the matrix							   " << endl;
  //  outfile2 << "      if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] &&			   " << endl;
  //  outfile2 << " 	   leftBr [iSigmaPhiSigmaEta]  <= sigmaPhiSigmaEta && sigmaPhiSigmaEta <rightBr [iSigmaPhiSigmaEta] ){	   " << endl;
  //  outfile2 << "	       Double_t tmp = h_corr[iEta]->GetBinContent(2*iSigmaPhiSigmaEta+1);				   " << endl;	   
  //  outfile2 << "	       for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  " << endl;
  //  outfile2 << "	       return  tmp;                                                                                        " << endl;
  //  outfile2 << "      }													   " << endl;
  //  outfile2 << "    } 													   " << endl;
  //  outfile2 << "  }													   " << endl;
  //  outfile2 << "	       for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  " << endl;
  //  outfile2 << "  return 1;												   " << endl;
  //  outfile2 << "}													   " << endl;
  
  outfile2.close();                                                                                                      
  
  // save also the histograms to be used as correction
  TFile *histFile= new TFile("./histCorrections_brEta.root", "RECREATE" );
  TH1F *h_corr[nBinsEta];
  for (Int_t i = 0; i<nBinsEta ; ++i){
    h_corr[i] = (TH1F*)h_CBetabr[i]->Clone(Form("h_corr_%d",i));
    h_corr[i]->Write();
  }
  histFile->Close();
    
  TCanvas *cc_BadClustersEtaBr = new TCanvas("cc_BadClustersEtaBr" ,"cc_BadClustersEtaBr", 600,600);
  cc_BadClustersEtaBr->cd();
  gPad->SetPhi(270+30);
  hh_EgenFracEtaBr->GetXaxis()->SetTitle("#sigma_{#phi} / #sigma_{#eta}");
  hh_EgenFracEtaBr->GetXaxis()->SetTitleOffset(1.6);
  hh_EgenFracEtaBr->GetYaxis()->SetTitle("#eta");
  hh_EgenFracEtaBr->GetYaxis()->SetTitleOffset(1.6);
  hh_EgenFracEtaBr->Draw("lego2");
  TCanvas *c_BadClustersET_EB     = new TCanvas("c_BadClustersET_EB"     ,"c_BadClustersET_EB"    , 600,600);
  c_BadClustersET_EB->cd();
  h_EgenFracET_EB ->GetXaxis()->SetTitle("ET [GeV]");
  h_EgenFracET_EB ->Draw("");
  TCanvas *c_BadClustersET_EE     = new TCanvas("c_BadClustersET_EE"     ,"c_BadClustersET_EE"    , 600,600);
  c_BadClustersET_EE->cd();
  h_EgenFracET_EE ->GetXaxis()->SetTitle("ET [GeV]");
  h_EgenFracET_EE ->Draw("");

  return;
}

void scCorrections::deriveCorrections_BrEtaET()
{
  // Add to the existing applyCorrection the ET dependence
  ofstream outfile2;  
  TString filename = "./applyScCorrections.C";
  outfile2.open(filename,ios::app);

  outfile2 << "                                                                         " << endl;
  outfile2 << "Double_t applyScCorrectionsET_EB(Double_t ET){  							   " << endl;
  outfile2 << "  const Int_t    nBinsET             = " << nBinsET<< ";             " << endl; 
  outfile2 << "  Double_t       leftET  [nBinsET];                    " << endl;
  outfile2 << "  Double_t       rightET [nBinsET];                    " << endl;
  outfile2 << "  Double_t       ETBins  [nBinsET*2];                  " << endl;
  // ET
  for (Int_t i = 0; i< nBinsET; ++i){ 
    outfile2 << "  leftET[" << i << "] =  " << leftET[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsET; ++i){ 
    outfile2 << "  rightET[" << i << "] =  " << rightET[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< 2*nBinsET; ++i){ 
    outfile2 << "  ETBins[" << i << "] =  " << ETBins[i] << " ; " << endl;
  }
  outfile2 << "  TH1F *h_CBET_EB    = new TH1F(\"h_CBET_EB\"    ,\"h_CBET_EB\"    ,nBinsET*2-1, ETBins); " << endl;

  // this is the residual ET dependence after BrEta correction
  for (Int_t iET = 0; iET< nBinsET; ++iET){ 
    outfile2 << "  h_CBET_EB->SetBinContent(" << 2*iET+1 << ", " << h_CBET_EB->GetBinContent(2*iET+1) << "); " << endl;
  }

  outfile2 << "  Int_t iET      = -1;                                                                   " << endl;              

  outfile2 << "  for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 				      " << endl;
  outfile2 << "    if ( leftET [iiET]  <= (ET)      &&       (ET) < rightET[iiET] ) {		      " << endl;
  outfile2 << "      iET  = iiET;  								      " << endl;
  outfile2 << "    }										      " << endl;
  outfile2 << "  }										      " << endl;

  // Extra protections										       
  outfile2 << "  if (ET < leftET  [0] )         { iET = 0;              if (DBG) std::cout << \" WARNING [applyScCorrections]: ET < leftET  [0] )       \" << std::endl;}   " << endl;
  outfile2 << "  if (ET > rightET [nBinsET-1] ) { iET = nBinsET-1;      if (DBG) std::cout << \" WARNING [applyScCorrections]: ET > rightET [nBinsET-1] \" << std::endl;}   " << endl;

  outfile2 << "  if (iET == -1) {delete  h_CBET_EB; return 1;}					      " << endl;
  
  outfile2 << "											      " << endl;
  outfile2 << "  Int_t binET =  2*iET+1 ; // h_CBET_EB->FindBin(ET);                                    " << endl;       
  outfile2 << "  Double_t tmp = h_CBET_EB->GetBinContent(binET);                                        " << endl;
  outfile2 << "  delete h_CBET_EB;                                                             	      " << endl;
  outfile2 << "  if ( 0< binET && binET < 2*nBinsET+1) return tmp;                          	      " << endl;
  outfile2 << "  else return 1.;                                                               	      " << endl;
  outfile2 << "                                                                                         " << endl;
  outfile2 << "}                                                                               " << endl;                                                            
  outfile2 << "                                                                                  " << endl; 
    
  outfile2 << "                                                                         " << endl;
  outfile2 << "Double_t applyScCorrectionsET_EE(Double_t ET){  							   " << endl;
  outfile2 << "  const Int_t    nBinsET             = " << nBinsET<< ";             " << endl; 
  outfile2 << "  Double_t       leftET  [nBinsET];                    " << endl;
  outfile2 << "  Double_t       rightET [nBinsET];                    " << endl;
  outfile2 << "  Double_t       ETBins  [nBinsET*2];                  " << endl;
  // ET
  for (Int_t i = 0; i< nBinsET; ++i){ 
    outfile2 << "  leftET[" << i << "] =  " << leftET[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsET; ++i){ 
    outfile2 << "  rightET[" << i << "] =  " << rightET[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< 2*nBinsET; ++i){ 
    outfile2 << "  ETBins[" << i << "] =  " << ETBins[i] << " ; " << endl;
  }
  outfile2 << "  TH1F *h_CBET_EE    = new TH1F(\"h_CBET_EE\"    ,\"h_CBET_EE\"    ,nBinsET*2-1, ETBins); " << endl;

  // this is the residual ET dependence after BrEta correction
  for (Int_t iET = 0; iET< nBinsET; ++iET){ 
    outfile2 << "  h_CBET_EE->SetBinContent(" << 2*iET+1 << ", " << h_CBET_EE->GetBinContent(2*iET+1) << "); " << endl;
  }

  outfile2 << "  Int_t iET      = -1;                                                                   " << endl;              

  outfile2 << "  for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 				      " << endl;
  outfile2 << "    if ( leftET [iiET]  <= (ET)      &&       (ET) < rightET[iiET] ) {		      " << endl;
  outfile2 << "      iET  = iiET;  								      " << endl;
  outfile2 << "    }										      " << endl;
  outfile2 << "  }										      " << endl;

  // Extra protections										       
  outfile2 << "  if (ET < leftET  [0] )         { iET = 0;              if (DBG) std::cout << \" WARNING [applyScCorrections]: ET < leftET  [0] )       \" << std::endl;}   " << endl;
  outfile2 << "  if (ET > rightET [nBinsET-1] ) { iET = nBinsET-1;      if (DBG) std::cout << \" WARNING [applyScCorrections]: ET > rightET [nBinsET-1] \" << std::endl;}   " << endl;

  outfile2 << "  if (iET == -1) {delete  h_CBET_EE; return 1;}					      " << endl;

  outfile2 << "											      " << endl;
  outfile2 << "  Int_t binET =  2*iET+1 ; // h_CBET_EE->FindBin(ET);                                    " << endl;       
  outfile2 << "  Double_t tmp = h_CBET_EE->GetBinContent(binET);                                        " << endl;
  outfile2 << "  delete h_CBET_EE;                                                             	      " << endl;
  outfile2 << "  if ( 0< binET && binET < 2*nBinsET+1) return tmp;                          	      " << endl;
  outfile2 << "  else return 1.;                                                               	      " << endl;
  outfile2 << "                                                                                         " << endl;
  outfile2 << "}                                                                               " << endl;                                                            
  outfile2 << "                                                                                  " << endl; 

  outfile2 << "#endif                                                 " << endl;
  
  outfile2.close();                                                                                                      
  
  // save also the histograms to be used as correction
  TFile *histFile= new TFile("./histCorrections_ET.root", "RECREATE" );
  h_CBET_EB->Write();  
  h_CBET_EE->Write();  
  histFile->Close();    
  return;
}



void scCorrections::getErecEgenFits_BrEta(correctionType applyCorrections, bool apply)
{
  Int_t minEntries      = 100;
  Int_t minChi2         = 3;

  // Fit the CB for each histogram of the ErecEgen matrix
  //--------------------------------------------------------------------------------
  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){
    for (Int_t iBr = 0; iBr < nBinsBr; ++iBr){
      Double_t mean    = 1;
      Double_t emean   = 0;
      Double_t chi2    = 0;      
      Double_t fracBad = 0;
      if (h_ErecEgen[iEta][iBr]->GetEntries() >minEntries) 
	{
	  fitCB(h_ErecEgen[iEta][iBr], mean, emean, chi2, fracBad);
	  if (!apply){
	    if (chi2 > minChi2){ // if the fit is bad
	      if ( 2*iEta+1 -2 > 1){ // set it to the previous eta bin if possible
		mean    = hh_CBetabr      ->GetBinContent(2*iBr+1  , 2*iEta+1  -2);
		emean   = hh_CBetabr      ->GetBinError  (2*iBr+1  , 2*iEta+1  -2);
		cout << "WARNING: BAD FIT set mean and error to the values of the previous eta bin " << mean << " " << emean << endl;
		//getchar();	      
	      }
	      else if (2*iBr+1 -2 > 1 ){ // otherwise set it to the previous Br bin if possible
		mean    = hh_CBetabr      ->GetBinContent(2*iBr+1 -2  , 2*iEta+1 );
		emean   = hh_CBetabr      ->GetBinError  (2*iBr+1 -2  , 2*iEta+1 );	    
		cout << "WARNING: BAD FIT set mean and error to the values of the previous br bin " << mean << " " << emean << endl;
		//getchar();
	      }
	      else{  // otherwise set it to 1. Don't know how to treat them
		mean  = 1;
		emean = 0;
		cout << "WARNING: BAD FIT set mean and error to " << mean << " " << emean << endl;
		//getchar();
	      }
	    }
	  }
	}
      else {// if not enough statistics
	if (!apply){
	  if  ( 2*iEta+1 -2 > 1){ // set it to the previous eta bin if possible
	    mean    = hh_CBetabr      ->GetBinContent(2*iBr+1   , 2*iEta+1 -2 );
	    emean   = hh_CBetabr      ->GetBinError  (2*iBr+1   , 2*iEta+1 -2 );
	    cout << "WARNING: NOT ENOUGH STAT set mean and error to the values of the previous eta bin " << mean << " " << emean << endl;
	    //getchar();
	  }
	  else if (2*iBr+1 -2 > 1 ){ // otherwise set it to the previous Br bin if possible
	    mean    = hh_CBetabr      ->GetBinContent(2*iBr+1 -2  , 2*iEta+1 );
	    emean   = hh_CBetabr      ->GetBinError  (2*iBr+1 -2  , 2*iEta+1 );
	    cout << "WARNING: NOT ENOUGH STAT set mean and error to the values of the previous br bin " << mean << " " << emean << endl;
	    //getchar();
	  }
	  else {  // otherwise set it to 1. Don't know how to treat them
	    mean  = 1;
	    emean = 0;
	    //getchar();
	    cout << "WARNING: NOT ENOUGH STAT set mean and error to " << mean << " " << emean << endl;
	  }
	}
      }
      hh_CBetabr      ->SetBinContent(2*iBr+1,2*iEta+1,mean  );
      hh_CBetabr      ->SetBinError  (2*iBr+1,2*iEta+1,emean );
      hh_Chi2etabr    ->SetBinContent(2*iBr+1,2*iEta+1,chi2  ); // dof = 4 = number dof floating parametes
      hh_EgenFracEtaBr->SetBinContent(2*iBr+1,2*iEta+1, fracBad);
      cout << " FRACTION BELOW " << EgenFrac << " = " << fracBad << endl;
    }
  }
  
  // 2D plot CB most probable value in bins (eta,br)
  TCanvas *cc_CBetabr = new TCanvas("cc_CBetabr","cc_CBetabr",600,600);
  cc_CBetabr->cd();
  //  hh_CBetabr->GetZaxis()->SetRangeUser(0.89,1.01);
  hh_CBetabr->GetZaxis()->SetRangeUser(0.80,1.01);
  hh_CBetabr->GetXaxis()->SetTitle("#sigma_{phi}/#sigma_{eta}");
  hh_CBetabr->GetXaxis()->SetTitleOffset(1.5);
  hh_CBetabr->GetYaxis()->SetTitle("#eta");
  hh_CBetabr->GetYaxis()->SetTitleOffset(1.5);
  hh_CBetabr->GetZaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
  hh_CBetabr->GetZaxis()->SetTitleOffset(1.2);
  hh_CBetabr->SetTitle(0);
  gPad->SetPhi(180+30);
  hh_CBetabr->Draw("lego2");  
  if      (applyCorrections == brEta  )  cc_CBetabr->Print("./plots_tmp/cc_CBetabr_corrBrEta.png");
  else if (applyCorrections == brEtaET)  cc_CBetabr->Print("./plots_tmp/cc_CBetabr_corrBrEtaET.png");
  else                                   cc_CBetabr->Print("./plots_tmp/cc_CBetabr.png");
  
  TCanvas *cc_Chi2etabr = new TCanvas("cc_Chi2etabr","cc_Chi2etabr",600,600);
  cc_Chi2etabr->cd();
  gPad->SetPhi(180+30);
  hh_Chi2etabr->GetZaxis()->SetRangeUser(0.,5);
  hh_Chi2etabr->Draw("lego2");
  
  // 1D projection over br
  //
  // transfer each eta bins slice in a TH1F
  for (Int_t i = 0; i< nBinsEta; ++i){
    for (Int_t j = 0; j< nBinsBr; ++j){
      h_CBetabr[i]->SetBinContent(2*j+1,hh_CBetabr->GetBinContent(2*j+1,2*i+1));
      h_CBetabr[i]->SetBinError  (2*j+1,hh_CBetabr->GetBinError  (2*j+1,2*i+1));
    }
  }  
  // overlap the eta slices in a plot
  TCanvas *c_CBetabr = new TCanvas("c_CBetabr","c_CBetabr",600,600);
  c_CBetabr->cd();
  c_CBetabr->SetLeftMargin(0.2);
  for (Int_t i = 0; i<nBinsEta ; ++i){
    h_CBetabr[i]-> SetMarkerStyle(20);
    h_CBetabr[i]-> SetMarkerColor(i+22);
    h_CBetabr[i]-> SetTitle(0);
    h_CBetabr[i]-> GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");
    h_CBetabr[i]-> GetYaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
    h_CBetabr[i]-> GetYaxis()->SetTitleOffset(2);
    if(i == 0) {
      h_CBetabr[i]-> Draw();
      h_CBetabr[i]-> GetYaxis()->SetRangeUser(0.9,1.05);
    }
    else h_CBetabr[i]-> Draw("same");
  }
  TBox *bandBr = new TBox(leftBr[0],0.999,rightBr[nBinsBr-1],1.001);
  bandBr->SetFillColor(1);
  bandBr->SetFillStyle(3004);
  bandBr->Draw();
  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->SetFillColor(0);
  TString txt[nBinsEta];
  for (Int_t i = 0; i<nBinsEta ; ++i){
    txt[i] = floatToString(leftEta[i])+" < #eta < "+ floatToString(rightEta[i]);
    leg->AddEntry(h_CBetabr[i],txt[i]);
  }
  leg->Draw();
  if      (applyCorrections == brEta  ) c_CBetabr->Print("./plots_tmp/c_CBetabr_corrBrEta.png");
  else if (applyCorrections == brEtaET) c_CBetabr->Print("./plots_tmp/c_CBetabr_corrBrEtaET.png");
  else                                  c_CBetabr->Print("./plots_tmp/c_CBetabr.png");


  // Fit the CB for each histogram of the ErecEgen ET-vector  *** BARREL ***
  //--------------------------------------------------------------------------------
  TCanvas *c_CBET_EB = new TCanvas("c_CBET_EB","c_CBET_EB",600,600);
  c_CBET_EB->cd();
  c_CBET_EB->SetLeftMargin(0.2);
  for (Int_t iET = 0; iET < nBinsET; ++iET){
    Double_t mean  = 1;
    Double_t emean = 0;
    Double_t chi2  = 0;
    Double_t fracBad = 0;
    if (h_ETcorrBrEta_EB[iET]->GetEntries() > minEntries) fitCB(h_ETcorrBrEta_EB[iET], mean, emean, chi2, fracBad);
    h_CBET_EB  ->SetBinContent(2*iET+1,mean  );
    h_CBET_EB  ->SetBinError  (2*iET+1,emean );
    h_Chi2ET_EB->SetBinContent(2*iET+1,chi2  ); // dof = 4 = number fof floating parametes
    h_EgenFracET_EB->SetBinContent(2*iET+1, fracBad);
  }
  h_CBET_EB  ->SetMarkerStyle(20);
  h_CBET_EB  ->SetMarkerColor(1);
  h_CBET_EB  ->SetTitle(0);
  h_CBET_EB  ->GetYaxis()->SetRangeUser(0.99,1.01);
  h_CBET_EB  ->GetXaxis()->SetTitle("ET [GeV]");
  h_CBET_EB  ->GetYaxis()->SetTitleOffset(2);
  h_CBET_EB  ->GetYaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
  h_CBET_EB  ->Draw();
  TBox *bandET_EB = new TBox(leftET[0],0.999,rightET[nBinsET-1],1.001);
  bandET_EB->SetFillColor(1);
  bandET_EB->SetFillStyle(3004);
  bandET_EB->Draw();
  if      (applyCorrections == brEta  ) c_CBET_EB->Print("./plots_tmp/c_CBET_EB_corrBrEta.png");
  else if (applyCorrections == brEtaET) c_CBET_EB->Print("./plots_tmp/c_CBET_EB_corrBrEtaET.png");
  else                                  c_CBET_EB->Print("./plots_tmp/c_CBET_EB.png");

  TCanvas *c_Chi2ET_EB = new TCanvas("c_Chi2ET_EB","c_Chi2ET_EB",600,600);
  c_Chi2ET_EB->cd();
  gPad->SetPhi(180+30);
  h_Chi2ET_EB->GetYaxis()->SetRangeUser(0.,5);
  h_Chi2ET_EB->Draw("");
  

  // Fit the CB for each histogram of the ErecEgen ET-vector  *** ENDCAP ***
  //--------------------------------------------------------------------------------
  TCanvas *c_CBET_EE = new TCanvas("c_CBET_EE","c_CBET_EE",600,600);
  c_CBET_EE->cd();
  c_CBET_EE->SetLeftMargin(0.2);
  for (Int_t iET = 0; iET < nBinsET; ++iET){
    Double_t mean  = 1;
    Double_t emean = 0;
    Double_t chi2  = 0;
    Double_t fracBad = 0;
    if (h_ETcorrBrEta_EE[iET]->GetEntries() > minEntries) fitCB(h_ETcorrBrEta_EE[iET], mean, emean, chi2, fracBad);
    h_CBET_EE  ->SetBinContent(2*iET+1,mean  );
    h_CBET_EE  ->SetBinError  (2*iET+1,emean );
    h_Chi2ET_EE->SetBinContent(2*iET+1,chi2  ); // dof = 4 = number fof floating parametes
    h_EgenFracET_EE->SetBinContent(2*iET+1, fracBad);
  }
  h_CBET_EE  ->SetMarkerStyle(20);
  h_CBET_EE  ->SetMarkerColor(1);
  h_CBET_EE  ->SetTitle(0);
  h_CBET_EE  ->GetYaxis()->SetRangeUser(0.99,1.01);
  h_CBET_EE  ->GetXaxis()->SetTitle("ET [GeV]");
  h_CBET_EE  ->GetYaxis()->SetTitleOffset(2);
  h_CBET_EE  ->GetYaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
  h_CBET_EE  ->Draw();
  TBox *bandET_EE = new TBox(leftET[0],0.999,rightET[nBinsET-1],1.001);
  bandET_EE->SetFillColor(1);
  bandET_EE->SetFillStyle(3004);
  bandET_EE->Draw();
  if      (applyCorrections == brEta  ) c_CBET_EE->Print("./plots_tmp/c_CBET_EE_corrBrEta.png");
  else if (applyCorrections == brEtaET) c_CBET_EE->Print("./plots_tmp/c_CBET_EE_corrBrEtaET.png");
  else                                  c_CBET_EE->Print("./plots_tmp/c_CBET_EE.png");

  TCanvas *c_Chi2ET_EE = new TCanvas("c_Chi2ET_EE","c_Chi2ET_EE",600,600);
  c_Chi2ET_EE->cd();
  gPad->SetPhi(180+30);
  h_Chi2ET_EE->GetYaxis()->SetRangeUser(0.,5);
  h_Chi2ET_EE->Draw("");
  

  return;
}




















/*

//================================================================================
//
// test R9 instead of sigma_phi/sigma_eta
//
//================================================================================

// matrix (eta,R9) of histograms used for the Erec/Egen CB fits
TH1F *h_ErecEgen_R9[nBinsEta][nBinsR9]; 

// 2D histogram containing the results of the CB fits
TH2F *hh_CBetaR9    = new TH2F("hh_CBetaR9"    ,"hh_CBetaR9"    ,nBinsR9*2-1, R9bins, nBinsEta*2-1, etabins);

// 1D to fit the corrections
TH1F *h_CBetaR9[nBinsEta];

// 2D histogram containing the chi2 of the CB fits
TH2F *hh_Chi2etaR9  = new TH2F("hh_Chi2etaR9"  ,"hh_Chi2etaR9"  ,nBinsR9*2-1, R9bins, nBinsEta*2-1, etabins);

// vector ET dependence 
TH1F *h_ETcorrR9Eta_EB[nBinsET];
TH1F *h_ETcorrR9Eta_EE[nBinsET];

TH2F *hh_EgenFracEtaR9  = new TH2F("hh_EgenFracEtaR9"  ,"hh_EgenFracEtaR9"  ,nBinsR9*2-1, R9bins, nBinsEta*2-1, etabins);



// Extract f(R9, eta)
void scCorrections::run_R9Eta()
{    
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  initHistograms_R9();
  correctionType Corr  = none;
  fillHistograms_R9(Corr);
  bool apply = false;
  getErecEgenFits_R9Eta(Corr, apply);
  deriveCorrections_R9Eta();
}

// Apply f(R9, eta) (closure test)
void scCorrections::run_apply_R9Eta()
{  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  initHistograms_R9();
  correctionType Corr  = R9Eta;
  fillHistograms_R9(Corr);
  bool apply = true;
  getErecEgenFits_R9Eta(Corr, apply);
}

 // Extract the ET correction after applying f(R9, eta)
void scCorrections::run_R9EtaET()
{    
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  // extract the BrEta corrections
  run_R9Eta();

  // apply the BrEta correction
  run_apply_R9Eta();

  // write out the correction macro
  deriveCorrections_R9EtaET();  
}

// Apply f(sigma_phi/sigma_eta, eta)*F(ET) (closure test)
void scCorrections::run_apply_R9EtaET()
{  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  //
  initHistograms_R9();
  correctionType Corr  = R9EtaET;
  fillHistograms_R9(Corr);
  bool apply = true;
  getErecEgenFits_R9Eta(Corr, apply);
}

void scCorrections::initHistograms_R9()
{
  // R9
  for (Int_t i = 0; i<nBinsEta; ++i){
    for (Int_t j = 0; j<nBinsR9; ++j){
      h_ErecEgen_R9[i][j] = new TH1F(Form("h_ErecEgen_R9_%d_%d",i,j),Form("h_ErecEgen_R9_%d_%d",i,j),nBinsErecEgen,ErecEgenMin, ErecEgenMax);
    }
  }
  for (Int_t i = 0; i<nBinsEta; ++i){
    h_CBetaR9[i] = new TH1F(Form("h_CBetaR9_%d",i), Form("h_CBetaR9_%d",i),nBinsR9*2-1, R9bins);
  }

  // ET - R9
  for (Int_t i = 0; i<nBinsET; ++i){
    h_ETcorrR9Eta_EB[i] = new TH1F(Form("h_ETcorrR9Eta_EB_%d",i), Form("h_ETcorrR9Eta_EB_%d",i), nBinsErecEgen,ErecEgenMin, ErecEgenMax);
    h_ETcorrR9Eta_EE[i] = new TH1F(Form("h_ETcorrR9Eta_EE_%d",i), Form("h_ETcorrR9Eta_EE_%d",i), nBinsErecEgen,ErecEgenMin, ErecEgenMax);
  }
  cout << "initHistograms_R9(): done" << endl;
}

//
void scCorrections::fillHistograms_R9(correctionType applyCorrections )
{
  // read the Correction parameters
  //
  // f(R9, eta ) from histogram 
  //
  TFile *histFile_R9Eta = 0;
  TH1F *h_corr[nBinsEta];
  if (applyCorrections == R9Eta || applyCorrections == R9EtaET){
    histFile_R9Eta = new TFile("./histCorrections_R9Eta.root","READ");
    if (histFile_R9Eta->IsZombie()) {cout << "ERROR: ./histCorrections_R9Eta.root does not exist " << endl; exit(-1);}
    //
    for (Int_t i = 0; i<nBinsEta; ++i){
      h_corr[i] = (TH1F*)histFile_R9Eta->Get(Form("h_corr_%d",i));
    }
  }  
  //
  // F(ET) correction from histogram ( *** BARREL / ENDCAP *** )
  //
  TFile *histFile_ET = 0;
  TH1F *h_corrET_EB = 0;
  TH1F *h_corrET_EE = 0;
  if (applyCorrections == R9EtaET){
    histFile_ET = new TFile("./histCorrections_ET_R9.root","READ");
    if (histFile_ET->IsZombie()) {cout << "ERROR: ./histCorrections_ET_R9.root does not exist " << endl; exit(-1);}
    //    
    h_corrET_EB = (TH1F*)histFile_ET->Get("h_CBET_EB");
    h_corrET_EE = (TH1F*)histFile_ET->Get("h_CBET_EE");
  }

  // To speed up the loop read only the used variables from the NTuple
  // Add ALL variables used, else they are set to 0 and the compiler won't notice it
  //
  if (fChain == 0) return;
  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("phoSCEta"      ,1); 
  fChain->SetBranchStatus("phoSCPhi"      ,1); 
  fChain->SetBranchStatus("phoCetaCorrE"  ,1); 
  fChain->SetBranchStatus("phoCetaCorrEt" ,1); 
  fChain->SetBranchStatus("phoSCBrem"     ,1); 
  fChain->SetBranchStatus("phoR9"         ,1);       
  fChain->SetBranchStatus("mcE"           ,1);
  fChain->SetBranchStatus("mcEta"         ,1);
  fChain->SetBranchStatus("mcPhi"         ,1);
  fChain->SetBranchStatus("eleSCEta"      ,1);        
  fChain->SetBranchStatus("eleSCPhi"      ,1);        
  fChain->SetBranchStatus("eleCetaCorrE"  ,1);  
  fChain->SetBranchStatus("eleCetaCorrEt" ,1); 
  fChain->SetBranchStatus("eleBrLinear"   ,1);   
  fChain->SetBranchStatus("eleE3x3"       ,1);
  fChain->SetBranchStatus("eleCetaCorrE"  ,1);
  fChain->SetBranchStatus("eleSCRawEn"    ,1);

  // Loop over all entries:
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry % 1000 == 0) cout << jentry << endl;

    // Loop over the reconstructed objects (electrons or photons)
    Int_t nobj = 0;
    if      (partType == "electron") nobj = nEle;
    else if (partType == "photon"  ) nobj = nPho;

    // skip pathological events
    if (nobj < 1 || nobj > 2) continue;
    
    for (Int_t iobj = 0; iobj < nobj; ++iobj) {

      Double_t nt_em_eta        = 0;     
      Double_t nt_em_phi        = 0;
      Double_t nt_emCorrEta_e   = 0;
      Double_t nt_emCorrEta_et  = 0;
      Double_t nt_em_br1        = 0;
      Double_t nt_em_r9         = 0;
      Double_t nt_mc_e	        = 0; 
      
      // Use the correct naming of your ntuple
      if (partType == "photon"){
	nt_em_eta        = phoSCEta[iobj];        
	nt_em_phi        = phoSCPhi[iobj];        
	nt_emCorrEta_e   = phoCetaCorrE[iobj];  
	nt_emCorrEta_et  = phoCetaCorrEt[iobj]; 
	nt_em_br1        = phoSCBrem[iobj];     
	nt_em_r9         = phoR9[iobj]; 
	nt_mc_e          = mcE[0]; // mcE[0] = mcE[1] by contruction of the di-X ntuples
	//
	// Cheap truth match:
	// the photons/electrons are back to back, just need to compare one reco-photon to the two mc-photons
	// 	Float_t dr0 = TMath::Sqrt(pow(phoSCEta[0]-mcEta[0],2) + pow(phoSCPhi[0]-mcPhi[0],2));
	// 	Float_t dr1 = TMath::Sqrt(pow(phoSCEta[0]-mcEta[1],2) + pow(phoSCPhi[0]-mcPhi[1],2));
	// 	if (dr0 < dr1 ) nt_mc_e = mcE[0];            
	// 	else nt_mc_e = mcE[1];  
	// 	cout << phoSCEta[0] << " " << mcEta[0] <<  " " << phoSCPhi[0] << " " << mcPhi[0] <<  " " << mcE[0]  << endl;
	// 	cout << phoSCEta[0] << " " << mcEta[1] <<  " " << phoSCPhi[0] << " " << mcPhi[1] <<  " " << mcE[1]  << endl;
	// 	cout << nt_mc_e << endl;
      }
      else if (partType == "electron"){
	nt_em_eta        = eleSCEta[iobj];        
	nt_em_phi        = eleSCPhi[iobj];        
	nt_emCorrEta_e   = eleCetaCorrE[iobj];  
	nt_emCorrEta_et  = eleCetaCorrEt[iobj]; 
	nt_em_br1        = eleBrLinear[iobj];   
	if ( eleCetaCorrE[iobj] != 0) nt_em_r9 =  eleE3x3[iobj]/eleSCRawEn[iobj]; 
	nt_mc_e          = mcE[0]; // mcE[0] = mcE[1] by contruction of the di-X ntuples
	//
	// Cheap truth match:
	// the photons/electrons are back to back, just need to compare one reco-photon to the two mc-photons
	// 	Float_t dr0 = TMath::Sqrt(pow(eleSCEta[0]-mcEta[0],2) + pow(eleSCPhi[0]-mcPhi[0],2));
	// 	Float_t dr1 = TMath::Sqrt(pow(eleSCEta[0]-mcEta[1],2) + pow(eleSCPhi[0]-mcPhi[1],2));
	// 	if (dr0 < dr1 ) nt_mc_e = mcE[0];            
	// 	else nt_mc_e = mcE[1];  
	// 	cout << eleSCEta[0] << " " << mcEta[0] <<  " " << eleSCPhi[0] << " " << mcPhi[0] <<  " " << mcE[0]  << endl;
	// 	cout << eleSCEta[0] << " " << mcEta[1] <<  " " << eleSCPhi[0] << " " << mcPhi[1] <<  " " << mcE[1]  << endl;
	// 	cout << nt_mc_e << endl;
      }

      // main event selection
      //---------------------------------
      if ( 
	  // select showering objects: 
	  ( (TMath::Abs(nt_em_eta) < etaCrackMin && (nt_em_r9 < minR9_EB)) || (TMath::Abs(nt_em_eta) > etaCrackMax && (nt_em_r9 < minR9_EE)) ) 
	  
	  &&

	  // remove the lowest ET objects
	  nt_emCorrEta_et > minETCut
	  
	  && 
	  
	  // remove the phi gaps between the super modules (only for the barrel):: local corrections will take care of them
	  // for the endcap there are no local corrections, so don't exclude any object
	  !isInPhiCracks(nt_em_phi, nt_em_eta,  fiducialPhiCut) 	  
	  
	   ){ 
	
	// Find the correct (R9,eta) bin for the h_ErecEgen_R9[nBinsEta][nBinsR9] and h_ETcorrR9Eta_EB[nBinsET];
	Int_t iEta    = -1;
	Int_t iR9     = -1;
	Int_t iET     = -1;
	Int_t iETcorr = -1;	
	//
	for (Int_t iiEta = 0; iiEta < nBinsEta; ++iiEta){ 
	  if ( leftEta[iiEta] <= TMath::Abs(nt_em_eta) && TMath::Abs(nt_em_eta) <rightEta[iiEta]) {
	    iEta = iiEta; 
	    // cout << "ETA " << leftEta[iiEta]  << " " << TMath::Abs(nt_em_eta) << " " << rightEta[iiEta] << " " << iEta << endl;
	  }
	}
	for (Int_t iiR9  = 0; iiR9  < nBinsR9;  ++iiR9 ){ 
	  if ( leftR9 [iiR9]  <= nt_em_r9             &&             nt_em_r9 <rightR9 [iiR9] ) {
	    iR9  = iiR9;  
	    // cout << "R9 " << leftR9[iiR9]  << " " << nt_em_r9 << " " << rightR9[iiR9] << " " << iR9 << endl;
	  }
	}
	// Find the correct (ET) bin for the h_ETcorrR9Eta_EX[nBinsET]
	for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 
	  if ( leftET [iiET]  <= nt_emCorrEta_et      &&       nt_emCorrEta_et < rightET[iiET] ) 
	    iET  = iiET;  
	  // cout << "ET " << leftR9[iiET]  << " " << nt_emCorrEta_et << " " << rightR9[iiET] << " " << nt_emCorrEta_et << endl;
	}		
	if (iEta == -1 || iR9 == -1 || iET == -1) continue;


	// Initialize the correction factors
	Double_t corrR9Eta = 1.; 
	Double_t corrEt    = 1.; 	
	
	// get the corrections f(R9,eta) and F(ET) from histogram
	//
	if (applyCorrections == R9Eta || applyCorrections == R9EtaET) {
	  corrR9Eta  = h_corr[iEta]->GetBinContent(2*iR9+1);                               
	  //
	  // Find the correct (ET corrected with f(R9,Eta) ) bin for the correction histogram 
	  for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 
	    if ( leftET [iiET]  <= (nt_emCorrEta_et / corrR9Eta)      &&       (nt_emCorrEta_et / corrR9Eta) < rightET[iiET] ) {
	      iETcorr  = iiET;  
	      // cout << "ET " << leftET [iiET] << " " <<  nt_emCorrEta_et / corrR9Eta  << " " <<  rightET[iiET] <<  " " << iETcorr <<endl;
	    }
	  }
	  if (iETcorr == -1) continue; // always skip the event, there is something patological	    
	  if (applyCorrections == R9EtaET) {
	    //
	    // histogram correction ET *** BARREL, ENDCAP ***		
	    if (TMath::Abs(nt_em_eta) < etaCrackMin ){                    
	      corrEt = h_corrET_EB    ->GetBinContent(2*iETcorr+1);		
	      //	    cout << " corr " << iETcorr << " " << h_corrET_EB    ->GetBinContent(2*iETcorr+1) <<endl;		
	    }								
	    if (TMath::Abs(nt_em_eta) > etaCrackMax ) {			
	      corrEt = h_corrET_EE    ->GetBinContent(2*iETcorr+1);		
	      //	    cout << "corr " << iETcorr << " " << h_corrET_EE    ->GetBinContent(2*iETcorr+1) <<endl;		
	    }		
	  }
	}

//MDDB	// get the corrections f(R9,eta) and F(ET) from File
//MDDB	//
//MDDB	if (applyCorrections == R9Eta || applyCorrections == R9EtaET) {
//MDDB	  corrR9Eta = applyScCorrectionsR9Eta(nt_em_eta,  nt_em_r9);                           //MDDB
//MDDB	  
//MDDB	  // Find the correct (ET corrected with f(R9,Eta) ) bin for the correction histogram 
//MDDB	  for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 
//MDDB	    if ( leftET [iiET]  <= (nt_emCorrEta_et / corrR9Eta)      &&       (nt_emCorrEta_et / corrR9Eta) < rightET[iiET] ) {
//MDDB	      iETcorr  = iiET;  
//MDDB	      // cout << "ET " << leftET [iiET] << " " <<  nt_emCorrEta_et / corrR9Eta  << " " <<  rightET[iiET] <<  " " << iETcorr <<endl;
//MDDB	    }
//MDDB	  }	
//MDDB	  if (iETcorr == -1) continue; // always skip the event, there is something patological	   
//MDDB	  if (applyCorrections == R9EtaET) {
//MDDB	    if (TMath::Abs(nt_em_eta) < etaCrackMin ) corrEt = applyScCorrectionsET_EB( nt_emCorrEta_et / corrR9Eta );   //MDDB
//MDDB	    if (TMath::Abs(nt_em_eta) > etaCrackMax ) corrEt = applyScCorrectionsET_EE( nt_emCorrEta_et / corrR9Eta );   //MDDB
//MDDB	  }
//MDDB	}

	// fill the matrix of histograms h_ErecEgen_R9[nBinsEta][nBinsR9] 
	// ===========================================================
	//
	// *** f(R9,eta) ***
	if (applyCorrections == R9Eta) {
	  // fill
	  if (nt_mc_e!=0 && corrR9Eta!=0) h_ErecEgen_R9[iEta][iR9]->Fill(nt_emCorrEta_e / nt_mc_e / corrR9Eta);	  
	}	
	// *** f(R9,eta)*F(ET) ***
	else if (applyCorrections == R9EtaET) {
	  // fill
	  if (nt_mc_e!=0 && corrEt!=0) h_ErecEgen_R9[iEta][iR9]->Fill(nt_emCorrEta_e / nt_mc_e / corrR9Eta /corrEt);	  
	}	
	// *** no corrections ***
	else { 
	  // fill
	  if (nt_mc_e!=0 )  h_ErecEgen_R9[iEta][iR9]->Fill(nt_emCorrEta_e / nt_mc_e );     
	}
	
	// fill the vector of ET dependence histograms
	// ===========================================
	//	
	// *** f(R9,eta) ***
	if (applyCorrections == R9Eta) {
	  if (TMath::Abs(nt_em_eta) < etaCrackMin ){
	    // fill
	    if (nt_mc_e!=0 && corrR9Eta!=0)    h_ETcorrR9Eta_EB[iET]->Fill(nt_emCorrEta_e / nt_mc_e / corrR9Eta);
	  }
	  else if (TMath::Abs(nt_em_eta) > etaCrackMax ) {
	    // fill
	    if (nt_mc_e!=0 && corrR9Eta!=0)    h_ETcorrR9Eta_EE[iET]->Fill(nt_emCorrEta_e / nt_mc_e / corrR9Eta);		
	  }
	}
	// *** f(R9,eta)*F(ET) ***
	else if (applyCorrections == R9EtaET) {		
	  // bin position for the R9,eta corrected ET
	  if (TMath::Abs(nt_em_eta) < etaCrackMin ){
	    
	    // fill
	    if (nt_mc_e!=0 && corrEt!=0) h_ETcorrR9Eta_EB[iETcorr]->Fill(nt_emCorrEta_e / nt_mc_e / corrR9Eta / corrEt );	   
	  }
	  if (TMath::Abs(nt_em_eta) > etaCrackMax ) {
	    // fill
	    if (nt_mc_e!=0 && corrEt!=0) h_ETcorrR9Eta_EE[iETcorr]->Fill(nt_emCorrEta_e / nt_mc_e / corrR9Eta / corrEt);	   
	  }	      
	}	    
	// *** no corrections ***
	else{
	  if (TMath::Abs(nt_em_eta) < etaCrackMin ){
	    h_ETcorrR9Eta_EB[iET]->Fill(nt_emCorrEta_e / nt_mc_e );
	  }
	  if (TMath::Abs(nt_em_eta) > etaCrackMax ) { 
	    h_ETcorrR9Eta_EE[iET]->Fill(nt_emCorrEta_e / nt_mc_e );
	  }
	}

      } // if selected event
    } // loop over obj
  } // loop over entries	
  
  if (applyCorrections == R9Eta || applyCorrections == R9EtaET) histFile_R9Eta->Close();    
  if (applyCorrections == R9EtaET) histFile_ET->Close();
  
  return;
  
}

void scCorrections::deriveCorrections_R9Eta()
{
  // Write out the correction macro
  ofstream outfile2;  
  TString filename = "./applyScCorrections.C";
  outfile2.open(filename);
  outfile2 << "#ifndef scCorrectionsR9ETAET_h                                              " << endl;
  outfile2 << "#define scCorrectionsR9ETAET_h                                              " << endl;
  outfile2 << "#include \"TFile.h\"                                                      " << endl;
  outfile2 << "#include \"TH1F.h\"				                       " << endl;
  outfile2 << "  	                                                               " << endl;
  outfile2 << "Double_t applyScCorrectionsR9Eta(Double_t eta, Double_t R9){  " << endl;
  outfile2 << "  	                                                               " << endl;
  outfile2 << "  const Int_t    nBinsEta = " << nBinsEta << ";                    " << endl;                                  
  outfile2 << "  Double_t       leftEta  ["  << nBinsEta << "];                   " << endl;
  outfile2 << "  Double_t       rightEta ["  << nBinsEta << "];                   " << endl;
  outfile2 << "  const Int_t    nBinsR9  = " << nBinsR9 << ";                     " << endl;                                    
  outfile2 << "  Double_t       leftR9  [" << nBinsR9 << "];                      " << endl;
  outfile2 << "  Double_t       rightR9 [" << nBinsR9 << "];                      " << endl;
  outfile2 << "  Double_t R9bins  [" << nBinsR9*2 << "];                    " << endl;
  outfile2 << "  TH1F *h_corr["<< nBinsEta << "];                                 " << endl;
  outfile2 << "  	                                                               " << endl;
  for (Int_t i = 0; i< nBinsEta; ++i){ 
    outfile2 << "  leftEta[" << i << "] =  " << leftEta[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsEta; ++i){ 
    outfile2 << "  rightEta[" << i<< "] =  " << rightEta[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsR9; ++i){ 
    outfile2 << "  leftR9[" << i << "] =  " << leftR9[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsR9; ++i){ 
    outfile2 << "  rightR9["<< i << "] =  " << rightR9[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< 2*nBinsR9; ++i){ 
    outfile2 << "  R9bins[" << i << "] =  " << R9bins[i] << " ; " << endl;
  }
  outfile2 << "  for (Int_t i = 0; i<nBinsEta; ++i){                                                   " << endl;
  outfile2 << "    h_corr[i] = new TH1F(Form(\"h_corr_%d\",i),Form(\"h_corr_%d\",i),nBinsR9*2-1, R9bins);  " << endl;
  outfile2 << "  }                                                                                     " << endl;
  for (Int_t i = 0; i<nBinsEta ; ++i){
    for (Int_t j = 0; j< nBinsR9; ++j){
      outfile2 << "  h_corr[" << i << "]->SetBinContent(" <<  2*j+1 << "," << h_CBetaR9[i]->GetBinContent(2*j+1) << " );"<< endl;
    }
    outfile2 << endl;
  }
  outfile2 << "   Int_t tmpEta = -1;                                                                                             " << endl;         
  outfile2 << "   for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								       " << endl;
  outfile2 << "     if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] ){				       " << endl;
  outfile2 << "       tmpEta = iEta;											       " << endl;
  outfile2 << "     }													       " << endl;
  outfile2 << "   }													       " << endl;
  outfile2 << "   													       " << endl;
  outfile2 << "   Int_t tmpR9 = -1;											       " << endl;
  outfile2 << "   for (Int_t iR9 = 0; iR9 < nBinsR9; ++iR9){			       " << endl;
  outfile2 << "     if (leftR9 [iR9]  <= R9 && R9 <rightR9 [iR9] ){      " << endl;	   
  outfile2 << "       tmpR9 = iR9;										       " << endl;
  outfile2 << "     }													       " << endl;
  outfile2 << "   }													       " << endl;
  outfile2 << "   if (tmpEta == -1 || tmpR9 == -1){									       " << endl;
  outfile2 << "     for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  		       " << endl;
  outfile2 << "     return  1;												       " << endl;
  outfile2 << "   }													       " << endl;
  outfile2 << "   													       " << endl;
  outfile2 << "   Double_t tmp = h_corr[tmpEta]->GetBinContent(2*tmpR9+1);				   		       " << endl;
  outfile2 << "   for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  		       " << endl;
  outfile2 << "   return  tmp;                                                                                                   " << endl;                         
  outfile2 << "}													   " << endl;  

  //  outfile2 << "  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								   " << endl;
  //  outfile2 << "    for (Int_t iR9 = 0; iR9 < nBinsR9; ++iR9){			   " << endl;
  //  outfile2 << "      													   " << endl;
  //  outfile2 << "      // select the right (eta,R9) bin in the matrix							   " << endl;
  //  outfile2 << "      if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] &&			   " << endl;
  //  outfile2 << " 	   leftR9 [iR9]  <= R9 && R9 <rightR9 [iR9] ){	   " << endl;
  //  outfile2 << "	       Double_t tmp = h_corr[iEta]->GetBinContent(2*iR9+1);				   " << endl;	   
  //  outfile2 << "	       for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  " << endl;
  //  outfile2 << "	       return  tmp;                                                                                        " << endl;
  //  outfile2 << "      }													   " << endl;
  //  outfile2 << "    } 													   " << endl;
  //  outfile2 << "  }													   " << endl;
  //  outfile2 << "	       for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  " << endl;
  //  outfile2 << "  return 1;												   " << endl;
  //  outfile2 << "}													   " << endl;
  
  outfile2.close();                                                                                                      
  
  // save also the histograms to be used as correction
  TFile *histFile= new TFile("./histCorrections_R9Eta.root", "RECREATE" );
  TH1F *h_corr[nBinsEta];
  for (Int_t i = 0; i<nBinsEta ; ++i){
    h_corr[i] = (TH1F*)h_CBetaR9[i]->Clone(Form("h_corr_%d",i));
    h_corr[i]->Write();
  }
  histFile->Close();
    
  TCanvas *cc_BadClustersEtaR9 = new TCanvas("cc_BadClustersEtaR9" ,"cc_BadClustersEtaR9", 600,600);
  cc_BadClustersEtaR9->cd();
  gPad->SetPhi(270+30);
  hh_EgenFracEtaR9->GetXaxis()->SetTitle("R9");
  hh_EgenFracEtaR9->GetXaxis()->SetTitleOffset(1.6);
  hh_EgenFracEtaR9->GetYaxis()->SetTitle("#eta");
  hh_EgenFracEtaR9->GetYaxis()->SetTitleOffset(1.6);
  hh_EgenFracEtaR9->Draw("lego2");
  TCanvas *c_BadClustersET_EB     = new TCanvas("c_BadClustersET_EB"     ,"c_BadClustersET_EB"    , 600,600);
  c_BadClustersET_EB->cd();
  h_EgenFracET_EB ->GetXaxis()->SetTitle("ET [GeV]");
  h_EgenFracET_EB ->Draw("");
  TCanvas *c_BadClustersET_EE     = new TCanvas("c_BadClustersET_EE"     ,"c_BadClustersET_EE"    , 600,600);
  c_BadClustersET_EE->cd();
  h_EgenFracET_EE ->GetXaxis()->SetTitle("ET [GeV]");
  h_EgenFracET_EE ->Draw("");

  return;
}


void scCorrections::deriveCorrections_R9EtaET()
{
  // Add to the existing applyCorrection the ET dependence
  ofstream outfile2;  
  TString filename = "./applyScCorrections.C";
  outfile2.open(filename,ios::app);

  outfile2 << "                                                                         " << endl;
  outfile2 << "Double_t applyScCorrectionsET_EB(Double_t ET){  							   " << endl;
  outfile2 << "  const Int_t    nBinsET             = " << nBinsET<< ";             " << endl; 
  outfile2 << "  Double_t       leftET  [nBinsET];                    " << endl;
  outfile2 << "  Double_t       rightET [nBinsET];                    " << endl;
  outfile2 << "  Double_t       ETBins  [nBinsET*2];                  " << endl;
  // ET
  for (Int_t i = 0; i< nBinsET; ++i){ 
    outfile2 << "  leftET[" << i << "] =  " << leftET[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsET; ++i){ 
    outfile2 << "  rightET[" << i << "] =  " << rightET[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< 2*nBinsET; ++i){ 
    outfile2 << "  ETBins[" << i << "] =  " << ETBins[i] << " ; " << endl;
  }
  outfile2 << "  TH1F *h_CBET_EB    = new TH1F(\"h_CBET_EB\"    ,\"h_CBET_EB\"    ,nBinsET*2-1, ETBins); " << endl;

  // this is the residual ET dependence after R9Eta correction
  for (Int_t iET = 0; iET< nBinsET; ++iET){ 
    outfile2 << "  h_CBET_EB->SetBinContent(" << 2*iET+1 << ", " << h_CBET_EB->GetBinContent(2*iET+1) << "); " << endl;
  }

  outfile2 << "  Int_t iET      = -1;                                                                   " << endl;              
  outfile2 << "  for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 				      " << endl;
  outfile2 << "    if ( leftET [iiET]  <= (ET)      &&       (ET) < rightET[iiET] ) {		      " << endl;
  outfile2 << "      iET  = iiET;  								      " << endl;
  outfile2 << "    }										      " << endl;
  outfile2 << "  }										      " << endl;
  outfile2 << "  if (iET == -1) {delete  h_CBET_EB; return 1;}					      " << endl;
  outfile2 << "											      " << endl;
  outfile2 << "  Int_t binET =  2*iET+1 ; // h_CBET_EB->FindBin(ET);                                    " << endl;       
  outfile2 << "  Double_t tmp = h_CBET_EB->GetBinContent(binET);                                        " << endl;
  outfile2 << "  delete h_CBET_EB;                                                             	      " << endl;
  outfile2 << "  if ( 0< binET && binET < 2*nBinsET+1) return tmp;                          	      " << endl;
  outfile2 << "  else return 1.;                                                               	      " << endl;
  outfile2 << "                                                                                         " << endl;
  outfile2 << "}                                                                               " << endl;                                                            
  outfile2 << "                                                                                  " << endl; 
    
  outfile2 << "                                                                         " << endl;
  outfile2 << "Double_t applyScCorrectionsET_EE(Double_t ET){  							   " << endl;
  outfile2 << "  const Int_t    nBinsET             = " << nBinsET<< ";             " << endl; 
  outfile2 << "  Double_t       leftET  [nBinsET];                    " << endl;
  outfile2 << "  Double_t       rightET [nBinsET];                    " << endl;
  outfile2 << "  Double_t       ETBins  [nBinsET*2];                  " << endl;
  // ET
  for (Int_t i = 0; i< nBinsET; ++i){ 
    outfile2 << "  leftET[" << i << "] =  " << leftET[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsET; ++i){ 
    outfile2 << "  rightET[" << i << "] =  " << rightET[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< 2*nBinsET; ++i){ 
    outfile2 << "  ETBins[" << i << "] =  " << ETBins[i] << " ; " << endl;
  }
  outfile2 << "  TH1F *h_CBET_EE    = new TH1F(\"h_CBET_EE\"    ,\"h_CBET_EE\"    ,nBinsET*2-1, ETBins); " << endl;

  // this is the residual ET dependence after R9Eta correction
  for (Int_t iET = 0; iET< nBinsET; ++iET){ 
    outfile2 << "  h_CBET_EE->SetBinContent(" << 2*iET+1 << ", " << h_CBET_EE->GetBinContent(2*iET+1) << "); " << endl;
  }

  outfile2 << "  Int_t iET      = -1;                                                                   " << endl;              
  outfile2 << "  for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 				      " << endl;
  outfile2 << "    if ( leftET [iiET]  <= (ET)      &&       (ET) < rightET[iiET] ) {		      " << endl;
  outfile2 << "      iET  = iiET;  								      " << endl;
  outfile2 << "    }										      " << endl;
  outfile2 << "  }										      " << endl;
  outfile2 << "  if (iET == -1) {delete  h_CBET_EE; return 1;}					      " << endl;
  outfile2 << "											      " << endl;
  outfile2 << "  Int_t binET =  2*iET+1 ; // h_CBET_EE->FindBin(ET);                                    " << endl;       
  outfile2 << "  Double_t tmp = h_CBET_EE->GetBinContent(binET);                                        " << endl;
  outfile2 << "  delete h_CBET_EE;                                                             	      " << endl;
  outfile2 << "  if ( 0< binET && binET < 2*nBinsET+1) return tmp;                          	      " << endl;
  outfile2 << "  else return 1.;                                                               	      " << endl;
  outfile2 << "                                                                                         " << endl;
  outfile2 << "}                                                                               " << endl;                                                            
  outfile2 << "                                                                                  " << endl; 

  outfile2 << "#endif                                                 " << endl;
  
  outfile2.close();                                                                                                      
  
  // save also the histograms to be used as correction
  TFile *histFile= new TFile("./histCorrections_ET_R9.root", "RECREATE" );
  h_CBET_EB->Write();  
  h_CBET_EE->Write();  
  histFile->Close();    
  return;
}


void scCorrections::getErecEgenFits_R9Eta(correctionType applyCorrections, bool apply)
{
  Int_t minEntries      = 100;
  Int_t minChi2         = 3;

  // Fit the CB for each histogram of the ErecEgen matrix
  //--------------------------------------------------------------------------------
  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){
    for (Int_t iR9 = 0; iR9 < nBinsR9; ++iR9){
      Double_t mean    = 1;
      Double_t emean   = 0;
      Double_t chi2    = 0;      
      Double_t fracBad = 0;
      if (h_ErecEgen_R9[iEta][iR9]->GetEntries() >minEntries) 
	{
	  fitCB(h_ErecEgen_R9[iEta][iR9], mean, emean, chi2, fracBad);
	  if (!apply){	    
	    if (chi2 > minChi2){ // if the fit is bad
	      if ( 2*iEta+1 -2 > 1){ // set it to the previous eta bin if possible
		mean    = hh_CBetaR9      ->GetBinContent(2*iR9+1  , 2*iEta+1  -2);
		emean   = hh_CBetaR9      ->GetBinError  (2*iR9+1  , 2*iEta+1  -2);
		cout << "WARNING: BAD FIT set mean and error to the values of the previous eta bin " << mean << " " << emean << endl;
		//getchar();	      
	      }
	      else if (2*iR9+1 -2 > 1 ){ // otherwise set it to the previous R9 bin if possible
		mean    = hh_CBetaR9      ->GetBinContent(2*iR9+1 -2  , 2*iEta+1 );
		emean   = hh_CBetaR9      ->GetBinError  (2*iR9+1 -2  , 2*iEta+1 );	    
		cout << "WARNING: BAD FIT set mean and error to the values of the previous R9 bin " << mean << " " << emean << endl;
		//getchar();
	      }
	      else{  // otherwise set it to 1. Don't know how to treat them
		mean  = 1;
		emean = 0;
		cout << "WARNING: BAD FIT set mean and error to " << mean << " " << emean << endl;
		//getchar();
	      }
	    }
	  }
	}
      else { // if not enough statistics
	if (!apply){
	  if  ( 2*iEta+1 -2 > 1){ // set it to the previous eta bin if possible
	    mean    = hh_CBetaR9      ->GetBinContent(2*iR9+1   , 2*iEta+1 -2 );
	    emean   = hh_CBetaR9      ->GetBinError  (2*iR9+1   , 2*iEta+1 -2 );
	    cout << "WARNING: NOT ENOUGH STAT set mean and error to the values of the previous eta bin " << mean << " " << emean << endl;
	    //getchar();
	  }
	  else if (2*iR9+1 -2 > 1 ){ // otherwise set it to the previous R9 bin if possible
	    mean    = hh_CBetaR9      ->GetBinContent(2*iR9+1 -2  , 2*iEta+1 );
	    emean   = hh_CBetaR9      ->GetBinError  (2*iR9+1 -2  , 2*iEta+1 );
	    cout << "WARNING: NOT ENOUGH STAT set mean and error to the values of the previous R9 bin " << mean << " " << emean << endl;
	    //getchar();
	  }
	  else {  // otherwise set it to 1. Don't know how to treat them
	    mean  = 1;
	    emean = 0;
	    //getchar();
	    cout << "WARNING: NOT ENOUGH STAT set mean and error to " << mean << " " << emean << endl;
	  }
	}
      }
      hh_CBetaR9      ->SetBinContent(2*iR9+1,2*iEta+1,mean  );
      hh_CBetaR9      ->SetBinError  (2*iR9+1,2*iEta+1,emean );
      hh_Chi2etaR9    ->SetBinContent(2*iR9+1,2*iEta+1,chi2  ); // dof = 4 = number dof floating parametes
      hh_EgenFracEtaR9->SetBinContent(2*iR9+1,2*iEta+1, fracBad);
      cout << " FRACTION BELOW " << EgenFrac << " = " << fracBad << endl;
    }
  }
  
  // 2D plot CB most probable value in bins (eta,R9)
  TCanvas *cc_CBetaR9 = new TCanvas("cc_CBetaR9","cc_CBetaR9",600,600);
  cc_CBetaR9->cd();
  hh_CBetaR9->GetZaxis()->SetRangeUser(0.89,1.01);
  hh_CBetaR9->GetXaxis()->SetTitle("R9");
  hh_CBetaR9->GetXaxis()->SetTitleOffset(1.5);
  hh_CBetaR9->GetYaxis()->SetTitle("#eta");
  hh_CBetaR9->GetYaxis()->SetTitleOffset(1.5);
  hh_CBetaR9->GetZaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
  hh_CBetaR9->GetZaxis()->SetTitleOffset(1.2);
  hh_CBetaR9->SetTitle(0);
  gPad->SetPhi(180+30);
  hh_CBetaR9->Draw("lego2");  
  if      (applyCorrections == R9Eta  )  cc_CBetaR9->Print("./plots_tmp/cc_CBetaR9_corrR9Eta.png");
  else if (applyCorrections == R9EtaET)  cc_CBetaR9->Print("./plots_tmp/cc_CBetaR9_corrR9EtaET.png");
  else                                   cc_CBetaR9->Print("./plots_tmp/cc_CBetaR9.png");
  
  TCanvas *cc_Chi2etaR9 = new TCanvas("cc_Chi2etaR9","cc_Chi2etaR9",600,600);
  cc_Chi2etaR9->cd();
  gPad->SetPhi(180+30);
  hh_Chi2etaR9->GetZaxis()->SetRangeUser(0.,5);
  hh_Chi2etaR9->Draw("lego2");
  
  // 1D projection over R9
  //
  // transfer each eta bins slice in a TH1F
  for (Int_t i = 0; i< nBinsEta; ++i){
    for (Int_t j = 0; j< nBinsR9; ++j){
      h_CBetaR9[i]->SetBinContent(2*j+1,hh_CBetaR9->GetBinContent(2*j+1,2*i+1));
      h_CBetaR9[i]->SetBinError  (2*j+1,hh_CBetaR9->GetBinError  (2*j+1,2*i+1));
    }
  }  
  // overlap the eta slices in a plot
  TCanvas *c_CBetaR9 = new TCanvas("c_CBetaR9","c_CBetaR9",600,600);
  c_CBetaR9->cd();
  c_CBetaR9->SetLeftMargin(0.2);
  for (Int_t i = 0; i<nBinsEta ; ++i){
    h_CBetaR9[i]-> SetMarkerStyle(20);
    h_CBetaR9[i]-> SetMarkerColor(i+22);
    h_CBetaR9[i]-> SetTitle(0);
    h_CBetaR9[i]-> GetXaxis()->SetTitle("R9");
    h_CBetaR9[i]-> GetYaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
    h_CBetaR9[i]-> GetYaxis()->SetTitleOffset(2);
    if(i == 0) {
      h_CBetaR9[i]-> Draw();
      h_CBetaR9[i]-> GetYaxis()->SetRangeUser(0.9,1.05);
    }
    else h_CBetaR9[i]-> Draw("same");
  }
  TBox *bandR9 = new TBox(leftR9[0],0.999,rightR9[nBinsR9-1],1.001);
  bandR9->SetFillColor(1);
  bandR9->SetFillStyle(3004);
  bandR9->Draw();
  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->SetFillColor(0);
  TString txt[nBinsEta];
  for (Int_t i = 0; i<nBinsEta ; ++i){
    txt[i] = floatToString(leftEta[i])+" < #eta < "+ floatToString(rightEta[i]);
    leg->AddEntry(h_CBetaR9[i],txt[i]);
  }
  leg->Draw();
  if      (applyCorrections == R9Eta  ) c_CBetaR9->Print("./plots_tmp/c_CBetaR9_corrR9Eta.png");
  else if (applyCorrections == R9EtaET) c_CBetaR9->Print("./plots_tmp/c_CBetaR9_corrR9EtaET.png");
  else                                  c_CBetaR9->Print("./plots_tmp/c_CBetaR9.png");


  // Fit the CB for each histogram of the ErecEgen ET-vector  *** BARREL ***
  //--------------------------------------------------------------------------------
  TCanvas *c_CBET_EB = new TCanvas("c_CBET_EB","c_CBET_EB",600,600);
  c_CBET_EB->cd();
  c_CBET_EB->SetLeftMargin(0.2);
  for (Int_t iET = 0; iET < nBinsET; ++iET){
    Double_t mean  = 1;
    Double_t emean = 0;
    Double_t chi2  = 0;
    Double_t fracBad = 0;
    if (h_ETcorrR9Eta_EB[iET]->GetEntries() > minEntries) fitCB(h_ETcorrR9Eta_EB[iET], mean, emean, chi2, fracBad);
    h_CBET_EB  ->SetBinContent(2*iET+1,mean  );
    h_CBET_EB  ->SetBinError  (2*iET+1,emean );
    h_Chi2ET_EB->SetBinContent(2*iET+1,chi2  ); // dof = 4 = number fof floating parametes
    h_EgenFracET_EB->SetBinContent(2*iET+1, fracBad);
  }
  h_CBET_EB  ->SetMarkerStyle(20);
  h_CBET_EB  ->SetMarkerColor(1);
  h_CBET_EB  ->SetTitle(0);
  h_CBET_EB  ->GetYaxis()->SetRangeUser(0.99,1.01);
  h_CBET_EB  ->GetXaxis()->SetTitle("ET [GeV]");
  h_CBET_EB  ->GetYaxis()->SetTitleOffset(2);
  h_CBET_EB  ->GetYaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
  h_CBET_EB  ->Draw();
  TBox *bandET_EB = new TBox(leftET[0],0.999,rightET[nBinsET-1],1.001);
  bandET_EB->SetFillColor(1);
  bandET_EB->SetFillStyle(3004);
  bandET_EB->Draw();
  if      (applyCorrections == R9Eta  ) c_CBET_EB->Print("./plots_tmp/c_CBET_EB_corrR9Eta.png");
  else if (applyCorrections == R9EtaET) c_CBET_EB->Print("./plots_tmp/c_CBET_EB_corrR9EtaET.png");
  else                                  c_CBET_EB->Print("./plots_tmp/c_CBET_EB.png");

  TCanvas *c_Chi2ET_EB = new TCanvas("c_Chi2ET_EB","c_Chi2ET_EB",600,600);
  c_Chi2ET_EB->cd();
  gPad->SetPhi(180+30);
  h_Chi2ET_EB->GetYaxis()->SetRangeUser(0.,5);
  h_Chi2ET_EB->Draw("");
  

  // Fit the CB for each histogram of the ErecEgen ET-vector  *** ENDCAP ***
  //--------------------------------------------------------------------------------
  TCanvas *c_CBET_EE = new TCanvas("c_CBET_EE","c_CBET_EE",600,600);
  c_CBET_EE->cd();
  c_CBET_EE->SetLeftMargin(0.2);
  for (Int_t iET = 0; iET < nBinsET; ++iET){
    Double_t mean  = 1;
    Double_t emean = 0;
    Double_t chi2  = 0;
    Double_t fracBad = 0;
    if (h_ETcorrR9Eta_EE[iET]->GetEntries() > minEntries) fitCB(h_ETcorrR9Eta_EE[iET], mean, emean, chi2, fracBad);
    h_CBET_EE  ->SetBinContent(2*iET+1,mean  );
    h_CBET_EE  ->SetBinError  (2*iET+1,emean );
    h_Chi2ET_EE->SetBinContent(2*iET+1,chi2  ); // dof = 4 = number fof floating parametes
    h_EgenFracET_EE->SetBinContent(2*iET+1, fracBad);
  }
  h_CBET_EE  ->SetMarkerStyle(20);
  h_CBET_EE  ->SetMarkerColor(1);
  h_CBET_EE  ->SetTitle(0);
  h_CBET_EE  ->GetYaxis()->SetRangeUser(0.99,1.01);
  h_CBET_EE  ->GetXaxis()->SetTitle("ET [GeV]");
  h_CBET_EE  ->GetYaxis()->SetTitleOffset(2);
  h_CBET_EE  ->GetYaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
  h_CBET_EE  ->Draw();
  TBox *bandET_EE = new TBox(leftET[0],0.999,rightET[nBinsET-1],1.001);
  bandET_EE->SetFillColor(1);
  bandET_EE->SetFillStyle(3004);
  bandET_EE->Draw();
  if      (applyCorrections == R9Eta  ) c_CBET_EE->Print("./plots_tmp/c_CBET_EE_corrR9Eta.png");
  else if (applyCorrections == R9EtaET) c_CBET_EE->Print("./plots_tmp/c_CBET_EE_corrR9EtaET.png");
  else                                  c_CBET_EE->Print("./plots_tmp/c_CBET_EE.png");

  TCanvas *c_Chi2ET_EE = new TCanvas("c_Chi2ET_EE","c_Chi2ET_EE",600,600);
  c_Chi2ET_EE->cd();
  gPad->SetPhi(180+30);
  h_Chi2ET_EE->GetYaxis()->SetRangeUser(0.,5);
  h_Chi2ET_EE->Draw("");
  

  return;
}


*/


















/*



//================================================================================
//
// test 3D : sigma_phi/sigma_eta, eta, ET in one step
//
//================================================================================

// matrix (br,eta,ET) of histograms used for the Erec/Egen CB fits
TH1F *h_ErecEgen_3D[nBinsBr][nBinsEta][nBinsET]; 

// 2D histogram containing the results of the CB fits
TH2F *hh_CBeta3D_etaBr[nBinsET];
//TH2F *hh_CBeta3D_BrET     = new TH2F("hh_CBetabr_BrET"     ,"hh_CBetabr_BrET"     ,nBinsBr *2-1, brbins , nBinsET*2-1, ETBins);
//TH2F *hh_CBeta3D_EtaET    = new TH2F("hh_CBetabr_EtaET"    ,"hh_CBetabr_EtaET"    ,nBinsEta*2-1, etabins, nBinsET*2-1, ETBins);

// 1D to fit the corrections
TH3F *hhh_CB3D =  new TH3F("hhh_CB3D"    ,"hhh_CB3D"    ,nBinsBr*2-1, brbins, nBinsEta*2-1, etabins,   nBinsET*2-1, ETBins);

// 2D histogram containing the chi2 of the CB fits
TH3F *hhh_Chi23D  = new TH3F("hhh_Chi23D"  ,"hhh_Chi23D"  ,nBinsBr*2-1, brbins, nBinsEta*2-1, etabins,   nBinsET*2-1, ETBins);

TH3F *hhh_EgenFrac3D  = new TH3F("hhh_EgenFrac3D"  ,"hhh_EgenFrac3D"  ,nBinsBr*2-1, brbins, nBinsEta*2-1, etabins,   nBinsET*2-1, ETBins);


// Extract f 3D
void scCorrections::run_3D()
{    
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  initHistograms_3D();
  correctionType Corr  = none;
  fillHistograms_3D(Corr);
  bool apply = false;
  getErecEgenFits_3D(Corr, apply);
  deriveCorrections_3D();
}

// Apply f 3D
void scCorrections::run_apply_3D()
{  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  initHistograms_3D();
  correctionType Corr  = C3D;
  fillHistograms_3D(Corr);
  bool apply = true;
  getErecEgenFits_3D(Corr, apply);
}

void scCorrections::initHistograms_3D()
{
  // 3D
  for (Int_t i = 0; i<nBinsBr; ++i){
    for (Int_t j = 0; j<nBinsEta; ++j){
      for (Int_t k = 0; k<nBinsET; ++k){      
	h_ErecEgen_3D[i][j][k] = new TH1F(Form("h_ErecEgen_3D_%d_%d_%d",i,j,k),Form("h_ErecEgen_3D_%d_%d_%d",i,j,k),nBinsErecEgen,ErecEgenMin, ErecEgenMax);
      }
    }
  }
  for (Int_t k = 0; k<nBinsET; ++k){          
    hh_CBeta3D_etaBr[k] = new TH2F(Form("hh_CBetabr_etaBr_%d",k),Form("hh_CBetabr_etaBr_%d",k) ,nBinsEta*2-1, etabins, nBinsBr*2-1, brbins);
  }
  cout << "initHistograms_3D(): done" << endl;
}
//
void scCorrections::fillHistograms_3D(correctionType applyCorrections )
{
  // read the Correction parameters
  //
  // 3D f(sigma_phi/sigma_eta, eta, ET ) from histogram 
  //
  TFile *histFile_brEtaET = 0;
  TH3F *h_corr = 0;
  if (applyCorrections == C3D){
    histFile_brEtaET = new TFile("./histCorrections_3D.root","READ");
    if (histFile_brEtaET->IsZombie()) {cout << "ERROR: ./histCorrections_3D.root does not exist " << endl; exit(-1);}
    //
    h_corr = (TH3F*)histFile_brEtaET->Get("hhh_CB3D");
  }


  // To speed up the loop read only the used variables from the NTuple
  // Add ALL variables used, else they are set to 0 and the compiler won't notice it
  //
  if (fChain == 0) return;
  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("phoSCEta"      ,1); 
  fChain->SetBranchStatus("phoSCPhi"      ,1); 
  fChain->SetBranchStatus("phoCetaCorrE"  ,1); 
  fChain->SetBranchStatus("phoCetaCorrEt" ,1); 
  fChain->SetBranchStatus("phoSCBrem"     ,1); 
  fChain->SetBranchStatus("phoR9"         ,1);       
  fChain->SetBranchStatus("mcE"           ,1);
  fChain->SetBranchStatus("mcEta"         ,1);
  fChain->SetBranchStatus("mcPhi"         ,1);
  fChain->SetBranchStatus("eleSCEta"      ,1);        
  fChain->SetBranchStatus("eleSCPhi"      ,1);        
  fChain->SetBranchStatus("eleCetaCorrE"  ,1);  
  fChain->SetBranchStatus("eleCetaCorrEt" ,1); 
  fChain->SetBranchStatus("eleBrLinear"   ,1);   
  fChain->SetBranchStatus("eleE3x3"       ,1);
  fChain->SetBranchStatus("eleCetaCorrE"  ,1);
  fChain->SetBranchStatus("eleSCRawEn"    ,1);

  // Loop over all entries:
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry % 1000 == 0) cout << jentry << endl;

    // Loop over the reconstructed objects (electrons or photons)
    Int_t nobj = 0;
    if      (partType == "electron") nobj = nEle;
    else if (partType == "photon"  ) nobj = nPho;
    
    // skip pathological events
    if (nobj < 1 || nobj > 2) continue;
    
    for (Int_t iobj = 0; iobj < nobj; ++iobj) {
      
      Double_t nt_em_eta        = 0;     
      Double_t nt_em_phi        = 0;
      Double_t nt_emCorrEta_e   = 0;
      Double_t nt_emCorrEta_et  = 0;
      Double_t nt_em_br1        = 0;
      Double_t nt_em_r9         = 0;
      Double_t nt_mc_e	        = 0; 
      
      // Use the correct naming of your ntuple
      if (partType == "photon"){
	nt_em_eta        = phoSCEta[iobj];        
	nt_em_phi        = phoSCPhi[iobj];        
	nt_emCorrEta_e   = phoCetaCorrE[iobj];  
	nt_emCorrEta_et  = phoCetaCorrEt[iobj]; 
	nt_em_br1        = phoSCBrem[iobj];     
	nt_em_r9         = phoR9[iobj]; 
	nt_mc_e          = mcE[0]; // mcE[0] = mcE[1] by contruction of the di-X ntuples
	//
	// Cheap truth match:
	// the photons/electrons are back to back, just need to compare one reco-photon to the two mc-photons
	// 	Float_t dr0 = TMath::Sqrt(pow(phoSCEta[0]-mcEta[0],2) + pow(phoSCPhi[0]-mcPhi[0],2));
	// 	Float_t dr1 = TMath::Sqrt(pow(phoSCEta[0]-mcEta[1],2) + pow(phoSCPhi[0]-mcPhi[1],2));
	// 	if (dr0 < dr1 ) nt_mc_e = mcE[0];            
	// 	else nt_mc_e = mcE[1];  
	// 	cout << phoSCEta[0] << " " << mcEta[0] <<  " " << phoSCPhi[0] << " " << mcPhi[0] <<  " " << mcE[0]  << endl;
	// 	cout << phoSCEta[0] << " " << mcEta[1] <<  " " << phoSCPhi[0] << " " << mcPhi[1] <<  " " << mcE[1]  << endl;
	// 	cout << nt_mc_e << endl;
      }
      else if (partType == "electron"){
	nt_em_eta        = eleSCEta[iobj];        
	nt_em_phi        = eleSCPhi[iobj];        
	nt_emCorrEta_e   = eleCetaCorrE[iobj];  
	nt_emCorrEta_et  = eleCetaCorrEt[iobj]; 
	nt_em_br1        = eleBrLinear[iobj];   
	if ( eleCetaCorrE[iobj] != 0) nt_em_r9 = eleE3x3[iobj]/eleSCRawEn[iobj]; 
	nt_mc_e          = mcE[0]; // mcE[0] = mcE[1] by contruction of the di-X ntuples
	//
	// Cheap truth match:
	// the photons/electrons are back to back, just need to compare one reco-photon to the two mc-photons
	// 	Float_t dr0 = TMath::Sqrt(pow(eleSCEta[0]-mcEta[0],2) + pow(eleSCPhi[0]-mcPhi[0],2));
	// 	Float_t dr1 = TMath::Sqrt(pow(eleSCEta[0]-mcEta[1],2) + pow(eleSCPhi[0]-mcPhi[1],2));
	// 	if (dr0 < dr1 ) nt_mc_e = mcE[0];            
	// 	else nt_mc_e = mcE[1];  
	// 	cout << eleSCEta[0] << " " << mcEta[0] <<  " " << eleSCPhi[0] << " " << mcPhi[0] <<  " " << mcE[0]  << endl;
	// 	cout << eleSCEta[0] << " " << mcEta[1] <<  " " << eleSCPhi[0] << " " << mcPhi[1] <<  " " << mcE[1]  << endl;
	// 	cout << nt_mc_e << endl;
      }

      // main event selection
      //---------------------------------
      if ( 
	  // select showering objects: 
	  ( (TMath::Abs(nt_em_eta) < etaCrackMin && (nt_em_r9 < minR9_EB)) || (TMath::Abs(nt_em_eta) > etaCrackMax && (nt_em_r9 < minR9_EE)) ) 
	  
	  &&
	  
	  // remove the lowest ET objects
	  nt_emCorrEta_et > minETCut    
	  
	  && 
	  
	  // remove the phi gaps between the super modules (only for the barrel):: local corrections will take care of them
	  // for the endcap there are no local corrections, so don't exclude any object
	  !isInPhiCracks(nt_em_phi, nt_em_eta,  fiducialPhiCut) 	  
	  
	   ){ 
	
	// Find the correct (br,eta,ET) bin for the h_ErecEgen_3D[nBinsBr][nBinsEta][nBinsET]
	Int_t iEta    = -1;
	Int_t iBr     = -1;
	Int_t iET     = -1;
	//
	for (Int_t iiEta = 0; iiEta < nBinsEta; ++iiEta){ 
	  if ( leftEta[iiEta] <= TMath::Abs(nt_em_eta) && TMath::Abs(nt_em_eta) <rightEta[iiEta]) {
	    iEta = iiEta; 
	    // cout << "ETA " << leftEta[iiEta]  << " " << TMath::Abs(nt_em_eta) << " " << rightEta[iiEta] << " " << iEta << endl;
	  }
	}
	for (Int_t iiBr  = 0; iiBr  < nBinsBr;  ++iiBr ){ 
	  if ( leftBr [iiBr]  <= nt_em_br1             &&             nt_em_br1 <rightBr [iiBr] ) {
	    iBr  = iiBr;  
	    // cout << "BR " << leftBr[iiBr]  << " " << nt_em_br1 << " " << rightBr[iiBr] << " " << iBr << endl;
	  }
	}
	// Find the correct (ET) bin for the h_ETcorrBrEta_EX[nBinsET]
	for (Int_t iiET  = 0; iiET  < nBinsET;  ++iiET ){ 
	  if ( leftET [iiET]  <= nt_emCorrEta_et      &&       nt_emCorrEta_et < rightET[iiET] ) 
	    iET  = iiET;  
	  // cout << "ET " << leftBr[iiET]  << " " << nt_emCorrEta_et << " " << rightBr[iiET] << " " << nt_emCorrEta_et << endl;
	}		
	if (iEta == -1 || iBr == -1 || iET == -1) continue;
	
	// Initialize the correction factors
	Double_t corr3D = 1.; 
	
	// get the corrections f(br,eta) and F(ET) from histogram
	//
	if (applyCorrections == C3D){
	  corr3D  = h_corr->GetBinContent(2*iBr+1, 2*iEta+1, 2*iET+1);
	}
	
	// fill the matrix of histograms h_ErecEgen_3D[nBinsBr][nBinsEta][nBinET]
	// ===========================================================
	//
	// *** f(sigma_phi/sigma_eta,eta) ***
	if (applyCorrections == C3D) {
	  // fill
	  if (nt_mc_e!=0 && corr3D!=0) h_ErecEgen_3D[iBr][iEta][iET]->Fill(nt_emCorrEta_e / nt_mc_e / corr3D);	  
	}	
	// *** no corrections ***
	else { 
	  // fill
	  if (nt_mc_e!=0 )  h_ErecEgen_3D[iBr][iEta][iET]->Fill(nt_emCorrEta_e / nt_mc_e );     
	}	
      } // if selected event
    } // loop over obj
  } // loop over entries	
  
  if (applyCorrections == C3D ) histFile_brEtaET->Close();    
  
  return;
  
}

void scCorrections::deriveCorrections_3D()
{
  TFile *histFile= new TFile("./histCorrections_3D.root", "RECREATE" );
  hhh_CB3D->Write();
  histFile->Close();

  return;
}

void scCorrections::getErecEgenFits_3D(correctionType applyCorrections, bool apply)
{
  Int_t minEntries      = 100;
  Int_t minChi2         = 3;

  // Fit the CB for each histogram of the ErecEgen matrix
  //--------------------------------------------------------------------------------
  for (Int_t iBr = 0; iBr < nBinsBr; ++iBr){
    for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){
      for (Int_t iET = 0; iET < nBinsET; ++iET){
	
	Double_t mean    = 1;
	Double_t emean   = 0;
	Double_t chi2    = 0;      
	Double_t fracBad = 0;
	if (h_ErecEgen_3D[iBr][iEta][iET]->GetEntries() >minEntries) 
	  {
	    fitCB(h_ErecEgen_3D[iBr][iEta][iET], mean, emean, chi2, fracBad);
	    if (!apply){	    
	      if (chi2 > minChi2){ // if the fit is bad
		if ( 2*iEta+1 -2 > 1){ // set it to the previous eta bin if possible
		  mean    = hhh_CB3D      ->GetBinContent(2*iBr+1  , 2*iEta+1  -2, 2*iET+1 );
		  emean   = hhh_CB3D      ->GetBinError  (2*iBr+1  , 2*iEta+1  -2, 2*iET+1 );
		  cout << "WARNING: BAD FIT set mean and error to the values of the previous eta bin " << mean << " " << emean << endl;
		  //getchar();	      
		}
		else if (2*iBr+1 -2 > 1 ){ // otherwise set it to the previous Br bin if possible
		  mean    = hhh_CB3D      ->GetBinContent(2*iBr+1 -2  , 2*iEta+1, 2*iET+1  );
		  emean   = hhh_CB3D      ->GetBinError  (2*iBr+1 -2  , 2*iEta+1, 2*iET+1  );	    
		  cout << "WARNING: BAD FIT set mean and error to the values of the previous br bin " << mean << " " << emean << endl;
		  //getchar();
		}
		else{  // otherwise set it to 1. Don't know how to treat them
		  mean  = 1;
		  emean = 0;
		  cout << "WARNING: BAD FIT set mean and error to " << mean << " " << emean << endl;
		  //getchar();
		}
	      }
	    }
	  }
	else { // if not enough statistics
	  if (!apply){	    
	    if  ( 2*iEta+1 -2 > 1){ // set it to the previous eta bin if possible
	      mean    = hhh_CB3D      ->GetBinContent(2*iBr+1   , 2*iEta+1 -2, 2*iET+1  );
	      emean   = hhh_CB3D      ->GetBinError  (2*iBr+1   , 2*iEta+1 -2, 2*iET+1  );
	      cout << "WARNING: NOT ENOUGH STAT set mean and error to the values of the previous eta bin " << mean << " " << emean << endl;
	      //getchar();
	    }
	    else if (2*iBr+1 -2 > 1 ){ // otherwise set it to the previous Br bin if possible
	      mean    = hhh_CB3D      ->GetBinContent(2*iBr+1 -2  , 2*iEta+1, 2*iET+1  );
	      emean   = hhh_CB3D      ->GetBinError  (2*iBr+1 -2  , 2*iEta+1, 2*iET+1  );
	      cout << "WARNING: NOT ENOUGH STAT set mean and error to the values of the previous br bin " << mean << " " << emean << endl;
	      //getchar();
	    }
	    else {  // otherwise set it to 1. Don't know how to treat them
	      mean  = 1;
	      emean = 0;
	      //getchar();
	      cout << "WARNING: NOT ENOUGH STAT set mean and error to " << mean << " " << emean << endl;
	    }
	  }
	}
	hhh_CB3D         ->SetBinContent(2*iBr+1, 2*iEta+1, 2*iET+1, mean   );
	hhh_CB3D         ->SetBinError  (2*iBr+1, 2*iEta+1, 2*iET+1, emean  );
	hhh_Chi23D       ->SetBinContent(2*iBr+1, 2*iEta+1, 2*iET+1, chi2   ); // dof = 4 = number dof floating parametes
	hhh_EgenFrac3D   ->SetBinContent(2*iBr+1, 2*iEta+1, 2*iET+1, fracBad);
	cout << " FRACTION BELOW " << EgenFrac << " = " << fracBad << endl;
      }
    }
  }

  // 2D plot CB most probable value in bins (eta,br)
  TCanvas *cc_CBetabr = new TCanvas("cc_CBetabr","cc_CBetabr",600,600);
  cc_CBetabr->cd();
  for (Int_t iET = 0; iET < nBinsET; ++iET) {    
    for (Int_t iEta = 0; iEta<nBinsEta; ++iEta){
      for (Int_t iBr = 0; iBr<nBinsBr; ++iBr){
	hh_CBeta3D_etaBr[iET]->SetBinContent(2*iEta+1,2*iBr+1,hhh_CB3D->GetBinContent(2*iBr+1, 2*iEta+1, 2*iET+1));
      }
    }
  }
  TCanvas *c = new TCanvas("c","c",600,600);
  c->Divide(3,3);
  for (Int_t iET = 0; iET < nBinsET; ++iET) {    
    hh_CBeta3D_etaBr[iET]->GetZaxis()->SetRangeUser(0.99,1.01);
    c->cd(iET+1);
    gPad->SetPhi(180+30);
    hh_CBeta3D_etaBr[iET]->Draw("lego2");    
  }

  return;
}


*/






















//================================================================================
//
// test f(sigma_phi/sigma_eta, eta) * F(E) : energy instead of ET
//
//================================================================================

// vector E dependence 
TH1F *h_EcorrBrEta_EB[nBinsE];
TH1F *h_CBE_EB    = new TH1F("h_CBE_EB"    ,"h_CBE_EB"    ,nBinsE*2-1, EBins);
TH1F *h_Chi2E_EB  = new TH1F("h_Chi2E_EB"  ,"h_Chi2E_EB"  ,nBinsE*2-1, EBins);
TH1F *h_EcorrBrEta_EE[nBinsE];
TH1F *h_CBE_EE    = new TH1F("h_CBE_EE"    ,"h_CBE_EE"    ,nBinsE*2-1, EBins);
TH1F *h_Chi2E_EE  = new TH1F("h_Chi2E_EE"  ,"h_Chi2E_EE"  ,nBinsE*2-1, EBins);



// Extract f(sigma_phi,sigma_eta, E)
void scCorrections::run_brEta_E() // bad name it's only brEta
{    
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  initHistograms_E();
  correctionType Corr  = none;
  fillHistograms_E(Corr);
  bool apply = false;
  getErecEgenFits_BrEta_E(Corr, apply);
  deriveCorrections_BrEta_E();

}

// Apply f(sigma_phi/sigma_eta, eta) (closure test)
void scCorrections::run_apply_brEta_E()
{  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  initHistograms_E();
  correctionType Corr  = brEta;
  fillHistograms_E(Corr);
  bool apply = true;
  getErecEgenFits_BrEta_E(Corr, apply);
}

// Extract the ET correction after applying f(sigma_phi,sigma_eta, eta)
void scCorrections::run_brEtaE()
{  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  // extract the BrEta corrections
  run_brEta_E();

  // apply the BrEta correction
  run_apply_brEta_E();
  
  // write out the correction macro
  deriveCorrections_BrEtaE();  
}

// Apply f(sigma_phi/sigma_eta, eta)*F(E) (closure test)
void scCorrections::run_apply_brEtaE()
{  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  //
  initHistograms_E();
  correctionType Corr  = brEtaE;
  fillHistograms_E(Corr);
  bool apply = true;
  getErecEgenFits_BrEta_E(Corr, apply);
}

//
void scCorrections::initHistograms_E()
{
  // Br
  for (Int_t i = 0; i<nBinsEta; ++i){
    for (Int_t j = 0; j<nBinsBr; ++j){
      h_ErecEgen[i][j] = new TH1F(Form("h_ErecEgen_%d_%d",i,j),Form("h_ErecEgen_%d_%d",i,j),nBinsErecEgen,ErecEgenMin, ErecEgenMax);
    }
  }
  for (Int_t i = 0; i<nBinsEta; ++i){
    h_CBetabr[i] = new TH1F(Form("h_CBetabr_%d",i), Form("h_CBetabr_%d",i),nBinsBr*2-1, brbins);
  }
  // E - Br
  for (Int_t i = 0; i<nBinsE; ++i){
    h_EcorrBrEta_EB[i] = new TH1F(Form("h_EcorrBrEta_EB_%d",i), Form("h_EcorrBrEta_EB_%d",i), nBinsErecEgen,ErecEgenMin, ErecEgenMax);
    h_EcorrBrEta_EE[i] = new TH1F(Form("h_EcorrBrEta_EE_%d",i), Form("h_EcorrBrEta_EE_%d",i), nBinsErecEgen,ErecEgenMin, ErecEgenMax);
  }
  cout << "initHistograms(): done" << endl;
}

//
void scCorrections::fillHistograms_E(correctionType applyCorrections )
{
  // read the Correction parameters
  //
  // f(sigma_phi/sigma_eta, eta ) from histogram 
  //
  TFile *histFile_brEta = 0;
  TH1F *h_corr[nBinsEta];
  if (applyCorrections == brEta || applyCorrections == brEtaE){
    histFile_brEta = new TFile("./histCorrections_brEta_E.root","READ");
    if (histFile_brEta->IsZombie()) {cout << "ERROR: ./histCorrections_brEta_E.root does not exist " << endl; exit(-1);}
    //
    for (Int_t i = 0; i<nBinsEta; ++i){
      h_corr[i] = (TH1F*)histFile_brEta->Get(Form("h_corr_%d",i));
    }
  }  
  //
  // F(E) correction from histogram ( *** BARREL / ENDCAP *** )
  //
  TFile *histFile_E = 0;
  TH1F *h_corrE_EB = 0;
  TH1F *h_corrE_EE = 0;
  if (applyCorrections == brEtaE){
    histFile_E = new TFile("./histCorrections_E.root","READ");
    if (histFile_E->IsZombie()) {cout << "ERROR: ./histCorrections_E.root does not exist " << endl; exit(-1);}
    //    
    h_corrE_EB = (TH1F*)histFile_E->Get("h_CBE_EB");
    h_corrE_EE = (TH1F*)histFile_E->Get("h_CBE_EE");
  }

  // To speed up the loop read only the used variables from the NTuple
  // Add ALL variables used, else they are set to 0 and the compiler won't notice it
  //
  if (fChain == 0) return;
  fChain->SetBranchStatus("*",1);
 

  // Loop over all entries:
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry % 1000 == 0) cout << jentry << endl;

     // Loop over the reconstructed objects (electrons or photons)
    Int_t nobj = 0;
    if      (partType == "electron") nobj = NEles;
    else if (partType == "photon"  ) nobj = NPhotons;

    // skip pathological events
    if (nobj < 1 || nobj > 2) continue;

    {
      bool goodevt = true;
      if      (partType == "electron") { 
	if (ElSCindex[0]==-1 || ElSCindex[1]==-1) goodevt=false;
	if (ElGenE[0]<0 || ElGenE[1]<0) goodevt=false;
      }
      else if (partType == "photon"  ) { 
	if (PhotSCindex[0]==-1 || PhotSCindex[1]==-1) goodevt=false;
	if (PhoMCmatchexitcode[0]<=0 || PhoMCmatchexitcode[1]<=0) goodevt=false;
      }
      if (!goodevt) continue;
    }
    
    for (Int_t iobj = 0; iobj < nobj; ++iobj) {
      
      Double_t nt_em_eta        = 0;     
      Double_t nt_em_phi        = 0;
      Double_t nt_emCorrEta_e   = 0;
      Double_t nt_emCorrEta_et  = 0;
      Double_t nt_em_br1        = 0;
      Double_t nt_em_r9         = 0;
      Double_t nt_mc_e	        = 0; 
      
      // Use the correct naming of your ntuple
      if (partType == "photon"){
	int scindex = PhotSCindex[iobj];
	nt_em_eta        = SCEta[scindex];        
	nt_em_phi        = SCPhi[scindex];        
	nt_emCorrEta_e   = (fabs(nt_em_eta)<1.4442) ? SCRaw[scindex]*getEtaCorrectionBarrel(SCEta[scindex]) : SCRaw[scindex]+SCPre[scindex]; 
	nt_emCorrEta_et  = nt_emCorrEta_e/cosh(nt_em_eta);
	nt_em_br1        = SCBrem[scindex];     
	nt_em_r9         = PhoR9[iobj]; 
	nt_mc_e          = GenPhotonPt[PhoMCmatchindex[iobj]]*cosh(GenPhotonEta[PhoMCmatchindex[iobj]]);
	//
	// Cheap truth match:
	// the photons/electrons are back to back, just need to compare one reco-photon to the two mc-photons
	// 	Float_t dr0 = TMath::Sqrt(pow(phoSCEta[0]-mcEta[0],2) + pow(phoSCPhi[0]-mcPhi[0],2));
	// 	Float_t dr1 = TMath::Sqrt(pow(phoSCEta[0]-mcEta[1],2) + pow(phoSCPhi[0]-mcPhi[1],2));
	// 	if (dr0 < dr1 ) nt_mc_e = mcE[0];            
	// 	else nt_mc_e = mcE[1];  
	// 	cout << phoSCEta[0] << " " << mcEta[0] <<  " " << phoSCPhi[0] << " " << mcPhi[0] <<  " " << mcE[0]  << endl;
	// 	cout << phoSCEta[0] << " " << mcEta[1] <<  " " << phoSCPhi[0] << " " << mcPhi[1] <<  " " << mcE[1]  << endl;
	// 	cout << nt_mc_e << endl;
      }
      else if (partType == "electron"){
	int scindex = ElSCindex[iobj];
	nt_em_eta        = SCEta[scindex];        
	nt_em_phi        = SCPhi[scindex];        
	nt_emCorrEta_e   = (fabs(nt_em_eta)<1.4442) ? SCRaw[scindex]*getEtaCorrectionBarrel(SCEta[scindex]) : SCRaw[scindex]+SCPre[scindex]; 
	nt_emCorrEta_et  = nt_emCorrEta_e/cosh(nt_em_eta);
	nt_em_br1        = SCBrem[scindex];     
	nt_em_r9         = SCR9[scindex]; 
	nt_mc_e          = ElGenE[iobj];
	//
	// Cheap truth match:
	// the photons/electrons are back to back, just need to compare one reco-photon to the two mc-photons
	// 	Float_t dr0 = TMath::Sqrt(pow(eleSCEta[0]-mcEta[0],2) + pow(eleSCPhi[0]-mcPhi[0],2));
	// 	Float_t dr1 = TMath::Sqrt(pow(eleSCEta[0]-mcEta[1],2) + pow(eleSCPhi[0]-mcPhi[1],2));
	// 	if (dr0 < dr1 ) nt_mc_e = mcE[0];            
	// 	else nt_mc_e = mcE[1];  
	// 	cout << eleSCEta[0] << " " << mcEta[0] <<  " " << eleSCPhi[0] << " " << mcPhi[0] <<  " " << mcE[0]  << endl;
	// 	cout << eleSCEta[0] << " " << mcEta[1] <<  " " << eleSCPhi[0] << " " << mcPhi[1] <<  " " << mcE[1]  << endl;
	// 	cout << nt_mc_e << endl;
      }


      // main event selection
      //---------------------------------
      if ( 
	  // select showering objects: 
	  ( (TMath::Abs(nt_em_eta) < etaCrackMin && (nt_em_r9 < minR9_EB)) || (TMath::Abs(nt_em_eta) > etaCrackMax && (nt_em_r9 < minR9_EE)) ) 
	  
	  &&
	  
	  // remove the lowest ET objects
	  nt_emCorrEta_et > minETCut    
	  
	  && 
	  
	  // remove the phi gaps between the super modules (only for the barrel):: local corrections will take care of them
	  // for the endcap there are no local corrections, so don't exclude any object
	  !isInPhiCracks(nt_em_phi, nt_em_eta,  fiducialPhiCut) 	  
	  
	  &&
	  // this is to get exactly the same br,eta corrections as in the ET case
	  (leftET [0]  <= nt_emCorrEta_et      && nt_emCorrEta_et < rightET[nBinsET-1] )
	  
	   ){ 
	
	// Find the correct (br,eta) bin for the h_ErecEgen[nBinsEta][nBinsBr] and h_EcorrBrEta_EB[nBinsE];
	Int_t iEta    = -1;
	Int_t iBr     = -1;
	Int_t iE      = -1;
	Int_t iEcorr  = -1;	
	//
	for (Int_t iiEta = 0; iiEta < nBinsEta; ++iiEta){ 
	  if ( leftEta[iiEta] <= TMath::Abs(nt_em_eta) && TMath::Abs(nt_em_eta) <rightEta[iiEta]) {
	    iEta = iiEta; 
	    // cout << "ETA " << leftEta[iiEta]  << " " << TMath::Abs(nt_em_eta) << " " << rightEta[iiEta] << " " << iEta << endl;
	  }
	}
	for (Int_t iiBr  = 0; iiBr  < nBinsBr;  ++iiBr ){ 
	  if ( leftBr [iiBr]  <= nt_em_br1             &&             nt_em_br1 <rightBr [iiBr] ) {
	    iBr  = iiBr;  
	    // cout << "BR " << leftBr[iiBr]  << " " << nt_em_br1 << " " << rightBr[iiBr] << " " << iBr << endl;
	  }
	}
	// Find the correct (E) bin for the h_EcorrBrEta_EX[nBinsE]
	for (Int_t iiE  = 0; iiE  < nBinsE;  ++iiE ){ 
	  if ( leftE [iiE]  <= nt_emCorrEta_e      &&       nt_emCorrEta_e < rightE[iiE] ) 
	    iE  = iiE;  
	  // cout << "E " << leftBr[iiE]  << " " << nt_emCorrEta_e << " " << rightBr[iiE] << " " << nt_emCorrEta_e << endl;
	}		
	if (iEta == -1 || iBr == -1 || iE == -1) continue;


	// Initialize the correction factors
	Double_t corrBrEta = 1.; 
	Double_t corrEt    = 1.; 	
	
 	// get the corrections f(br,eta) and F(E) from histogram                                                                               //MDDB
 	//																       //MDDB
 	if (applyCorrections == brEta || applyCorrections == brEtaE) {									       //MDDB
 	  corrBrEta  = h_corr[iEta]->GetBinContent(2*iBr+1);                               						       //MDDB
 	  //																       //MDDB
 	  // Find the correct (E corrected with f(br,Eta) ) bin for the correction histogram 						       //MDDB
 	  for (Int_t iiE  = 0; iiE  < nBinsE;  ++iiE ){ 										       //MDDB
 	    if ( leftE [iiE]  <= (nt_emCorrEta_e / corrBrEta)      &&       (nt_emCorrEta_e / corrBrEta) < rightE[iiE] ) {		       //MDDB
 	      iEcorr  = iiE;  														       //MDDB
 	      // cout << "E " << leftE [iiE] << " " <<  nt_emCorrEta_e / corrBrEta  << " " <<  rightE[iiE] <<  " " << iEcorr <<endl;	       //MDDB
 	    }																       //MDDB
 	  }																       //MDDB
 	  if (iEcorr == -1) continue; // always skip the event, there is something patological	    					       //MDDB
 	  if (applyCorrections == brEtaE) {												       //MDDB
 	    //																       //MDDB
 	    // histogram correction E *** BARREL, ENDCAP ***										       //MDDB
 	    if (TMath::Abs(nt_em_eta) < etaCrackMin ){                    								       //MDDB
 	      corrEt = h_corrE_EB    ->GetBinContent(2*iEcorr+1);									       //MDDB
 	      //	    cout << " corr " << iEcorr << " " << h_corrE_EB    ->GetBinContent(2*iEcorr+1) <<endl;			       //MDDB
 	    }																       //MDDB
 	    if (TMath::Abs(nt_em_eta) > etaCrackMax ) {											       //MDDB
 	      corrEt = h_corrE_EE    ->GetBinContent(2*iEcorr+1);									       //MDDB
 	      //	    cout << "corr " << iEcorr << " " << h_corrE_EE    ->GetBinContent(2*iEcorr+1) <<endl;			       //MDDB
 	    }																       //MDDB
 	  }																       //MDDB
 	}																       //MDDB

//MDDB 	// get the corrections f(br,eta) and F(E) from File                                                                                    
//MDDB 	//
//MDDB 	if (applyCorrections == brEta || applyCorrections == brEtaE) {
//MDDB 	  corrBrEta = applyScCorrectionsBrEta(nt_em_eta,  nt_em_br1);                           
//MDDB 	  
//MDDB 	  // Find the correct (E corrected with f(br,Eta) ) bin for the correction histogram 
//MDDB 	  for (Int_t iiE  = 0; iiE  < nBinsE;  ++iiE ){ 
//MDDB 	    if ( leftE [iiE]  <= (nt_emCorrEta_e / corrBrEta)      &&       (nt_emCorrEta_e / corrBrEta) < rightE[iiE] ) {
//MDDB 	      iEcorr  = iiE;  
//MDDB 	      // cout << "E " << leftE [iiE] << " " <<  nt_emCorrEta_e / corrBrEta  << " " <<  rightE[iiE] <<  " " << iEcorr <<endl;
//MDDB 	    }
//MDDB 	  }	
//MDDB 	  if (iEcorr == -1) continue; // always skip the event, there is something patological	   
//MDDB 	  if (applyCorrections == brEtaE) {
//MDDB 	    if (TMath::Abs(nt_em_eta) < etaCrackMin ) corrEt = applyScCorrectionsE_EB( nt_emCorrEta_e / corrBrEta );   
//MDDB 	    if (TMath::Abs(nt_em_eta) > etaCrackMax ) corrEt = applyScCorrectionsE_EE( nt_emCorrEta_e / corrBrEta );   
//MDDB 	  }
//MDDB 	}

	// fill the matrix of histograms h_ErecEgen[nBinsEta][nBinsBr] 
	// ===========================================================
	//
	// *** f(sigma_phi/sigma_eta,eta) ***
	if (applyCorrections == brEta) {
	  // fill
	  if (nt_mc_e!=0 && corrBrEta!=0) h_ErecEgen[iEta][iBr]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta);	  
	}	
	// *** f(sigma_phi/sigma_eta,eta)*F(E) ***
	else if (applyCorrections == brEtaE) {
	  // fill
	  if (nt_mc_e!=0 && corrEt!=0) h_ErecEgen[iEta][iBr]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta /corrEt);	  
	}	
	// *** no corrections ***
	else { 
	  // fill
	  if (nt_mc_e!=0 )  h_ErecEgen[iEta][iBr]->Fill(nt_emCorrEta_e / nt_mc_e );     
	}
	
	// fill the vector of E dependence histograms
	// ===========================================
	//	
	// *** f(sigma_phi/sigma_eta,eta) ***
	if (applyCorrections == brEta) {
	  if (TMath::Abs(nt_em_eta) < etaCrackMin ){
	    // fill
	    if (nt_mc_e!=0 && corrBrEta!=0)    h_EcorrBrEta_EB[iE]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta);
	  }
	  else if (TMath::Abs(nt_em_eta) > etaCrackMax ) {
	    // fill
	    if (nt_mc_e!=0 && corrBrEta!=0)    h_EcorrBrEta_EE[iE]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta);		
	  }
	}
	// *** f(sigma_phi/sigma_eta,eta)*F(E) ***
	else if (applyCorrections == brEtaE) {		
	  // bin position for the br,eta corrected E
	  if (TMath::Abs(nt_em_eta) < etaCrackMin ){
	    
	    // fill
	    if (nt_mc_e!=0 && corrEt!=0) h_EcorrBrEta_EB[iEcorr]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta / corrEt );	   
	  }
	  if (TMath::Abs(nt_em_eta) > etaCrackMax ) {
	    // fill
	    if (nt_mc_e!=0 && corrEt!=0) h_EcorrBrEta_EE[iEcorr]->Fill(nt_emCorrEta_e / nt_mc_e / corrBrEta / corrEt);	   
	  }	      
	}	    
	// *** no corrections ***
	else{
	  if (TMath::Abs(nt_em_eta) < etaCrackMin ){
	    h_EcorrBrEta_EB[iE]->Fill(nt_emCorrEta_e / nt_mc_e );
	  }
	  if (TMath::Abs(nt_em_eta) > etaCrackMax ) { 
	    h_EcorrBrEta_EE[iE]->Fill(nt_emCorrEta_e / nt_mc_e );
	  }
	}

      } // if selected event
    } // loop over obj
  } // loop over entries	
  
  if (applyCorrections == brEta || applyCorrections == brEtaE) histFile_brEta->Close();    
  if (applyCorrections == brEtaE) histFile_E->Close();
  
  return;
  
}

void scCorrections::deriveCorrections_BrEta_E()
{
  // Write out the correction macro
  ofstream outfile2;  
  TString filename = "./applyScCorrections.C";
  outfile2.open(filename);
  outfile2 << "#ifndef scCorrectionsBrETAE_h                                              " << endl;
  outfile2 << "#define scCorrectionsBrETAE_h                                              " << endl;
  outfile2 << "#include \"TFile.h\"                                                      " << endl;
  outfile2 << "#include \"TH1F.h\"				                       " << endl;
  outfile2 << "#include <iostream>  			                                " << endl;
  outfile2 << "  			                                                " << endl;
  outfile2 << "bool DBG = false;			                                " << endl;
  outfile2 << "  	                                                                " << endl;
  outfile2 << "Double_t applyScCorrectionsBrEta(Double_t eta, Double_t sigmaPhiSigmaEta){  " << endl;
  outfile2 << "  	                                                               " << endl;
  outfile2 << "  const Int_t    nBinsEta = " << nBinsEta << ";                    " << endl;                                  
  outfile2 << "  Double_t       leftEta  ["  << nBinsEta << "];                   " << endl;
  outfile2 << "  Double_t       rightEta ["  << nBinsEta << "];                   " << endl;
  outfile2 << "  const Int_t    nBinsBr  = " << nBinsBr << ";                     " << endl;                                    
  outfile2 << "  Double_t       leftBr  [" << nBinsBr << "];                      " << endl;
  outfile2 << "  Double_t       rightBr [" << nBinsBr << "];                      " << endl;
  outfile2 << "  Double_t brbins  [" << nBinsBr*2 << "];                    " << endl;
  outfile2 << "  TH1F *h_corr["<< nBinsEta << "];                                 " << endl;
  outfile2 << "  	                                                               " << endl;
  for (Int_t i = 0; i< nBinsEta; ++i){ 
    outfile2 << "  leftEta[" << i << "] =  " << leftEta[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsEta; ++i){ 
    outfile2 << "  rightEta[" << i<< "] =  " << rightEta[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsBr; ++i){ 
    outfile2 << "  leftBr[" << i << "] =  " << leftBr[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsBr; ++i){ 
    outfile2 << "  rightBr["<< i << "] =  " << rightBr[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< 2*nBinsBr; ++i){ 
    outfile2 << "  brbins[" << i << "] =  " << brbins[i] << " ; " << endl;
  }
  outfile2 << "  for (Int_t i = 0; i<nBinsEta; ++i){                                                   " << endl;
  outfile2 << "    h_corr[i] = new TH1F(Form(\"h_corr_%d\",i),Form(\"h_corr_%d\",i),nBinsBr*2-1, brbins);  " << endl;
  outfile2 << "  }                                                                                     " << endl;
  for (Int_t i = 0; i<nBinsEta ; ++i){
    for (Int_t j = 0; j< nBinsBr; ++j){
      outfile2 << "  h_corr[" << i << "]->SetBinContent(" <<  2*j+1 << "," << h_CBetabr[i]->GetBinContent(2*j+1) << " );"<< endl;
    }
    outfile2 << endl;
  }

  outfile2 << "     Int_t tmpEta = -1;                                                                                                                                                                        " << endl;
  outfile2 << "     Int_t tmpBr  = -1;																					      " << endl;
  outfile2 << "   																							      " << endl;
  outfile2 << "     // Extra protections										       										      " << endl;
  outfile2 << "      																							      " << endl;
  outfile2 << "     if (TMath::Abs(eta)  <   leftEta[0]           ) { tmpEta = 0;          if (DBG) std::cout << \" WARNING [applyScCorrections]: (TMath::Abs(eta)  <   leftEta[0]          \" << std::endl;}   " << endl;
  outfile2 << "     if (TMath::Abs(eta)  >=  rightEta[nBinsEta-1] ) { tmpEta = nBinsEta-1; if (DBG) std::cout << \" WARNING [applyScCorrections]: TMath::Abs(eta)  >=  rightEta[nBinsEta-1] \" << std::endl;}   " << endl;
  outfile2 << "   																							      " << endl;
  outfile2 << "     if (sigmaPhiSigmaEta <  leftBr [0]            ) {tmpBr = 0;            if (DBG) std::cout << \" WARNING [applyScCorrections]: sigmaPhiSigmaEta <  leftBr [0]            \" << std::endl;}   " << endl;
  outfile2 << "     if (sigmaPhiSigmaEta >= rightBr[nBinsBr]      ) {tmpBr = nBinsBr  -1;  if (DBG) std::cout << \" WARNING [applyScCorrections]: sigmaPhiSigmaEta >= rightBr[nBinsBr]      \" << std::endl;}   " << endl;    
  outfile2 << "                                                                                                                                                                                               " << endl;
  outfile2 << "     for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								              								      " << endl;
  outfile2 << "       if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] ){				       									      " << endl;
  outfile2 << "         tmpEta = iEta;											       										      " << endl;
  outfile2 << "       }													       										      " << endl;
  outfile2 << "     }													         									      " << endl;
  outfile2 << "   																							      " << endl;
  outfile2 << "     for (Int_t iSigmaPhiSigmaEta = 0; iSigmaPhiSigmaEta < nBinsBr; ++iSigmaPhiSigmaEta){			       									      " << endl;
  outfile2 << "       if (leftBr [iSigmaPhiSigmaEta]  <= sigmaPhiSigmaEta && sigmaPhiSigmaEta <rightBr [iSigmaPhiSigmaEta] ){      									      " << endl;
  outfile2 << "         tmpBr = iSigmaPhiSigmaEta;										       									      " << endl;
  outfile2 << "       }													       										      " << endl;
  outfile2 << "     }													       										      " << endl;
  outfile2 << "     																							      " << endl;
  outfile2 << "     // Interpolation																					      " << endl;
  outfile2 << "     Double_t tmpInter = 1;																				      " << endl;
  outfile2 << "     // In eta cracks/gaps 																				      " << endl;
  outfile2 << "     if (tmpEta == -1 && tmpBr != -1 ) { // need to interpolate only eta, if br is out of range skip this										      " << endl;
  outfile2 << "       																							      " << endl;
  outfile2 << "       if (TMath::Abs(eta) >=rightEta[nBinsEta-1] ) { // out of ECAL boundary														      " << endl;
  outfile2 << "         for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];    															      " << endl;
  outfile2 << "         return 1; // don't correct																			      " << endl;
  outfile2 << "       }																							      " << endl;
  outfile2 << "       																							      " << endl;
  outfile2 << "       // central bin at eta = 0																				      " << endl;
  outfile2 << "       if (0 <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[0] ) {															      " << endl;
  outfile2 << "         																						      " << endl;
  outfile2 << "         tmpInter = h_corr[0]->GetBinContent(2*tmpBr+1);																	      " << endl;
  outfile2 << "         																						      " << endl;
  outfile2 << "       }																							      " << endl;
  outfile2 << "       else { // all other crack-bins																			      " << endl;
  outfile2 << "         																						      " << endl;
  outfile2 << "         for (Int_t iEta = 0; iEta < nBinsEta-1; ++iEta){								       								      " << endl;
  outfile2 << "   	if (rightEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <leftEta[iEta+1] ){													      " << endl;
  outfile2 << "   	  tmpInter = ( h_corr[iEta]  ->GetBinContent(2*tmpBr+1) + 															      " << endl;
  outfile2 << "   		       h_corr[iEta+1]->GetBinContent(2*tmpBr+1) ) / 2. ;														      " << endl;
  outfile2 << "   	}																						      " << endl;
  outfile2 << "         }																						      " << endl;
  outfile2 << "       }																							      " << endl;
  outfile2 << "       for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];    															      " << endl;
  outfile2 << "       return tmpInter;																					      " << endl;
  outfile2 << "     }  																							      " << endl;
  outfile2 << "     // end interpolation                                                                                                                                                                      " << endl;

//   outfile2 << "   Int_t tmpEta = -1;                                                                                             " << endl;         
//   outfile2 << "   for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								       " << endl;
//   outfile2 << "     if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] ){				       " << endl;
//   outfile2 << "       tmpEta = iEta;											       " << endl;
//   outfile2 << "     }													       " << endl;
//   outfile2 << "   }													       " << endl;
//   outfile2 << "   													       " << endl;
//   outfile2 << "   Int_t tmpBr = -1;											       " << endl;
//   outfile2 << "   for (Int_t iSigmaPhiSigmaEta = 0; iSigmaPhiSigmaEta < nBinsBr; ++iSigmaPhiSigmaEta){			       " << endl;
//   outfile2 << "     if (leftBr [iSigmaPhiSigmaEta]  <= sigmaPhiSigmaEta && sigmaPhiSigmaEta <rightBr [iSigmaPhiSigmaEta] ){      " << endl;	   
//   outfile2 << "       tmpBr = iSigmaPhiSigmaEta;										       " << endl;
//   outfile2 << "     }													       " << endl;
//   outfile2 << "   }	                											       " << endl;

  outfile2 << "   if (tmpEta == -1 || tmpBr == -1){									       " << endl;
  outfile2 << "     for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  		       " << endl;
  outfile2 << "     return  1; // don't correct											       " << endl;
  outfile2 << "   }													       " << endl;
  outfile2 << "   													       " << endl;
  outfile2 << "   Double_t tmp = h_corr[tmpEta]->GetBinContent(2*tmpBr+1);				   		       " << endl;
  outfile2 << "   for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  		       " << endl;
  outfile2 << "   return  tmp;                                                                                                   " << endl;                         
  outfile2 << "}													   " << endl;  

  //  outfile2 << "  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){								   " << endl;
  //  outfile2 << "    for (Int_t iSigmaPhiSigmaEta = 0; iSigmaPhiSigmaEta < nBinsBr; ++iSigmaPhiSigmaEta){			   " << endl;
  //  outfile2 << "      													   " << endl;
  //  outfile2 << "      // select the right (eta,br) bin in the matrix							   " << endl;
  //  outfile2 << "      if ( leftEta[iEta] <= TMath::Abs(eta) && TMath::Abs(eta) <rightEta[iEta] &&			   " << endl;
  //  outfile2 << " 	   leftBr [iSigmaPhiSigmaEta]  <= sigmaPhiSigmaEta && sigmaPhiSigmaEta <rightBr [iSigmaPhiSigmaEta] ){	   " << endl;
  //  outfile2 << "	       Double_t tmp = h_corr[iEta]->GetBinContent(2*iSigmaPhiSigmaEta+1);				   " << endl;	   
  //  outfile2 << "	       for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  " << endl;
  //  outfile2 << "	       return  tmp;                                                                                        " << endl;
  //  outfile2 << "      }													   " << endl;
  //  outfile2 << "    } 													   " << endl;
  //  outfile2 << "  }													   " << endl;
  //  outfile2 << "	       for (Int_t i = 0; i < nBinsEta; ++i) delete h_corr[i];                                  " << endl;
  //  outfile2 << "  return 1;												   " << endl;
  //  outfile2 << "}													   " << endl;
  
  outfile2.close();                                                                                                      
  
  // save also the histograms to be used as correction
  TFile *histFile= new TFile("./histCorrections_brEta_E.root", "RECREATE" );
  TH1F *h_corr[nBinsEta];
  for (Int_t i = 0; i<nBinsEta ; ++i){
    h_corr[i] = (TH1F*)h_CBetabr[i]->Clone(Form("h_corr_%d",i));
    h_corr[i]->Write();
  }
  histFile->Close();
    
  TCanvas *cc_BadClustersEtaBr = new TCanvas("cc_BadClustersEtaBr" ,"cc_BadClustersEtaBr", 600,600);
  cc_BadClustersEtaBr->cd();
  gPad->SetPhi(270+30);
  hh_EgenFracEtaBr->GetXaxis()->SetTitle("#sigma_{#phi} / #sigma_{#eta}");
  hh_EgenFracEtaBr->GetXaxis()->SetTitleOffset(1.6);
  hh_EgenFracEtaBr->GetYaxis()->SetTitle("#eta");
  hh_EgenFracEtaBr->GetYaxis()->SetTitleOffset(1.6);
  hh_EgenFracEtaBr->Draw("lego2");
  TCanvas *c_BadClustersET_EB     = new TCanvas("c_BadClustersET_EB"     ,"c_BadClustersET_EB"    , 600,600);
  c_BadClustersET_EB->cd();
  h_EgenFracET_EB ->GetXaxis()->SetTitle("ET [GeV]");
  h_EgenFracET_EB ->Draw("");
  TCanvas *c_BadClustersET_EE     = new TCanvas("c_BadClustersET_EE"     ,"c_BadClustersET_EE"    , 600,600);
  c_BadClustersET_EE->cd();
  h_EgenFracET_EE ->GetXaxis()->SetTitle("ET [GeV]");
  h_EgenFracET_EE ->Draw("");

  return;
}

void scCorrections::deriveCorrections_BrEtaE()
{
  // Add to the existing applyCorrection the E dependence
  ofstream outfile2;  
  TString filename = "./applyScCorrections.C";
  outfile2.open(filename,ios::app);

  outfile2 << "                                                                         " << endl;
  outfile2 << "Double_t applyScCorrectionsE_EB(Double_t E){  							   " << endl;
  outfile2 << "  const Int_t    nBinsE             = " << nBinsE<< ";             " << endl; 
  outfile2 << "  Double_t       leftE  [nBinsE];                    " << endl;
  outfile2 << "  Double_t       rightE [nBinsE];                    " << endl;
  outfile2 << "  Double_t       EBins  [nBinsE*2];                  " << endl;
  // E
  for (Int_t i = 0; i< nBinsE; ++i){ 
    outfile2 << "  leftE[" << i << "] =  " << leftE[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsE; ++i){ 
    outfile2 << "  rightE[" << i << "] =  " << rightE[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< 2*nBinsE; ++i){ 
    outfile2 << "  EBins[" << i << "] =  " << EBins[i] << " ; " << endl;
  }
  outfile2 << "  TH1F *h_CBE_EB    = new TH1F(\"h_CBE_EB\"    ,\"h_CBE_EB\"    ,nBinsE*2-1, EBins); " << endl;

  // this is the residual E dependence after BrEta correction
  for (Int_t iE = 0; iE< nBinsE; ++iE){ 
    outfile2 << "  h_CBE_EB->SetBinContent(" << 2*iE+1 << ", " << h_CBE_EB->GetBinContent(2*iE+1) << "); " << endl;
  }

  outfile2 << "  Int_t iE      = -1;                                                                   " << endl;              

  outfile2 << "  for (Int_t iiE  = 0; iiE  < nBinsE;  ++iiE ){ 				      " << endl;
  outfile2 << "    if ( leftE [iiE]  <= (E)      &&       (E) < rightE[iiE] ) {		      " << endl;
  outfile2 << "      iE  = iiE;  								      " << endl;
  outfile2 << "    }										      " << endl;
  outfile2 << "  }										      " << endl;

  // Extra protections										       
  outfile2 << "  if (E < leftE  [0] )         { iE = 0;              if (DBG) std::cout << \" WARNING [applyScCorrections]: E < leftE  [0] )       \" << std::endl;}   " << endl;
  outfile2 << "  if (E > rightE [nBinsE-1] ) { iE = nBinsE-1;      if (DBG) std::cout << \" WARNING [applyScCorrections]: E > rightE [nBinsE-1] \" << std::endl;}   " << endl;

  outfile2 << "  if (iE == -1) {delete  h_CBE_EB; return 1;}					      " << endl;

  outfile2 << "											      " << endl;
  outfile2 << "  Int_t binE =  2*iE+1 ; // h_CBE_EB->FindBin(E);                                    " << endl;       
  outfile2 << "  Double_t tmp = h_CBE_EB->GetBinContent(binE);                                        " << endl;
  outfile2 << "  delete h_CBE_EB;                                                             	      " << endl;
  outfile2 << "  if ( 0< binE && binE < 2*nBinsE+1) return tmp;                          	      " << endl;
  outfile2 << "  else return 1.;                                                               	      " << endl;
  outfile2 << "                                                                                         " << endl;
  outfile2 << "}                                                                               " << endl;                                                            
  outfile2 << "                                                                                  " << endl; 
    
  outfile2 << "                                                                         " << endl;
  outfile2 << "Double_t applyScCorrectionsE_EE(Double_t E){  							   " << endl;
  outfile2 << "  const Int_t    nBinsE             = " << nBinsE<< ";             " << endl; 
  outfile2 << "  Double_t       leftE  [nBinsE];                    " << endl;
  outfile2 << "  Double_t       rightE [nBinsE];                    " << endl;
  outfile2 << "  Double_t       EBins  [nBinsE*2];                  " << endl;
  // E
  for (Int_t i = 0; i< nBinsE; ++i){ 
    outfile2 << "  leftE[" << i << "] =  " << leftE[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< nBinsE; ++i){ 
    outfile2 << "  rightE[" << i << "] =  " << rightE[i] << " ; " << endl;
  }
  for (Int_t i = 0; i< 2*nBinsE; ++i){ 
    outfile2 << "  EBins[" << i << "] =  " << EBins[i] << " ; " << endl;
  }
  outfile2 << "  TH1F *h_CBE_EE    = new TH1F(\"h_CBE_EE\"    ,\"h_CBE_EE\"    ,nBinsE*2-1, EBins); " << endl;

  // this is the residual E dependence after BrEta correction
  for (Int_t iE = 0; iE< nBinsE; ++iE){ 
    outfile2 << "  h_CBE_EE->SetBinContent(" << 2*iE+1 << ", " << h_CBE_EE->GetBinContent(2*iE+1) << "); " << endl;
  }

  outfile2 << "  Int_t iE      = -1;                                                                   " << endl;              
  outfile2 << "  for (Int_t iiE  = 0; iiE  < nBinsE;  ++iiE ){ 				      " << endl;
  outfile2 << "    if ( leftE [iiE]  <= (E)      &&       (E) < rightE[iiE] ) {		      " << endl;
  outfile2 << "      iE  = iiE;  								      " << endl;
  outfile2 << "    }										      " << endl;
  outfile2 << "  }										      " << endl;
  
  // Extra protections										       
  outfile2 << "  if (E < leftE  [0] )         { iE = 0;              if (DBG) std::cout << \" WARNING [applyScCorrections]: E < leftE  [0] )       \" << std::endl;}   " << endl;
  outfile2 << "  if (E > rightE [nBinsE-1] ) { iE = nBinsE-1;      if (DBG) std::cout << \" WARNING [applyScCorrections]: E > rightE [nBinsE-1] \" << std::endl;}   " << endl;

  outfile2 << "  if (iE == -1) {delete  h_CBE_EE; return 1;}					      " << endl;

  outfile2 << "											      " << endl;
  outfile2 << "  Int_t binE =  2*iE+1 ; // h_CBE_EE->FindBin(E);                                    " << endl;       
  outfile2 << "  Double_t tmp = h_CBE_EE->GetBinContent(binE);                                        " << endl;
  outfile2 << "  delete h_CBE_EE;                                                             	      " << endl;
  outfile2 << "  if ( 0< binE && binE < 2*nBinsE+1) return tmp;                          	      " << endl;
  outfile2 << "  else return 1.;                                                               	      " << endl;
  outfile2 << "                                                                                         " << endl;
  outfile2 << "}                                                                               " << endl;                                                            
  outfile2 << "                                                                                  " << endl; 

  outfile2 << "#endif                                                 " << endl;
  
  outfile2.close();                                                                                                      
  
  // save also the histograms to be used as correction
  TFile *histFile= new TFile("./histCorrections_E.root", "RECREATE" );
  h_CBE_EB->Write();  
  h_CBE_EE->Write();  
  histFile->Close();    
  return;
}

void scCorrections::getErecEgenFits_BrEta_E(correctionType applyCorrections, bool apply)
{
  Int_t minEntries      = 100;
  Int_t minChi2         = 3;

  // Fit the CB for each histogram of the ErecEgen matrix
  //--------------------------------------------------------------------------------
  for (Int_t iEta = 0; iEta < nBinsEta; ++iEta){
    for (Int_t iBr = 0; iBr < nBinsBr; ++iBr){
      Double_t mean    = 1;
      Double_t emean   = 0;
      Double_t chi2    = 0;      
      Double_t fracBad = 0;
      if (h_ErecEgen[iEta][iBr]->GetEntries() >minEntries) 
	{
	  fitCB(h_ErecEgen[iEta][iBr], mean, emean, chi2, fracBad);
	  if (!apply){	    	    
	    if (chi2 > minChi2){ // if the fit is bad
	      if ( 2*iEta+1 -2 > 1){ // set it to the previous eta bin if possible
		mean    = hh_CBetabr      ->GetBinContent(2*iBr+1  , 2*iEta+1  -2);
		emean   = hh_CBetabr      ->GetBinError  (2*iBr+1  , 2*iEta+1  -2);
		cout << "WARNING: BAD FIT set mean and error to the values of the previous eta bin " << mean << " " << emean << endl;
		//getchar();	      
	      }
	      else if (2*iBr+1 -2 > 1 ){ // otherwise set it to the previous Br bin if possible
		mean    = hh_CBetabr      ->GetBinContent(2*iBr+1 -2  , 2*iEta+1 );
		emean   = hh_CBetabr      ->GetBinError  (2*iBr+1 -2  , 2*iEta+1 );	    
		cout << "WARNING: BAD FIT set mean and error to the values of the previous br bin " << mean << " " << emean << endl;
		//getchar();
	      }
	      else{  // otherwise set it to 1. Don't know how to treat them
		mean  = 1;
		emean = 0;
		cout << "WARNING: BAD FIT set mean and error to " << mean << " " << emean << endl;
		//getchar();
	      }
	    }
	  }
	}
      else { // if not enough statistics
	if (!apply){	    	  
	  if  ( 2*iEta+1 -2 > 1){ // set it to the previous eta bin if possible
	    mean    = hh_CBetabr      ->GetBinContent(2*iBr+1   , 2*iEta+1 -2 );
	    emean   = hh_CBetabr      ->GetBinError  (2*iBr+1   , 2*iEta+1 -2 );
	    cout << "WARNING: NOT ENOUGH STAT set mean and error to the values of the previous eta bin " << mean << " " << emean << endl;
	    //getchar();
	  }
	  else if (2*iBr+1 -2 > 1 ){ // otherwise set it to the previous Br bin if possible
	    mean    = hh_CBetabr      ->GetBinContent(2*iBr+1 -2  , 2*iEta+1 );
	    emean   = hh_CBetabr      ->GetBinError  (2*iBr+1 -2  , 2*iEta+1 );
	    cout << "WARNING: NOT ENOUGH STAT set mean and error to the values of the previous br bin " << mean << " " << emean << endl;
	    //getchar();
	  }
	  else {  // otherwise set it to 1. Don't know how to treat them
	    mean  = 1;
	    emean = 0;
	    //getchar();
	    cout << "WARNING: NOT ENOUGH STAT set mean and error to " << mean << " " << emean << endl;
	  }
	}
      }
      hh_CBetabr      ->SetBinContent(2*iBr+1,2*iEta+1,mean  );
      hh_CBetabr      ->SetBinError  (2*iBr+1,2*iEta+1,emean );
      hh_Chi2etabr    ->SetBinContent(2*iBr+1,2*iEta+1,chi2  ); // dof = 4 = number dof floating parametes
      cout << " FRACTION BELOW " << EgenFrac << " = " << fracBad << endl;
    }
  }
  // 2D plot CB most probable value in bins (eta,br)
  TCanvas *cc_CBetabr = new TCanvas("cc_CBetabr","cc_CBetabr",600,600);
  cc_CBetabr->cd();
  hh_CBetabr->GetZaxis()->SetRangeUser(0.89,1.01);
  hh_CBetabr->GetXaxis()->SetTitle("#sigma_{phi}/#sigma_{eta}");
  hh_CBetabr->GetXaxis()->SetTitleOffset(1.5);
  hh_CBetabr->GetYaxis()->SetTitle("#eta");
  hh_CBetabr->GetYaxis()->SetTitleOffset(1.5);
  hh_CBetabr->GetZaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
  hh_CBetabr->GetZaxis()->SetTitleOffset(1.2);
  hh_CBetabr->SetTitle(0);
  gPad->SetPhi(180+30);
  hh_CBetabr->Draw("lego2");  
  if      (applyCorrections == brEta  )  cc_CBetabr->Print("./plots_tmp/cc_CBetabr_corrBrEta.png");
  else if (applyCorrections == brEtaE)  cc_CBetabr->Print("./plots_tmp/cc_CBetabr_corrBrEtaE.png");
  else                                   cc_CBetabr->Print("./plots_tmp/cc_CBetabr.png");
  
  TCanvas *cc_Chi2etabr = new TCanvas("cc_Chi2etabr","cc_Chi2etabr",600,600);
  cc_Chi2etabr->cd();
  gPad->SetPhi(180+30);
  hh_Chi2etabr->GetZaxis()->SetRangeUser(0.,5);
  hh_Chi2etabr->Draw("lego2");
  
  // 1D projection over br
  //
  // transfer each eta bins slice in a TH1F
  for (Int_t i = 0; i< nBinsEta; ++i){
    for (Int_t j = 0; j< nBinsBr; ++j){
      h_CBetabr[i]->SetBinContent(2*j+1,hh_CBetabr->GetBinContent(2*j+1,2*i+1));
      h_CBetabr[i]->SetBinError  (2*j+1,hh_CBetabr->GetBinError  (2*j+1,2*i+1));
    }
  }  
  // overlap the eta slices in a plot
  TCanvas *c_CBetabr = new TCanvas("c_CBetabr","c_CBetabr",600,600);
  c_CBetabr->cd();
  c_CBetabr->SetLeftMargin(0.2);
  for (Int_t i = 0; i<nBinsEta ; ++i){
    h_CBetabr[i]-> SetMarkerStyle(20);
    h_CBetabr[i]-> SetMarkerColor(i+22);
    h_CBetabr[i]-> SetTitle(0);
    h_CBetabr[i]-> GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");
    h_CBetabr[i]-> GetYaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
    h_CBetabr[i]-> GetYaxis()->SetTitleOffset(2);
    if(i == 0) {
      h_CBetabr[i]-> Draw();
      h_CBetabr[i]-> GetYaxis()->SetRangeUser(0.9,1.05);
    }
    else h_CBetabr[i]-> Draw("same");
  }
  TBox *bandBr = new TBox(leftBr[0],0.999,rightBr[nBinsBr-1],1.001);
  bandBr->SetFillColor(1);
  bandBr->SetFillStyle(3004);
  bandBr->Draw();
  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->SetFillColor(0);
  TString txt[nBinsEta];
  for (Int_t i = 0; i<nBinsEta ; ++i){
    txt[i] = floatToString(leftEta[i])+" < #eta < "+ floatToString(rightEta[i]);
    leg->AddEntry(h_CBetabr[i],txt[i]);
  }
  leg->Draw();
  if      (applyCorrections == brEta  ) c_CBetabr->Print("./plots_tmp/c_CBetabr_corrBrEta.png");
  else if (applyCorrections == brEtaE) c_CBetabr->Print("./plots_tmp/c_CBetabr_corrBrEtaE.png");
  else                                  c_CBetabr->Print("./plots_tmp/c_CBetabr.png");


  // Fit the CB for each histogram of the ErecEgen E-vector  *** BARREL ***
  //--------------------------------------------------------------------------------
  TCanvas *c_CBE_EB = new TCanvas("c_CBE_EB","c_CBE_EB",600,600);
  c_CBE_EB->cd();
  c_CBE_EB->SetLeftMargin(0.2);
  for (Int_t iE = 0; iE < nBinsE; ++iE){
    Double_t mean  = 1;
    Double_t emean = 0;
    Double_t chi2  = 0;
    Double_t fracBad = 0;
    if (h_EcorrBrEta_EB[iE]->GetEntries() > minEntries) fitCB(h_EcorrBrEta_EB[iE], mean, emean, chi2, fracBad);
    h_CBE_EB  ->SetBinContent(2*iE+1,mean  );
    h_CBE_EB  ->SetBinError  (2*iE+1,emean );
  }
  h_CBE_EB  ->SetMarkerStyle(20);
  h_CBE_EB  ->SetMarkerColor(1);
  h_CBE_EB  ->SetTitle(0);
  h_CBE_EB  ->GetYaxis()->SetRangeUser(0.99,1.01);
  h_CBE_EB  ->GetXaxis()->SetTitle("E [GeV]");
  h_CBE_EB  ->GetYaxis()->SetTitleOffset(2);
  h_CBE_EB  ->GetYaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
  h_CBE_EB  ->Draw();
  TBox *bandE_EB = new TBox(leftE[0],0.999,rightE[nBinsE-1],1.001);
  bandE_EB->SetFillColor(1);
  bandE_EB->SetFillStyle(3004);
  bandE_EB->Draw();
  if      (applyCorrections == brEta  ) c_CBE_EB->Print("./plots_tmp/c_CBE_EB_corrBrEta.png");
  else if (applyCorrections == brEtaE) c_CBE_EB->Print("./plots_tmp/c_CBE_EB_corrBrEtaE.png");
  else                                  c_CBE_EB->Print("./plots_tmp/c_CBE_EB.png");

  // Fit the CB for each histogram of the ErecEgen E-vector  *** ENDCAP ***
  //--------------------------------------------------------------------------------
  TCanvas *c_CBE_EE = new TCanvas("c_CBE_EE","c_CBE_EE",600,600);
  c_CBE_EE->cd();
  c_CBE_EE->SetLeftMargin(0.2);
  for (Int_t iE = 0; iE < nBinsE; ++iE){
    Double_t mean  = 1;
    Double_t emean = 0;
    Double_t chi2  = 0;
    Double_t fracBad = 0;
    if (h_EcorrBrEta_EE[iE]->GetEntries() > minEntries) fitCB(h_EcorrBrEta_EE[iE], mean, emean, chi2, fracBad);
    h_CBE_EE  ->SetBinContent(2*iE+1,mean  );
    h_CBE_EE  ->SetBinError  (2*iE+1,emean );
    h_Chi2E_EE->SetBinContent(2*iE+1,chi2  ); // dof = 4 = number fof floating parametes
  }
  h_CBE_EE  ->SetMarkerStyle(20);
  h_CBE_EE  ->SetMarkerColor(1);
  h_CBE_EE  ->SetTitle(0);
  h_CBE_EE  ->GetYaxis()->SetRangeUser(0.99,1.01);
  h_CBE_EE  ->GetXaxis()->SetTitle("E [GeV]");
  h_CBE_EE  ->GetYaxis()->SetTitleOffset(2);
  h_CBE_EE  ->GetYaxis()->SetTitle("CrystalBall     mp(E_{rec}/E_{gen})");
  h_CBE_EE  ->Draw();
  TBox *bandE_EE = new TBox(leftE[0],0.999,rightE[nBinsE-1],1.001);
  bandE_EE->SetFillColor(1);
  bandE_EE->SetFillStyle(3004);
  bandE_EE->Draw();
  if      (applyCorrections == brEta  ) c_CBE_EE->Print("./plots_tmp/c_CBE_EE_corrBrEta.png");
  else if (applyCorrections == brEtaE) c_CBE_EE->Print("./plots_tmp/c_CBE_EE_corrBrEtaE.png");
  else                                  c_CBE_EE->Print("./plots_tmp/c_CBE_EE.png");

  return;
}




