//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 30 14:21:13 2011 by ROOT version 5.30/02
// from TTree Tree/Tree
// found on file: mc_inclusive.root
//////////////////////////////////////////////////////////

#ifndef template_production_h
#define template_production_h

#include "binsdef.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "TRandom3.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TF1.h"
#include "RooBinning.h"
#include "TString.h"
#include "TH1.h"
#include "TProfile.h"
#include "TMath.h"
#include "TRandom3.h"
#include <map>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include <vector>
#include <algorithm> 

using namespace std;
using namespace RooFit;

class template_production {
public :

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

//   vector<float> invmass_vector;
//   vector<float> diphotonpt_vector;

   // Declaration of leaf types
   Float_t         event_luminormfactor;
   Float_t         event_Kfactor;
   Float_t         event_weight;
   Float_t         event_rho;
   Float_t         event_sigma;
   Int_t           event_nPU;
   Int_t           event_PUOOTnumInteractionsEarly;
   Int_t           event_PUOOTnumInteractionsLate;
   Int_t           event_nRecVtx;
   Int_t           event_pass12whoissiglike;
   Int_t           event_CSCTightHaloID;
   Int_t           event_NMuons;
   Int_t           event_NMuonsTot;
   Float_t         dipho_mgg_photon;
   Float_t         dipho_mgg_newCorr;
   Float_t         dipho_mgg_newCorrLocal;
   Float_t         pholead_eta;
   Float_t         photrail_eta;
   Float_t         pholead_phi;
   Float_t         photrail_phi;
   Float_t         pholead_px;
   Float_t         photrail_px;
   Float_t         pholead_py;
   Float_t         photrail_py;
   Float_t         pholead_pt;
   Float_t         photrail_pt;
   Float_t         pholead_pz;
   Float_t         photrail_pz;
   Float_t         pholead_energy;
   Float_t         photrail_energy;
   Float_t         pholead_energySCdefault;
   Float_t         photrail_energySCdefault;
   Float_t         pholead_energyNewCorr;
   Float_t         photrail_energyNewCorr;
   Float_t         pholead_energyNewCorrLocal;
   Float_t         photrail_energyNewCorrLocal;
   Float_t         pholead_SCeta;
   Float_t         photrail_SCeta;
   Float_t         pholead_SCphi;
   Float_t         photrail_SCphi;
   Int_t           pholead_PhoHasPixSeed;
   Int_t           pholead_PhoHasConvTrks;
   Int_t           pholead_PhoScSeedSeverity;
   Int_t           photrail_PhoHasPixSeed;
   Int_t           photrail_PhoHasConvTrks;
   Int_t           photrail_PhoScSeedSeverity;
   Float_t         pholead_r9;
   Float_t         photrail_r9;
   Float_t         pholead_sieie;
   Float_t         photrail_sieie;
   Float_t         pholead_hoe;
   Float_t         photrail_hoe;
   Float_t         pholead_brem;
   Float_t         photrail_brem;
   Float_t         pholead_sigmaPhi;
   Float_t         photrail_sigmaPhi;
   Float_t         pholead_sigmaEta;
   Float_t         photrail_sigmaEta;
   Float_t         pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx;
   Float_t         photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx;
   Float_t         pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx;
   Float_t         photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx;
   Float_t         pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx;
   Float_t         photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx;
   Float_t         pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;
   Float_t         photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;
   Float_t         pholead_pho_Cone01NeutralHadronIso_mvVtx;
   Float_t         photrail_pho_Cone01NeutralHadronIso_mvVtx;
   Float_t         pholead_pho_Cone02NeutralHadronIso_mvVtx;
   Float_t         photrail_pho_Cone02NeutralHadronIso_mvVtx;
   Float_t         pholead_pho_Cone03NeutralHadronIso_mvVtx;
   Float_t         photrail_pho_Cone03NeutralHadronIso_mvVtx;
   Float_t         pholead_pho_Cone04NeutralHadronIso_mvVtx;
   Float_t         photrail_pho_Cone04NeutralHadronIso_mvVtx;
   Float_t         pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01;
   Float_t         photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01;
   Float_t         pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01;
   Float_t         photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01;
   Float_t         pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01;
   Float_t         photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01;
   Float_t         pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
   Float_t         photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
   Float_t         pholead_pho_Cone03PFCombinedIso;
   Float_t         photrail_pho_Cone03PFCombinedIso;
   Float_t         pholead_pho_Cone04PFCombinedIso;
   Float_t         photrail_pho_Cone04PFCombinedIso;
   Int_t           pholead_PhoPassConvSafeElectronVeto;
   Int_t           photrail_PhoPassConvSafeElectronVeto;
   Float_t         pholead_GenPhotonIsoDR04;
   Float_t         photrail_GenPhotonIsoDR04;
   Float_t         pholead_PhoIso03Ecal;
   Float_t         pholead_PhoIso03Hcal;
   Float_t         pholead_PhoIso03TrkSolid;
   Float_t         pholead_PhoIso03TrkHollow;
   Float_t         pholead_PhoIso03;
   Float_t         pholead_PhoIso04Ecal;
   Float_t         pholead_PhoIso04Hcal;
   Float_t         pholead_PhoIso04TrkSolid;
   Float_t         pholead_PhoIso04TrkHollow;
   Float_t         pholead_PhoIso04;
   Float_t         photrail_PhoIso03Ecal;
   Float_t         photrail_PhoIso03Hcal;
   Float_t         photrail_PhoIso03TrkSolid;
   Float_t         photrail_PhoIso03TrkHollow;
   Float_t         photrail_PhoIso03;
   Float_t         photrail_PhoIso04Ecal;
   Float_t         photrail_PhoIso04Hcal;
   Float_t         photrail_PhoIso04TrkSolid;
   Float_t         photrail_PhoIso04TrkHollow;
   Float_t         photrail_PhoIso04;
   Float_t         pholead_PhoS4OverS1;
   Float_t         pholead_PhoSigmaEtaEta;
   Float_t         pholead_PhoE1x5;
   Float_t         pholead_PhoE2x5;
   Float_t         pholead_PhoE3x3;
   Float_t         pholead_PhoE5x5;
   Float_t         pholead_PhomaxEnergyXtal;
   Float_t         pholead_PhoIso03HcalDepth1;
   Float_t         pholead_PhoIso03HcalDepth2;
   Float_t         pholead_PhoIso04HcalDepth1;
   Float_t         pholead_PhoIso04HcalDepth2;
   Int_t           pholead_PhoIso03nTrksSolid;
   Int_t           pholead_PhoIso03nTrksHollow;
   Int_t           pholead_PhoIso04nTrksSolid;
   Int_t           pholead_PhoIso04nTrksHollow;
   Float_t         pholead_Pho_ChargedHadronIso;
   Float_t         pholead_Pho_NeutralHadronIso;
   Float_t         pholead_Pho_PhotonIso;
   Int_t           pholead_Pho_isPFPhoton;
   Int_t           pholead_Pho_isPFElectron;
   Float_t         photrail_PhoS4OverS1;
   Float_t         photrail_PhoSigmaEtaEta;
   Float_t         photrail_PhoE1x5;
   Float_t         photrail_PhoE2x5;
   Float_t         photrail_PhoE3x3;
   Float_t         photrail_PhoE5x5;
   Float_t         photrail_PhomaxEnergyXtal;
   Float_t         photrail_PhoIso03HcalDepth1;
   Float_t         photrail_PhoIso03HcalDepth2;
   Float_t         photrail_PhoIso04HcalDepth1;
   Float_t         photrail_PhoIso04HcalDepth2;
   Int_t           photrail_PhoIso03nTrksSolid;
   Int_t           photrail_PhoIso03nTrksHollow;
   Int_t           photrail_PhoIso04nTrksSolid;
   Int_t           photrail_PhoIso04nTrksHollow;
   Float_t         photrail_Pho_ChargedHadronIso;
   Float_t         photrail_Pho_NeutralHadronIso;
   Float_t         photrail_Pho_PhotonIso;
   Int_t           photrail_Pho_isPFPhoton;
   Int_t           photrail_Pho_isPFElectron;
   Int_t           pholead_PhoMCmatchindex;
   Int_t           pholead_PhoMCmatchexitcode;
   Int_t           photrail_PhoMCmatchindex;
   Int_t           photrail_PhoMCmatchexitcode;
   Int_t           pholead_Npfcandphotonincone;
   Int_t           pholead_Npfcandchargedincone;
   Int_t           pholead_Npfcandneutralincone;
   Int_t           photrail_Npfcandphotonincone;
   Int_t           photrail_Npfcandchargedincone;
   Int_t           photrail_Npfcandneutralincone;
   Float_t         pholead_scareaSF;
   Float_t         photrail_scareaSF;
   Float_t         pholead_photonpfcandenergies[30];
   Float_t         pholead_photonpfcandets[30];
   Float_t         pholead_photonpfcanddetas[30];
   Float_t         pholead_photonpfcanddphis[30];
   Float_t         photrail_photonpfcandenergies[30];
   Float_t         photrail_photonpfcandets[30];
   Float_t         photrail_photonpfcanddetas[30];
   Float_t         photrail_photonpfcanddphis[30];
   
   // List of branches
   TBranch        *b_event_luminormfactor;   //!
   TBranch        *b_event_Kfactor;   //!
   TBranch        *b_event_weight;   //!
   TBranch        *b_event_rho;   //!
   TBranch        *b_event_sigma;   //!
   TBranch        *b_event_nPU;   //!
   TBranch        *b_event_PUOOTnumInteractionsEarly;   //!
   TBranch        *b_event_PUOOTnumInteractionsLate;   //!
   TBranch        *b_event_nRecVtx;   //!
   TBranch        *b_event_pass12whoissiglike;   //!
   TBranch        *b_event_CSCTightHaloID; //!
   TBranch        *b_event_NMuons; //!
   TBranch        *b_event_NMuonsTot; //!
   TBranch        *b_dipho_mgg_photon;   //!
   TBranch        *b_dipho_mgg_newCorr;   //!
   TBranch        *b_dipho_mgg_newCorrLocal;   //!
   TBranch        *b_pholead_eta;   //!
   TBranch        *b_photrail_eta;   //!
   TBranch        *b_pholead_phi;   //!
   TBranch        *b_photrail_phi;   //!
   TBranch        *b_pholead_px;   //!
   TBranch        *b_photrail_px;   //!
   TBranch        *b_pholead_py;   //!
   TBranch        *b_photrail_py;   //!
   TBranch        *b_pholead_pt;   //!
   TBranch        *b_photrail_pt;   //!
   TBranch        *b_pholead_pz;   //!
   TBranch        *b_photrail_pz;   //!
   TBranch        *b_pholead_energy;   //!
   TBranch        *b_photrail_energy;   //!
   TBranch        *b_pholead_energySCdefault;   //!
   TBranch        *b_photrail_energySCdefault;   //!
   TBranch        *b_pholead_energyNewCorr;   //!
   TBranch        *b_photrail_energyNewCorr;   //!
   TBranch        *b_pholead_energyNewCorrLocal;   //!
   TBranch        *b_photrail_energyNewCorrLocal;   //!
   TBranch        *b_pholead_SCeta;   //!
   TBranch        *b_photrail_SCeta;   //!
   TBranch        *b_pholead_SCphi;   //!
   TBranch        *b_photrail_SCphi;   //!
   TBranch        *b_pholead_PhoHasPixSeed;   //!
   TBranch        *b_pholead_PhoHasConvTrks;   //!
   TBranch        *b_pholead_PhoScSeedSeverity;   //!
   TBranch        *b_photrail_PhoHasPixSeed;   //!
   TBranch        *b_photrail_PhoHasConvTrks;   //!
   TBranch        *b_photrail_PhoScSeedSeverity;   //!
   TBranch        *b_pholead_r9;   //!
   TBranch        *b_photrail_r9;   //!
   TBranch        *b_pholead_sieie;   //!
   TBranch        *b_photrail_sieie;   //!
   TBranch        *b_pholead_hoe;   //!
   TBranch        *b_photrail_hoe;   //!
   TBranch        *b_pholead_brem;   //!
   TBranch        *b_photrail_brem;   //!
   TBranch        *b_pholead_sigmaPhi;   //!
   TBranch        *b_photrail_sigmaPhi;   //!
   TBranch        *b_pholead_sigmaEta;   //!
   TBranch        *b_photrail_sigmaEta;   //!
   TBranch        *b_pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx;   //!
   TBranch        *b_photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx;   //!
   TBranch        *b_pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx;   //!
   TBranch        *b_photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx;   //!
   TBranch        *b_pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx;   //!
   TBranch        *b_photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx;   //!
   TBranch        *b_pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;   //!
   TBranch        *b_photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;   //!
   TBranch        *b_pholead_pho_Cone01NeutralHadronIso_mvVtx;   //!
   TBranch        *b_photrail_pho_Cone01NeutralHadronIso_mvVtx;   //!
   TBranch        *b_pholead_pho_Cone02NeutralHadronIso_mvVtx;   //!
   TBranch        *b_photrail_pho_Cone02NeutralHadronIso_mvVtx;   //!
   TBranch        *b_pholead_pho_Cone03NeutralHadronIso_mvVtx;   //!
   TBranch        *b_photrail_pho_Cone03NeutralHadronIso_mvVtx;   //!
   TBranch        *b_pholead_pho_Cone04NeutralHadronIso_mvVtx;   //!
   TBranch        *b_photrail_pho_Cone04NeutralHadronIso_mvVtx;   //!
   TBranch        *b_pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01;   //!
   TBranch        *b_photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01;   //!
   TBranch        *b_pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01;   //!
   TBranch        *b_photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01;   //!
   TBranch        *b_pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01;   //!
   TBranch        *b_photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01;   //!
   TBranch        *b_pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;   //!
   TBranch        *b_photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;   //!
   TBranch        *b_pholead_pho_Cone03PFCombinedIso;   //!
   TBranch        *b_photrail_pho_Cone03PFCombinedIso;   //!
   TBranch        *b_pholead_pho_Cone04PFCombinedIso;   //!
   TBranch        *b_photrail_pho_Cone04PFCombinedIso;   //!
   TBranch        *b_pholead_PhoPassConvSafeElectronVeto;   //!
   TBranch        *b_photrail_PhoPassConvSafeElectronVeto;   //!
   TBranch        *b_pholead_GenPhotonIsoDR04;   //!
   TBranch        *b_photrail_GenPhotonIsoDR04;   //!
   TBranch        *b_pholead_PhoIso03Ecal;   //!
   TBranch        *b_pholead_PhoIso03Hcal;   //!
   TBranch        *b_pholead_PhoIso03TrkSolid;   //!
   TBranch        *b_pholead_PhoIso03TrkHollow;   //!
   TBranch        *b_pholead_PhoIso03;   //!
   TBranch        *b_pholead_PhoIso04Ecal;   //!
   TBranch        *b_pholead_PhoIso04Hcal;   //!
   TBranch        *b_pholead_PhoIso04TrkSolid;   //!
   TBranch        *b_pholead_PhoIso04TrkHollow;   //!
   TBranch        *b_pholead_PhoIso04;   //!
   TBranch        *b_photrail_PhoIso03Ecal;   //!
   TBranch        *b_photrail_PhoIso03Hcal;   //!
   TBranch        *b_photrail_PhoIso03TrkSolid;   //!
   TBranch        *b_photrail_PhoIso03TrkHollow;   //!
   TBranch        *b_photrail_PhoIso03;   //!
   TBranch        *b_photrail_PhoIso04Ecal;   //!
   TBranch        *b_photrail_PhoIso04Hcal;   //!
   TBranch        *b_photrail_PhoIso04TrkSolid;   //!
   TBranch        *b_photrail_PhoIso04TrkHollow;   //!
   TBranch        *b_photrail_PhoIso04;   //!
   TBranch        *b_pholead_PhoS4OverS1;   //!
   TBranch        *b_pholead_PhoSigmaEtaEta;   //!
   TBranch        *b_pholead_PhoE1x5;   //!
   TBranch        *b_pholead_PhoE2x5;   //!
   TBranch        *b_pholead_PhoE3x3;   //!
   TBranch        *b_pholead_PhoE5x5;   //!
   TBranch        *b_pholead_PhomaxEnergyXtal;   //!
   TBranch        *b_pholead_PhoIso03HcalDepth1;   //!
   TBranch        *b_pholead_PhoIso03HcalDepth2;   //!
   TBranch        *b_pholead_PhoIso04HcalDepth1;   //!
   TBranch        *b_pholead_PhoIso04HcalDepth2;   //!
   TBranch        *b_pholead_PhoIso03nTrksSolid;   //!
   TBranch        *b_pholead_PhoIso03nTrksHollow;   //!
   TBranch        *b_pholead_PhoIso04nTrksSolid;   //!
   TBranch        *b_pholead_PhoIso04nTrksHollow;   //!
   TBranch        *b_pholead_Pho_ChargedHadronIso;   //!
   TBranch        *b_pholead_Pho_NeutralHadronIso;   //!
   TBranch        *b_pholead_Pho_PhotonIso;   //!
   TBranch        *b_pholead_Pho_isPFPhoton;   //!
   TBranch        *b_pholead_Pho_isPFElectron;   //!
   TBranch        *b_photrail_PhoS4OverS1;   //!
   TBranch        *b_photrail_PhoSigmaEtaEta;   //!
   TBranch        *b_photrail_PhoE1x5;   //!
   TBranch        *b_photrail_PhoE2x5;   //!
   TBranch        *b_photrail_PhoE3x3;   //!
   TBranch        *b_photrail_PhoE5x5;   //!
   TBranch        *b_photrail_PhomaxEnergyXtal;   //!
   TBranch        *b_photrail_PhoIso03HcalDepth1;   //!
   TBranch        *b_photrail_PhoIso03HcalDepth2;   //!
   TBranch        *b_photrail_PhoIso04HcalDepth1;   //!
   TBranch        *b_photrail_PhoIso04HcalDepth2;   //!
   TBranch        *b_photrail_PhoIso03nTrksSolid;   //!
   TBranch        *b_photrail_PhoIso03nTrksHollow;   //!
   TBranch        *b_photrail_PhoIso04nTrksSolid;   //!
   TBranch        *b_photrail_PhoIso04nTrksHollow;   //!
   TBranch        *b_photrail_Pho_ChargedHadronIso;   //!
   TBranch        *b_photrail_Pho_NeutralHadronIso;   //!
   TBranch        *b_photrail_Pho_PhotonIso;   //!
   TBranch        *b_photrail_Pho_isPFPhoton;   //!
   TBranch        *b_photrail_Pho_isPFElectron;   //!
   TBranch        *b_pholead_PhoMCmatchindex;   //!
   TBranch        *b_pholead_PhoMCmatchexitcode;   //!
   TBranch        *b_photrail_PhoMCmatchindex;   //!
   TBranch        *b_photrail_PhoMCmatchexitcode;   //!
   TBranch        *b_pholead_Npfcandphotonincone;
   TBranch        *b_pholead_Npfcandchargedincone;
   TBranch        *b_pholead_Npfcandneutralincone;
   TBranch        *b_photrail_Npfcandphotonincone;
   TBranch        *b_photrail_Npfcandchargedincone;
   TBranch        *b_photrail_Npfcandneutralincone;
   TBranch        *b_pholead_scareaSF;
   TBranch        *b_photrail_scareaSF;
   TBranch        *b_pholead_photonpfcandets;
   TBranch        *b_pholead_photonpfcandenergies;
   TBranch        *b_pholead_photonpfcanddetas;
   TBranch        *b_pholead_photonpfcanddphis;
   TBranch        *b_photrail_photonpfcandets;
   TBranch        *b_photrail_photonpfcandenergies;
   TBranch        *b_photrail_photonpfcanddetas;
   TBranch        *b_photrail_photonpfcanddphis;


   template_production(TTree *tree=0);
   virtual ~template_production();
   /* virtual Int_t    Cut(Long64_t entry); */
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   virtual void     Loop(int maxevents = -1);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void     WriteOutput(const char* filename, const TString dirname);

   void Setup(Bool_t _isdata, TString _mode, TString _differentialvariable);

   TString differentialvariable;

   std::vector<std::vector<TProfile*> > GetPUScaling(bool doEB, TString diffvar);

   TRandom3 *randomgen;

   float AbsDeltaPhi(double phi1, double phi2);

   static const int n_templates=n_bins;

   bool dosignal;

  RooRealVar *roovar1;
  RooRealVar *roovar2;
  RooRealVar *roopt1;
  RooRealVar *roopt2;
  RooRealVar *roosieie1;
  RooRealVar *roosieie2;
  RooRealVar *rooeta1;
  RooRealVar *rooeta2;
  RooRealVar *roorho;
  RooRealVar *roosigma;
  RooRealVar *rooweight;

  RooRealVar *roovar_invmass;
  RooRealVar *roovar_diphotonpt;
  RooRealVar *roovar_costhetastar;
  RooRealVar *roovar_dphi;
  RooRealVar *roovar_dR;

  RooArgSet *rooargset_diffvariables;

  RooDataSet *roodset_signal[2][2];
  RooDataSet *roodset_background[2][2];


   std::map<TString, RooDataSet*> obs_roodset;
   std::map<TString, RooDataSet*> template2d_roodset;

   std::map<TString, TProfile*> true_purity;
   std::map<TString, TH1F*> true_purity_isppevent;
   std::map<TString, TH1F*> true_purity_isnotppevent;

   TString get_name_obs_roodset(int region, TString diffvariable, int bin);
   TString get_name_true_purity(int region, TString diffvariable);
   TString get_name_true_purity_ispp(int region, TString diffvariable);
   TString get_name_true_purity_isnotpp(int region, TString diffvariable);
   TString get_name_template2d_roodset(int region, TString sigorbkg);

   Float_t pholead_outvar;
   Float_t photrail_outvar;
   Int_t candcounter;

   Bool_t initialized;

   Bool_t isdata;

   TString mode;

   void FillDiffVariables();

   Bool_t dosignaltemplate;
   Bool_t dobackgroundtemplate;
   Bool_t dodistribution;
   Bool_t do2dtemplate;
   Bool_t do2ptemplate;
   Bool_t do1p1ftemplate;
   Bool_t do2ftemplate;

   Int_t Choose_bin_invmass(float invmass, int region);
   Int_t Choose_bin_diphotonpt(float diphotonpt, int region);
   Int_t Choose_bin_costhetastar(float costhetastar, int region);
   Int_t Choose_bin_dphi(float dphi, int region);
   Int_t Choose_bin_dR(float dR, int region);

   Int_t Choose_bin_pt(float pt, int region);
   Int_t Choose_bin_eta(float eta, int region);
   Int_t Choose_bin_sieie(float sieie, int region);

   float getpuenergy(int reg, float eta);

};

#endif

#ifdef template_production_cxx
template_production::template_production(TTree *tree)
{

  TH1F::SetDefaultSumw2(kTRUE);

   if (tree==0) std::cout << "Tree not ready!" << std::endl;
   if (!tree) return;
   fChain = tree;

   initialized=false;
   dosignaltemplate=false;
   dobackgroundtemplate=false;
   dodistribution=false;
   do2dtemplate = false;
   do2ptemplate = false;
   do1p1ftemplate = false;
   do2ftemplate = false;

}

void template_production::Setup(Bool_t _isdata, TString _mode, TString _differentialvariable){

  isdata=_isdata;
  mode=_mode;
  differentialvariable=_differentialvariable;

  Init();

  if (mode=="standard" || mode=="preselection_diphoton" || mode=="standard_2frag" || mode=="standard_pixelrev") dodistribution=true;
  if (mode=="standard_2pgen" || mode=="standard_1p1fbothgen" || mode=="standard_2fgen") dodistribution=true;
  if (mode=="signal" || mode=="fragmentation" || mode=="nofragmentation" || mode=="signal_2frag" || mode=="randomcone" || mode=="cutPFchargediso_signal" || mode=="cutPFchargediso_randomcone") dosignaltemplate=true;
  if (mode=="background" || mode=="sieiesideband" || mode=="cutPFchargediso_background" || mode=="cutPFchargediso_sieiesideband") dobackgroundtemplate=true;
  if (mode=="sigsig" || mode=="2pgen" || mode=="zmumu" || mode=="zee" || mode=="2pgen_2frag") do2ptemplate=true; 
  if (mode=="sigbkg" || mode=="1p1fbothgen" || mode=="1prcone1fgen" || mode=="1pgen1fside" || mode=="1p1fbothgen_2frag" || mode=="1pgen1fside_2frag") do1p1ftemplate=true; 
  if (mode=="bkgbkg" || mode=="2fgen") do2ftemplate=true; 
  do2dtemplate = (do2ptemplate || do1p1ftemplate || do2ftemplate);


  randomgen = new TRandom3(0);

  roovar1 = new RooRealVar("roovar1","roovar1",leftrange,rightrange);
  roovar2 = new RooRealVar("roovar2","roovar2",leftrange,rightrange);
  rooeta1 = new RooRealVar("rooeta1","rooeta1",0,2.5);
  rooeta2 = new RooRealVar("rooeta2","rooeta2",0,2.5);
  roopt1 = new RooRealVar("roopt1","roopt1",25,1000);
  roopt2 = new RooRealVar("roopt2","roopt2",25,1000);
  roosieie1 = new RooRealVar("roosieie1","roosieie1",0,0.045);
  roosieie2 = new RooRealVar("roosieie2","roosieie2",0,0.045);
  roorho = new RooRealVar("roorho","roorho",0,50);
  roosigma = new RooRealVar("roosigma","roosigma",0,50);
  rooweight = new RooRealVar("rooweight","rooweight",0,5);

  roovar_invmass = new RooRealVar("roovar_invmass","roovar_invmass",binsdef_diphoton_invmass_EBEB[0],binsdef_diphoton_invmass_EBEB[n_templates_invmass_EBEB]);
  roovar_diphotonpt = new RooRealVar("roovar_diphotonpt","roovar_diphotonpt",binsdef_diphoton_diphotonpt_EBEB[0],binsdef_diphoton_diphotonpt_EBEB[n_templates_diphotonpt_EBEB]);
  roovar_costhetastar = new RooRealVar("roovar_costhetastar","roovar_costhetastar",binsdef_diphoton_costhetastar_EBEB[0],binsdef_diphoton_costhetastar_EBEB[n_templates_costhetastar_EBEB]);
  roovar_dphi = new RooRealVar("roovar_dphi","roovar_dphi",binsdef_diphoton_dphi_EBEB[0],binsdef_diphoton_dphi_EBEB[n_templates_dphi_EBEB]);
  roovar_dR = new RooRealVar("roovar_dR","roovar_dR",binsdef_diphoton_dR_EBEB[0],binsdef_diphoton_dR_EBEB[n_templates_dR_EBEB]);

  rooargset_diffvariables = new RooArgSet(*roovar_invmass,*roovar_diphotonpt,*roovar_costhetastar,*roovar_dphi,*roovar_dR);

  for (int i=0; i<2; i++){
      TString name_signal="signal";
      TString reg;
      if (i==0) reg="EB"; else if (i==1) reg="EE";
      TString t2=Form("roodset_%s_%s_rv%d",name_signal.Data(),reg.Data(),1);
      roodset_signal[i][0] = new RooDataSet(t2.Data(),t2.Data(),RooArgSet(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*rooweight),WeightVar(*rooweight));
      t2=Form("roodset_%s_%s_rv%d",name_signal.Data(),reg.Data(),2);
      roodset_signal[i][1] = new RooDataSet(t2.Data(),t2.Data(),RooArgSet(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma,*rooweight),WeightVar(*rooweight));
    }
  for (int i=0; i<2; i++){
      TString name_background="background";
      TString reg;
      if (i==0) reg="EB"; else if (i==1) reg="EE";
      TString t2=Form("roodset_%s_%s_rv%d",name_background.Data(),reg.Data(),1);
      roodset_background[i][0] = new RooDataSet(t2.Data(),t2.Data(),RooArgSet(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*rooweight),WeightVar(*rooweight));
      t2=Form("roodset_%s_%s_rv%d",name_background.Data(),reg.Data(),2);
      roodset_background[i][1] = new RooDataSet(t2.Data(),t2.Data(),RooArgSet(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma,*rooweight),WeightVar(*rooweight));
    }


  for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){

    for (int i=0; i<3; i++) {
      TString reg;
      if (i==0) reg="EBEB"; else if (i==1) reg="EBEE"; else if (i==2) reg="EEEE"; else if (i==3) reg="EEEB";
      for (int j=0; j<n_templates+1; j++) {
	TString t2=Form("obs_roodset_%s_%s_b%d",reg.Data(),diffvariable->Data(),j);
	RooArgSet args(*roovar1,*roovar2,*roopt1,*roosieie1,*rooeta1,*roopt2,*roosieie2,*rooeta2);
	args.add(RooArgSet(*roorho,*roosigma,*rooweight));
	args.add(*rooargset_diffvariables);
	obs_roodset[t2] = new RooDataSet(t2.Data(),t2.Data(),args,WeightVar(*rooweight));
      }

      int bins_to_run=-1; 
      float *binsdef=NULL;

      TString splitting = reg;
      if (splitting=="EEEB") splitting="EBEE";
      if (*diffvariable=="invmass"){
	if (splitting=="EBEB")      bins_to_run+=n_templates_invmass_EBEB;
	else if (splitting=="EBEE") bins_to_run+=n_templates_invmass_EBEE;
	else if (splitting=="EEEE") bins_to_run+=n_templates_invmass_EEEE; 
	if (splitting=="EBEB")      binsdef=binsdef_diphoton_invmass_EBEB;
	else if (splitting=="EBEE") binsdef=binsdef_diphoton_invmass_EBEE;
	else if (splitting=="EEEE") binsdef=binsdef_diphoton_invmass_EEEE;
      }
      if (*diffvariable=="diphotonpt"){
	if (splitting=="EBEB")      bins_to_run+=n_templates_diphotonpt_EBEB;
	else if (splitting=="EBEE") bins_to_run+=n_templates_diphotonpt_EBEE;
	else if (splitting=="EEEE") bins_to_run+=n_templates_diphotonpt_EEEE; 
	if (splitting=="EBEB")      binsdef=binsdef_diphoton_diphotonpt_EBEB;
	else if (splitting=="EBEE") binsdef=binsdef_diphoton_diphotonpt_EBEE;
	else if (splitting=="EEEE") binsdef=binsdef_diphoton_diphotonpt_EEEE;
      }
      if (*diffvariable=="costhetastar"){
	if (splitting=="EBEB")      bins_to_run+=n_templates_costhetastar_EBEB;
	else if (splitting=="EBEE") bins_to_run+=n_templates_costhetastar_EBEE;
	else if (splitting=="EEEE") bins_to_run+=n_templates_costhetastar_EEEE; 
	if (splitting=="EBEB")      binsdef=binsdef_diphoton_costhetastar_EBEB;
	else if (splitting=="EBEE") binsdef=binsdef_diphoton_costhetastar_EBEE;
	else if (splitting=="EEEE") binsdef=binsdef_diphoton_costhetastar_EEEE;
      }
      if (*diffvariable=="dphi"){
	if (splitting=="EBEB")      bins_to_run+=n_templates_dphi_EBEB;
	else if (splitting=="EBEE") bins_to_run+=n_templates_dphi_EBEE;
	else if (splitting=="EEEE") bins_to_run+=n_templates_dphi_EEEE; 
	if (splitting=="EBEB")      binsdef=binsdef_diphoton_dphi_EBEB;
	else if (splitting=="EBEE") binsdef=binsdef_diphoton_dphi_EBEE;
	else if (splitting=="EEEE") binsdef=binsdef_diphoton_dphi_EEEE;
      }
      if (*diffvariable=="dR"){
	if (splitting=="EBEB")      bins_to_run+=n_templates_dR_EBEB;
	else if (splitting=="EBEE") bins_to_run+=n_templates_dR_EBEE;
	else if (splitting=="EEEE") bins_to_run+=n_templates_dR_EEEE; 
	if (splitting=="EBEB")      binsdef=binsdef_diphoton_dR_EBEB;
	else if (splitting=="EBEE") binsdef=binsdef_diphoton_dR_EBEE;
	else if (splitting=="EEEE") binsdef=binsdef_diphoton_dR_EEEE;
      }


      TString t3=Form("true_purity_%s_%s",reg.Data(),diffvariable->Data());
      true_purity[t3] = new TProfile(t3.Data(),t3.Data(),bins_to_run,binsdef);
      TString t3_ispp=Form("ispp_true_purity_%s_%s",reg.Data(),diffvariable->Data());
      true_purity_isppevent[t3_ispp] = new TH1F(t3_ispp.Data(),t3_ispp.Data(),bins_to_run,binsdef);
      TString t3_isnotpp=Form("isnotpp_true_purity_%s_%s",reg.Data(),diffvariable->Data());
      true_purity_isnotppevent[t3_isnotpp] = new TH1F(t3_isnotpp.Data(),t3_isnotpp.Data(),bins_to_run,binsdef);

    }

  }

  std::vector<TString> tobuild;
  if (do2ptemplate) tobuild.push_back(TString("sigsig"));
  else if (do2ftemplate) tobuild.push_back(TString("bkgbkg"));
  else if (do1p1ftemplate) {tobuild.push_back(TString("sigbkg")); tobuild.push_back(TString("bkgsig"));}

  for (int i=0; i<3; i++)
    for (unsigned int j=0; j<tobuild.size(); j++) {
      TString t2 = get_name_template2d_roodset(i,tobuild[j]);
      RooArgSet args(*roovar1,*roovar2,*roopt1,*roosieie1,*rooeta1,*roopt2,*roosieie2,*rooeta2);
      args.add(RooArgSet(*roorho,*roosigma,*rooweight));
      args.add(*rooargset_diffvariables);
      template2d_roodset[t2]= new RooDataSet(t2.Data(),t2.Data(),args,WeightVar(*rooweight));
    }
  


  
  initialized=true;
  
};


template_production::~template_production()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   delete randomgen;


}

Int_t template_production::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t template_production::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void template_production::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers

   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_luminormfactor", &event_luminormfactor, &b_event_luminormfactor);
   fChain->SetBranchAddress("event_Kfactor", &event_Kfactor, &b_event_Kfactor);
   fChain->SetBranchAddress("event_weight", &event_weight, &b_event_weight);
   fChain->SetBranchAddress("event_rho", &event_rho, &b_event_rho);
   fChain->SetBranchAddress("event_sigma", &event_sigma, &b_event_sigma);
   fChain->SetBranchAddress("event_nPU", &event_nPU, &b_event_nPU);
   fChain->SetBranchAddress("event_PUOOTnumInteractionsEarly", &event_PUOOTnumInteractionsEarly, &b_event_PUOOTnumInteractionsEarly);
   fChain->SetBranchAddress("event_PUOOTnumInteractionsLate", &event_PUOOTnumInteractionsLate, &b_event_PUOOTnumInteractionsLate);
   fChain->SetBranchAddress("event_nRecVtx", &event_nRecVtx, &b_event_nRecVtx);
   fChain->SetBranchAddress("event_pass12whoissiglike", &event_pass12whoissiglike, &b_event_pass12whoissiglike);
   fChain->SetBranchAddress("event_CSCTightHaloID", &event_CSCTightHaloID, &b_event_CSCTightHaloID);
   fChain->SetBranchAddress("event_NMuons", &event_NMuons, &b_event_NMuons);
   fChain->SetBranchAddress("event_NMuonsTot", &event_NMuonsTot, &b_event_NMuonsTot);
   fChain->SetBranchAddress("dipho_mgg_photon", &dipho_mgg_photon, &b_dipho_mgg_photon);
   fChain->SetBranchAddress("dipho_mgg_newCorr", &dipho_mgg_newCorr, &b_dipho_mgg_newCorr);
   fChain->SetBranchAddress("dipho_mgg_newCorrLocal", &dipho_mgg_newCorrLocal, &b_dipho_mgg_newCorrLocal);
   fChain->SetBranchAddress("pholead_eta", &pholead_eta, &b_pholead_eta);
   fChain->SetBranchAddress("photrail_eta", &photrail_eta, &b_photrail_eta);
   fChain->SetBranchAddress("pholead_phi", &pholead_phi, &b_pholead_phi);
   fChain->SetBranchAddress("photrail_phi", &photrail_phi, &b_photrail_phi);
   fChain->SetBranchAddress("pholead_px", &pholead_px, &b_pholead_px);
   fChain->SetBranchAddress("photrail_px", &photrail_px, &b_photrail_px);
   fChain->SetBranchAddress("pholead_py", &pholead_py, &b_pholead_py);
   fChain->SetBranchAddress("photrail_py", &photrail_py, &b_photrail_py);
   fChain->SetBranchAddress("pholead_pt", &pholead_pt, &b_pholead_pt);
   fChain->SetBranchAddress("photrail_pt", &photrail_pt, &b_photrail_pt);
   fChain->SetBranchAddress("pholead_pz", &pholead_pz, &b_pholead_pz);
   fChain->SetBranchAddress("photrail_pz", &photrail_pz, &b_photrail_pz);
   fChain->SetBranchAddress("pholead_energy", &pholead_energy, &b_pholead_energy);
   fChain->SetBranchAddress("photrail_energy", &photrail_energy, &b_photrail_energy);
   fChain->SetBranchAddress("pholead_energySCdefault", &pholead_energySCdefault, &b_pholead_energySCdefault);
   fChain->SetBranchAddress("photrail_energySCdefault", &photrail_energySCdefault, &b_photrail_energySCdefault);
   fChain->SetBranchAddress("pholead_energyNewCorr", &pholead_energyNewCorr, &b_pholead_energyNewCorr);
   fChain->SetBranchAddress("photrail_energyNewCorr", &photrail_energyNewCorr, &b_photrail_energyNewCorr);
   fChain->SetBranchAddress("pholead_energyNewCorrLocal", &pholead_energyNewCorrLocal, &b_pholead_energyNewCorrLocal);
   fChain->SetBranchAddress("photrail_energyNewCorrLocal", &photrail_energyNewCorrLocal, &b_photrail_energyNewCorrLocal);
   fChain->SetBranchAddress("pholead_SCeta", &pholead_SCeta, &b_pholead_SCeta);
   fChain->SetBranchAddress("photrail_SCeta", &photrail_SCeta, &b_photrail_SCeta);
   fChain->SetBranchAddress("pholead_SCphi", &pholead_SCphi, &b_pholead_SCphi);
   fChain->SetBranchAddress("photrail_SCphi", &photrail_SCphi, &b_photrail_SCphi);
   fChain->SetBranchAddress("pholead_PhoHasPixSeed", &pholead_PhoHasPixSeed, &b_pholead_PhoHasPixSeed);
   fChain->SetBranchAddress("pholead_PhoHasConvTrks", &pholead_PhoHasConvTrks, &b_pholead_PhoHasConvTrks);
   fChain->SetBranchAddress("pholead_PhoScSeedSeverity", &pholead_PhoScSeedSeverity, &b_pholead_PhoScSeedSeverity);
   fChain->SetBranchAddress("photrail_PhoHasPixSeed", &photrail_PhoHasPixSeed, &b_photrail_PhoHasPixSeed);
   fChain->SetBranchAddress("photrail_PhoHasConvTrks", &photrail_PhoHasConvTrks, &b_photrail_PhoHasConvTrks);
   fChain->SetBranchAddress("photrail_PhoScSeedSeverity", &photrail_PhoScSeedSeverity, &b_photrail_PhoScSeedSeverity);
   fChain->SetBranchAddress("pholead_r9", &pholead_r9, &b_pholead_r9);
   fChain->SetBranchAddress("photrail_r9", &photrail_r9, &b_photrail_r9);
   fChain->SetBranchAddress("pholead_sieie", &pholead_sieie, &b_pholead_sieie);
   fChain->SetBranchAddress("photrail_sieie", &photrail_sieie, &b_photrail_sieie);
   fChain->SetBranchAddress("pholead_hoe", &pholead_hoe, &b_pholead_hoe);
   fChain->SetBranchAddress("photrail_hoe", &photrail_hoe, &b_photrail_hoe);
   fChain->SetBranchAddress("pholead_brem", &pholead_brem, &b_pholead_brem);
   fChain->SetBranchAddress("photrail_brem", &photrail_brem, &b_photrail_brem);
   fChain->SetBranchAddress("pholead_sigmaPhi", &pholead_sigmaPhi, &b_pholead_sigmaPhi);
   fChain->SetBranchAddress("photrail_sigmaPhi", &photrail_sigmaPhi, &b_photrail_sigmaPhi);
   fChain->SetBranchAddress("pholead_sigmaEta", &pholead_sigmaEta, &b_pholead_sigmaEta);
   fChain->SetBranchAddress("photrail_sigmaEta", &photrail_sigmaEta, &b_photrail_sigmaEta);
   fChain->SetBranchAddress("pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx", &pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx, &b_pholead_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx);
   fChain->SetBranchAddress("photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx", &photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx, &b_photrail_pho_Cone01PhotonIso_dEta015EB_dR070EE_mvVtx);
   fChain->SetBranchAddress("pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx", &pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx, &b_pholead_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx);
   fChain->SetBranchAddress("photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx", &photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx, &b_photrail_pho_Cone02PhotonIso_dEta015EB_dR070EE_mvVtx);
   fChain->SetBranchAddress("pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx", &pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx, &b_pholead_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx);
   fChain->SetBranchAddress("photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx", &photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx, &b_photrail_pho_Cone03PhotonIso_dEta015EB_dR070EE_mvVtx);
   fChain->SetBranchAddress("pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx", &pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx, &b_pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx);
   fChain->SetBranchAddress("photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx", &photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx, &b_photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx);
   fChain->SetBranchAddress("pholead_pho_Cone01NeutralHadronIso_mvVtx", &pholead_pho_Cone01NeutralHadronIso_mvVtx, &b_pholead_pho_Cone01NeutralHadronIso_mvVtx);
   fChain->SetBranchAddress("photrail_pho_Cone01NeutralHadronIso_mvVtx", &photrail_pho_Cone01NeutralHadronIso_mvVtx, &b_photrail_pho_Cone01NeutralHadronIso_mvVtx);
   fChain->SetBranchAddress("pholead_pho_Cone02NeutralHadronIso_mvVtx", &pholead_pho_Cone02NeutralHadronIso_mvVtx, &b_pholead_pho_Cone02NeutralHadronIso_mvVtx);
   fChain->SetBranchAddress("photrail_pho_Cone02NeutralHadronIso_mvVtx", &photrail_pho_Cone02NeutralHadronIso_mvVtx, &b_photrail_pho_Cone02NeutralHadronIso_mvVtx);
   fChain->SetBranchAddress("pholead_pho_Cone03NeutralHadronIso_mvVtx", &pholead_pho_Cone03NeutralHadronIso_mvVtx, &b_pholead_pho_Cone03NeutralHadronIso_mvVtx);
   fChain->SetBranchAddress("photrail_pho_Cone03NeutralHadronIso_mvVtx", &photrail_pho_Cone03NeutralHadronIso_mvVtx, &b_photrail_pho_Cone03NeutralHadronIso_mvVtx);
   fChain->SetBranchAddress("pholead_pho_Cone04NeutralHadronIso_mvVtx", &pholead_pho_Cone04NeutralHadronIso_mvVtx, &b_pholead_pho_Cone04NeutralHadronIso_mvVtx);
   fChain->SetBranchAddress("photrail_pho_Cone04NeutralHadronIso_mvVtx", &photrail_pho_Cone04NeutralHadronIso_mvVtx, &b_photrail_pho_Cone04NeutralHadronIso_mvVtx);
   fChain->SetBranchAddress("pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01", &pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01, &b_pholead_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01);
   fChain->SetBranchAddress("photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01", &photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01, &b_photrail_pho_Cone01ChargedHadronIso_dR02_dz02_dxy01);
   fChain->SetBranchAddress("pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01", &pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01, &b_pholead_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01);
   fChain->SetBranchAddress("photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01", &photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01, &b_photrail_pho_Cone02ChargedHadronIso_dR02_dz02_dxy01);
   fChain->SetBranchAddress("pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01", &pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01, &b_pholead_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01);
   fChain->SetBranchAddress("photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01", &photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01, &b_photrail_pho_Cone03ChargedHadronIso_dR02_dz02_dxy01);
   fChain->SetBranchAddress("pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01", &pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01, &b_pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01);
   fChain->SetBranchAddress("photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01", &photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01, &b_photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01);
   fChain->SetBranchAddress("pholead_pho_Cone03PFCombinedIso", &pholead_pho_Cone03PFCombinedIso, &b_pholead_pho_Cone03PFCombinedIso);
   fChain->SetBranchAddress("photrail_pho_Cone03PFCombinedIso", &photrail_pho_Cone03PFCombinedIso, &b_photrail_pho_Cone03PFCombinedIso);
   fChain->SetBranchAddress("pholead_pho_Cone04PFCombinedIso", &pholead_pho_Cone04PFCombinedIso, &b_pholead_pho_Cone04PFCombinedIso);
   fChain->SetBranchAddress("photrail_pho_Cone04PFCombinedIso", &photrail_pho_Cone04PFCombinedIso, &b_photrail_pho_Cone04PFCombinedIso);
   fChain->SetBranchAddress("pholead_PhoPassConvSafeElectronVeto", &pholead_PhoPassConvSafeElectronVeto, &b_pholead_PhoPassConvSafeElectronVeto);
   fChain->SetBranchAddress("photrail_PhoPassConvSafeElectronVeto", &photrail_PhoPassConvSafeElectronVeto, &b_photrail_PhoPassConvSafeElectronVeto);
   fChain->SetBranchAddress("pholead_GenPhotonIsoDR04", &pholead_GenPhotonIsoDR04, &b_pholead_GenPhotonIsoDR04);
   fChain->SetBranchAddress("photrail_GenPhotonIsoDR04", &photrail_GenPhotonIsoDR04, &b_photrail_GenPhotonIsoDR04);
   fChain->SetBranchAddress("pholead_PhoIso03Ecal", &pholead_PhoIso03Ecal, &b_pholead_PhoIso03Ecal);
   fChain->SetBranchAddress("pholead_PhoIso03Hcal", &pholead_PhoIso03Hcal, &b_pholead_PhoIso03Hcal);
   fChain->SetBranchAddress("pholead_PhoIso03TrkSolid", &pholead_PhoIso03TrkSolid, &b_pholead_PhoIso03TrkSolid);
   fChain->SetBranchAddress("pholead_PhoIso03TrkHollow", &pholead_PhoIso03TrkHollow, &b_pholead_PhoIso03TrkHollow);
   fChain->SetBranchAddress("pholead_PhoIso03", &pholead_PhoIso03, &b_pholead_PhoIso03);
   fChain->SetBranchAddress("pholead_PhoIso04Ecal", &pholead_PhoIso04Ecal, &b_pholead_PhoIso04Ecal);
   fChain->SetBranchAddress("pholead_PhoIso04Hcal", &pholead_PhoIso04Hcal, &b_pholead_PhoIso04Hcal);
   fChain->SetBranchAddress("pholead_PhoIso04TrkSolid", &pholead_PhoIso04TrkSolid, &b_pholead_PhoIso04TrkSolid);
   fChain->SetBranchAddress("pholead_PhoIso04TrkHollow", &pholead_PhoIso04TrkHollow, &b_pholead_PhoIso04TrkHollow);
   fChain->SetBranchAddress("pholead_PhoIso04", &pholead_PhoIso04, &b_pholead_PhoIso04);
   fChain->SetBranchAddress("photrail_PhoIso03Ecal", &photrail_PhoIso03Ecal, &b_photrail_PhoIso03Ecal);
   fChain->SetBranchAddress("photrail_PhoIso03Hcal", &photrail_PhoIso03Hcal, &b_photrail_PhoIso03Hcal);
   fChain->SetBranchAddress("photrail_PhoIso03TrkSolid", &photrail_PhoIso03TrkSolid, &b_photrail_PhoIso03TrkSolid);
   fChain->SetBranchAddress("photrail_PhoIso03TrkHollow", &photrail_PhoIso03TrkHollow, &b_photrail_PhoIso03TrkHollow);
   fChain->SetBranchAddress("photrail_PhoIso03", &photrail_PhoIso03, &b_photrail_PhoIso03);
   fChain->SetBranchAddress("photrail_PhoIso04Ecal", &photrail_PhoIso04Ecal, &b_photrail_PhoIso04Ecal);
   fChain->SetBranchAddress("photrail_PhoIso04Hcal", &photrail_PhoIso04Hcal, &b_photrail_PhoIso04Hcal);
   fChain->SetBranchAddress("photrail_PhoIso04TrkSolid", &photrail_PhoIso04TrkSolid, &b_photrail_PhoIso04TrkSolid);
   fChain->SetBranchAddress("photrail_PhoIso04TrkHollow", &photrail_PhoIso04TrkHollow, &b_photrail_PhoIso04TrkHollow);
   fChain->SetBranchAddress("photrail_PhoIso04", &photrail_PhoIso04, &b_photrail_PhoIso04);
   fChain->SetBranchAddress("pholead_PhoS4OverS1", &pholead_PhoS4OverS1, &b_pholead_PhoS4OverS1);
   fChain->SetBranchAddress("pholead_PhoSigmaEtaEta", &pholead_PhoSigmaEtaEta, &b_pholead_PhoSigmaEtaEta);
   fChain->SetBranchAddress("pholead_PhoE1x5", &pholead_PhoE1x5, &b_pholead_PhoE1x5);
   fChain->SetBranchAddress("pholead_PhoE2x5", &pholead_PhoE2x5, &b_pholead_PhoE2x5);
   fChain->SetBranchAddress("pholead_PhoE3x3", &pholead_PhoE3x3, &b_pholead_PhoE3x3);
   fChain->SetBranchAddress("pholead_PhoE5x5", &pholead_PhoE5x5, &b_pholead_PhoE5x5);
   fChain->SetBranchAddress("pholead_PhomaxEnergyXtal", &pholead_PhomaxEnergyXtal, &b_pholead_PhomaxEnergyXtal);
   fChain->SetBranchAddress("pholead_PhoIso03HcalDepth1", &pholead_PhoIso03HcalDepth1, &b_pholead_PhoIso03HcalDepth1);
   fChain->SetBranchAddress("pholead_PhoIso03HcalDepth2", &pholead_PhoIso03HcalDepth2, &b_pholead_PhoIso03HcalDepth2);
   fChain->SetBranchAddress("pholead_PhoIso04HcalDepth1", &pholead_PhoIso04HcalDepth1, &b_pholead_PhoIso04HcalDepth1);
   fChain->SetBranchAddress("pholead_PhoIso04HcalDepth2", &pholead_PhoIso04HcalDepth2, &b_pholead_PhoIso04HcalDepth2);
   fChain->SetBranchAddress("pholead_PhoIso03nTrksSolid", &pholead_PhoIso03nTrksSolid, &b_pholead_PhoIso03nTrksSolid);
   fChain->SetBranchAddress("pholead_PhoIso03nTrksHollow", &pholead_PhoIso03nTrksHollow, &b_pholead_PhoIso03nTrksHollow);
   fChain->SetBranchAddress("pholead_PhoIso04nTrksSolid", &pholead_PhoIso04nTrksSolid, &b_pholead_PhoIso04nTrksSolid);
   fChain->SetBranchAddress("pholead_PhoIso04nTrksHollow", &pholead_PhoIso04nTrksHollow, &b_pholead_PhoIso04nTrksHollow);
   fChain->SetBranchAddress("pholead_Pho_ChargedHadronIso", &pholead_Pho_ChargedHadronIso, &b_pholead_Pho_ChargedHadronIso);
   fChain->SetBranchAddress("pholead_Pho_NeutralHadronIso", &pholead_Pho_NeutralHadronIso, &b_pholead_Pho_NeutralHadronIso);
   fChain->SetBranchAddress("pholead_Pho_PhotonIso", &pholead_Pho_PhotonIso, &b_pholead_Pho_PhotonIso);
   fChain->SetBranchAddress("pholead_Pho_isPFPhoton", &pholead_Pho_isPFPhoton, &b_pholead_Pho_isPFPhoton);
   fChain->SetBranchAddress("pholead_Pho_isPFElectron", &pholead_Pho_isPFElectron, &b_pholead_Pho_isPFElectron);
   fChain->SetBranchAddress("photrail_PhoS4OverS1", &photrail_PhoS4OverS1, &b_photrail_PhoS4OverS1);
   fChain->SetBranchAddress("photrail_PhoSigmaEtaEta", &photrail_PhoSigmaEtaEta, &b_photrail_PhoSigmaEtaEta);
   fChain->SetBranchAddress("photrail_PhoE1x5", &photrail_PhoE1x5, &b_photrail_PhoE1x5);
   fChain->SetBranchAddress("photrail_PhoE2x5", &photrail_PhoE2x5, &b_photrail_PhoE2x5);
   fChain->SetBranchAddress("photrail_PhoE3x3", &photrail_PhoE3x3, &b_photrail_PhoE3x3);
   fChain->SetBranchAddress("photrail_PhoE5x5", &photrail_PhoE5x5, &b_photrail_PhoE5x5);
   fChain->SetBranchAddress("photrail_PhomaxEnergyXtal", &photrail_PhomaxEnergyXtal, &b_photrail_PhomaxEnergyXtal);
   fChain->SetBranchAddress("photrail_PhoIso03HcalDepth1", &photrail_PhoIso03HcalDepth1, &b_photrail_PhoIso03HcalDepth1);
   fChain->SetBranchAddress("photrail_PhoIso03HcalDepth2", &photrail_PhoIso03HcalDepth2, &b_photrail_PhoIso03HcalDepth2);
   fChain->SetBranchAddress("photrail_PhoIso04HcalDepth1", &photrail_PhoIso04HcalDepth1, &b_photrail_PhoIso04HcalDepth1);
   fChain->SetBranchAddress("photrail_PhoIso04HcalDepth2", &photrail_PhoIso04HcalDepth2, &b_photrail_PhoIso04HcalDepth2);
   fChain->SetBranchAddress("photrail_PhoIso03nTrksSolid", &photrail_PhoIso03nTrksSolid, &b_photrail_PhoIso03nTrksSolid);
   fChain->SetBranchAddress("photrail_PhoIso03nTrksHollow", &photrail_PhoIso03nTrksHollow, &b_photrail_PhoIso03nTrksHollow);
   fChain->SetBranchAddress("photrail_PhoIso04nTrksSolid", &photrail_PhoIso04nTrksSolid, &b_photrail_PhoIso04nTrksSolid);
   fChain->SetBranchAddress("photrail_PhoIso04nTrksHollow", &photrail_PhoIso04nTrksHollow, &b_photrail_PhoIso04nTrksHollow);
   fChain->SetBranchAddress("photrail_Pho_ChargedHadronIso", &photrail_Pho_ChargedHadronIso, &b_photrail_Pho_ChargedHadronIso);
   fChain->SetBranchAddress("photrail_Pho_NeutralHadronIso", &photrail_Pho_NeutralHadronIso, &b_photrail_Pho_NeutralHadronIso);
   fChain->SetBranchAddress("photrail_Pho_PhotonIso", &photrail_Pho_PhotonIso, &b_photrail_Pho_PhotonIso);
   fChain->SetBranchAddress("photrail_Pho_isPFPhoton", &photrail_Pho_isPFPhoton, &b_photrail_Pho_isPFPhoton);
   fChain->SetBranchAddress("photrail_Pho_isPFElectron", &photrail_Pho_isPFElectron, &b_photrail_Pho_isPFElectron);
   fChain->SetBranchAddress("pholead_PhoMCmatchindex", &pholead_PhoMCmatchindex, &b_pholead_PhoMCmatchindex);
   fChain->SetBranchAddress("pholead_PhoMCmatchexitcode", &pholead_PhoMCmatchexitcode, &b_pholead_PhoMCmatchexitcode);
   fChain->SetBranchAddress("photrail_PhoMCmatchindex", &photrail_PhoMCmatchindex, &b_photrail_PhoMCmatchindex);
   fChain->SetBranchAddress("photrail_PhoMCmatchexitcode", &photrail_PhoMCmatchexitcode, &b_photrail_PhoMCmatchexitcode);
   fChain->SetBranchAddress("pholead_Npfcandphotonincone",&pholead_Npfcandphotonincone, &b_pholead_Npfcandphotonincone);
   fChain->SetBranchAddress("pholead_Npfcandchargedincone",&pholead_Npfcandchargedincone, &b_pholead_Npfcandchargedincone);
   fChain->SetBranchAddress("pholead_Npfcandneutralincone",&pholead_Npfcandneutralincone, &b_pholead_Npfcandneutralincone);
   fChain->SetBranchAddress("photrail_Npfcandphotonincone",&photrail_Npfcandphotonincone, &b_photrail_Npfcandphotonincone);
   fChain->SetBranchAddress("photrail_Npfcandchargedincone",&photrail_Npfcandchargedincone, &b_photrail_Npfcandchargedincone);
   fChain->SetBranchAddress("photrail_Npfcandneutralincone",&photrail_Npfcandneutralincone, &b_photrail_Npfcandneutralincone);
   fChain->SetBranchAddress("pholead_scareaSF",&pholead_scareaSF, &b_pholead_scareaSF);
   fChain->SetBranchAddress("photrail_scareaSF",&photrail_scareaSF, &b_photrail_scareaSF);
   fChain->SetBranchAddress("pholead_photonpfcandenergies",&pholead_photonpfcandenergies, &b_pholead_photonpfcandenergies);
   fChain->SetBranchAddress("pholead_photonpfcandets",&pholead_photonpfcandets, &b_pholead_photonpfcandets);
   fChain->SetBranchAddress("pholead_photonpfcanddetas",&pholead_photonpfcanddetas, &b_pholead_photonpfcanddetas);
   fChain->SetBranchAddress("pholead_photonpfcanddphis",&pholead_photonpfcanddphis, &b_pholead_photonpfcanddphis);
   fChain->SetBranchAddress("photrail_photonpfcandenergies",&photrail_photonpfcandenergies, &b_photrail_photonpfcandenergies);
   fChain->SetBranchAddress("photrail_photonpfcandets",&photrail_photonpfcandets, &b_photrail_photonpfcandets);
   fChain->SetBranchAddress("photrail_photonpfcanddetas",&photrail_photonpfcanddetas, &b_photrail_photonpfcanddetas);
   fChain->SetBranchAddress("photrail_photonpfcanddphis",&photrail_photonpfcanddphis, &b_photrail_photonpfcanddphis);

   Notify();
}

Bool_t template_production::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void template_production::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

/* Int_t template_production::Cut(Long64_t entry) */
/* { */
/* // This function may be called from Loop. */
/* // returns  1 if entry is accepted. */
/* // returns -1 otherwise. */
/*    return 1; */
/* } */

void template_production::WriteOutput(const char* filename, const TString _dirname){
  TFile *out = TFile::Open(filename,"update");
  TString dirname("");
  if (isdata) dirname.Append("data_"); else dirname.Append("mc_");
  dirname.Append(_dirname.Data());
  out->mkdir(dirname.Data());
  out->cd(dirname.Data());

  out->mkdir(Form("%s/roofit",dirname.Data()));
  out->cd(Form("%s/roofit",dirname.Data()));

  roovar1->Write();
  roovar2->Write();
  roopt1->Write();
  roopt2->Write();
  roosieie1->Write();
  roosieie2->Write();
  rooeta1->Write();
  rooeta2->Write();
  roorho->Write();
  roosigma->Write();
  rooweight->Write();

  roovar_invmass->Write();
  roovar_diphotonpt->Write();
  roovar_costhetastar->Write();
  roovar_dphi->Write();
  roovar_dR->Write();
  

  if (dosignaltemplate || dobackgroundtemplate) {


    for (int k=0; k<2; k++) for (int i=0; i<2; i++) (roodset_signal[i][k])->Write();
    for (int k=0; k<2; k++) for (int i=0; i<2; i++) (roodset_background[i][k])->Write();

  }

  if (do2dtemplate){ 

    for (std::map<TString, RooDataSet*>::const_iterator it = template2d_roodset.begin(); it!=template2d_roodset.end(); it++) (it->second)->Write();

  }

  if (dodistribution) {

    for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){
      for (int i=0; i<3; i++)
	for (int j=0; j<n_templates; j++) {
	  obs_roodset[get_name_obs_roodset(i,*diffvariable,n_templates)]->append(*(obs_roodset[get_name_obs_roodset(i,*diffvariable,j)]));
	}
    }

    for (std::map<TString, RooDataSet*>::const_iterator it = obs_roodset.begin(); it!=obs_roodset.end(); it++) (it->second)->Write();

    out->mkdir(Form("%s/gen_purity",dirname.Data()));
    out->cd(Form("%s/gen_purity",dirname.Data()));

    for (std::map<TString, TProfile*>::const_iterator it = true_purity.begin(); it!=true_purity.end(); it++) (it->second)->Write();
    for (std::map<TString, TH1F*>::const_iterator it = true_purity_isppevent.begin(); it!=true_purity_isppevent.end(); it++) (it->second)->Write();
    for (std::map<TString, TH1F*>::const_iterator it = true_purity_isnotppevent.begin(); it!=true_purity_isnotppevent.end(); it++) (it->second)->Write();

  }

  out->mkdir(Form("%s/plots",dirname.Data()));
  out->cd(Form("%s/plots",dirname.Data()));

  {

    RooArgList toplot;
    toplot.add(*rooargset_diffvariables);
    toplot.add(RooArgSet(*roovar1,*roovar2,*roopt1,*roosieie1,*rooeta1,*roopt2,*roosieie2,*rooeta2));
    toplot.add(RooArgSet(*roorho,*roosigma));
    
    for (int k=0; k<toplot.getSize(); k++){
      
      if (do2dtemplate){
	for (std::map<TString, RooDataSet*>::const_iterator it = template2d_roodset.begin(); it!=template2d_roodset.end(); it++) {
	  RooRealVar *thisvar = (RooRealVar*)(toplot.at(k));
	  RooPlot *thisplot = thisvar->frame(Name(Form("%s_%s",thisvar->GetName(),(it->second)->GetName())));
	  (it->second)->plotOn(thisplot);
	  thisplot->Write();
	}
      }
      if (dodistribution){
	for (std::map<TString, RooDataSet*>::const_iterator it = obs_roodset.begin(); it!=obs_roodset.end(); it++) {
	  RooRealVar *thisvar = (RooRealVar*)(toplot.at(k));
	  RooPlot *thisplot = thisvar->frame(Name(Form("%s_%s",thisvar->GetName(),(it->second)->GetName())));
	  (it->second)->plotOn(thisplot);
	  thisplot->Write();
	}
      }

    }

  }

  std::cout << "Writing output..." << std::endl;


  out->Close();
};



TString template_production::get_name_obs_roodset(int region, TString diffvariable, int bin){
  TString name_signal="obs_roodset";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s_b%d",name_signal.Data(),reg.Data(),diffvariable.Data(),bin);
  return t;
};

TString template_production::get_name_true_purity(int region, TString diffvariable){
  TString name_signal="true_purity";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s",name_signal.Data(),reg.Data(),diffvariable.Data());
  return t;
};

TString template_production::get_name_true_purity_ispp(int region, TString diffvariable){
  TString name_signal="ispp_true_purity";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s",name_signal.Data(),reg.Data(),diffvariable.Data());
  return t;
};

TString template_production::get_name_true_purity_isnotpp(int region, TString diffvariable){
  TString name_signal="isnotpp_true_purity";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s",name_signal.Data(),reg.Data(),diffvariable.Data());
  return t;
};

TString template_production::get_name_template2d_roodset(int region, TString sigorbkg){
  TString name_signal="template_roodset";
  TString reg;
  if (region==0) reg="EBEB"; else if (region==1) reg="EBEE"; else if (region==2) reg="EEEE"; else if (region==3) reg="EEEB";
  TString t=Form("%s_%s_%s",name_signal.Data(),reg.Data(),sigorbkg.Data());
  return t;
};

Int_t template_production::Choose_bin_invmass(float invmass, int region){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_diphoton_invmass_EBEB; index=n_templates_invmass_EBEB;}
  if (region==2) {cuts=binsdef_diphoton_invmass_EEEE; index=n_templates_invmass_EEEE;}
  if (region==3) {cuts=binsdef_diphoton_invmass_EBEE; index=n_templates_invmass_EBEE;}
  if (region==4) {cuts=binsdef_diphoton_invmass_EBEE; index=n_templates_invmass_EBEE;}
  if (region==1) {cuts=binsdef_diphoton_invmass_EBEE; index=n_templates_invmass_EBEE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]=9999;

  if (invmass<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value mass " << invmass << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((invmass>=cuts[i]) && (invmass<cuts[i+1])) return i;
  
  std::cout << "WARNING: called bin choice for out-of-range value mass " << invmass << std::endl;
  return -999;


};

Int_t template_production::Choose_bin_diphotonpt(float diphotonpt, int region){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_diphoton_diphotonpt_EBEB; index=n_templates_diphotonpt_EBEB;}
  if (region==2) {cuts=binsdef_diphoton_diphotonpt_EEEE; index=n_templates_diphotonpt_EEEE;}
  if (region==3) {cuts=binsdef_diphoton_diphotonpt_EBEE; index=n_templates_diphotonpt_EBEE;}
  if (region==4) {cuts=binsdef_diphoton_diphotonpt_EBEE; index=n_templates_diphotonpt_EBEE;}
  if (region==1) {cuts=binsdef_diphoton_diphotonpt_EBEE; index=n_templates_diphotonpt_EBEE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]=9999;

  if (diphotonpt<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value diphotonpt " << diphotonpt << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((diphotonpt>=cuts[i]) && (diphotonpt<cuts[i+1])) return i;
  
  std::cout << "WARNING: called bin choice for out-of-range value diphotonpt " << diphotonpt << std::endl;
  return -999;


};

Int_t template_production::Choose_bin_costhetastar(float costhetastar, int region){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_diphoton_costhetastar_EBEB; index=n_templates_costhetastar_EBEB;}
  if (region==2) {cuts=binsdef_diphoton_costhetastar_EEEE; index=n_templates_costhetastar_EEEE;}
  if (region==3) {cuts=binsdef_diphoton_costhetastar_EBEE; index=n_templates_costhetastar_EBEE;}
  if (region==4) {cuts=binsdef_diphoton_costhetastar_EBEE; index=n_templates_costhetastar_EBEE;}
  if (region==1) {cuts=binsdef_diphoton_costhetastar_EBEE; index=n_templates_costhetastar_EBEE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]=9999;

  if (costhetastar<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << costhetastar << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((costhetastar>=cuts[i]) && (costhetastar<cuts[i+1])) return i;
  
  std::cout << "WARNING: called bin choice for out-of-range value " << costhetastar << std::endl;
  return -999;


};

Int_t template_production::Choose_bin_dphi(float dphi, int region){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_diphoton_dphi_EBEB; index=n_templates_dphi_EBEB;}
  if (region==2) {cuts=binsdef_diphoton_dphi_EEEE; index=n_templates_dphi_EEEE;}
  if (region==3) {cuts=binsdef_diphoton_dphi_EBEE; index=n_templates_dphi_EBEE;}
  if (region==4) {cuts=binsdef_diphoton_dphi_EBEE; index=n_templates_dphi_EBEE;}
  if (region==1) {cuts=binsdef_diphoton_dphi_EBEE; index=n_templates_dphi_EBEE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]=9999;

  if (dphi<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << dphi << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((dphi>=cuts[i]) && (dphi<cuts[i+1])) return i;
  
  std::cout << "WARNING: called bin choice for out-of-range value " << dphi << std::endl;
  return -999;


};

Int_t template_production::Choose_bin_dR(float dR, int region){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_diphoton_dR_EBEB; index=n_templates_dR_EBEB;}
  if (region==2) {cuts=binsdef_diphoton_dR_EEEE; index=n_templates_dR_EEEE;}
  if (region==3) {cuts=binsdef_diphoton_dR_EBEE; index=n_templates_dR_EBEE;}
  if (region==4) {cuts=binsdef_diphoton_dR_EBEE; index=n_templates_dR_EBEE;}
  if (region==1) {cuts=binsdef_diphoton_dR_EBEE; index=n_templates_dR_EBEE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]=9999;

  if (dR<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << dR << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((dR>=cuts[i]) && (dR<cuts[i+1])) return i;
  
  std::cout << "WARNING: called bin choice for out-of-range value " << dR << std::endl;
  return -999;


};


Int_t template_production::Choose_bin_pt(float pt, int region){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_single_gamma_EB; index=n_templates_EB;}
  if (region==1) {cuts=binsdef_single_gamma_EE; index=n_templates_EE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]=9999;

  if (pt<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << pt << " cuts[0]=" << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((pt>=cuts[i]) && (pt<cuts[i+1])) return i;
  
  std::cout << "WARNING: called bin choice for out-of-range value " << pt << std::endl;
  return -999;


};

Int_t template_production::Choose_bin_eta(float eta, int region){

  eta=fabs(eta);

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_single_gamma_EB_eta; index=n_templates_EB;}
  if (region==1) {cuts=binsdef_single_gamma_EE_eta; index=n_templates_EE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]=9999;

  if (eta<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << eta << " cuts[0]=" << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((eta>=cuts[i]) && (eta<cuts[i+1])) return i;

  std::cout << "WARNING: called bin choice for out-of-range value " << eta << std::endl;
  return -999;

};


Int_t template_production::Choose_bin_sieie(float sieie, int region){

  if (region==1) return 0;

  if (sieie<0.008) return 0;
  if (sieie<0.009) return 1;
  if (sieie<0.010) return 2;
  if (sieie<0.011) return 3;
  if (sieie<0.012) return 4;
  if (sieie<0.013) return 5;
  if (sieie<0.014) return 6;
  return 7;

};

float template_production::AbsDeltaPhi(double phi1, double phi2){
  // From cmssw reco::deltaPhi()
  double result = phi1 - phi2;
  while( result >   TMath::Pi() ) result -= TMath::TwoPi();
  while( result <= -TMath::Pi() ) result += TMath::TwoPi();
  return TMath::Abs(result);
}

void template_production::FillDiffVariables(){

  roovar_invmass->setVal(dipho_mgg_photon);
  {
    float px = pholead_px+photrail_px;
    float py = pholead_py+photrail_py;
    float pt = sqrt(px*px+py*py);    
    roovar_diphotonpt->setVal(pt);
  }
  {
    TLorentzVector pho1(pholead_px,pholead_py,pholead_pz,pholead_energy);
    TLorentzVector pho2(photrail_px,photrail_py,photrail_pz,photrail_energy);
    
    // COS THETASTAR HX
    //	  TVector3 boost = (pho1+pho2).BoostVector();
    //	  TLorentzVector boostedpho1 = pho1;
    //	  boostedpho1.Boost(-boost);
    //	  float thetastar1 = boostedpho1.Angle(boost);
    //	  bin_couple = Choose_bin_costhetastar(fabs(TMath::Cos(thetastar1)),event_ok_for_dataset_local);
    //	  value_diffvariable=fabs(TMath::Cos(thetastar1));
    
    // DELTA ETA
    //	  value_diffvariable = fabs(TMath::TanH((pho1.Rapidity()-pho2.Rapidity())/2));
    //	  bin_couple = Choose_bin_costhetastar(value_diffvariable,event_ok_for_dataset_local);
    
    // COS THETASTAR CS
    TLorentzVector b1,b2,diphoton;
    b1.SetPx(0); b1.SetPy(0); b1.SetPz( 3500); b1.SetE(3500);
    b2.SetPx(0); b2.SetPy(0); b2.SetPz(-3500); b2.SetE(3500);
    TLorentzVector boostedpho1 = pho1; 
    TLorentzVector boostedpho2 = pho2; 
    TLorentzVector boostedb1 = b1; 
    TLorentzVector boostedb2 = b2; 
    TVector3 boost = (pho1+pho2).BoostVector();
    boostedpho1.Boost(-boost);
    boostedpho2.Boost(-boost);
    boostedb1.Boost(-boost);
    boostedb2.Boost(-boost);
    TVector3 direction_cs = (boostedb1.Vect().Unit()-boostedb2.Vect().Unit()).Unit();
    roovar_costhetastar->setVal( fabs(TMath::Cos(direction_cs.Angle(boostedpho1.Vect()))) );
  }
  {
    float phi1 = pholead_phi;
    float phi2 = photrail_phi;
    float dphi = AbsDeltaPhi(phi1,phi2);
    roovar_dphi->setVal(dphi);
  }
  {
    float phi1 = pholead_phi;
    float phi2 = photrail_phi;
    float dphi = AbsDeltaPhi(phi1,phi2);
    float deta = pholead_eta-photrail_eta;
    float dR = sqrt(deta*deta+dphi*dphi);
    roovar_dR->setVal(dR);
  }
  
  return;

};

#endif // #ifdef template_production_cxx


