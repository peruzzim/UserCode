//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 30 14:21:13 2011 by ROOT version 5.30/02
// from TTree Tree/Tree
// found on file: mc_inclusive.root
//////////////////////////////////////////////////////////

#ifndef template_production_h
#define template_production_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "TRandom3.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TF1.h"

using namespace std;
using namespace RooFit;

class template_production {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         event_luminormfactor;
   Float_t         event_Kfactor;
   Float_t         event_weight;
   Float_t         event_rho;
   Int_t           event_nPU;
   Int_t           event_nRecVtx;
   Float_t         dipho_mgg_photon;
   Float_t         dipho_mgg_newCorr;
   Float_t         dipho_mgg_newCorrLocal;
   Float_t         pholead_eta;
   Float_t         photrail_eta;
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
   Float_t         pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt5;
   Float_t         pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt0;
   Float_t         pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt5;
   Float_t         pholead_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
   Float_t         pholead_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
   Float_t         pholead_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
   Float_t         pholead_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
   Float_t         pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0;
   Float_t         pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5;
   Float_t         pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks;
   Float_t         pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks;
   Float_t         pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0;
   Float_t         pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5;
   Float_t         pholead_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
   Float_t         pholead_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
   Float_t         pholead_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
   Float_t         pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
   Float_t         pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old;
   Float_t         pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old;
   Float_t         pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old;
   Float_t         pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old;
   Float_t         pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0;
   Float_t         pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
   Float_t         pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
   Float_t         pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0;
   Float_t         pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
   Float_t         pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;
   Float_t         pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0;
   Float_t         pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
   Float_t         pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
   Float_t         pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0;
   Float_t         pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
   Float_t         pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;
   Float_t         pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0;
   Float_t         pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
   Float_t         pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
   Float_t         pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0;
   Float_t         pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
   Float_t         pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;
   Float_t         pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0;
   Float_t         pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
   Float_t         pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
   Float_t         pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0;
   Float_t         pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
   Float_t         pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;
   Float_t         photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt5;
   Float_t         photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt0;
   Float_t         photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt5;
   Float_t         photrail_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
   Float_t         photrail_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
   Float_t         photrail_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
   Float_t         photrail_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;
   Float_t         photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0;
   Float_t         photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5;
   Float_t         photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks;
   Float_t         photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks;
   Float_t         photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0;
   Float_t         photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5;
   Float_t         photrail_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
   Float_t         photrail_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
   Float_t         photrail_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
   Float_t         photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx;
   Float_t         photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old;
   Float_t         photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old;
   Float_t         photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old;
   Float_t         photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old;
   Float_t         photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0;
   Float_t         photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
   Float_t         photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
   Float_t         photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0;
   Float_t         photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
   Float_t         photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;
   Float_t         photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0;
   Float_t         photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
   Float_t         photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
   Float_t         photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0;
   Float_t         photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
   Float_t         photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;
   Float_t         photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0;
   Float_t         photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
   Float_t         photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
   Float_t         photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0;
   Float_t         photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
   Float_t         photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;
   Float_t         photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0;
   Float_t         photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;
   Float_t         photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;
   Float_t         photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0;
   Float_t         photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;
   Float_t         photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;
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

   // List of branches
   TBranch        *b_event_luminormfactor;   //!
   TBranch        *b_event_Kfactor;   //!
   TBranch        *b_event_weight;   //!
   TBranch        *b_event_rho;   //!
   TBranch        *b_event_nPU;   //!
   TBranch        *b_event_nRecVtx;   //!
   TBranch        *b_dipho_mgg_photon;   //!
   TBranch        *b_dipho_mgg_newCorr;   //!
   TBranch        *b_dipho_mgg_newCorrLocal;   //!
   TBranch        *b_pholead_eta;   //!
   TBranch        *b_photrail_eta;   //!
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
   TBranch        *b_pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt5;   //!
   TBranch        *b_pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt0;   //!
   TBranch        *b_pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt5;   //!
   TBranch        *b_pholead_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_pholead_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_pholead_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_pholead_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0;   //!
   TBranch        *b_pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5;   //!
   TBranch        *b_pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks;   //!
   TBranch        *b_pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks;   //!
   TBranch        *b_pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0;   //!
   TBranch        *b_pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5;   //!
   TBranch        *b_pholead_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_pholead_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_pholead_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old;   //!
   TBranch        *b_pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old;   //!
   TBranch        *b_pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old;   //!
   TBranch        *b_pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old;   //!
   TBranch        *b_pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt5;   //!
   TBranch        *b_photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt0;   //!
   TBranch        *b_photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt5;   //!
   TBranch        *b_photrail_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_photrail_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_photrail_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_photrail_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx;   //!
   TBranch        *b_photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0;   //!
   TBranch        *b_photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5;   //!
   TBranch        *b_photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks;   //!
   TBranch        *b_photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks;   //!
   TBranch        *b_photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0;   //!
   TBranch        *b_photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5;   //!
   TBranch        *b_photrail_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_photrail_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_photrail_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx;   //!
   TBranch        *b_photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old;   //!
   TBranch        *b_photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old;   //!
   TBranch        *b_photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old;   //!
   TBranch        *b_photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old;   //!
   TBranch        *b_photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0;   //!
   TBranch        *b_photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU;   //!
   TBranch        *b_photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0;   //!
   TBranch        *b_photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01;   //!
   TBranch        *b_photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU;   //!
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

   template_production(TTree *tree=0);
   virtual ~template_production();
   /* virtual Int_t    Cut(Long64_t entry); */
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void     WriteOutput(const char* filename, const bool isdata, const TString dirname);
   void     ShowHistogram(int i, int j, int k, int l, int m);

   void Setup(TString _varname, Float_t _leftrange, Float_t _rightrange, Int_t _nbins, Bool_t _isdata, Bool_t _dopucorr);

   TProfile** GetPUScaling();

   TRandom3 *randomgen;

   static const int n_templates=3;

   RooRealVar *roovar[2][3][2];
   RooDataHist *roohist[3][2][3][n_templates][2];
   TH1F *templatehist[3][2][n_templates];
   RooDataSet *roodset[4][n_templates][2];
   RooDataSet *roodset_single[2][n_templates];

   TBranch *b_pholead_outvar;
   Float_t pholead_outvar;

   TBranch *b_photrail_outvar;
   Float_t photrail_outvar;

   Bool_t initialized;

   Bool_t isdata;

   Bool_t doabsolute;
   Bool_t dopucorr;

   Float_t leftrange, rightrange;
   Int_t nbins;

   TString varname;
   TString sidebandoption;

   TString get_roohist_name(TString varname, TString st, TString reg, TString count, int _binnumber,TString _ord);
   TString get_roodset_name(TString varname, TString reg, int _binnumber, TString _ord);
   TString get_roodset_name_single(TString _varname, TString _reg, int _binnumber);
   TString get_roovar_name(TString _varname, int i, int j, TString _ord);

   Int_t Choose_bin_invmass(float invmass);

};

#endif

#ifdef template_production_cxx
template_production::template_production(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mc_inclusive.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("mc_inclusive.root");
//      }
//      f->GetObject("Tree",tree);
//
//   }

   if (tree==0) std::cout << "Tree not ready!" << std::endl;
   if (!tree) return;
   fChain = tree;

   initialized=false;

}

void template_production::Setup(TString _varname, Float_t _leftrange, Float_t _rightrange, Int_t _nbins, Bool_t _isdata, Bool_t _dopucorr){

  varname=_varname;
  leftrange=_leftrange;
  rightrange=_rightrange;
  nbins=_nbins;
  isdata=_isdata;

  doabsolute=false;
  if (varname=="PhoIso04"){
     std::cout << "Building templates for absolute isolation" << std::endl;
     doabsolute=true;
  }

  dopucorr=_dopucorr;
  if (varname!="PhoIso04") {
    dopucorr=false;
    std::cout << "WARNING: NOT applying PU correction for this unknown variable" << std::endl;
  }
  if (dopucorr) std::cout << "Applying PU correction" << std::endl;


  Init();

  // roovar[EB,EE][1,2,both][noorder,order]
  for (int i=0; i<2; i++) for (int j=0; j<3; j++) for (int k=0; k<2; k++) {
	TString t=varname;
	if (i==0) t.Append("_EB"); else t.Append("_EE");
	if (j==0) t.Append("_1"); else if (j==1) t.Append("_2"); else t.Append("_both");
	if (k==0) t.Append("_noorder"); else t.Append("_order");
	roovar[i][j][k] = new RooRealVar(t.Data(),t.Data(),leftrange,rightrange);
	roovar[i][j][k]->setBins(nbins);
	std::cout << t.Data() << std::endl;
      }

   
  randomgen = new TRandom3(0);
  
  
  // roohist[sig,bkg,all][EB,EE][1,2,both][n_templates][noorder,order]
  for (int i=0; i<3; i++) {
    for (int j=0; j<2; j++) {
      for (int k=0; k<3; k++) {
	for (int l=0; l<n_templates; l++){
	  for (int m=0; m<2; m++){
	    TString st;
	    if (i==0) st="sig"; else if (i==1) st="bkg"; else if (i==2) st="all";
	    TString reg;
	    if (j==0) reg="EB"; else if (j==1) reg="EE";
	    TString count;
	    if (k==0) count="1"; else if (k==1) count="2"; else count="both";
	    TString ord;
	    if (m==0) ord="noorder"; else ord="order";
	    TString t=Form("roohist_%s_%s_%s_%s_b%d_%s",varname.Data(),st.Data(),reg.Data(),count.Data(),l,ord.Data());
	    std::cout << t.Data() << " " << roovar[j][k][m]->getTitle().Data() << std::endl;
	    roohist[i][j][k][l][m] = new RooDataHist(t.Data(),t.Data(),RooArgSet(*(roovar[j][k][m])));
	  }
	}
      }
    }
  }
  
  // templatehist[sig,bkg,all][EB,EE][n_templates]
  for (int i=0; i<3; i++) {
    for (int j=0;j<2;j++) {
      for (int l=0; l<n_templates; l++){
	TString st;
	if (i==0) st="sig"; else if (i==1) st="bkg"; else if (i==2) st="all";
	TString reg;
	if (j==0) reg="EB"; else if (j==1) reg="EE";
	TString t=Form("templatehist_%s_%s_%s_b%d",varname.Data(),st.Data(),reg.Data(),l);
	templatehist[i][j][l] = new TH1F(t.Data(),t.Data(),nbins,leftrange,rightrange);
      }
    }
  }
  
  
  // roodataset[EBEB,EBEE,EEEE,EEEB][n_templates][noorder,order]
  for (int i=0; i<4; i++){
    for (int l=0; l<n_templates; l++){
      for (int m=0; m<2; m++){
	TString reg;
	if (i==0) reg="EBEB"; else if (i==1) reg="EBEE"; else if (i==2) reg="EEEE"; else if (i==3) reg="EEEB";
	TString ord;
	if (m==0) ord="noorder"; else ord="order";
	TString t=Form("roodataset_%s_%s_b%d_%s",varname.Data(),reg.Data(),l,ord.Data());
	RooArgSet vars;
	if (i==0) { vars.add(*(roovar[0][0][m])); vars.add(*(roovar[0][1][m])); }
	if (i==2) { vars.add(*(roovar[1][0][m])); vars.add(*(roovar[1][1][m])); }
	if ((i==1 || i==3) && m==0) { vars.add(*(roovar[0][2][0])); vars.add(*(roovar[1][2][0])); }
	if (i==1 && m==1) { vars.add(*(roovar[0][0][1])); vars.add(*(roovar[1][1][1])); }
	if (i==3 && m==1) { vars.add(*(roovar[1][0][1])); vars.add(*(roovar[0][1][1])); }
	roodset[i][l][m] = new RooDataSet(t.Data(),t.Data(),vars);
      }
    }
  }
  
  // roodataset_single[EB,EE][bin]
  for (int i=0; i<2; i++) for (int j=0; j<n_templates; j++) {
      TString reg; if (i==0) reg="EB"; else reg="EE";
      TString t=Form("roodataset_single_%s_%s_b%d",varname.Data(),reg.Data(),j);
      roodset_single[i][j] = new RooDataSet(t.Data(),t.Data(),RooArgList(*roovar[i][2][0]));
    }

  initialized=true;
  
};

template_production::~template_production()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   delete randomgen;

   for (int m=0; m<2; m++) for (int j=0;j<2;j++) for (int k=0;k<3;k++) delete roovar[j][k][m];
   for (int m=0; m<2; m++) for (int l=0; l<n_templates; l++) for (int i=0; i<3; i++) for (int j=0;j<2;j++) for (int k=0;k<3;k++) delete roohist[i][j][k][l][m];
   for (int l=0; l<n_templates; l++) for (int i=0; i<3; i++) for (int j=0;j<2;j++) delete templatehist[i][j][l];
   for (int m=0; m<2; m++) for (int l=0; l<n_templates; l++) for (int i=0; i<4; i++) delete roodset[i][l][m];
   for (int l=0; l<n_templates; l++) for (int i=0; i<2; i++) delete roodset_single[i][l];

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
   fChain->SetBranchAddress("event_nPU", &event_nPU, &b_event_nPU);
   fChain->SetBranchAddress("event_nRecVtx", &event_nRecVtx, &b_event_nRecVtx);
   fChain->SetBranchAddress("dipho_mgg_photon", &dipho_mgg_photon, &b_dipho_mgg_photon);
   fChain->SetBranchAddress("dipho_mgg_newCorr", &dipho_mgg_newCorr, &b_dipho_mgg_newCorr);
   fChain->SetBranchAddress("dipho_mgg_newCorrLocal", &dipho_mgg_newCorrLocal, &b_dipho_mgg_newCorrLocal);
   fChain->SetBranchAddress("pholead_eta", &pholead_eta, &b_pholead_eta);
   fChain->SetBranchAddress("photrail_eta", &photrail_eta, &b_photrail_eta);
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
   fChain->SetBranchAddress("pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt5", &pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt5, &b_pholead_Pho_Cone04PhotonIso_dR0_dEta0_pt5);
   fChain->SetBranchAddress("pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt0", &pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt0, &b_pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt0);
   fChain->SetBranchAddress("pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt5", &pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt5, &b_pholead_Pho_Cone04PhotonIso_dR8_dEta0_pt5);
   fChain->SetBranchAddress("pholead_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", &pholead_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_pholead_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("pholead_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", &pholead_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_pholead_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("pholead_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", &pholead_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_pholead_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("pholead_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", &pholead_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_pholead_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0", &pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0, &b_pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0);
   fChain->SetBranchAddress("pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5", &pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5, &b_pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5);
   fChain->SetBranchAddress("pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks", &pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks, &b_pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks);
   fChain->SetBranchAddress("pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks", &pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks, &b_pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks);
   fChain->SetBranchAddress("pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0", &pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0, &b_pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0);
   fChain->SetBranchAddress("pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5", &pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5, &b_pholead_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5);
   fChain->SetBranchAddress("pholead_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx", &pholead_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_pholead_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("pholead_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx", &pholead_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_pholead_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("pholead_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx", &pholead_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_pholead_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx", &pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_pholead_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old", &pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old, &b_pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old);
   fChain->SetBranchAddress("pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old", &pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old, &b_pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old);
   fChain->SetBranchAddress("pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old", &pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old, &b_pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old);
   fChain->SetBranchAddress("pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old", &pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old, &b_pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old);
   fChain->SetBranchAddress("pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0", &pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", &pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", &pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_pholead_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0", &pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", &pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", &pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_pholead_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0", &pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", &pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", &pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_pholead_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0", &pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", &pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", &pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_pholead_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0", &pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", &pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", &pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_pholead_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0", &pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", &pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", &pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_pholead_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0", &pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", &pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", &pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_pholead_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0", &pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", &pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", &pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_pholead_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt5", &photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt5, &b_photrail_Pho_Cone04PhotonIso_dR0_dEta0_pt5);
   fChain->SetBranchAddress("photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt0", &photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt0, &b_photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt0);
   fChain->SetBranchAddress("photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt5", &photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt5, &b_photrail_Pho_Cone04PhotonIso_dR8_dEta0_pt5);
   fChain->SetBranchAddress("photrail_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", &photrail_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_photrail_Pho_Cone01PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("photrail_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", &photrail_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_photrail_Pho_Cone02PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("photrail_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", &photrail_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_photrail_Pho_Cone03PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("photrail_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx", &photrail_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx, &b_photrail_Pho_Cone04PhotonIso_dR045EB070EE_dEta015_pt08EB1EE_mvVtx);
   fChain->SetBranchAddress("photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0", &photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0, &b_photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0);
   fChain->SetBranchAddress("photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5", &photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5, &b_photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5);
   fChain->SetBranchAddress("photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks", &photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks, &b_photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_nocracks);
   fChain->SetBranchAddress("photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks", &photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks, &b_photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt5_nocracks);
   fChain->SetBranchAddress("photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0", &photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0, &b_photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt0);
   fChain->SetBranchAddress("photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5", &photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5, &b_photrail_Pho_Cone04NeutralHadronIso_dR7_dEta0_pt5);
   fChain->SetBranchAddress("photrail_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx", &photrail_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_photrail_Pho_Cone01NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("photrail_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx", &photrail_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_photrail_Pho_Cone02NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("photrail_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx", &photrail_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_photrail_Pho_Cone03NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx", &photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx, &b_photrail_Pho_Cone04NeutralHadronIso_dR0_dEta0_pt0_mvVtx);
   fChain->SetBranchAddress("photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old", &photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old, &b_photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0_old);
   fChain->SetBranchAddress("photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old", &photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old, &b_photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU_old);
   fChain->SetBranchAddress("photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old", &photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old, &b_photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0_old);
   fChain->SetBranchAddress("photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old", &photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old, &b_photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU_old);
   fChain->SetBranchAddress("photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0", &photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", &photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", &photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_photrail_Pho_Cone01ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0", &photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", &photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", &photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_photrail_Pho_Cone01ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0", &photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", &photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", &photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_photrail_Pho_Cone02ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0", &photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", &photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", &photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_photrail_Pho_Cone02ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0", &photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", &photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", &photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_photrail_Pho_Cone03ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0", &photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", &photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", &photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_photrail_Pho_Cone03ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0", &photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0, &b_photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz0);
   fChain->SetBranchAddress("photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01", &photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01, &b_photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU", &photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU, &b_photrail_Pho_Cone04ChargedHadronIso_dR0_dEta0_pt0_PFnoPU);
   fChain->SetBranchAddress("photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0", &photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0, &b_photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz0);
   fChain->SetBranchAddress("photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01", &photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01, &b_photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_dz1_dxy01);
   fChain->SetBranchAddress("photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU", &photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU, &b_photrail_Pho_Cone04ChargedHadronIso_dR015_dEta0_pt0_PFnoPU);
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
void template_production::WriteOutput(const char* filename, const bool _isdata, const TString _dirname){
  TFile *out = TFile::Open(filename,"update");
  TString dirname("");
  if (_isdata) dirname.Append("data_"); else dirname.Append("mc_");
  dirname.Append(_dirname.Data());
  out->mkdir(dirname.Data());
  out->cd(dirname.Data());

  for (int m=0; m<2; m++) for (int j=0;j<2;j++) for (int k=0;k<3;k++)  roovar[j][k][m]->Write();
  for (int m=0; m<2; m++) for (int l=0; l<n_templates; l++) for (int i=0; i<3; i++) for (int j=0;j<2;j++) for (int k=0;k<3;k++)  roohist[i][j][k][l][m]->Write();
  for (int l=0; l<n_templates; l++) for (int i=0; i<3; i++) for (int j=0;j<2;j++)  templatehist[i][j][l]->Write();
  for (int m=0; m<2; m++) for (int l=0; l<n_templates; l++) for (int i=0; i<4; i++)  roodset[i][l][m]->Write();
  for (int l=0; l<n_templates; l++) for (int i=0; i<2; i++)  roodset_single[i][l]->Write();

  std::cout << "output written" << std::endl;

  out->Close();
};

void template_production::ShowHistogram(int i, int j, int k, int l, int m){
  RooPlot *outroovarframe = roovar[j][k][m]->frame(Name(varname.Data()),Title(varname.Data()));
  roohist[i][j][k][l][m]->plotOn(outroovarframe);
  outroovarframe->Draw();
};

TString template_production::get_roohist_name(TString _varname, TString _st, TString _reg, TString _count, int _binnumber, TString _ord){
  TString a(Form("roohist_%s_%s_%s_%s_b%d_%s",_varname.Data(),_st.Data(),_reg.Data(),_count.Data(),_binnumber,_ord.Data()));
  return a;
};

TString template_production::get_roodset_name(TString _varname, TString _reg, int _binnumber, TString _ord){
  TString a(Form("roodataset_%s_%s_b%d_%s",_varname.Data(),_reg.Data(),_binnumber,_ord.Data()));
  return a;
};

TString template_production::get_roodset_name_single(TString _varname, TString _reg, int _binnumber){
  TString a(Form("roodataset_single_%s_%s_b%d",_varname.Data(),_reg.Data(),_binnumber));
  return a;
};

TString template_production::get_roovar_name(TString _varname, int i, int j, TString _ord){
      TString t=_varname;
      if (i==0) t.Append("_EB"); else t.Append("_EE");
      if (j==0) t.Append("_1"); else if (j==1) t.Append("_2"); else t.Append("_both");
      t.Append("_");
      t.Append(_ord);
      return t;
};

Int_t template_production::Choose_bin_invmass(float invmass){

  const float cuts[n_templates+1] = {80,120,160,9999};

  if (invmass<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << invmass << std::endl;
    return -999;
  }

  for (int i=0; i<n_templates; i++) if ((invmass>=cuts[i]) && (invmass<cuts[i+1])) return i;
  
  std::cout << "WARNING: called bin choice for out-of-range value " << invmass << std::endl;
  return -999;


};

#endif // #ifdef template_production_cxx


