//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 20 11:30:23 2012 by ROOT version 5.32/01
// from TTree Tree_standard_sel/Tree_standard_sel
// found on file: signal_input.root
//////////////////////////////////////////////////////////

#ifndef efficiency_measure_h
#define efficiency_measure_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "binsdef.h"

using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class efficiency_measure {
public :

   TH1F *w_eff_gg[3];
   TH1F *w_eff_sg[2];
   TH1F *w_tot_gg[3];
   TH1F *w_tot_sg[2];
   TH1F *w_passing_gg[3];
   TH1F *w_passing_sg[2];

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

   efficiency_measure(const char* filename, const char* outname);
   virtual ~efficiency_measure();
   //   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual void     LoopOne(TString diffvariable, TFile *outf);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


   TString outname;

   //   std::vector<TString> diffvariables_list;

   float getpuenergy(int reg, float eta);
   Int_t Choose_bin_eta(float eta, int region);

};

#endif

#ifdef efficiency_measure_cxx
efficiency_measure::efficiency_measure(const char* filename, const char* _outname) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  TTree *tree=NULL;

   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f || !f->IsOpen()) {
         f = new TFile(filename);
      }
      f->GetObject("Tree_standard_sel",tree);

   }
   Init(tree);

   outname = TString(_outname);

//   diffvariables_list.push_back(TString("invmass"));
//   diffvariables_list.push_back(TString("diphotonpt"));
//   diffvariables_list.push_back(TString("costhetastar"));
//   diffvariables_list.push_back(TString("dphi"));


}

efficiency_measure::~efficiency_measure()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t efficiency_measure::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t efficiency_measure::LoadTree(Long64_t entry)
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

void efficiency_measure::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
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

Bool_t efficiency_measure::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void efficiency_measure::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
//Int_t efficiency_measure::Cut(Long64_t entry)
//{
//// This function may be called from Loop.
//// returns  1 if entry is accepted.
//// returns -1 otherwise.
//   return 1;
//}

Int_t efficiency_measure::Choose_bin_eta(float eta, int region){

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



#endif // #ifdef efficiency_measure_cxx
