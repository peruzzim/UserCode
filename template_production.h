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
#include "RooRealVar.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TF1.h"
#include "RooBinning.h"
#include "TString.h"
#include "TH1.h"

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

   template_production(TTree *tree=0);
   virtual ~template_production();
   /* virtual Int_t    Cut(Long64_t entry); */
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void     WriteOutput(const char* filename, const TString dirname);

   void Setup(Bool_t _isdata, TString _mode, TString _differentialvariable);

   TString differentialvariable;

   TProfile** GetPUScaling(bool doEB, TString diffvar);

   TRandom3 *randomgen;

   static const int n_templates=n_bins;

   bool dosignal;

   TH1F *histo_pt[2];

   TH1F *template_signal[2][n_templates+1];
   TH1F *template_background[2][n_templates+1];

   TH1F *obs_hist_single[2][n_templates];
   TH2F *obs_hist[3][n_templates];

   bool pt_reweighting_initialized;
   bool do_pt_reweighting;

   Float_t pholead_outvar;
   Float_t photrail_outvar;

   Bool_t initialized;

   Bool_t isdata;

   TString mode;

   Int_t n_histobins;
   Float_t leftrange;
   Float_t rightrange;

   Bool_t dosignaltemplate;
   Bool_t dobackgroundtemplate;
   Bool_t dodistribution;

   Int_t Choose_bin_invmass(float invmass, int region);
   Int_t Choose_bin_pt(float pt, int region);
   Int_t Choose_bin_eta(float eta, int region);
   Int_t Choose_bin_sieie(float sieie, int region);

   float FindPtWeight(float pt, float eta);
   void Initialize_Pt_Reweighting(TString dset1, TString dset2, TString temp1, TString temp2);
   void SetNoPtReweighting();

   TH1F *histo_pt_reweighting[2];

};

#endif

#ifdef template_production_cxx
template_production::template_production(TTree *tree)
{

  TH1F::SetDefaultSumw2(kTRUE);

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
   dosignaltemplate=false;
   dobackgroundtemplate=false;
   dodistribution=false;

   pt_reweighting_initialized = 0;
   do_pt_reweighting = 0;

   n_histobins = 400;
   leftrange = 0.0;
   rightrange = 5.0;

}

void template_production::Setup(Bool_t _isdata, TString _mode, TString _differentialvariable){

  for (int i=0; i<5; i++) std::cout << "WARNING: NO PT REWEIGHTING FOR 2D HISTOS" << std::endl;

  isdata=_isdata;
  mode=_mode;
  differentialvariable=_differentialvariable;

  Init();

  if (mode=="standard") dodistribution=true;
  if (mode=="signal" || mode=="randomcone") dosignaltemplate=true;
  if (mode=="background" || mode=="impinging" || mode=="sieiesideband" || mode=="combisosideband") dobackgroundtemplate=true;
   
  if (mode=="sieiesideband") for (int i=0; i<5; i++) std::cout << "Warning: large sieie sideband (0.014/0.031)" << std::endl;

  randomgen = new TRandom3(0);

  for (int i=0; i<2; i++){
    TString name="histo_pt_";
    TString reg;
    if (i==0) reg="EB"; else if (i==1) reg="EE";
    name.Append(reg);
    histo_pt[i] = new TH1F(name.Data(),name.Data(),172,40,900);
  }

  
  // template_{signal,background}[EB,EE][n_templates]
  for (int i=0; i<2; i++)
    for (int j=0; j<n_templates+1; j++) {
      TString name_signal="template_signal";
      TString reg;
      if (i==0) reg="EB"; else if (i==1) reg="EE";
      TString t=Form("%s_%s_b%d",name_signal.Data(),reg.Data(),j);
      template_signal[i][j] = new TH1F(t.Data(),t.Data(),n_histobins,leftrange,rightrange);
      template_signal[i][j]->Sumw2();
    }
  for (int i=0; i<2; i++)
    for (int j=0; j<n_templates+1; j++) {
      TString name_background="template_background";
      TString reg;
      if (i==0) reg="EB"; else if (i==1) reg="EE";
      TString t=Form("%s_%s_b%d",name_background.Data(),reg.Data(),j);
      template_background[i][j] = new TH1F(t.Data(),t.Data(),n_histobins,leftrange,rightrange);
      template_background[i][j]->Sumw2();
    }
      

  // obs_hist{,_single}

  for (int i=0; i<2; i++)
    for (int j=0; j<n_templates; j++) {
      TString name_signal="obs_hist_single";
      TString reg;
      if (i==0) reg="EB"; else if (i==1) reg="EE";
      TString t=Form("%s_%s_b%d",name_signal.Data(),reg.Data(),j);
      obs_hist_single[i][j] = new TH1F(t.Data(),t.Data(),n_histobins,leftrange,rightrange);
      obs_hist_single[i][j]->Sumw2();
    }
  for (int i=0; i<3; i++)
    for (int j=0; j<n_templates; j++) {
      TString name_signal="obs_hist";
      TString reg;
      if (i==0) reg="EBEB"; else if (i==1) reg="EBEE"; else if (i==2) reg="EEEE"; else if (i==3) reg="EEEB";
      TString t=Form("%s_%s_b%d",name_signal.Data(),reg.Data(),j);
      obs_hist[i][j] = new TH2F(t.Data(),t.Data(),n_histobins,leftrange,rightrange,n_histobins,leftrange,rightrange);
      obs_hist[i][j]->Sumw2();
    }


  
  initialized=true;
  
};

void template_production::Initialize_Pt_Reweighting(TString dset1, TString dset2, TString temp1, TString temp2){

  TString regions[2];
  regions[0]="EB";
  regions[1]="EE";

  for (int i=0; i<2; i++){

    TString reg=regions[i];

    TString file1="out_";
    file1.Append(dset1);
    file1.Append("_");
    file1.Append(temp1);
    file1.Append(".root");

    TString file2="out_";
    file2.Append(dset2);
    file2.Append("_");
    file2.Append(temp2);
    file2.Append(".root");


    TFile *f1 = new TFile(file1.Data(),"read");
    TFile *f2 = new TFile(file2.Data(),"read");


    TString name1="";
    if (dset1=="data") name1.Append("data_Tree_"); else name1.Append("mc_Tree_");
    if (temp1=="bkg") name1.Append("background_template/");
    if (temp1=="sig") name1.Append("signal_template/");
    if (temp1=="rcone") name1.Append("randomcone_signal_template/");
    if (temp1=="impinging") name1.Append("impinging_track_template/");
    if (temp1=="sieiesideband") name1.Append("sieiesideband_sel/");
    if (temp1=="combisosideband") name1.Append("combisosideband_sel/");
    name1.Append("histo_pt_");
    name1.Append(reg);

    TString name2="";
    if (dset2=="data") name2.Append("data_Tree_"); else name2.Append("mc_Tree_");
    if (temp2=="bkg") name2.Append("background_template/");
    if (temp2=="sig") name2.Append("signal_template/");
    if (temp2=="rcone") name2.Append("randomcone_signal_template/");
    if (temp2=="impinging") name2.Append("impinging_track_template/");
    if (temp2=="sieiesideband") name2.Append("sieiesideband_sel/");
    if (temp2=="combisosideband") name2.Append("combisosideband_sel/");
    name2.Append("histo_pt_");
    name2.Append(reg);

    TH1F *h[2];
    f1->GetObject(name1,h[0]);
    f2->GetObject(name2,h[1]);
    assert(h[0]!=NULL);
    assert(h[1]!=NULL);

    h[0]->Print();
    h[1]->Print();

    TH1F *newhist = (TH1F*)(h[1]->Clone(Form("reweight_%s",reg.Data())));
    assert(newhist!=NULL);
    newhist->Print();

    newhist->Scale(1.0/newhist->Integral());
    h[0]->Scale(1.0/h[0]->Integral());

    newhist->Divide(h[0]);

    pt_reweighting_initialized = 1;
    do_pt_reweighting = 1;

    histo_pt_reweighting[i] = newhist;

  }

};

void template_production::SetNoPtReweighting(){
  pt_reweighting_initialized = 1;
  do_pt_reweighting = 0;
};

float template_production::FindPtWeight(float pt, float eta){

  if (!pt_reweighting_initialized){
    std::cout << "PT REWEIGHTING NOT INITIALIZED" << std::endl;
    return -999;
  }

  if (!do_pt_reweighting) return 1;

  TH1F *h;
  if (fabs(eta)<1.5) h=histo_pt_reweighting[0]; else h=histo_pt_reweighting[1];

  if (pt>h->GetXaxis()->GetXmax() || pt<h->GetXaxis()->GetXmin()) return 1;
  float res = h->GetBinContent(h->FindBin(pt));
  if (res==0) res=1;
  //  std::cout << res << std::endl;
  return res;
  
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

  if (dosignaltemplate || dobackgroundtemplate) {

    for (int i=0; i<2; i++) histo_pt[i]->Write();

    for (int i=0; i<2; i++) for (int l=0; l<n_templates; l++) template_signal[i][n_templates]->Add(template_signal[i][l]);
    for (int i=0; i<2; i++) for (int l=0; l<n_templates; l++) template_background[i][n_templates]->Add(template_background[i][l]);

    for (int i=0; i<2; i++) for (int l=0; l<n_templates+1; l++) template_signal[i][l]->Write();
    for (int i=0; i<2; i++) for (int l=0; l<n_templates+1; l++) template_background[i][l]->Write();

  }

  if (dodistribution) {
  for (int i=0; i<2; i++) for (int l=0; l<n_templates; l++) obs_hist_single[i][l]->Write();
  for (int i=0; i<3; i++) for (int l=0; l<n_templates; l++) obs_hist[i][l]->Write();
  }

  std::cout << "output written" << std::endl;

  out->Close();
};


Int_t template_production::Choose_bin_invmass(float invmass, int region){

  int index;

  float *cuts=NULL;

  if (region==0) {cuts=binsdef_diphoton_EBEB; index=n_templates_EBEB;}
  if (region==2) {cuts=binsdef_diphoton_EEEE; index=n_templates_EEEE;}
  if (region==3) {cuts=binsdef_diphoton_EBEE; index=n_templates_EBEE;}
  if (region==4) {cuts=binsdef_diphoton_EBEE; index=n_templates_EBEE;}
  if (region==1) {cuts=binsdef_diphoton_EBEE; index=n_templates_EBEE;}

  assert (cuts!=NULL);
  assert (index!=0);

  cuts[index]=9999;

  if (invmass<cuts[0]){
    std::cout << "WARNING: called bin choice for out-of-range value " << invmass << " cuts[0]= " << cuts[0] << std::endl;
    return -999;
  }

  for (int i=0; i<index; i++) if ((invmass>=cuts[i]) && (invmass<cuts[i+1])) return i;
  
  std::cout << "WARNING: called bin choice for out-of-range value " << invmass << std::endl;
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

#endif // #ifdef template_production_cxx


