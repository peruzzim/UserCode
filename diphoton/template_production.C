#ifndef template_production_cxx
#define template_production_cxx
#include "template_production.h"

#include <algorithm>

using namespace std;

bool do_scan_cone = false;
bool do_event_mixing = false; // 2events

void template_production::Loop(int maxevents)
{

  if (fChain == 0) return;
  if (!initialized){
    std::cout << "Not initialized! Call Setup() first." << std::endl;
    return;
  }


    TKDTreeID *kdtree[2];
    std::vector<float> match_pho_pt[2];
    std::vector<float> match_pho_eta[2];
    std::vector<float> match_evt_rho[2];

    TFile *matchingfile;
    TTree *matchingtree;
    UInt_t matchingtree_event_fileuuid;
    Int_t matchingtree_event_run;
    Int_t matchingtree_event_lumi;
    Int_t matchingtree_event_number;
    Int_t matchingtree_index_1event_sigsig_1[nclosestmore];
    Int_t matchingtree_index_1event_sigsig_2[nclosestmore];
    Int_t matchingtree_index_1event_sigbkg_1[nclosestmore];
    Int_t matchingtree_index_1event_sigbkg_2[nclosestmore];
    Int_t matchingtree_index_1event_bkgsig_1[nclosestmore];
    Int_t matchingtree_index_1event_bkgsig_2[nclosestmore];
    Int_t matchingtree_index_1event_bkgbkg_1[nclosestmore];
    Int_t matchingtree_index_1event_bkgbkg_2[nclosestmore];
    Int_t matchingtree_index_2events_sigsig_1[nclosestmore];
    Int_t matchingtree_index_2events_sigsig_2[nclosestmore];
    Int_t matchingtree_index_2events_sigbkg_1[nclosestmore];
    Int_t matchingtree_index_2events_sigbkg_2[nclosestmore];
    Int_t matchingtree_index_2events_bkgsig_1[nclosestmore];
    Int_t matchingtree_index_2events_bkgsig_2[nclosestmore];
    Int_t matchingtree_index_2events_bkgbkg_1[nclosestmore];
    Int_t matchingtree_index_2events_bkgbkg_2[nclosestmore];
    Float_t matchingtree_deta_1event_sigsig_1[nclosestmore];
    Float_t matchingtree_deta_1event_sigsig_2[nclosestmore];
    Float_t matchingtree_deta_1event_sigbkg_1[nclosestmore];
    Float_t matchingtree_deta_1event_sigbkg_2[nclosestmore];
    Float_t matchingtree_deta_1event_bkgsig_1[nclosestmore];
    Float_t matchingtree_deta_1event_bkgsig_2[nclosestmore];
    Float_t matchingtree_deta_1event_bkgbkg_1[nclosestmore];
    Float_t matchingtree_deta_1event_bkgbkg_2[nclosestmore];
    Float_t matchingtree_deta_2events_sigsig_1[nclosestmore];
    Float_t matchingtree_deta_2events_sigsig_2[nclosestmore];
    Float_t matchingtree_deta_2events_sigbkg_1[nclosestmore];
    Float_t matchingtree_deta_2events_sigbkg_2[nclosestmore];
    Float_t matchingtree_deta_2events_bkgsig_1[nclosestmore];
    Float_t matchingtree_deta_2events_bkgsig_2[nclosestmore];
    Float_t matchingtree_deta_2events_bkgbkg_1[nclosestmore];
    Float_t matchingtree_deta_2events_bkgbkg_2[nclosestmore];
    Float_t matchingtree_drho_1event_sigsig_1[nclosestmore];
    Float_t matchingtree_drho_1event_sigsig_2[nclosestmore];
    Float_t matchingtree_drho_1event_sigbkg_1[nclosestmore];
    Float_t matchingtree_drho_1event_sigbkg_2[nclosestmore];
    Float_t matchingtree_drho_1event_bkgsig_1[nclosestmore];
    Float_t matchingtree_drho_1event_bkgsig_2[nclosestmore];
    Float_t matchingtree_drho_1event_bkgbkg_1[nclosestmore];
    Float_t matchingtree_drho_1event_bkgbkg_2[nclosestmore];
    Float_t matchingtree_drho_2events_sigsig_1[nclosestmore];
    Float_t matchingtree_drho_2events_sigsig_2[nclosestmore];
    Float_t matchingtree_drho_2events_sigbkg_1[nclosestmore];
    Float_t matchingtree_drho_2events_sigbkg_2[nclosestmore];
    Float_t matchingtree_drho_2events_bkgsig_1[nclosestmore];
    Float_t matchingtree_drho_2events_bkgsig_2[nclosestmore];
    Float_t matchingtree_drho_2events_bkgbkg_1[nclosestmore];
    Float_t matchingtree_drho_2events_bkgbkg_2[nclosestmore];
    Float_t matchingtree_dpt_1event_sigsig_1[nclosestmore];
    Float_t matchingtree_dpt_1event_sigsig_2[nclosestmore];
    Float_t matchingtree_dpt_1event_sigbkg_1[nclosestmore];
    Float_t matchingtree_dpt_1event_sigbkg_2[nclosestmore];
    Float_t matchingtree_dpt_1event_bkgsig_1[nclosestmore];
    Float_t matchingtree_dpt_1event_bkgsig_2[nclosestmore];
    Float_t matchingtree_dpt_1event_bkgbkg_1[nclosestmore];
    Float_t matchingtree_dpt_1event_bkgbkg_2[nclosestmore];
    Float_t matchingtree_dpt_2events_sigsig_1[nclosestmore];
    Float_t matchingtree_dpt_2events_sigsig_2[nclosestmore];
    Float_t matchingtree_dpt_2events_sigbkg_1[nclosestmore];
    Float_t matchingtree_dpt_2events_sigbkg_2[nclosestmore];
    Float_t matchingtree_dpt_2events_bkgsig_1[nclosestmore];
    Float_t matchingtree_dpt_2events_bkgsig_2[nclosestmore];
    Float_t matchingtree_dpt_2events_bkgbkg_1[nclosestmore];
    Float_t matchingtree_dpt_2events_bkgbkg_2[nclosestmore];

    if (mode=="standard_domatching"){
    
      TTree *mytree[2];
      TFile *f = new TFile(inputfilename.Data());
      f->GetObject("Tree_1Drandomcone_template",mytree[0]);
      f->GetObject("Tree_2Drandomconesideband_template",mytree[1]);
      assert(mytree[0]); assert(mytree[1]);
    
      for (int i=0; i<2; i++){
      
	int mynentries = mytree[i]->GetEntriesFast();

	cout << mynentries << " entries found" << endl;

	float pho1_pt;
	float pho1_eta;
	float pho1_phi;
	float pho2_pt;
	float pho2_eta;
	float pho2_phi;
	float evt_rho;
	int pass12whoissiglike;
	TBranch *b_pho1_pt;
	TBranch *b_pho1_eta;
	TBranch *b_pho1_phi;
	TBranch *b_pho2_pt;
	TBranch *b_pho2_eta;
	TBranch *b_pho2_phi;
	TBranch *b_evt_rho;
	TBranch *b_pass12whoissiglike;
	Double_t *array_pho_pt = new Double_t[mynentries];
	Double_t *array_pho_eta = new Double_t[mynentries];
	Double_t *array_pho_phi = new Double_t[mynentries];
	Double_t *array_evt_rho = new Double_t[mynentries];
	Double_t *array_otherpho_eta = new Double_t[mynentries];
	mytree[i]->SetBranchStatus("*",0);
	mytree[i]->SetBranchStatus("pholead_pt",1);
	mytree[i]->SetBranchStatus("pholead_SCeta",1);
	mytree[i]->SetBranchStatus("pholead_SCphi",1);
	mytree[i]->SetBranchStatus("photrail_pt",1);
	mytree[i]->SetBranchStatus("photrail_SCeta",1);
	mytree[i]->SetBranchStatus("photrail_SCphi",1);
	mytree[i]->SetBranchStatus("event_rho",1);
	mytree[i]->SetBranchStatus("event_pass12whoissiglike");
	mytree[i]->SetBranchAddress("pholead_pt",&pho1_pt,&b_pho1_pt);
	mytree[i]->SetBranchAddress("pholead_SCeta",&pho1_eta,&b_pho1_eta);
	mytree[i]->SetBranchAddress("pholead_SCphi",&pho1_phi,&b_pho1_phi);
	mytree[i]->SetBranchAddress("photrail_pt",&pho2_pt,&b_pho2_pt);
	mytree[i]->SetBranchAddress("photrail_SCeta",&pho2_eta,&b_pho2_eta);
	mytree[i]->SetBranchAddress("photrail_SCphi",&pho2_phi,&b_pho2_phi);
	mytree[i]->SetBranchAddress("event_rho",&evt_rho,&b_evt_rho);
	mytree[i]->SetBranchAddress("event_pass12whoissiglike",&pass12whoissiglike, &b_pass12whoissiglike);

	for (int k=0; k<mynentries; k++){
	  mytree[i]->GetEntry(k);
	  if (i==0) pass12whoissiglike=1;
	  float pho_pt = (pass12whoissiglike==1) ? pho1_pt : pho2_pt;
	  float pho_eta = (pass12whoissiglike==1) ? pho1_eta : pho2_eta;
	  float pho_phi = (pass12whoissiglike==1) ? pho1_phi : pho2_phi;
	  const float pi = TMath::Pi();
	  pho_phi+=pi/2; // match preferentially to rotated pi/2
	  while (pho_phi > pi) pho_phi -= 2*pi;
	  while (pho_phi <= -pi) pho_phi += 2*pi;
	  float otherpho_eta = (pass12whoissiglike==0) ? pho1_eta : pho2_eta;
	  array_pho_eta[k] = (pho_eta+(fabs(pho_eta)>1.4442)*1000*pho_eta/fabs(pho_eta))/0.1; 
	  array_pho_phi[k] = pho_phi/0.4;
	  array_evt_rho[k] = evt_rho/2.0;
	  array_pho_pt[k] = TMath::Log(pho_pt)/0.2;
	  array_otherpho_eta[k] = (fabs(otherpho_eta)<1.4442) ? 0 : 1000;
	  match_pho_eta[i].push_back(pho_eta); 
	  match_evt_rho[i].push_back(evt_rho);
	  match_pho_pt[i].push_back(pho_pt);
	}

	cout << "arrays filled" << endl;

	kdtree[i] = new TKDTreeID(mynentries,(i==0) ? 2 : 5,1);
	kdtree[i]->SetData(0,array_pho_eta);
	kdtree[i]->SetData(1,array_evt_rho);
	if (i==1) kdtree[i]->SetData(2,array_pho_pt);
	if (i==1) kdtree[i]->SetData(3,array_otherpho_eta);
	if (i==1) kdtree[i]->SetData(4,array_pho_phi);
	cout << "SetData completed" << endl;
	kdtree[i]->Build();
	cout << "KD-tree built" << endl;
      
      }


    matchingfile = new TFile("matchingfile.root","recreate");
    matchingfile->cd();
    matchingtree = new TTree("matchingtree","matchingtree");
    matchingtree->Branch("matchingtree_event_fileuuid",&matchingtree_event_fileuuid,"matchingtree_event_fileuuid/i");
    matchingtree->Branch("matchingtree_event_run",&matchingtree_event_run,"matchingtree_event_run/I");
    matchingtree->Branch("matchingtree_event_lumi",&matchingtree_event_lumi,"matchingtree_event_lumi/I");
    matchingtree->Branch("matchingtree_event_number",&matchingtree_event_number,"matchingtree_event_number/I");
    matchingtree->Branch("matchingtree_index_1event_sigsig_1",&matchingtree_index_1event_sigsig_1,Form("matchingtree_index_1event_sigsig_1[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_1event_sigsig_2",&matchingtree_index_1event_sigsig_2,Form("matchingtree_index_1event_sigsig_2[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_1event_sigbkg_1",&matchingtree_index_1event_sigbkg_1,Form("matchingtree_index_1event_sigbkg_1[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_1event_sigbkg_2",&matchingtree_index_1event_sigbkg_2,Form("matchingtree_index_1event_sigbkg_2[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_1event_bkgsig_1",&matchingtree_index_1event_bkgsig_1,Form("matchingtree_index_1event_bkgsig_1[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_1event_bkgsig_2",&matchingtree_index_1event_bkgsig_2,Form("matchingtree_index_1event_bkgsig_2[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_1event_bkgbkg_1",&matchingtree_index_1event_bkgbkg_1,Form("matchingtree_index_1event_bkgbkg_1[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_1event_bkgbkg_2",&matchingtree_index_1event_bkgbkg_2,Form("matchingtree_index_1event_bkgbkg_2[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_2events_sigsig_1",&matchingtree_index_2events_sigsig_1,Form("matchingtree_index_2events_sigsig_1[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_2events_sigsig_2",&matchingtree_index_2events_sigsig_2,Form("matchingtree_index_2events_sigsig_2[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_2events_sigbkg_1",&matchingtree_index_2events_sigbkg_1,Form("matchingtree_index_2events_sigbkg_1[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_2events_sigbkg_2",&matchingtree_index_2events_sigbkg_2,Form("matchingtree_index_2events_sigbkg_2[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_2events_bkgsig_1",&matchingtree_index_2events_bkgsig_1,Form("matchingtree_index_2events_bkgsig_1[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_2events_bkgsig_2",&matchingtree_index_2events_bkgsig_2,Form("matchingtree_index_2events_bkgsig_2[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_2events_bkgbkg_1",&matchingtree_index_2events_bkgbkg_1,Form("matchingtree_index_2events_bkgbkg_1[%d]/I",nclosestmore));
    matchingtree->Branch("matchingtree_index_2events_bkgbkg_2",&matchingtree_index_2events_bkgbkg_2,Form("matchingtree_index_2events_bkgbkg_2[%d]/I",nclosestmore));

    matchingtree->Branch("matchingtree_deta_1event_sigsig_1",&matchingtree_deta_1event_sigsig_1,Form("matchingtree_deta_1event_sigsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_1event_sigsig_2",&matchingtree_deta_1event_sigsig_2,Form("matchingtree_deta_1event_sigsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_1event_sigbkg_1",&matchingtree_deta_1event_sigbkg_1,Form("matchingtree_deta_1event_sigbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_1event_sigbkg_2",&matchingtree_deta_1event_sigbkg_2,Form("matchingtree_deta_1event_sigbkg_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_1event_bkgsig_1",&matchingtree_deta_1event_bkgsig_1,Form("matchingtree_deta_1event_bkgsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_1event_bkgsig_2",&matchingtree_deta_1event_bkgsig_2,Form("matchingtree_deta_1event_bkgsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_1event_bkgbkg_1",&matchingtree_deta_1event_bkgbkg_1,Form("matchingtree_deta_1event_bkgbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_1event_bkgbkg_2",&matchingtree_deta_1event_bkgbkg_2,Form("matchingtree_deta_1event_bkgbkg_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_2events_sigsig_1",&matchingtree_deta_2events_sigsig_1,Form("matchingtree_deta_2events_sigsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_2events_sigsig_2",&matchingtree_deta_2events_sigsig_2,Form("matchingtree_deta_2events_sigsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_2events_sigbkg_1",&matchingtree_deta_2events_sigbkg_1,Form("matchingtree_deta_2events_sigbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_2events_sigbkg_2",&matchingtree_deta_2events_sigbkg_2,Form("matchingtree_deta_2events_sigbkg_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_2events_bkgsig_1",&matchingtree_deta_2events_bkgsig_1,Form("matchingtree_deta_2events_bkgsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_2events_bkgsig_2",&matchingtree_deta_2events_bkgsig_2,Form("matchingtree_deta_2events_bkgsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_2events_bkgbkg_1",&matchingtree_deta_2events_bkgbkg_1,Form("matchingtree_deta_2events_bkgbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_deta_2events_bkgbkg_2",&matchingtree_deta_2events_bkgbkg_2,Form("matchingtree_deta_2events_bkgbkg_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_1event_sigsig_1",&matchingtree_drho_1event_sigsig_1,Form("matchingtree_drho_1event_sigsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_1event_sigsig_2",&matchingtree_drho_1event_sigsig_2,Form("matchingtree_drho_1event_sigsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_1event_sigbkg_1",&matchingtree_drho_1event_sigbkg_1,Form("matchingtree_drho_1event_sigbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_1event_sigbkg_2",&matchingtree_drho_1event_sigbkg_2,Form("matchingtree_drho_1event_sigbkg_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_1event_bkgsig_1",&matchingtree_drho_1event_bkgsig_1,Form("matchingtree_drho_1event_bkgsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_1event_bkgsig_2",&matchingtree_drho_1event_bkgsig_2,Form("matchingtree_drho_1event_bkgsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_1event_bkgbkg_1",&matchingtree_drho_1event_bkgbkg_1,Form("matchingtree_drho_1event_bkgbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_1event_bkgbkg_2",&matchingtree_drho_1event_bkgbkg_2,Form("matchingtree_drho_1event_bkgbkg_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_2events_sigsig_1",&matchingtree_drho_2events_sigsig_1,Form("matchingtree_drho_2events_sigsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_2events_sigsig_2",&matchingtree_drho_2events_sigsig_2,Form("matchingtree_drho_2events_sigsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_2events_sigbkg_1",&matchingtree_drho_2events_sigbkg_1,Form("matchingtree_drho_2events_sigbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_2events_sigbkg_2",&matchingtree_drho_2events_sigbkg_2,Form("matchingtree_drho_2events_sigbkg_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_2events_bkgsig_1",&matchingtree_drho_2events_bkgsig_1,Form("matchingtree_drho_2events_bkgsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_2events_bkgsig_2",&matchingtree_drho_2events_bkgsig_2,Form("matchingtree_drho_2events_bkgsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_2events_bkgbkg_1",&matchingtree_drho_2events_bkgbkg_1,Form("matchingtree_drho_2events_bkgbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_drho_2events_bkgbkg_2",&matchingtree_drho_2events_bkgbkg_2,Form("matchingtree_drho_2events_bkgbkg_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_1event_sigsig_1",&matchingtree_dpt_1event_sigsig_1,Form("matchingtree_dpt_1event_sigsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_1event_sigsig_2",&matchingtree_dpt_1event_sigsig_2,Form("matchingtree_dpt_1event_sigsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_1event_sigbkg_1",&matchingtree_dpt_1event_sigbkg_1,Form("matchingtree_dpt_1event_sigbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_1event_sigbkg_2",&matchingtree_dpt_1event_sigbkg_2,Form("matchingtree_dpt_1event_sigbkg_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_1event_bkgsig_1",&matchingtree_dpt_1event_bkgsig_1,Form("matchingtree_dpt_1event_bkgsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_1event_bkgsig_2",&matchingtree_dpt_1event_bkgsig_2,Form("matchingtree_dpt_1event_bkgsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_1event_bkgbkg_1",&matchingtree_dpt_1event_bkgbkg_1,Form("matchingtree_dpt_1event_bkgbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_1event_bkgbkg_2",&matchingtree_dpt_1event_bkgbkg_2,Form("matchingtree_dpt_1event_bkgbkg_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_2events_sigsig_1",&matchingtree_dpt_2events_sigsig_1,Form("matchingtree_dpt_2events_sigsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_2events_sigsig_2",&matchingtree_dpt_2events_sigsig_2,Form("matchingtree_dpt_2events_sigsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_2events_sigbkg_1",&matchingtree_dpt_2events_sigbkg_1,Form("matchingtree_dpt_2events_sigbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_2events_sigbkg_2",&matchingtree_dpt_2events_sigbkg_2,Form("matchingtree_dpt_2events_sigbkg_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_2events_bkgsig_1",&matchingtree_dpt_2events_bkgsig_1,Form("matchingtree_dpt_2events_bkgsig_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_2events_bkgsig_2",&matchingtree_dpt_2events_bkgsig_2,Form("matchingtree_dpt_2events_bkgsig_2[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_2events_bkgbkg_1",&matchingtree_dpt_2events_bkgbkg_1,Form("matchingtree_dpt_2events_bkgbkg_1[%d]/F",nclosestmore));
    matchingtree->Branch("matchingtree_dpt_2events_bkgbkg_2",&matchingtree_dpt_2events_bkgbkg_2,Form("matchingtree_dpt_2events_bkgbkg_2[%d]/F",nclosestmore));

    }


  Long64_t nentries = fChain->GetEntriesFast();
  int limit_entries = maxevents;
  //int limit_entries = 1e+4;
  //int limit_entries = -1;
  const float thr = float(limit_entries)/float(nentries);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if (jentry%100000==0) std::cout << "Processing entry " << jentry << std::endl;

    if (limit_entries>0){
      if (randomgen->Uniform(0,1) > thr) continue;
    }

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //    if (event_nRecVtx>5) continue;

    //    if (jentry==1e+4) break;

    // initial kinematic selection
    //    if (dodistribution) if (pholead_pt<40 || photrail_pt<30 || dipho_mgg_photon<80) continue;
    if (dodistribution) if (pholead_pt<40 || photrail_pt<25) continue;
    if (do2dtemplate) if (pholead_pt<40 || photrail_pt<25) continue;
    float cutpt = (mode=="muon") ? 10 : 25;
    if (dosignaltemplate || dobackgroundtemplate) if (pholead_pt<cutpt) continue;

    if (mode=="zee") {
      if (fabs(dipho_mgg_photon-91.2)>10) continue;
    }

    if (mode=="standard_pixelrev"){
      if (dipho_mgg_photon<60 || dipho_mgg_photon>150) continue;
    }

    //    if (mode=="background") if (pholead_PhoMCmatchexitcode!=3) continue;

    //    if (mode=="background") if (pholead_PhoMCmatchexitcode!=1 && pholead_PhoMCmatchexitcode!=2) continue;

    //    if (mode=="background") if (pholead_PhoMCmatchexitcode!=0) continue;  

    //    if (isdata && event_CSCTightHaloID>0) continue;
    //if (mode!="muon" && event_NMuons>0) continue;
    //if (mode!="muon" && event_NMuonsTot>0) continue;
    //    if (!isdata && (event_PUOOTnumInteractionsEarly<20 && event_PUOOTnumInteractionsLate<20)) continue;
    //    if (!isdata && event_PUOOTnumInteractionsEarly>3) continue;

    if (mode=="fragmentation") if (pholead_PhoMCmatchexitcode!=1) continue;
    if (mode=="nofragmentation") if (pholead_PhoMCmatchexitcode!=2) continue;

    if (mode=="cutPFchargediso_signal" || mode=="cutPFchargediso_background" || mode=="cutPFchargediso_randomcone" || mode=="cutPFchargediso_sieiesideband") if (pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01>0.1) continue;





    Int_t event_ok_for_dataset=-1;

    Int_t reg_lead;
    Int_t reg_trail;

    // dataset 0:EBEB 3/4->1:EBEE 2:EEEE

    if (fabs(pholead_SCeta)<1.4442 && fabs(photrail_SCeta)<1.4442) {
      event_ok_for_dataset=0;
      reg_lead=0;
      reg_trail=0;
    }
    else if (fabs(pholead_SCeta)>1.56 && fabs(photrail_SCeta)>1.56) {
      event_ok_for_dataset=2;
      reg_lead=1;
      reg_trail=1;
    }
    else if (fabs(pholead_SCeta)<1.4442 && fabs(photrail_SCeta)>1.56) {
      event_ok_for_dataset=3;
      reg_lead=0;
      reg_trail=1;
    }
    else if (fabs(pholead_SCeta)>1.56 && fabs(photrail_SCeta)<1.4442) {
      event_ok_for_dataset=4;
      reg_lead=1;
      reg_trail=0;
    }
    else std::cout << "We have a problem here!!!" << std::endl;

//    if (differentialvariable=="photoniso"){
//      pholead_outvar=pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;
//      photrail_outvar=photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;
//      candcounter=pholead_Npfcandphotonincone;
//    }
//    else if (differentialvariable=="combiso"){
//      pholead_outvar=pholead_pho_Cone04PFCombinedIso;
//      photrail_outvar=photrail_pho_Cone04PFCombinedIso;
//      candcounter=pholead_Npfcandphotonincone+pholead_Npfcandchargedincone+pholead_Npfcandneutralincone;
//    }
//    else if (differentialvariable=="chargediso"){
//      pholead_outvar=pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
//      photrail_outvar=photrail_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
//      candcounter=pholead_Npfcandchargedincone;
//    }
//    else if (differentialvariable=="neutraliso"){
//      pholead_outvar=pholead_pho_Cone04NeutralHadronIso_mvVtx;
//      photrail_outvar=photrail_pho_Cone04NeutralHadronIso_mvVtx;
//      candcounter=pholead_Npfcandneutralincone;
//    }


    pholead_outvar = -999;
    photrail_outvar = -999;
    candcounter = -999;

    bool recalc_lead = false;
    bool recalc_trail = false;

    if (dosignaltemplate||dobackgroundtemplate) recalc_lead=true;
    if (do2ptemplate || do1p1ftemplate || do2ftemplate) {recalc_lead=true; recalc_trail=true;} 
    if (dodistribution) {recalc_lead=true; recalc_trail=true;}



    if (recalc_lead){
      float et_recalc = 0;
      float e_recalc = 0;
      int number_recalc = 0;

      bool printout=false;

      if (printout) std::cout << "---" << std::endl;

      for (int i=0; i<pholead_Npfcandphotonincone; i++) {

	float et=pholead_photonpfcandets[i];
	float e=pholead_photonpfcandenergies[i];
	float deta=pholead_photonpfcanddetas[i];
	float dphi=pholead_photonpfcanddphis[i];
	float dR=sqrt(deta*deta+dphi*dphi);
	float eta=fabs(TMath::ACosH(e/et));
	if (printout) {
	  std::cout << et << " " << e << " " << deta << " " << dphi << " " << dR << " " << eta << std::endl;
	}
	if (eta>1.4442 && eta<1.56) continue;
	if (eta>2.5) continue;

#include "cleaning.cc"

	if (fabs(pholead_SCeta)<1.4442 && eta>1.4442) continue;
	if (fabs(pholead_SCeta)>1.56 && eta<1.56) continue;
	et_recalc+=et;
	e_recalc+=e;
	number_recalc++;

      }

      pholead_outvar=et_recalc;

      if (printout) std::cout << "---" << std::endl;

    }

    if (recalc_trail){
      float et_recalc = 0;
      float e_recalc = 0;
      int number_recalc = 0;

      bool printout=false;

      if (printout) std::cout << "---" << std::endl;

      for (int i=0; i<photrail_Npfcandphotonincone; i++) {

	float et=photrail_photonpfcandets[i];
	float e=photrail_photonpfcandenergies[i];
	float deta=photrail_photonpfcanddetas[i];
	float dphi=photrail_photonpfcanddphis[i];
	float dR=sqrt(deta*deta+dphi*dphi);
	float eta=fabs(TMath::ACosH(e/et));
	if (printout) {
	  std::cout << et << " " << e << " " << deta << " " << dphi << " " << dR << " " << eta << std::endl;
	}
	if (eta>1.4442 && eta<1.56) continue;
	if (eta>2.5) continue;

#include "cleaning.cc"

	if (fabs(photrail_SCeta)<1.4442 && eta>1.4442) continue;
	if (fabs(photrail_SCeta)>1.56 && eta<1.56) continue;
	et_recalc+=et;
	e_recalc+=e;
	number_recalc++;

      }

      photrail_outvar=et_recalc;

      if (printout) std::cout << "---" << std::endl;

    }

//    if (fabs(pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx-pholead_outvar)>0) cout << pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx << " " << pholead_outvar << endl;
//    if (fabs(photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx-photrail_outvar)>0) cout << photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx << " " << photrail_outvar << endl;

    roorho->setVal(event_rho);
    roosigma->setVal(event_sigma);

    pholead_outvar-=getpuenergy(reg_lead,pholead_SCeta);
    if (dodistribution) photrail_outvar-=getpuenergy(reg_trail,photrail_SCeta);
    if (do2ptemplate || do1p1ftemplate || do2ftemplate) photrail_outvar-=getpuenergy(reg_trail,photrail_SCeta);
    

    if (recalc_lead){
      if (pholead_outvar<-100) std::cout << "PROBLEM WITH ISOLATION CALCULATION!!!" << std::endl;
      assert (pholead_outvar>=-100);
      //      if (pholead_outvar<leftrange) {/*std::cout << "Warning: fixing underflow " << pholead_outvar << std::endl;*/ pholead_outvar=leftrange+1e-5;}
      if (pholead_outvar<=leftrange) continue;
      if (pholead_outvar>=rightrange) continue;
      //      if (pholead_outvar>=rightrange) pholead_outvar=rightrange-1e-5; // overflow in last bin 
    }
    if (recalc_trail){
      if (photrail_outvar<-100) std::cout << "PROBLEM WITH ISOLATION CALCULATION!!!" << std::endl;
      assert (photrail_outvar>=-100);
      //      if (photrail_outvar<leftrange) {/*std::cout << "Warning: fixing underflow " << photrail_outvar << std::endl;*/ photrail_outvar=leftrange+1e-5;}
      if (photrail_outvar<=leftrange) continue;
      if (photrail_outvar>=rightrange) continue;
      //      if (photrail_outvar>=rightrange) photrail_outvar=rightrange-1e-5; // overflow in last bin 
    }

    Float_t weight=event_luminormfactor*event_Kfactor*event_weight;

    if (mode=="standard_2frag" || mode=="2pgen_2frag" || mode=="1p1fbothgen_2frag" || mode=="1pgen1fside_2frag") {
      if (pholead_PhoMCmatchexitcode==1 && pholead_GenPhotonIsoDR04<5) weight*=2;
      if (photrail_PhoMCmatchexitcode==1 && photrail_GenPhotonIsoDR04<5) weight*=2;
    }
    if (mode=="signal_2frag"){
      if (pholead_PhoMCmatchexitcode==1 && pholead_GenPhotonIsoDR04<5) weight*=2;
    }

    if (dosignaltemplate||dobackgroundtemplate){

      if (dosignaltemplate){
	roovar1->setVal(pholead_outvar);
	roovar2->setVal(pholead_outvar);
	roopt1->setVal(pholead_pt);
	roopt2->setVal(pholead_pt);
	roosieie1->setVal(pholead_sieie);
	roosieie2->setVal(pholead_sieie);
	rooeta1->setVal(fabs(pholead_SCeta));
	rooeta2->setVal(fabs(pholead_SCeta));
	rooweight->setVal(weight);
	roodset_signal[reg_lead][0]->add(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma),weight);
	roodset_signal[reg_lead][1]->add(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma),weight);
	if (do_scan_cone) for (int k=0; k<50; k++) scan_cone_histos[reg_lead][k]->Fill(pholead_test_rotatedphotoniso[k]-getpuenergy(reg_lead,pholead_SCeta),weight);
	if (do_scan_cone) for (int k=0; k<50; k++) if (pholead_test_rotatedwithcheckphotoniso[k]>-100) scan_conewithcheck_histos[reg_lead][k]->Fill(pholead_test_rotatedwithcheckphotoniso[k]-getpuenergy(reg_lead,pholead_SCeta),weight);
      }
      
      if (dobackgroundtemplate){
	  roovar1->setVal(pholead_outvar);
	  roovar2->setVal(pholead_outvar); 
	  roopt1->setVal(pholead_pt);
	  roopt2->setVal(pholead_pt);
	  roosieie1->setVal(pholead_sieie);
	  roosieie2->setVal(pholead_sieie);
	  rooeta1->setVal(fabs(pholead_SCeta));
	  rooeta2->setVal(fabs(pholead_SCeta));
	  rooweight->setVal(weight);
	  roodset_background[reg_lead][0]->add(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma),weight);
	  roodset_background[reg_lead][1]->add(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma),weight);
	  if (do_scan_cone) for (int k=0; k<50; k++) scan_cone_histos[reg_lead][k]->Fill(pholead_test_rotatedphotoniso[k]-getpuenergy(reg_lead,pholead_SCeta),weight);
	  if (do_scan_cone) for (int k=0; k<50; k++) if (pholead_test_rotatedwithcheckphotoniso[k]>-100) scan_conewithcheck_histos[reg_lead][k]->Fill(pholead_test_rotatedwithcheckphotoniso[k]-getpuenergy(reg_lead,pholead_SCeta),weight);
      }
      
    }

    if (do2dtemplate){

	int event_ok_for_dataset_local = event_ok_for_dataset;

	float in1=pholead_outvar;
	float in2=photrail_outvar;
	float ptin1=pholead_pt;
	float ptin2=photrail_pt;
	float sieiein1=pholead_sieie;
	float sieiein2=photrail_sieie;	
	float etain1=fabs(pholead_SCeta);
	float etain2=fabs(photrail_SCeta);

	bool doswap = false;
	if ((event_ok_for_dataset_local==0 || event_ok_for_dataset_local==2) && (randomgen->Uniform()>0.5)) doswap=true;
	if (event_ok_for_dataset_local==4) doswap=true;
	if (event_ok_for_dataset_local==3 || event_ok_for_dataset_local==4) event_ok_for_dataset_local=1;

	if (doswap){
	  float temp;
	  temp=in1; in1=in2; in2=temp;
	  temp=ptin1; ptin1=ptin2; ptin2=temp;
	  temp=sieiein1; sieiein1=sieiein2; sieiein2=temp;
	  temp=etain1; etain1=etain2; etain2=temp;
	  event_pass12whoissiglike=!event_pass12whoissiglike;
	}

	TString sigorbkg;
	if (do2ptemplate) sigorbkg=TString("sigsig");
	if (do2ftemplate) sigorbkg=TString("bkgbkg");
	if (do1p1ftemplate && event_pass12whoissiglike==0) sigorbkg=TString("sigbkg");
	if (do1p1ftemplate && event_pass12whoissiglike==1) sigorbkg=TString("bkgsig");
	
	roovar1->setVal(in1);
	roovar2->setVal(in2);
	roopt1->setVal(ptin1);
	roopt2->setVal(ptin2);
	roosieie1->setVal(sieiein1);
	roosieie2->setVal(sieiein2);
	rooeta1->setVal(etain1);
	rooeta2->setVal(etain2);
	rooweight->setVal(weight);
	RooArgSet args(*roovar1,*roovar2,*roopt1,*roopt2,*roosieie1,*roosieie2,*rooeta1,*rooeta2);
	args.add(RooArgSet(*roorho,*roosigma));
	FillDiffVariables(); // WARNING: WORKS ONLY IF DIFF VARIABLES ARE NOT SENSITIVE TO SWAPPING 1 WITH 2
	template2d_roodset[get_name_template2d_roodset(event_ok_for_dataset_local,sigorbkg)]->add(args,weight);

    } // end if do2dtemplate

    if (dodistribution && event_ok_for_dataset>-1){

	 int event_ok_for_dataset_local = event_ok_for_dataset;

	 FillDiffVariables(); // WARNING: WORKS ONLY IF DIFF VARIABLES ARE NOT SENSITIVE TO SWAPPING 1 WITH 2

	 float in1=pholead_outvar;
	 float in2=photrail_outvar;
	 float ptin1=pholead_pt;
	 float ptin2=photrail_pt;
	 float sieiein1=pholead_sieie;
	 float sieiein2=photrail_sieie;
	 float etain1=fabs(pholead_SCeta);
	 float etain2=fabs(photrail_SCeta);

	 bool doswap=false;

	 if ((event_ok_for_dataset_local==0 || event_ok_for_dataset_local==2) && (randomgen->Uniform()>0.5)) doswap=true;

	 if (event_ok_for_dataset_local==4) doswap=true;

	 if (event_ok_for_dataset_local==3 || event_ok_for_dataset_local==4) event_ok_for_dataset_local=1;

	 if (doswap){
	   float temp;
	   temp=in1; in1=in2; in2=temp;
	   temp=ptin1; ptin1=ptin2; ptin2=temp;
	   temp=sieiein1; sieiein1=sieiein2; sieiein2=temp;
	   temp=etain1; etain1=etain2; etain2=temp;
	 }

	 roovar1->setVal(in1);
	 roovar2->setVal(in2);
	 roopt1->setVal(ptin1);
	 roopt2->setVal(ptin2);
	 roosieie1->setVal(sieiein1);
	 roosieie2->setVal(sieiein2);
	 rooeta1->setVal(etain1);
	 rooeta2->setVal(etain2);
	 rooweight->setVal(weight);
	 RooArgSet args(*roovar1,*roovar2,*roopt1,*roopt2,*roosieie1,*roosieie2,*rooeta1,*rooeta2);
	 args.add(RooArgSet(*roorho,*roosigma));
	 args.add(*rooargset_diffvariables);
	 RooArgSet args2(*roovar1,*roovar2,*roopt1,*roopt2,*roosieie1,*roosieie2,*rooeta1,*rooeta2);
	 args2.add(RooArgSet(*roorho,*roosigma));

       for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){

	 Int_t bin_couple = -999;
	 float value_diffvariable;

	if (*diffvariable==TString("invmass")) {
	  value_diffvariable=roovar_invmass->getVal();
	  bin_couple = Choose_bin_invmass(value_diffvariable,event_ok_for_dataset_local);
	  //	  invmass_vector.push_back(value_diffvariable);
	}
	if (*diffvariable==TString("diphotonpt")){
	  value_diffvariable=roovar_diphotonpt->getVal();
	  bin_couple = Choose_bin_diphotonpt(value_diffvariable,event_ok_for_dataset_local);
	  //	  diphotonpt_vector.push_back(value_diffvariable);
	}
	if (*diffvariable==TString("costhetastar")){
	  value_diffvariable = roovar_costhetastar->getVal();
	  bin_couple = Choose_bin_costhetastar(value_diffvariable,event_ok_for_dataset_local); 
	}
	if (*diffvariable==TString("dphi")){
	  value_diffvariable=roovar_dphi->getVal();
	  bin_couple = Choose_bin_dphi(value_diffvariable,event_ok_for_dataset_local);
	}
	if (*diffvariable==TString("dR")){
	  value_diffvariable=roovar_dR->getVal();
	  bin_couple = Choose_bin_dR(value_diffvariable,event_ok_for_dataset_local);
	}
      
	if (bin_couple<0) continue;
	
	obs_roodset[get_name_obs_roodset(event_ok_for_dataset_local,*diffvariable,bin_couple)]->add(args,weight);

	if (donewtemplates) {
	  Float_t phoiso_1[2][2][nclosest];
	  Float_t phoiso_2[2][2][nclosest];
	  Float_t rewinfo_1[2][2][nclosest*6];
	  Float_t rewinfo_2[2][2][nclosest*6];
	  for (int l=0; l<nclosest; l++){

	    phoiso_1[0][0][l]=(do_event_mixing ? phoiso_template_2events_sigsig_1[l] : phoiso_template_1event_sigsig_1[l]);
	    phoiso_1[0][1][l]=(do_event_mixing ? phoiso_template_2events_sigbkg_1[l] : phoiso_template_1event_sigbkg_1[l]);
	    phoiso_1[1][0][l]=(do_event_mixing ? phoiso_template_2events_bkgsig_1[l] : phoiso_template_1event_bkgsig_1[l]);
	    phoiso_1[1][1][l]=(do_event_mixing ? phoiso_template_2events_bkgbkg_1[l] : phoiso_template_1event_bkgbkg_1[l]);
	    phoiso_2[0][0][l]=(do_event_mixing ? phoiso_template_2events_sigsig_2[l] : phoiso_template_1event_sigsig_2[l]);
	    phoiso_2[0][1][l]=(do_event_mixing ? phoiso_template_2events_sigbkg_2[l] : phoiso_template_1event_sigbkg_2[l]);
	    phoiso_2[1][0][l]=(do_event_mixing ? phoiso_template_2events_bkgsig_2[l] : phoiso_template_1event_bkgsig_2[l]);
	    phoiso_2[1][1][l]=(do_event_mixing ? phoiso_template_2events_bkgbkg_2[l] : phoiso_template_1event_bkgbkg_2[l]);


	    memcpy(&(rewinfo_1[0][0][l*6]),(do_event_mixing ? &(rewinfo_template_2events_sigsig_1[l*6]) : &(rewinfo_template_1event_sigsig_1[l*6])),6*sizeof(float));
	    memcpy(&(rewinfo_1[0][1][l*6]),(do_event_mixing ? &(rewinfo_template_2events_sigbkg_1[l*6]) : &(rewinfo_template_1event_sigbkg_1[l*6])),6*sizeof(float));
	    memcpy(&(rewinfo_1[1][0][l*6]),(do_event_mixing ? &(rewinfo_template_2events_bkgsig_1[l*6]) : &(rewinfo_template_1event_bkgsig_1[l*6])),6*sizeof(float));
	    memcpy(&(rewinfo_1[1][1][l*6]),(do_event_mixing ? &(rewinfo_template_2events_bkgbkg_1[l*6]) : &(rewinfo_template_1event_bkgbkg_1[l*6])),6*sizeof(float));
	    memcpy(&(rewinfo_2[0][0][l*6]),(do_event_mixing ? &(rewinfo_template_2events_sigsig_2[l*6]) : &(rewinfo_template_1event_sigsig_2[l*6])),6*sizeof(float));
	    memcpy(&(rewinfo_2[0][1][l*6]),(do_event_mixing ? &(rewinfo_template_2events_sigbkg_2[l*6]) : &(rewinfo_template_1event_sigbkg_2[l*6])),6*sizeof(float));
	    memcpy(&(rewinfo_2[1][0][l*6]),(do_event_mixing ? &(rewinfo_template_2events_bkgsig_2[l*6]) : &(rewinfo_template_1event_bkgsig_2[l*6])),6*sizeof(float));
	    memcpy(&(rewinfo_2[1][1][l*6]),(do_event_mixing ? &(rewinfo_template_2events_bkgbkg_2[l*6]) : &(rewinfo_template_1event_bkgbkg_2[l*6])),6*sizeof(float));

	  }

	  if (do_event_mixing) {cout << "WRONG!!!" << endl; assert(1==0);} // da implementare per il caso 2events

	  for (int n1=0; n1<2; n1++) for (int n2=0; n2<2; n2++) for (int l=0; l<nclosest; l++){
	    if (whichnewtemplate==0 && (n1!=0 || n2!=0)) continue;
	    if (whichnewtemplate==1 && (n1+n2!=1)) continue;
	    if (whichnewtemplate==2 && (n1!=1 || n2!=1)) continue;
	    float fill1=-999;
	    float fill2=-999;
	    float filleta1=-999;
	    float filleta2=-999;
	    float fillpt1=-999;
	    float fillpt2=-999;
	    float fillrho=-999;
	    float fillsigma=-999;

	    if (whichnewtemplate==0){
	      if (!doswap){
		fill1=(phoiso_1[n1][n2][l]);
		fill2=(phoiso_2[n1][n2][l]);
		filleta1=rewinfo_1[n1][n2][l*6+0];
                filleta2=rewinfo_1[n1][n2][l*6+0]-pholead_SCeta+photrail_SCeta;
		fillpt1=rewinfo_1[n1][n2][l*6+2];
                fillpt2=rewinfo_1[n1][n2][l*6+3];
		fillrho=rewinfo_1[n1][n2][l*6+4];
		fillsigma=rewinfo_1[n1][n2][l*6+5];
	      }
	      else {
		fill1=(phoiso_2[n1][n2][l]);
                fill2=(phoiso_1[n1][n2][l]);
		filleta1=rewinfo_1[n1][n2][l*6+0]-pholead_SCeta+photrail_SCeta;
		filleta2=rewinfo_1[n1][n2][l*6+0];
		fillpt1=rewinfo_1[n1][n2][l*6+3];
		fillpt2=rewinfo_1[n1][n2][l*6+2];
		fillrho=rewinfo_1[n1][n2][l*6+4];
		fillsigma=rewinfo_1[n1][n2][l*6+5];
	      }
	    }
	    else if (n1==0 && n2==1){
	      if (!doswap){
		fill1=(phoiso_1[n1][n2][l]);
		fill2=(phoiso_2[n1][n2][l]);
		filleta1=rewinfo_2[n1][n2][l*6+0]-photrail_SCeta+pholead_SCeta;
		filleta2=rewinfo_2[n1][n2][l*6+0];
		fillpt1=rewinfo_2[n1][n2][l*6+3];
		fillpt2=rewinfo_2[n1][n2][l*6+2];
		fillrho=rewinfo_2[n1][n2][l*6+4];
		fillsigma=rewinfo_2[n1][n2][l*6+5];
	      }
	      else {
		fill1=(phoiso_2[!n1][!n2][l]);
                fill2=(phoiso_1[!n1][!n2][l]);
		filleta1=rewinfo_2[n1][n2][l*6+0];
		filleta2=rewinfo_2[n1][n2][l*6+0]-photrail_SCeta+pholead_SCeta;
		fillpt1=rewinfo_2[n1][n2][l*6+2];
		fillpt2=rewinfo_2[n1][n2][l*6+3];
		fillrho=rewinfo_2[n1][n2][l*6+4];
		fillsigma=rewinfo_2[n1][n2][l*6+5];
	      }
	    }
	    else if (n1==1 && n2==0){
	      if (!doswap){
		fill1=(phoiso_1[n1][n2][l]);
		fill2=(phoiso_2[n1][n2][l]);
		filleta1=rewinfo_1[n1][n2][l*6+0];
		filleta2=rewinfo_1[n1][n2][l*6+0]-pholead_SCeta+photrail_SCeta;
		fillpt1=rewinfo_1[n1][n2][l*6+2];
		fillpt2=rewinfo_1[n1][n2][l*6+3];
		fillrho=rewinfo_1[n1][n2][l*6+4];
		fillsigma=rewinfo_1[n1][n2][l*6+5];
	      }
	      else {
		fill1=(phoiso_2[!n1][!n2][l]);
                fill2=(phoiso_1[!n1][!n2][l]);
		filleta1=rewinfo_1[n1][n2][l*6+0]-pholead_SCeta+photrail_SCeta;
		filleta2=rewinfo_1[n1][n2][l*6+0];
		fillpt1=rewinfo_1[n1][n2][l*6+3];
		fillpt2=rewinfo_1[n1][n2][l*6+2];
		fillrho=rewinfo_1[n1][n2][l*6+4];
		fillsigma=rewinfo_1[n1][n2][l*6+5];
	      }
	    }
	    else if (whichnewtemplate==2){
	      if (!doswap){
		fill1=(phoiso_1[n1][n2][l]);
		fill2=(phoiso_2[n1][n2][l]);
		filleta1=rewinfo_1[n1][n2][l*6+0];
		filleta2=rewinfo_2[n1][n2][l*6+0]; // not perfect, e' solo uno dei 2
		fillpt1=rewinfo_1[n1][n2][l*6+2];
		fillpt2=rewinfo_2[n1][n2][l*6+2];
		fillrho=(rewinfo_1[n1][n2][l*6+4]+rewinfo_2[n1][n2][l*6+4])/2;
		fillsigma=sqrt(pow(rewinfo_1[n1][n2][l*6+5],2)+pow(rewinfo_2[n1][n2][l*6+5],2))/sqrt(2);
	      }
	      else {
		fill1=(phoiso_2[n1][n2][l]);
                fill2=(phoiso_1[n1][n2][l]);
		filleta1=rewinfo_2[n1][n2][l*6+0];
		filleta2=rewinfo_1[n1][n2][l*6+0]; // not perfect, e' solo uno dei 2
		fillpt1=rewinfo_2[n1][n2][l*6+2];
		fillpt2=rewinfo_1[n1][n2][l*6+2];
		fillrho=(rewinfo_1[n1][n2][l*6+4]+rewinfo_2[n1][n2][l*6+4])/2;
		fillsigma=sqrt(pow(rewinfo_1[n1][n2][l*6+5],2)+pow(rewinfo_2[n1][n2][l*6+5],2))/sqrt(2);
	      }
	    }

	    if ((fabs(filleta1)>2.5) || (fabs(filleta1)>1.4442 && fabs(filleta1)<1.56)) continue;
	    if ((fabs(filleta2)>2.5) || (fabs(filleta2)>1.4442 && fabs(filleta2)<1.56)) continue;

	    fill1-=fillrho*geteffarea((fabs(filleta1)>1.4442),fabs(filleta1));
	    fill2-=fillrho*geteffarea((fabs(filleta2)>1.4442),fabs(filleta2));

//	    cout << n1 << n2 << l << " ";
//	    cout << pholead_SCeta << " " << (!doswap ? filleta1 : filleta2) << " " << (!doswap ? fillpt1 : fillpt2) << " "; 
//	    cout << photrail_SCeta << " " << (!doswap ? filleta2 : filleta1) << " " << (!doswap ? fillpt2 : fillpt1) << " "; 
//	    cout << endl;


	    if (fill1<=leftrange || fill2<=leftrange) continue;
	    if (fill1>=rightrange || fill2>=rightrange) continue;
	    roovar1->setVal(fill1);
	    roovar2->setVal(fill2);
	    rooeta1->setVal(fabs(filleta1));
	    rooeta2->setVal(fabs(filleta2));
	    roopt1->setVal(fillpt1);
	    roopt2->setVal(fillpt2);
	    roorho->setVal(fillrho);
	    roosigma->setVal(fillsigma);
	    if (n1==0 && n2==0) newtempl_roodset[get_name_newtempl_roodset(event_ok_for_dataset_local,*diffvariable,bin_couple,"sigsig")]->add(args2,weight);
	    else if (n1==0 && n2==1) newtempl_roodset[get_name_newtempl_roodset(event_ok_for_dataset_local,*diffvariable,bin_couple,"sigbkg")]->add(args2,weight);
	    else if (n1==1 && n2==0) newtempl_roodset[get_name_newtempl_roodset(event_ok_for_dataset_local,*diffvariable,bin_couple,"bkgsig")]->add(args2,weight);
	    else if (n1==1 && n2==1) newtempl_roodset[get_name_newtempl_roodset(event_ok_for_dataset_local,*diffvariable,bin_couple,"bkgbkg")]->add(args2,weight);
	  }
	}

	if (!isdata){
	  int isppevent = 0;
	  if ( ((pholead_PhoMCmatchexitcode==1 || pholead_PhoMCmatchexitcode==2) && pholead_GenPhotonIsoDR04<5) && ((photrail_PhoMCmatchexitcode==1 || photrail_PhoMCmatchexitcode==2) && photrail_GenPhotonIsoDR04<5) ) isppevent=1;
	  true_purity[get_name_true_purity(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,isppevent,weight);
	  if (isppevent) true_purity_isppevent[get_name_true_purity_ispp(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,weight);
	  else true_purity_isnotppevent[get_name_true_purity_isnotpp(event_ok_for_dataset_local,*diffvariable)]->Fill(value_diffvariable,weight);
	}

       }



       if (mode=="standard_domatching") { 

	const bool printout = false;

	if (printout){
	  std::cout << "MATCHING" << std::endl;
	  std::cout << pholead_SCeta << " " << event_rho << " " << pholead_pt << std::endl;
	  std::cout << photrail_SCeta << " " << event_rho << " " << photrail_pt << std::endl;
	  std::cout << "---" << std::endl;
	}

	matchingtree_event_fileuuid = event_fileuuid;
	matchingtree_event_run = event_run;
	matchingtree_event_lumi = event_lumi;
	matchingtree_event_number = event_number;

	int *matches1 = new int[nclosestmore];
	Double_t *dists1 = new Double_t[nclosestmore];
	int *matches2 = new int[nclosestmore];
	Double_t *dists2 = new Double_t[nclosestmore];
	Double_t p1[5];
	Double_t p2[5];
	p1[0]=(pholead_SCeta+(fabs(pholead_SCeta)>1.4442)*1000*pholead_SCeta/fabs(pholead_SCeta))/0.1;
	p1[1]=event_rho/2.0;
	p1[2]=TMath::Log(pholead_pt)/0.2;
	p1[3]=(fabs(photrail_SCeta)<1.4442) ? 0 : 1000;
	p1[4]=pholead_SCphi/0.4;
	p2[0]=(photrail_SCeta+(fabs(photrail_SCeta)>1.4442)*1000*photrail_SCeta/fabs(photrail_SCeta))/0.1;
	p2[1]=event_rho/2.0;
	p2[2]=TMath::Log(photrail_pt)/0.2;
	p2[3]=(fabs(pholead_SCeta)<1.4442) ? 0 : 1000;
	p2[4]=photrail_SCphi/0.4;

	for (int n1=0; n1<2; n1++) for (int n2=0; n2<2; n2++){
	    kdtree[n1]->FindNearestNeighbors(p1,nclosestmore,matches1,dists1);
	    kdtree[n2]->FindNearestNeighbors(p2,nclosestmore,matches2,dists2);

	    for (int l=0; l<nclosestmore; l++){

	      bool ismigrerror = false;

	      if (fabs(pholead_SCeta)<1.4442) if (fabs(match_pho_eta[n1].at(matches1[l]))>1.4442) ismigrerror = true;
	      if (fabs(pholead_SCeta)>1.4442) if (fabs(match_pho_eta[n1].at(matches1[l]))<1.4442) ismigrerror = true;
	      if (fabs(photrail_SCeta)<1.4442) if (fabs(match_pho_eta[n2].at(matches2[l]))>1.4442) ismigrerror = true;
	      if (fabs(photrail_SCeta)>1.4442) if (fabs(match_pho_eta[n2].at(matches2[l]))<1.4442) ismigrerror = true;

	      if (ismigrerror){

		cout << "-" << endl;
		cout << n1 << " " << n2 << " " << l << endl;
		cout << pholead_SCeta << " " << photrail_SCeta << endl;
		cout << matches1[l] << " " << matches2[l] << endl;
		cout << match_pho_eta[n1].at(matches1[l]) << " " << match_pho_eta[n2].at(matches2[l]) << endl;
		cout << "-" << endl;

	      }

	      if (n1==0 && n2==0){	      
		matchingtree_index_1event_sigsig_1[l] = matches1[l];
		matchingtree_index_1event_sigsig_2[l] = -999;
		matchingtree_deta_1event_sigsig_1[l] = match_pho_eta[n1].at(matches1[l])-pholead_SCeta;
		matchingtree_deta_1event_sigsig_2[l] = -999;
		matchingtree_drho_1event_sigsig_1[l] = (match_evt_rho[n1].at(matches1[l])-event_rho)/event_sigma;
		matchingtree_drho_1event_sigsig_2[l] = -999;
		matchingtree_dpt_1event_sigsig_1[l] = -999;
		matchingtree_dpt_1event_sigsig_2[l] = -999;
	      }
	      if (n1==0 && n2==1){	      
		matchingtree_index_1event_sigbkg_1[l] = -999;
		matchingtree_index_1event_sigbkg_2[l] = matches2[l];
		matchingtree_deta_1event_sigbkg_1[l] = -999;
		matchingtree_deta_1event_sigbkg_2[l] = match_pho_eta[n2].at(matches2[l])-photrail_SCeta;
		matchingtree_drho_1event_sigbkg_1[l] = -999;
		matchingtree_drho_1event_sigbkg_2[l] = (match_evt_rho[n2].at(matches2[l])-event_rho)/event_sigma;
		matchingtree_dpt_1event_sigbkg_1[l] = -999;
		matchingtree_dpt_1event_sigbkg_2[l] = match_pho_pt[n2].at(matches2[l])/photrail_pt;
	      }
	      if (n1==1 && n2==0){	      
		matchingtree_index_1event_bkgsig_1[l] = matches1[l];
		matchingtree_index_1event_bkgsig_2[l] = -999;
		matchingtree_deta_1event_bkgsig_1[l] = match_pho_eta[n1].at(matches1[l])-pholead_SCeta;
		matchingtree_deta_1event_bkgsig_2[l] = -999;
		matchingtree_drho_1event_bkgsig_1[l] = (match_evt_rho[n1].at(matches1[l])-event_rho)/event_sigma;
		matchingtree_drho_1event_bkgsig_2[l] = -999;
		matchingtree_dpt_1event_bkgsig_1[l] = match_pho_pt[n1].at(matches1[l])/pholead_pt;
		matchingtree_dpt_1event_bkgsig_2[l] = -999;
	      }
	      if (n1==1 && n2==1){	      
		matchingtree_index_1event_bkgbkg_1[l] = matches1[l];
		matchingtree_index_1event_bkgbkg_2[l] = matches2[l];
		matchingtree_deta_1event_bkgbkg_1[l] = match_pho_eta[n1].at(matches1[l])-pholead_SCeta;
		matchingtree_deta_1event_bkgbkg_2[l] = match_pho_eta[n2].at(matches2[l])-photrail_SCeta;
		matchingtree_drho_1event_bkgbkg_1[l] = (match_evt_rho[n1].at(matches1[l])-event_rho)/event_sigma;
		matchingtree_drho_1event_bkgbkg_2[l] = (match_evt_rho[n2].at(matches2[l])-event_rho)/event_sigma;
		matchingtree_dpt_1event_bkgbkg_1[l] = match_pho_pt[n1].at(matches1[l])/pholead_pt;
		matchingtree_dpt_1event_bkgbkg_2[l] = match_pho_pt[n2].at(matches2[l])/photrail_pt;
	      }
	    }
	}




       p1[1]=event_rho/2.0*randomgen->Uniform(0,1);
       p2[1]=event_rho/2.0-p1[1];

	for (int n1=0; n1<2; n1++) for (int n2=0; n2<2; n2++){
	    kdtree[n1]->FindNearestNeighbors(p1,nclosestmore,matches1,dists1);
	    kdtree[n2]->FindNearestNeighbors(p2,nclosestmore,matches2,dists2);

	    for (int l=0; l<nclosestmore; l++){

	      if (fabs(pholead_SCeta)<1.4442) if (fabs(match_pho_eta[n1].at(matches1[l]))>1.4442) cout << "MIGRATION ERROR" << endl;
	      if (fabs(pholead_SCeta)>1.4442) if (fabs(match_pho_eta[n1].at(matches1[l]))<1.4442) cout << "MIGRATION ERROR" << endl;
	      if (fabs(photrail_SCeta)<1.4442) if (fabs(match_pho_eta[n2].at(matches2[l]))>1.4442) cout << "MIGRATION ERROR" << endl;
	      if (fabs(photrail_SCeta)>1.4442) if (fabs(match_pho_eta[n2].at(matches2[l]))<1.4442) cout << "MIGRATION ERROR" << endl;

	      if (n1==0 && n2==0){	      
		matchingtree_index_2events_sigsig_1[l] = matches1[l];
		matchingtree_index_2events_sigsig_2[l] = matches2[l];
		matchingtree_deta_2events_sigsig_1[l] = match_pho_eta[n1].at(matches1[l])-pholead_SCeta;
		matchingtree_deta_2events_sigsig_2[l] = match_pho_eta[n2].at(matches2[l])-photrail_SCeta;
		matchingtree_drho_2events_sigsig_1[l] = (match_evt_rho[n1].at(matches1[l])+match_evt_rho[n2].at(matches2[l])-event_rho)/event_sigma;
		matchingtree_drho_2events_sigsig_2[l] = (match_evt_rho[n1].at(matches1[l])+match_evt_rho[n2].at(matches2[l])-event_rho)/event_sigma;
		matchingtree_dpt_2events_sigsig_1[l] = -999;
		matchingtree_dpt_2events_sigsig_2[l] = -999;
	      }
	      if (n1==0 && n2==1){	      
		matchingtree_index_2events_sigbkg_1[l] = matches1[l];
		matchingtree_index_2events_sigbkg_2[l] = matches2[l];
		matchingtree_deta_2events_sigbkg_1[l] = match_pho_eta[n1].at(matches1[l])-pholead_SCeta;
		matchingtree_deta_2events_sigbkg_2[l] = match_pho_eta[n2].at(matches2[l])-photrail_SCeta;
		matchingtree_drho_2events_sigbkg_1[l] = (match_evt_rho[n1].at(matches1[l])+match_evt_rho[n2].at(matches2[l])-event_rho)/event_sigma;
		matchingtree_drho_2events_sigbkg_2[l] = (match_evt_rho[n1].at(matches1[l])+match_evt_rho[n2].at(matches2[l])-event_rho)/event_sigma;
		matchingtree_dpt_2events_sigbkg_1[l] = -999;
		matchingtree_dpt_2events_sigbkg_2[l] = match_pho_pt[n2].at(matches2[l])/photrail_pt;
	      }
	      if (n1==1 && n2==0){	      
		matchingtree_index_2events_bkgsig_1[l] = matches1[l];
		matchingtree_index_2events_bkgsig_2[l] = matches2[l];
		matchingtree_deta_2events_bkgsig_1[l] = match_pho_eta[n1].at(matches1[l])-pholead_SCeta;
		matchingtree_deta_2events_bkgsig_2[l] = match_pho_eta[n2].at(matches2[l])-photrail_SCeta;
		matchingtree_drho_2events_bkgsig_1[l] = (match_evt_rho[n1].at(matches1[l])+match_evt_rho[n2].at(matches2[l])-event_rho)/event_sigma;
		matchingtree_drho_2events_bkgsig_2[l] = (match_evt_rho[n1].at(matches1[l])+match_evt_rho[n2].at(matches2[l])-event_rho)/event_sigma;
		matchingtree_dpt_2events_bkgsig_1[l] = match_pho_pt[n1].at(matches1[l])/pholead_pt;
		matchingtree_dpt_2events_bkgsig_2[l] = -999;
	      }
	      if (n1==1 && n2==1){	      
		matchingtree_index_2events_bkgbkg_1[l] = matches1[l];
		matchingtree_index_2events_bkgbkg_2[l] = matches2[l];
		matchingtree_deta_2events_bkgbkg_1[l] = match_pho_eta[n1].at(matches1[l])-pholead_SCeta;
		matchingtree_deta_2events_bkgbkg_2[l] = match_pho_eta[n2].at(matches2[l])-photrail_SCeta;
		matchingtree_drho_2events_bkgbkg_1[l] = (match_evt_rho[n1].at(matches1[l])+match_evt_rho[n2].at(matches2[l])-event_rho)/event_sigma;
		matchingtree_drho_2events_bkgbkg_2[l] = (match_evt_rho[n1].at(matches1[l])+match_evt_rho[n2].at(matches2[l])-event_rho)/event_sigma;
		matchingtree_dpt_2events_bkgbkg_1[l] = match_pho_pt[n1].at(matches1[l])/pholead_pt;
		matchingtree_dpt_2events_bkgbkg_2[l] = match_pho_pt[n2].at(matches2[l])/photrail_pt;
	      }
	    }
	  }

	matchingtree->Fill();
	delete[] matches1;
	delete[] dists1;
	delete[] matches2;
	delete[] dists2;
      }
     



 
    } // end if dodistribution
    



  } // end event loop
  std::cout << "Event loop finished" << std::endl;

  if (mode=="standard_domatching"){
    matchingfile->cd();
    matchingtree->Write();
    matchingfile->Close();
  }

//  if (invmass_vector.size()>0){
//    std::sort(invmass_vector.begin(),invmass_vector.end());
//    std::sort(diphotonpt_vector.begin(),diphotonpt_vector.end());
//    for (int i=1; i<10; i++) std::cout << "invmass" << invmass_vector.at(invmass_vector.size()-i) << std::endl;
//    for (int i=1; i<10; i++) std::cout << "diphotonpt" << diphotonpt_vector.at(diphotonpt_vector.size()-i) << std::endl;
//  }

};


#endif

void gen_templates(TString filename="input.root", TString mode="", bool isdata=1, const char* outfile="out.root", TString differentialvariable="photoniso", int maxevents=-1){
  
  TFile *outF = TFile::Open(outfile,"recreate");
  outF->Close();

  TFile *file;
  file = TFile::Open(filename.Data(),"read");

  TTree *t;

  TString treename[18];

  treename[0] = TString("Tree_2Dstandard_selection");
  treename[1] = TString("Tree_1Drandomcone_template");
  treename[2] = TString("Tree_1Dsideband_template");
  treename[3] = TString("Tree_2DZee_pixelvetoreversed_selection");
  treename[4] = TString("Tree_1Dpreselection");
  treename[5] = TString("Tree_1Dselection");
  treename[6] = TString("Tree_2Drandomcone_template");
  treename[7] = TString("Tree_2Drandomconesideband_template");
  treename[8] = TString("Tree_2Dsideband_template");
  treename[9] = TString("Tree_2Dstandard_preselection");
  treename[10] = TString("Tree_2DZmumu_selection");
  treename[11] = TString("Tree_1Dsignal_template");
  treename[12] = TString("Tree_1Dbackground_template");
  treename[13] = TString("Tree_2Dtruesigsig_template");
  treename[14] = TString("Tree_2Dtruesigbkg_template");
  treename[15] = TString("Tree_2Dtruebkgbkg_template");
  treename[16] = TString("Tree_2Drconeplusgenfake_template");
  treename[17] = TString("Tree_2Dgenpromptplussideband_template");


  TString treename_chosen="";
  if (mode=="standard") treename_chosen=treename[0];
  if (mode=="standard_domatching") treename_chosen=treename[0];
  if (mode=="standard_newtemplates_sigsig") treename_chosen=treename[0];
  if (mode=="standard_newtemplates_sigbkg") treename_chosen=treename[0];
  if (mode=="standard_newtemplates_bkgbkg") treename_chosen=treename[0];
  if (mode=="standard_2frag") treename_chosen=treename[0];
  if (mode=="standard_2pgen") treename_chosen=treename[13];
  if (mode=="standard_1p1fbothgen") treename_chosen=treename[14];
  if (mode=="standard_2fgen") treename_chosen=treename[15];
  if (mode=="signal") treename_chosen=treename[11];
  if (mode=="signal_2frag") treename_chosen=treename[11];
  if (mode=="fragmentation") treename_chosen=treename[11];
  if (mode=="nofragmentation") treename_chosen=treename[11];
  if (mode=="cutPFchargediso_signal") treename_chosen=treename[11];
  if (mode=="background") treename_chosen=treename[12];
  if (mode=="cutPFchargediso_background") treename_chosen=treename[12];
  if (mode=="randomcone") treename_chosen=treename[1];
  if (mode=="cutPFchargediso_randomcone") treename_chosen=treename[1];
  if (mode=="sieiesideband") treename_chosen=treename[2];
  if (mode=="cutPFchargediso_sieiesideband") treename_chosen=treename[2];
  if (mode=="zmumu") treename_chosen=treename[10];
  if (mode=="sigsig") treename_chosen=treename[6];
  if (mode=="sigbkg") treename_chosen=treename[7];
  if (mode=="bkgbkg") treename_chosen=treename[8];
  if (mode=="zee") treename_chosen=treename[3];
  if (mode=="standard_pixelrev") treename_chosen=treename[3];
  if (mode=="2pgen") treename_chosen=treename[13];
  if (mode=="2pgen_2frag") treename_chosen=treename[13];
  if (mode=="1p1fbothgen") treename_chosen=treename[14];
  if (mode=="1p1fbothgen_2frag") treename_chosen=treename[14];
  if (mode=="2fgen") treename_chosen=treename[15];
  if (mode=="1prcone1fgen") treename_chosen=treename[16];
  if (mode=="1pgen1fside") treename_chosen=treename[17];
  if (mode=="1pgen1fside_2frag") treename_chosen=treename[17];


  file->GetObject(treename_chosen.Data(),t);

  std::cout << "Processing selection " << treename_chosen.Data() << std::endl;
  
  template_production *temp = new template_production(t);
  temp->Setup(isdata,mode,differentialvariable);
  temp->inputfilename=filename;

  if (maxevents>0) temp->Loop(maxevents); else temp->Loop();
  std::cout << "Exited from event loop" << std::endl;
  temp->WriteOutput(outfile);
  std::cout << "Written output" << std::endl;

  file->Close();

};



void get_eff_area(TString filename, bool doEB, TString comp){

  const char* outfile=Form("puscaling_%s.root",(doEB) ? "EB" : "EE");
  
  TFile *outF = TFile::Open(outfile,"recreate");
  

  TFile *file[1];
  file[0] = TFile::Open(filename.Data());

  TTree *t;

  TString treename;
  treename = TString("Tree_1Drandomcone_template");


  std::vector<std::vector<TProfile*> > output;

  file[0]->GetObject(treename.Data(),t);
  template_production *temp = new template_production(t);
  output=temp->GetPUScaling(doEB,comp);

//  output[0]->Print();
//  output[1]->Print();
//  output[2]->Print();

  TF1 *f_iso[n_bins];
  TF1 *f_rho[n_bins];
  TF1 *f_iso_pu[n_bins];


  for (int i=0; i<n_bins; i++) f_iso[i] = new TF1(Form("f_iso_%d",i),"pol1(0)",5,20);
  for (int i=0; i<n_bins; i++) f_rho[i] = new TF1(Form("f_rho_%d",i),"pol1(0)",5,20);
  for (int i=0; i<n_bins; i++) f_iso_pu[i] = new TF1(Form("f_iso_pu_%d",i),"pol1(0)",5,20);

  for (int i=0; i<n_bins; i++){
    output[0][i]->Fit(f_iso[i],"R");
    output[1][i]->Fit(f_rho[i],"R");
    output[2][i]->Fit(f_iso_pu[i],"R");  
    std::cout << "RESULT (eff. area)" << std::endl << f_iso[i]->GetParameter(1)/f_rho[i]->GetParameter(1) << std::endl;
  }


  for (int i=0; i<n_bins; i++) gROOT->ProcessLine(Form(".! echo %s %s bin%d %e >> puscaling.txt",comp.Data(),(doEB) ? "EB" : "EE",i,f_iso[i]->GetParameter(1)/f_rho[i]->GetParameter(1)));

  outF->cd();

  for (int i=0; i<n_bins; i++){
    output[0][i]->Write();
    output[1][i]->Write();
    output[2][i]->Write();
    f_iso[i]->Write();
    f_rho[i]->Write();
    f_iso_pu[i]->Write();
  }

  file[0]->Close();

  outF->Close();

};

std::vector<std::vector<TProfile*> > template_production::GetPUScaling(bool doEB, TString diffvar){

  TProfile *prof_iso[n_bins];
  TProfile *prof_rho[n_bins];
  TProfile *prof_iso_pu[n_bins];

  for (int i=0; i<n_bins; i++) prof_iso[i] = new TProfile(Form("prof_iso_%d",i),Form("prof_iso_%d",i),30,0,30);
  for (int i=0; i<n_bins; i++) prof_rho[i] = new TProfile(Form("prof_rho_%d",i),Form("prof_rho_%d",i),30,0,30);
  for (int i=0; i<n_bins; i++) prof_iso_pu[i] = new TProfile(Form("prof_iso_pu_%d",i),Form("prof_iso_pu_%d",i),30,0,30);

  //  prof_rho->SetLineColor(kRed);
  //  prof_iso_pu->SetLineColor(kBlue);

  Init();

  if (fChain == 0){
    std::cout << "No chain!" << std::endl;
    return std::vector<std::vector<TProfile*> >();
  } 

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100000==0) { 
      std::cout << "Processing entry " << jentry << std::endl;
    }


    float pholead_phoiso=-999;
    { // recalc phoiso w/cleaning
      float et_recalc = 0;
      float e_recalc = 0;
      int number_recalc = 0;


      for (int i=0; i<pholead_Npfcandphotonincone; i++) {

	float et=pholead_photonpfcandets[i];
	float e=pholead_photonpfcandenergies[i];
	//	float deta=pholead_photonpfcanddetas[i];
	//	float dphi=pholead_photonpfcanddphis[i];
	//	float dR=sqrt(deta*deta+dphi*dphi);
	float eta=fabs(TMath::ACosH(e/et));

	if (eta>1.4442 && eta<1.56) continue;
	if (eta>2.5) continue;

#include "cleaning.cc"

	if (fabs(pholead_SCeta)<1.4442 && eta>1.4442) continue;
	if (fabs(pholead_SCeta)>1.56 && eta<1.56) continue;
	et_recalc+=et;
	e_recalc+=e;
	number_recalc++;

      }


      pholead_phoiso=et_recalc;
    }



    if (diffvar=="photoniso"){
      pholead_outvar=pholead_phoiso;
      assert (pholead_outvar>-100);
    }
    else if (diffvar=="combiso"){
      pholead_outvar=pholead_phoiso+pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01+pholead_pho_Cone04NeutralHadronIso_mvVtx;
    }
    else if (diffvar=="chargediso"){
      pholead_outvar=pholead_pho_Cone04ChargedHadronIso_dR02_dz02_dxy01;
    }
    else if (diffvar=="neutraliso"){
      pholead_outvar=pholead_pho_Cone04NeutralHadronIso_mvVtx;
    }

    if (pholead_pt<25) continue;

    Float_t weight=event_luminormfactor*event_Kfactor*event_weight;

    if (!doEB){
      if (fabs(pholead_SCeta)<1.56) continue;
    }
    else {
      if (fabs(pholead_SCeta)>1.4442) continue;
    }

    Int_t bin_lead = Choose_bin_eta(pholead_SCeta,(doEB) ? 0 : 1);

    prof_iso[bin_lead]->Fill(event_nRecVtx,pholead_outvar,weight);
    prof_iso_pu[bin_lead]->Fill(event_nRecVtx,pholead_outvar,weight);
    prof_rho[bin_lead]->Fill(event_nRecVtx,event_rho*3.14*0.4*0.4,weight);

  }

  std::vector<std::vector<TProfile*> > out;

  out.resize(3);
  out[0].resize(n_bins);
  out[1].resize(n_bins);
  out[2].resize(n_bins);

  for (int i=0; i<n_bins; i++) out[0][i]=prof_iso[i];
  for (int i=0; i<n_bins; i++) out[1][i]=prof_rho[i];
  for (int i=0; i<n_bins; i++) out[2][i]=prof_iso_pu[i];

//  prof_iso->Print();
//  prof_rho->Print();
//  prof_iso_pu->Print();

  return out;

};

float template_production::getpuenergy(int reg, float eta){

  int bin = Choose_bin_eta(fabs(eta),reg);
  float eff_area;
  if (isdata) eff_area = (reg==0) ? eff_areas_EB_data[bin] : eff_areas_EE_data[bin];
  else eff_area = (reg==0) ? eff_areas_EB_mc[bin] : eff_areas_EE_mc[bin];

  return 0.4*0.4*3.14*event_rho*eff_area;

};

float template_production::geteffarea(int reg, float eta){

  int bin = Choose_bin_eta(fabs(eta),reg);
  float eff_area;
  if (isdata) eff_area = (reg==0) ? eff_areas_EB_data[bin] : eff_areas_EE_data[bin];
  else eff_area = (reg==0) ? eff_areas_EB_mc[bin] : eff_areas_EE_mc[bin];

  return 0.4*0.4*3.14*eff_area;

};
