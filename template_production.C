#ifndef template_production_cxx
#define template_production_cxx
#include "template_production.h"

using namespace std;


void template_production::Loop()
{
  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;
  if (!initialized){
    std::cout << "Not initialized! Call Setup() first." << std::endl;
    return;
  }

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100000==0) std::cout << "Processing entry " << jentry << std::endl;

    pholead_outvar=pholead_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;
    photrail_outvar=photrail_pho_Cone04PhotonIso_dEta015EB_dR070EE_mvVtx;

    // initial kinematic selection
    if (dodistribution) if (pholead_pt<40 || photrail_pt<30 || dipho_mgg_photon<80) continue;
    if (dosignaltemplate || dobackgroundtemplate) if (pholead_pt<40) continue;

    Int_t event_ok_for_dataset=-1;

    Int_t reg_lead;
    Int_t reg_trail;

    // categorization:
    // templates 0:EB 1:EE (only EBEB and EEEE events used)
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

    Float_t weight=event_luminormfactor*event_Kfactor*event_weight;

    Int_t bin_lead = Choose_bin_pt(pholead_pt,reg_lead);

    Int_t bin_trail = -999;
    Int_t bin_couple = -999;
    if (dodistribution) {
      bin_trail = Choose_bin_pt(photrail_pt,reg_trail);
      bin_couple = Choose_bin_invmass(dipho_mgg_photon,event_ok_for_dataset);
    }
    

    // Int_t bin_lead = Choose_bin_eta(pholead_SCeta,reg_lead);
    // Int_t bin_trail = Choose_bin_eta(photrail_SCeta,reg_trail);

    if (pholead_outvar<leftrange) pholead_outvar=leftrange;
    if (pholead_outvar>rightrange) pholead_outvar=rightrange;
    if (photrail_outvar<leftrange) photrail_outvar=leftrange;
    if (photrail_outvar>rightrange) photrail_outvar=rightrange;

    if (dosignaltemplate){
      template_signal[reg_lead][bin_lead]->Fill(pholead_outvar,weight);
    }
    if (dobackgroundtemplate){
      template_background[reg_lead][bin_lead]->Fill(pholead_outvar,weight);
    }

    if (dodistribution && event_ok_for_dataset>-1){
      
      obs_hist_single[reg_lead][bin_lead]->Fill(pholead_outvar,weight);
      obs_hist_single[reg_trail][bin_trail]->Fill(photrail_outvar,weight);

      float in1=pholead_outvar;
      float in2=photrail_outvar;

      bool doswap=false;

      if ((event_ok_for_dataset==0 || event_ok_for_dataset==2) && (randomgen->Uniform()>0.5)) doswap=true;

      if (event_ok_for_dataset==4) doswap=true;

      if (event_ok_for_dataset==3 || event_ok_for_dataset==4) event_ok_for_dataset=1;

      if (doswap){
	  in2=pholead_outvar;
	  in1=photrail_outvar;
      }

      obs_hist[event_ok_for_dataset][bin_couple]->Fill(in1,in2,weight);
	
    }
    


  } // end event loop
  std::cout << "ended event loop" << std::endl;
};


#endif

void gen_templates(TString filename="input.root", TString mode="", bool isdata=1, const char* outfile="out.root"){
  
  TFile *outF = TFile::Open(outfile,"recreate");
  outF->Close();

  TFile *file;
  file = TFile::Open(filename.Data(),"read");

  TTree *t;

  TString treename[5];
  treename[0] = TString("Tree_standard_sel");
  treename[1] = TString("Tree_signal_template");
  treename[2] = TString("Tree_background_template");
  treename[3] = TString("Tree_DY_sel");
  treename[4] = TString("Tree_randomcone_signal_template");

  TString treename_chosen="";
  if (mode=="standard") treename_chosen=treename[0];
  if (mode=="signal") treename_chosen=treename[1];
  if (mode=="background") treename_chosen=treename[2];
  if (mode=="randomcone") treename_chosen=treename[4];

  file->GetObject(treename_chosen.Data(),t);

  std::cout << "Processing selection " << treename_chosen.Data() << std::endl;
  template_production *temp = new template_production(t);
  temp->Setup(isdata,mode);
  temp->Loop();
  std::cout << "exited from event loop" << std::endl;
  temp->WriteOutput(outfile,treename_chosen.Data());
  std::cout << "written output" << std::endl;

  file->Close();

};



void get_eff_area(bool doEB){

  const char* outfile=Form("puscaling_%s.root",(doEB) ? "EB" : "EE");
  
  TFile *outF = TFile::Open(outfile,"recreate");
  
  TString dir("/Users/peruzzi/nobackup/ntuples/gg_minitree_020612_purew2011B_noselection/");

  TString filenameMC=dir;
  filenameMC.Append("DiPhotonJets_7TeV_madgraph_Fall11_PU_S6_START42_V14B_v1_AODSIM.root");

  TFile *file[1];
  file[0] = TFile::Open(filenameMC.Data());

  TTree *t;

  TString treename[4];
  treename[3] = TString("Tree_signal_template");

  TProfile** output;

  file[0]->GetObject(treename[3].Data(),t);
  template_production *temp = new template_production(t);
  output=temp->GetPUScaling(doEB);

  output[0]->Print();
  output[1]->Print();
  output[2]->Print();

  TF1 *f_iso = new TF1("f_iso","pol1(0)",0,30);
  TF1 *f_rho = new TF1("f_rho","pol1(0)",0,30);
  TF1 *f_iso_pu = new TF1("f_iso_pu","pol1(0)",0,30);

  output[0]->Fit(f_iso);
  output[1]->Fit(f_rho);
  output[2]->Fit(f_iso_pu);

 
  outF->cd();
  output[0]->Write();
  output[1]->Write();
  output[2]->Write();
  f_iso->Write();
  f_rho->Write();
  f_iso_pu->Write();

  file[0]->Close();

  outF->Close();
};

TProfile** template_production::GetPUScaling(bool doEB){


  TProfile *prof_iso = new TProfile("prof_iso","prof_iso",30,0,30);
  TProfile *prof_rho = new TProfile("prof_rho","prof_rho",30,0,30);
  TProfile *prof_iso_pu = new TProfile("prof_iso_pu","prof_iso_pu",30,0,30);

  prof_rho->SetLineColor(kRed);
  prof_iso_pu->SetLineColor(kBlue);

  Init();

  if (fChain == 0){
    std::cout << "No chain!" << std::endl;
    return NULL;
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

    pholead_outvar=pholead_pho_Cone04PFCombinedIso*pholead_pt;
    
    if (pholead_pt<40) continue;

    Float_t weight=event_luminormfactor*event_Kfactor*event_weight;

    if (!doEB){
      if (fabs(pholead_SCeta)<1.4442) continue;
    }
    else {
      if (fabs(pholead_SCeta)>1.56) continue;
    }

    prof_iso->Fill(event_nRecVtx,pholead_outvar,weight/2);

     float eff_area_fraction_EB=0.406;
     float eff_area_fraction_EE=0.528;
     const float dR=0.4;
     pholead_outvar-=event_rho*3.14*dR*dR*((fabs(pholead_SCeta)<1.4442) ? eff_area_fraction_EB : eff_area_fraction_EE);
    
    prof_iso_pu->Fill(event_nRecVtx,pholead_outvar,weight/2);


    prof_rho->Fill(event_nRecVtx,event_rho*3.14*0.4*0.4,weight);

    

  }

  TProfile** out = new TProfile*[3];
  out[0]=prof_iso;
  out[1]=prof_rho;
  out[2]=prof_iso_pu;

  prof_iso->Print();
  prof_rho->Print();
  prof_iso_pu->Print();

  return out;

};

