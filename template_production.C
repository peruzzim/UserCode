#ifndef template_production_cxx
#define template_production_cxx
#include "template_production.h"


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

    if (pholead_pt<40 || photrail_pt<30 || dipho_mgg_photon<80) continue;

    Bool_t barrel;

    if (fabs(pholead_eta)<1.4442 && fabs(photrail_eta)<1.4442) barrel=true;
    else if (fabs(pholead_eta)>1.56 && fabs(photrail_eta)>1.56) barrel=false;
    else continue;

    bool badevent=false;
    
    if (varname=="sieie"){
      if (sidebandoption==TString("sideband")) {
	if (Cut_egm_10_006_sieierelaxed_sideband() < 0) badevent=true;
      }
      else {
	if  (Cut_egm_10_006_sieierelaxed() < 0) badevent=true;
      }
    }
    else {
      std::cout << "Variable name not known!!!" << std::endl;
      continue;
    }

    if (badevent) continue;
      
    Int_t issignal;

    outroovar->setVal(pholead_outvar);
    roohout[barrel][2]->add(*outroovar,event_luminormfactor*event_Kfactor*event_weight);
    if (!isdata) {
      if (pholead_PhoMCmatchexitcode==1 || pholead_PhoMCmatchexitcode==2) issignal=1;
      roohout[barrel][issignal]->add(*outroovar,event_luminormfactor*event_Kfactor*event_weight);
    }
      
    outroovar->setVal(photrail_outvar);
    roohout[barrel][2]->add(*outroovar,event_luminormfactor*event_Kfactor*event_weight);
    if (!isdata) {
      if (photrail_PhoMCmatchexitcode==1 || photrail_PhoMCmatchexitcode==2) issignal=1;
      roohout[barrel][issignal]->add(*outroovar,event_luminormfactor*event_Kfactor*event_weight);
    }
   

  }
};

Int_t template_production::Cut_egm_10_006_sieierelaxed(){

// egm10006 loose, sietaieta relaxed

  if (pholead_PhoIso04Ecal>4.2) return -1;
  if (pholead_PhoIso04Hcal>2.2) return -1;
  if (pholead_PhoIso04TrkHollow>2.0) return -1;
  if (pholead_hoe>0.05) return -1;  

  if (photrail_PhoIso04Ecal>4.2) return -1; 
  if (photrail_PhoIso04Hcal>2.2) return -1; 
  if (photrail_PhoIso04TrkHollow>2.0) return -1; 
  if (photrail_hoe>0.05) return -1; 
  
   return 1;
};

Int_t template_production::Cut_egm_10_006_sieierelaxed_sideband(){

// egm10006 loose, sietaieta relaxed, trkiso sideband (2<trkiso<5)
  
  if (pholead_PhoIso04Ecal>4.2) return -1; 
  if (pholead_PhoIso04Hcal>2.2) return -1; 
  if (pholead_PhoIso04TrkHollow<2.0) return -1; 
  if (pholead_PhoIso04TrkHollow>5.0) return -1; 
  if (pholead_hoe>0.05) return -1; 

  if (photrail_PhoIso04Ecal>4.2) return -1; 
  if (photrail_PhoIso04Hcal>2.2) return -1; 
  if (photrail_PhoIso04TrkHollow<2.0) return -1; 
  if (photrail_PhoIso04TrkHollow>5.0) return -1; 
  if (photrail_hoe>0.05) return -1; 
  
   return 1;
};

#endif

