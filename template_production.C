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

    // initial kinematic selection
    if (pholead_pt<40 || photrail_pt<30 || dipho_mgg_photon<80) continue;

    Int_t event_ok_for_templates=-1;
    Int_t event_ok_for_dataset=-1;

    // categorization:
    // templates 0:EB 1:EE
    // dataset 0:EBEB 3/4->1:EBEE 2:EEEE

    if (fabs(pholead_eta)<1.4442 && fabs(photrail_eta)<1.4442) {
      event_ok_for_templates=0;
      event_ok_for_dataset=0;
    }
    else if (fabs(pholead_eta)>1.56 && fabs(photrail_eta)>1.56) {
      event_ok_for_templates=1;
      event_ok_for_dataset=2;
    }
    else if (fabs(pholead_eta)<1.4442 && fabs(photrail_eta)>1.56) {
      event_ok_for_templates=-1;
      event_ok_for_dataset=3;
    }
    else if (fabs(pholead_eta)>1.56 && fabs(photrail_eta)<1.4442) {
      event_ok_for_templates=-1;
      event_ok_for_dataset=4;
    }

    
    if (varname=="sieie"){
      if (sidebandoption==TString("sideband")) {
	if (Cut_egm_10_006_sieierelaxed_sideband() < 0) {
	  event_ok_for_templates=-1;
	  event_ok_for_dataset=-1;
	}
      }
      else if  (Cut_egm_10_006_sieierelaxed() < 0) {
	  event_ok_for_templates=-1;
	  event_ok_for_dataset=-1;
	}
    }
    else {
      std::cout << "Variable name not known!!!" << std::endl;
      continue;
    }
      
    Float_t weight=event_luminormfactor*event_Kfactor*event_weight;

    if (event_ok_for_templates==0 || event_ok_for_templates==1){ // 0:EB, 1:EE

      Int_t indsignal;

      roovar[event_ok_for_templates][0]->setVal(pholead_outvar);
      roovar[event_ok_for_templates][1]->setVal(pholead_outvar);
      roohist[2][event_ok_for_templates][0]->add(*(roovar[event_ok_for_templates][0]),weight);
      roohist[2][event_ok_for_templates][1]->add(*(roovar[event_ok_for_templates][1]),weight);
      templatehist[2][event_ok_for_templates]->Fill(pholead_outvar,weight);
      if (!isdata) {
	if (pholead_PhoMCmatchexitcode==1 || pholead_PhoMCmatchexitcode==2) indsignal=0; else indsignal=1;
	roohist[indsignal][event_ok_for_templates][0]->add(*(roovar[event_ok_for_templates][0]),weight);
	roohist[indsignal][event_ok_for_templates][1]->add(*(roovar[event_ok_for_templates][1]),weight);
	templatehist[indsignal][event_ok_for_templates]->Fill(pholead_outvar,weight);
      }
    
   
      roovar[event_ok_for_templates][0]->setVal(photrail_outvar);
      roovar[event_ok_for_templates][1]->setVal(photrail_outvar);
      roohist[2][event_ok_for_templates][0]->add(*(roovar[event_ok_for_templates][0]),weight);
      roohist[2][event_ok_for_templates][1]->add(*(roovar[event_ok_for_templates][1]),weight);
      templatehist[2][event_ok_for_templates]->Fill(photrail_outvar,weight);
      if (!isdata) {
	if (photrail_PhoMCmatchexitcode==1 || photrail_PhoMCmatchexitcode==2) indsignal=0; else indsignal=1;
	roohist[indsignal][event_ok_for_templates][0]->add(*(roovar[event_ok_for_templates][0]),weight);
	roohist[indsignal][event_ok_for_templates][1]->add(*(roovar[event_ok_for_templates][1]),weight);
	templatehist[indsignal][event_ok_for_templates]->Fill(photrail_outvar,weight);
      }
    
    }

    if (event_ok_for_dataset>-1){

      RooArgSet list;
      Int_t filln=-1;
     
      if (event_ok_for_dataset==0){ // EBEB
	if (randomgen->Uniform()>0.5){
	  roovar[0][0]->setVal(pholead_outvar);
	  roovar[0][1]->setVal(photrail_outvar);
	}
	else {
	  roovar[0][1]->setVal(pholead_outvar);
	  roovar[0][0]->setVal(photrail_outvar);
	}
	list.add(*(roovar[0][0]));
	list.add(*(roovar[0][1]));
	filln=0;
      }
      else if (event_ok_for_dataset==2){ // EEEE
	if (randomgen->Uniform()>0.5){
	  roovar[1][0]->setVal(pholead_outvar);
	  roovar[1][1]->setVal(photrail_outvar);
	}
	else {
	  roovar[1][1]->setVal(pholead_outvar);
	  roovar[1][0]->setVal(photrail_outvar);
	}
	list.add(*(roovar[1][0]));
	list.add(*(roovar[1][1]));
	filln=2;
      }
      else if (event_ok_for_dataset==3){ // EBEE ordered
	roovar[0][0]->setVal(pholead_outvar);
	roovar[1][1]->setVal(photrail_outvar);
	list.add(*(roovar[0][0]));
	list.add(*(roovar[1][1]));
	filln=1;
      }
      else if (event_ok_for_dataset==4){ // EEEB ordered
	roovar[1][1]->setVal(pholead_outvar);
	roovar[0][0]->setVal(photrail_outvar);
	list.add(*(roovar[0][0]));
	list.add(*(roovar[1][1]));
	filln=1;
      }

      roodset[filln]->add(list,weight);

    }










   

  } // end event loop
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

void gen_templates(const char* filename, TString varname, Float_t leftrange, Float_t rightrange, Int_t nbins, Bool_t isdata, TString sidebandoption="", const char* outfile="out.root"){
  
  TFile *file = TFile::Open(filename);
  TTree *t;
  file->GetObject("Tree",t);
  
  template_production *temp = new template_production(t);
  temp->Setup(varname,leftrange,rightrange,nbins,isdata,sidebandoption);
  temp->Loop();
  temp->WriteOutput(outfile);
  
  file->Close();

};

