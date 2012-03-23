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

    // variable definition
    if (varname=="PhoIso04"){
      pholead_outvar=pholead_PhoIso04;
      photrail_outvar=photrail_PhoIso04;
    }
    else if (varname=="sieie"){
      pholead_outvar=pholead_sieie;
      photrail_outvar=photrail_sieie;
    }

    // absolute variables
   if (doabsolute){
     pholead_outvar*=pholead_pt;
     photrail_outvar*=photrail_pt;
   }

   // pu subtraction
   if (dopucorr){
     float eff_area_fraction=0.564;
     const float dR=0.4;
     pholead_outvar-=eff_area_fraction*event_rho*3.14*dR*dR;
     photrail_outvar-=eff_area_fraction*event_rho*3.14*dR*dR;
   }

   //   std::cout << pholead_outvar << " " << pholead_PhoIso04 << " " << pholead_PhoIso04*pholead_pt << " " << pholead_pt << std::endl;

    // initial kinematic selection
    if (pholead_pt<40 || photrail_pt<30 || dipho_mgg_photon<80) continue;

    if (pholead_outvar>rightrange || pholead_outvar<leftrange) continue;
    if (photrail_outvar>rightrange || photrail_outvar<leftrange) continue;

    Int_t event_ok_for_templates=-1;
    Int_t event_ok_for_dataset=-1;

    // categorization:
    // templates 0:EB 1:EE (only EBEB and EEEE events used)
    // dataset 0:EBEB 3/4->1:EBEE 2:EEEE

    if (fabs(pholead_SCeta)<1.4442 && fabs(photrail_SCeta)<1.4442) {
      event_ok_for_templates=0;
      event_ok_for_dataset=0;
    }
    else if (fabs(pholead_SCeta)>1.56 && fabs(photrail_SCeta)>1.56) {
      event_ok_for_templates=1;
      event_ok_for_dataset=2;
    }
    else if (fabs(pholead_SCeta)<1.4442 && fabs(photrail_SCeta)>1.56) {
      event_ok_for_templates=-1;
      event_ok_for_dataset=3;
    }
    else if (fabs(pholead_SCeta)>1.56 && fabs(photrail_SCeta)<1.4442) {
      event_ok_for_templates=-1;
      event_ok_for_dataset=4;
    }

      
    Float_t weight=event_luminormfactor*event_Kfactor*event_weight;




    if (event_ok_for_templates==0 || event_ok_for_templates==1){ // 0:EB, 1:EE

      Int_t indsignal;

      roovar[event_ok_for_templates][0]->setVal(pholead_outvar);
      roovar[event_ok_for_templates][1]->setVal(photrail_outvar);

      roohist[2][event_ok_for_templates][0]->add(*(roovar[event_ok_for_templates][0]),weight);
      roohist[2][event_ok_for_templates][1]->add(*(roovar[event_ok_for_templates][1]),weight);

      templatehist[2][event_ok_for_templates]->Fill(pholead_outvar,weight/2);
      templatehist[2][event_ok_for_templates]->Fill(photrail_outvar,weight/2);

      if (!isdata) {
	if (pholead_PhoMCmatchexitcode==1 || pholead_PhoMCmatchexitcode==2) indsignal=0; else indsignal=1;
	roohist[indsignal][event_ok_for_templates][0]->add(*(roovar[event_ok_for_templates][0]),weight);
	roohist[indsignal][event_ok_for_templates][1]->add(*(roovar[event_ok_for_templates][1]),weight);
	templatehist[indsignal][event_ok_for_templates]->Fill(pholead_outvar,weight/2);
	templatehist[indsignal][event_ok_for_templates]->Fill(photrail_outvar,weight/2);
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
	roovar[1][0]->setVal(pholead_outvar);
	roovar[0][1]->setVal(photrail_outvar);
	list.add(*(roovar[1][0]));
	list.add(*(roovar[0][1]));
	filln=1;
      }

      roodset[filln]->add(list,weight);

    }


  } // end event loop
};


#endif

void gen_templates(TString varname, Float_t leftrange, Float_t rightrange, Int_t nbins, const char* outfile="out.root",bool dopucorr=false){
  
  TFile *outF = TFile::Open(outfile,"recreate");
  outF->Close();

  TString dir("/Users/peruzzi/nobackup/vbox_sync/gg_minitree_020501_16mar12_");
  dir.Append(varname);
  dir.Append("/");

  TString filenameDATA=dir;
  filenameDATA.Append("data_16mar12_");
  filenameDATA.Append(varname);
  filenameDATA.Append(".root");

  TString filenameMC=dir;
  filenameMC.Append("mc_16mar12_");
  filenameMC.Append(varname);
  filenameMC.Append(".root");

  TFile *file[2];
  file[0] = TFile::Open(filenameMC.Data());
  file[1] = TFile::Open(filenameDATA.Data());


  TTree *t;

  TString treename[4];
  treename[0] = TString("Tree_standard_sel");
  treename[1] = TString("Tree_sideband_sel");
  treename[2] = TString("Tree_inclusive_sel");
  treename[3] = TString("Tree_DY_sel");

  for (int isdata=0; isdata<2; isdata++){
    for (int sel_cat=0; sel_cat<4; sel_cat++){
      std::cout << "Processing isdata=" << isdata << " selection " << treename[sel_cat].Data() << std::endl;
      file[isdata]->GetObject(treename[sel_cat].Data(),t);
      template_production *temp = new template_production(t);
      temp->Setup(varname,leftrange,rightrange,nbins,isdata,dopucorr);
      temp->Loop();
      temp->WriteOutput(outfile,isdata,treename[sel_cat].Data());
      //      delete temp;
    }
  }


  file[0]->Close();
  file[1]->Close();

};


void get_eff_area(){
  TString varname("PhoIso04");
  const char* outfile="puscaling.root";
  
  TFile *outF = TFile::Open(outfile,"recreate");
  
  TString dir("/Users/peruzzi/nobackup/vbox_sync/gg_minitree_020501_16mar12_");
  dir.Append(varname);
  dir.Append("/");

  TString filenameMC=dir;
  filenameMC.Append("DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_Summer11_PU_S4_START42_V11_v1.root");

  TFile *file[1];
  file[0] = TFile::Open(filenameMC.Data());

  TTree *t;

  TString treename[4];
  treename[3] = TString("Tree_DY_sel");

  TProfile** output;

  file[0]->GetObject(treename[3].Data(),t);
  template_production *temp = new template_production(t);
  output=temp->GetPUScaling();

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

TProfile** template_production::GetPUScaling(){


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

    pholead_outvar=pholead_PhoIso04;
    photrail_outvar=photrail_PhoIso04;
    
    pholead_outvar*=pholead_pt;
    photrail_outvar*=photrail_pt;
    
    if (pholead_pt<40 || photrail_pt<30 || dipho_mgg_photon<80) continue;

    Float_t weight=event_luminormfactor*event_Kfactor*event_weight;

    //    std::cout << "Filling" << pholead_outvar << " " << photrail_outvar << " " << event_rho << std::endl;

    prof_iso->Fill(event_nRecVtx,pholead_outvar,weight/2);
    prof_iso->Fill(event_nRecVtx,photrail_outvar,weight/2);

    const float eff_area_fraction=0.564;
    const float dR=0.4;
    pholead_outvar-=eff_area_fraction*event_rho*3.14*dR*dR;
    photrail_outvar-=eff_area_fraction*event_rho*3.14*dR*dR;
    
    prof_iso_pu->Fill(event_nRecVtx,pholead_outvar,weight/2);
    prof_iso_pu->Fill(event_nRecVtx,photrail_outvar,weight/2);


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
