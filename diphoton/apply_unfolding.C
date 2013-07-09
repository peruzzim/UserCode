#include "binsdef.h"
#include "RooUnfold-1.1.1/src/RooUnfold.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TRandom3.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TProfile.h"
#include "TF1.h"
#include "TString.h"
#include "TH1.h"
#include "TProfile.h"
#include "TMath.h"
#include "TRandom3.h"
#include <map>
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>
#include <algorithm> 

using namespace std;


void apply_unfolding(){

  TFile *unf_file = new TFile("plots/unfolding.root","read");

  TFile *outfile = new TFile("Unfolding_Results.root","recreate");  

  for (int i=0; i<3; i++){
  for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){

    if (*diffvariable=="dR") continue;

    TString reg;
    if (i==0) reg="EBEB"; else if (i==1) reg="EBEE"; else if (i==2) reg="EEEE";

    RooUnfoldResponse *resp = NULL;
    unf_file->GetObject(Form("response_%s_%s",reg.Data(),diffvariable->Data()),resp);
    resp->Print();

    TFile *histo_file = new TFile(Form("plots/histo_xsec_%s_%s.root",diffvariable->Data(),reg.Data()),"read");
    TH1F *histo_reco = NULL;
    histo_file->GetObject("xsec_ngammagammayield",histo_reco);
    histo_reco->Print();

    RooUnfoldBayes *algo = new RooUnfoldBayes(resp,histo_reco,4);
    TH1D *_unfolded = (TH1D*)(algo->Hreco());
    _unfolded->Print();

    TH1F unfolded;
    _unfolded->Copy(unfolded);
    
    std::map<TString,TString> translation2;
    translation2.insert(std::pair<TString,TString>(TString("invmass"),TString("mgg")));
    translation2.insert(std::pair<TString,TString>(TString("diphotonpt"),TString("pt")));
    translation2.insert(std::pair<TString,TString>(TString("costhetastar"),TString("costt")));
    translation2.insert(std::pair<TString,TString>(TString("dphi"),TString("phi")));
    translation2.insert(std::pair<TString,TString>(TString("dR"),TString("dR")));
    unfolded.SetName(Form("Unfolding_Nevt_%s_%s",translation2[TString(diffvariable->Data())].Data(),reg.Data()));
    unfolded.SetTitle(Form("Unfolding_Nevt_%s_%s",translation2[TString(diffvariable->Data())].Data(),reg.Data()));
    
    unfolded.Print();

    TH1F unfolded_relsyserr;
    unfolded.Copy(unfolded_relsyserr);
    unfolded_relsyserr.SetName(Form("Unfolding_RelativeSysErr_%s_%s",translation2[TString(diffvariable->Data())].Data(),reg.Data()));
    unfolded_relsyserr.SetTitle(Form("Unfolding_RelativeSysErr_%s_%s",translation2[TString(diffvariable->Data())].Data(),reg.Data()));
    unfolded_relsyserr.Reset();
    unfolded_relsyserr.Print();

    outfile->cd();
    unfolded.Write();
    unfolded_relsyserr.Write();

  }
  }

  outfile->Close();





}
