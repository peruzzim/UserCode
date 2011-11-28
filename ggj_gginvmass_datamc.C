#include <string>
#include <vector>
#include <stdlib>
#include <iostream>

using namespace std;

TCut pholead_cut_ecaliso =   "pholead_PhoIso04Ecal<4.2";
TCut pholead_cut_hcaliso =   "pholead_PhoIso04Hcal<2.2";
TCut pholead_cut_trkiso =    "pholead_PhoIso04TrkHollow<2.0";
TCut pholead_cut_hovere =    "pholead_hoe<0.05";
TCut pholead_cut_sietaieta = "(TMath::Abs(pholead_SCeta)<1.4442 && pholead_sieie<0.01) || (TMath::Abs(pholead_SCeta)>1.56 && pholead_sieie<0.03)";
TCut pholead_cut_egm_10_006_loose = pholead_cut_ecaliso * pholead_cut_hcaliso * pholead_cut_trkiso * pholead_cut_hovere * pholead_cut_sietaieta;

TCut photrail_cut_ecaliso =   "photrail_PhoIso04Ecal<4.2";
TCut photrail_cut_hcaliso =   "photrail_PhoIso04Hcal<2.2";
TCut photrail_cut_trkiso =    "photrail_PhoIso04TrkHollow<2.0";
TCut photrail_cut_hovere =    "photrail_hoe<0.05";
TCut photrail_cut_sietaieta = "(TMath::Abs(photrail_SCeta)<1.4442 && photrail_sieie<0.01) || (TMath::Abs(photrail_SCeta)>1.56 && photrail_sieie<0.03)";
TCut photrail_cut_egm_10_006_loose = photrail_cut_ecaliso * photrail_cut_hcaliso * photrail_cut_trkiso * photrail_cut_hovere * photrail_cut_sietaieta;


void ggj_gginvmass_datamc(const char* plot, TCut cut="", float scalefactormc=2.1){

  cut = cut*"event_luminormfactor*event_weight*event_Kfactor*(pholead_pt>40)*(photrail_pt>30)";

  const int ndata=3;
  const int nmc=6;

  TFile *fdata[ndata];
  TFile *fmc[nmc];

  TString basedir("/Users/peruzzi/work/trees/gg_minitree_020501/");

  fdata[0]= TFile::Open(basedir+"Photon_Run2011A_03Oct2011_v1.root","read");
  fdata[1]= TFile::Open(basedir+"Photon_Run2011A_05Aug2011_v1.root","read");
  fdata[2]= TFile::Open(basedir+"Photon_Run2011A_05Jul2011ReReco_ECAL_v1.root","read");

  fmc[0]= TFile::Open(basedir+"QCD_Pt_40_doubleEMEnriched_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2.root","read");
  fmc[1]= TFile::Open(basedir+"QCD_Pt_30to40_doubleEMEnriched_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2.root","read");
  fmc[2]= TFile::Open(basedir+"GJet_Pt_20_doubleEMEnriched_TuneZ2_7TeV_pythia6_Summer11_PU_S4_START42_V11_v1.root","read");
  fmc[3]= TFile::Open(basedir+"DiPhotonBox_Pt_25To250_7TeV_pythia6_Summer11_PU_S4_START42_V11_v2.root","read");
  fmc[4]= TFile::Open(basedir+"DiPhotonJets_7TeV_madgraph_Summer11_PU_S4_START42_V11_v1.root","read");
  fmc[5]= TFile::Open(basedir+"DYJetsToLL_TuneZ2_M_50_7TeV_madgraph_tauola_Summer11_PU_S4_START42_V11_v1.root","read");

  TTree *tdata[ndata];
  TTree *tmc[nmc];

  for (int i=0;i<ndata;i++) fdata[i]->GetObject("Tree",tdata[i]);
  for (int i=0;i<nmc;i++) fmc[i]->GetObject("Tree",tmc[i]);

  TH1F *hdata[ndata];
  TH1F *hmc[nmc];

  for (int i=0;i<ndata;i++){
    TString plotstring(plot);
    plotstring.Append(">>htemp");
    tdata[i]->Draw(plotstring.Data(),cut);
    hdata[i] = new TH1F(*htemp);
  }

 for (int i=0;i<nmc;i++){
    TString plotstring(plot);
    plotstring.Append(">>htemp");
    tmc[i]->Draw(plotstring.Data(),cut);
    hmc[i] = new TH1F(*htemp);
    hmc[i]->SetFillStyle(3001);
    hmc[i]->SetFillColor(i+1);
    hmc[i]->SetLineColor(kBlack);
    hmc[i]->Scale(scalefactormc);
  }

 TH1F *htotdata = new TH1F(*hdata[0]);
 for (int i=0; i<ndata; i++) htotdata->Add(hdata[i]);
 htotdata->SetFillStyle(0);
 htotdata->SetLineColor(1);
 htotdata->SetLineWidth(3);
 htotdata->SetMarkerStyle(20);

 THStack *htotmc =  new THStack("htotmc",plot);
 for (int i=0; i<nmc; i++) htotmc->Add(hmc[i]);
 
  htotmc->Draw();
  htotdata->Draw("sameE1");

  TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(htotdata,"data 2011A","f");
  leg->AddEntry(hmc[0],"QCD 40","f");
  leg->AddEntry(hmc[1],"QCD 30_40","f");
  leg->AddEntry(hmc[2],"GJet","f");
  leg->AddEntry(hmc[3],"DiPhotonBox","f");
  leg->AddEntry(hmc[4],"DiPhotonJets","f");
  leg->AddEntry(hmc[5],"DYJetsToLL","f");
  leg->Draw();
  
}
