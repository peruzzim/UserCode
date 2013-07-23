#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

const int nclosest = 10;

void prepare_matchingfile_forstep2(TString matching, TString data){

  TFile *m_file = new TFile(matching.Data(),"read");
  TTree *matchingtree = (TTree*)(m_file->Get("matchingtree"));

   UInt_t          matchingtree_event_fileuuid;
   Int_t           matchingtree_event_run;
   Int_t           matchingtree_event_lumi;
   Int_t           matchingtree_event_number;
   Int_t           matchingtree_index_sigsig_1[nclosest];
   Int_t           matchingtree_index_sigsig_2[nclosest];
   Int_t           matchingtree_index_sigbkg_1[nclosest];
   Int_t           matchingtree_index_sigbkg_2[nclosest];
   Int_t           matchingtree_index_bkgsig_1[nclosest];
   Int_t           matchingtree_index_bkgsig_2[nclosest];
   Int_t           matchingtree_index_bkgbkg_1[nclosest];
   Int_t           matchingtree_index_bkgbkg_2[nclosest];
   TBranch        *b_matchingtree_event_fileuuid;   //!
   TBranch        *b_matchingtree_event_run;   //!
   TBranch        *b_matchingtree_event_lumi;   //!
   TBranch        *b_matchingtree_event_number;   //!
   TBranch        *b_matchingtree_index_sigsig_1;   //!
   TBranch        *b_matchingtree_index_sigsig_2;   //!
   TBranch        *b_matchingtree_index_sigbkg_1;   //!
   TBranch        *b_matchingtree_index_sigbkg_2;   //!
   TBranch        *b_matchingtree_index_bkgsig_1;   //!
   TBranch        *b_matchingtree_index_bkgsig_2;   //!
   TBranch        *b_matchingtree_index_bkgbkg_1;   //!
   TBranch        *b_matchingtree_index_bkgbkg_2;   //!
   matchingtree->SetBranchAddress("matchingtree_event_fileuuid", &matchingtree_event_fileuuid, &b_matchingtree_event_fileuuid);
   matchingtree->SetBranchAddress("matchingtree_event_run", &matchingtree_event_run, &b_matchingtree_event_run);
   matchingtree->SetBranchAddress("matchingtree_event_lumi", &matchingtree_event_lumi, &b_matchingtree_event_lumi);
   matchingtree->SetBranchAddress("matchingtree_event_number", &matchingtree_event_number, &b_matchingtree_event_number);
   matchingtree->SetBranchAddress("matchingtree_index_sigsig_1", matchingtree_index_sigsig_1, &b_matchingtree_index_sigsig_1);
   matchingtree->SetBranchAddress("matchingtree_index_sigsig_2", matchingtree_index_sigsig_2, &b_matchingtree_index_sigsig_2);
   matchingtree->SetBranchAddress("matchingtree_index_sigbkg_1", matchingtree_index_sigbkg_1, &b_matchingtree_index_sigbkg_1);
   matchingtree->SetBranchAddress("matchingtree_index_sigbkg_2", matchingtree_index_sigbkg_2, &b_matchingtree_index_sigbkg_2);
   matchingtree->SetBranchAddress("matchingtree_index_bkgsig_1", matchingtree_index_bkgsig_1, &b_matchingtree_index_bkgsig_1);
   matchingtree->SetBranchAddress("matchingtree_index_bkgsig_2", matchingtree_index_bkgsig_2, &b_matchingtree_index_bkgsig_2);
   matchingtree->SetBranchAddress("matchingtree_index_bkgbkg_1", matchingtree_index_bkgbkg_1, &b_matchingtree_index_bkgbkg_1);
   matchingtree->SetBranchAddress("matchingtree_index_bkgbkg_2", matchingtree_index_bkgbkg_2, &b_matchingtree_index_bkgbkg_2);


  TFile *d_file = new TFile(data.Data(),"read");
  TTree *InputTree[2];
  InputTree[0] = (TTree*)(d_file->Get("Tree_1Drandomcone_template"));
  InputTree[1] = (TTree*)(d_file->Get("Tree_1Dsideband_template"));

  for (int i=0; i<2; i++){
    InputTree[i]->SetBranchStatus("*",0);
    InputTree[i]->SetBranchStatus("pholead_SCeta",1);
    InputTree[i]->SetBranchStatus("pholead_SCphi",1);
    InputTree[i]->SetBranchStatus("allphotonpfcand_count",1);
    InputTree[i]->SetBranchStatus("allphotonpfcand_pt",1);
    InputTree[i]->SetBranchStatus("allphotonpfcand_eta",1);
    InputTree[i]->SetBranchStatus("allphotonpfcand_phi",1);
    InputTree[i]->SetBranchStatus("allphotonpfcand_vx",1);
    InputTree[i]->SetBranchStatus("allphotonpfcand_vy",1);
    InputTree[i]->SetBranchStatus("allphotonpfcand_vz",1);
  }


   std::vector<UInt_t> lista_uuid;
   lista_uuid.push_back(-1);
   for (Long64_t i=0;i<matchingtree->GetEntriesFast(); i++) {matchingtree->GetEntry(i); lista_uuid.push_back(matchingtree_event_fileuuid);}
   std::sort(lista_uuid.begin(),lista_uuid.end());
   unsigned int size_beforeunique = lista_uuid.size();
   lista_uuid.erase(std::unique(lista_uuid.begin(),lista_uuid.end()),lista_uuid.end());
   assert(lista_uuid.size()==size_beforeunique);
   cout << lista_uuid.size() << " uuids found" << endl;

   std::map<long, TFile*> newmatchingfile;
   std::map<long, TTree*> newmatchingtree;
   std::map<long, TTree*> NewInputTree[2];

   Int_t matchingtree_index_1[2][2][nclosest];
   Int_t matchingtree_index_2[2][2][nclosest];

   UInt_t thisuuid=0;
   int index_uuid=0;

     for (Long64_t i=0;i<matchingtree->GetEntriesFast(); i++) {
       matchingtree->GetEntry(i);
       if (matchingtree_event_fileuuid!=thisuuid) {
	 cout << thisuuid << endl;
	 if (thisuuid!=0){
	   newmatchingtree[thisuuid]->Write();
	   for (int a=0; a<2; a++) NewInputTree[a][thisuuid]->Write(); 
	   newmatchingfile[thisuuid]->Close();
	 }
	 thisuuid=matchingtree_event_fileuuid;
	 cout << ++index_uuid << ". processing uuid " << thisuuid << endl;
	 TFile *newfile = new TFile(Form("matchingtree_%u.root",thisuuid),"recreate");
	 newmatchingfile.insert(std::pair<long,TFile*>(thisuuid,newfile));
	 newfile->cd();
	 newmatchingtree.insert(std::pair<long,TTree*>(thisuuid,matchingtree->CloneTree(0)));
	 for (int m=0; m<2; m++) assert( NewInputTree[m].find(thisuuid)==NewInputTree[m].end()  );
	 for (int m=0; m<2; m++) NewInputTree[m].insert(std::pair<long,TTree*>(thisuuid,InputTree[m]->CloneTree(0)));
       }
       for (int l=0; l<nclosest; l++){
	 matchingtree_index_1[0][0][l]=matchingtree_index_sigsig_1[l];
	 matchingtree_index_1[0][1][l]=matchingtree_index_sigbkg_1[l];
	 matchingtree_index_1[1][0][l]=matchingtree_index_bkgsig_1[l];
	 matchingtree_index_1[1][1][l]=matchingtree_index_bkgbkg_1[l];
	 matchingtree_index_2[0][0][l]=matchingtree_index_sigsig_2[l];
	 matchingtree_index_2[0][1][l]=matchingtree_index_sigbkg_2[l];
	 matchingtree_index_2[1][0][l]=matchingtree_index_bkgsig_2[l];
	 matchingtree_index_2[1][1][l]=matchingtree_index_bkgbkg_2[l];
       }
       for (int n1=0; n1<2; n1++) for (int n2=0; n2<2; n2++) for (int l=0; l<nclosest; l++){
	 InputTree[n1]->GetEntry(matchingtree_index_1[n1][n2][l]);
	 NewInputTree[n1][matchingtree_event_fileuuid]->Fill();
	 InputTree[n2]->GetEntry(matchingtree_index_2[n1][n2][l]);
	 NewInputTree[n2][matchingtree_event_fileuuid]->Fill();
	 matchingtree_index_1[n1][n2][l] = NewInputTree[n1][matchingtree_event_fileuuid]->GetEntriesFast();
	 matchingtree_index_2[n1][n2][l] = NewInputTree[n2][matchingtree_event_fileuuid]->GetEntriesFast();
       }
       newmatchingtree[matchingtree_event_fileuuid]->Fill();
     }



}
