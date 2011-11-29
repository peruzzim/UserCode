#define get_sieie_template_cxx
#include "get_sieie_template.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>



using namespace std;

void get_sieie_template::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L get_sieie_template.C
//      Root > get_sieie_template t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

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

   Long64_t nentries = fChain->GetEntriesFast();
   
   fChain->SetBranchStatus("*",0);

   fChain->SetBranchStatus("event_luminormfactor",1);
   fChain->SetBranchStatus("event_Kfactor",1);
   fChain->SetBranchStatus("event_weight",1);

   fChain->SetBranchStatus("pholead_sieie",1);
   fChain->SetBranchStatus("pholead_PhoHasPixSeed",1);
   fChain->SetBranchStatus("pholead_PhoIso04Ecal",1);
   fChain->SetBranchStatus("pholead_PhoIso04Hcal",1);
   fChain->SetBranchStatus("pholead_PhoIso04TrkHollow",1);
   fChain->SetBranchStatus("pholead_hoe",1);
   fChain->SetBranchStatus("pholead_PhoMCmatchexitcode",1);

   fChain->SetBranchStatus("photrail_sieie",1);
   fChain->SetBranchStatus("photrail_PhoHasPixSeed",1);
   fChain->SetBranchStatus("photrail_PhoIso04Ecal",1);
   fChain->SetBranchStatus("photrail_PhoIso04Hcal",1);
   fChain->SetBranchStatus("photrail_PhoIso04TrkHollow",1);
   fChain->SetBranchStatus("photrail_hoe",1);
   fChain->SetBranchStatus("photrail_PhoMCmatchexitcode",1);


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     if (jentry%100000==0) std::cout << "Processing entry " << jentry << std::endl;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

       if (Cut(ientry) < 0) continue;

       Int_t issignal=0;

       if (randomgen->Uniform()>0.5) { // work on lead
	 if (!isdata && (pholead_PhoMCmatchexitcode==1 || pholead_PhoMCmatchexitcode==2)) issignal=1;
	 sieievar->setVal(pholead_sieie);
       }
       else { // work on trail 
	 if (!isdata && (photrail_PhoMCmatchexitcode==1 || photrail_PhoMCmatchexitcode==2)) issignal=1;
	 sieievar->setVal(photrail_sieie);
       }

       hsieie[issignal]->add(RooArgSet(*sieievar),event_luminormfactor*event_Kfactor*event_weight);	 




   }

}
