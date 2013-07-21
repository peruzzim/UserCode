void prepare_matchingfile_forstep2(TString matching, TString data, TString outfile="matchingfile_forstep2.root"){

  TFile *m_file = new TFile(matching.Data(),"read");
  TTree *matchingtree = (TTree*)(m_file->Get("matchingtree"));

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

  TFile *newfile = new TFile(outfile.Data(),"recreate");
  newfile->cd();
  TTree *new_matchingtree = matchingtree->CloneTree();
  TTree *new_input0 = InputTree[0]->CloneTree();
  TTree *new_input1 = InputTree[1]->CloneTree();

  new_matchingtree->Write();
  new_input0->Write();
  new_input1->Write();
  newfile->Close();

}
