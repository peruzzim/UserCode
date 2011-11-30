#include "TFile.h"
#include "TTree.h"


void gen_templates(const char* filename, bool isdata, bool sideband, bool forbarrel, const char* outfile){

  TFile *file = TFile::Open(filename);
  TTree *t;
  file->GetObject("Tree",t);

  get_sieie_template *temp = new get_sieie_template(t);
  temp->SetParams(isdata,sideband,forbarrel);
  temp->Loop();
  temp->WriteOutput(outfile);

  
 file->Close();

 
}
