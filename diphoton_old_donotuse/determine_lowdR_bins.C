#include "binsdef.h"
#include "TFile.h"
#include "TString.h"
#include "RooDataSet.h"
#include <iostream>

using namespace RooFit;
using namespace std;

void determine_lowdR_bins(TString file){

  TFile *f = new TFile(file.Data());

  for (std::vector<TString>::const_iterator diffvariable = diffvariables_list.begin(); diffvariable!=diffvariables_list.end(); diffvariable++){

    TString region[3]={"EBEB","EBEE","EEEE"};

    for (int reg=0; reg<3; reg++){

    for (int bin=0; bin<=n_bins; bin++){
      
      RooDataSet *dset = NULL;
      f->GetObject(Form("roofit/obs_roodset_%s_%s_b%d",region[reg].Data(),diffvariable->Data(),bin),dset);
      if (!dset) continue;

      float x = dset->sumEntries("roovar_dR<1.0")/dset->sumEntries();
      if (x>0.1) cout << dset->GetName() << " " << x << " <-------------------------" << endl;


    }

    }





  }









};
