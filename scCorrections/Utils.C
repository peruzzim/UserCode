#include <iostream>
#include <sstream>
#include "TString.h"
#include "TMath.h"

TString IntToString(int number){
  ostringstream oss;
  oss << number;
  return oss.str();
}

TString floatToString(float number){
  ostringstream oss;
  oss << number;
  return oss.str();
}

bool isInPhiCracks(Double_t phi, Double_t eta, Double_t fiducialCut){
  
  if (TMath::Abs(eta) > 1.44) return false;

  // tranform radiants [-pi,pi] in degrees [0,360]
  phi = (phi+TMath::Pi()) *180/TMath::Pi();
  
  // each supermodule is 20 degrees wide
  Double_t moduleWidth = 20;
  
  // the first module is centered at phi=0, so the first cracks are +10 and -10  
  Double_t phi0 = 10.;
  
  bool OK = false;
  for (Int_t i = 0 ; i < 18; ++i){    
    if ((phi0 + moduleWidth*i -fiducialCut) <= phi && phi <= (phi0 + moduleWidth*i + fiducialCut)) OK = true;
    //    std::cout << " PHI " << (phi0 + moduleWidth*i -fiducialCut) << " " << phi << " " <<  (phi0 + moduleWidth*i + fiducialCut)  << " " << OK << std::endl; ;      
  }
  //  std::cout << "is in phi crack ? " << OK << std::endl;;

  return OK;
}

bool isInEtaCracks(Double_t eta){
  
  const Int_t nBinsEta = 8;
  Double_t leftEta [nBinsEta]       = {0.02, 0.46, 0.81, 1.16, 1.56};
  Double_t rightEta[nBinsEta]       = {0.42, 0.77, 1.13, 1.44, 2.5};

//   Double_t leftEta [nBinsEta]       = {0.02, 0.46, 0.81, 1.16, 1.5,   1.65,  1.9, 2.2 };
//   Double_t rightEta[nBinsEta]       = {0.42, 0.77, 1.13, 1.46, 1.65,  1.9 ,  2.2, 2.5 };

  bool OK = false;
  for (Int_t i = 0; i< nBinsEta; ++i){
    if (leftEta[i] < TMath::Abs(eta) && TMath::Abs(eta) < rightEta[i] ) OK = true;
    //    std::cout << leftEta[i] << " " << TMath::Abs(eta) << " " << rightEta[i] <<  " " << OK << std::endl;;
  }
  //  std::cout << "IS IN CRACK ? " << !OK << std::endl;;
  return !OK;
}

double f5x5( double iEta ) {
  if ( iEta < 40.2198 ) return 1;
  return 1 - 3.03103e-6*(iEta - 40.2198)*(iEta - 40.2198);
};

float getEtaCorrectionBarrel(float eta){
  return 1.0/f5x5((int)(TMath::Abs(eta)*(5/0.087)));
};
