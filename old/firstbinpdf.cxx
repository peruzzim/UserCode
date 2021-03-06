/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "firstbinpdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 
#include <algorithm>
#include <iostream>

using namespace std;

ClassImp(firstbinpdf) 

 firstbinpdf::firstbinpdf(const char *name, const char *title, 
                        RooAbsReal& _roovar) :
   RooAbsPdf(name,title), 
   roovar("roovar","roovar",this,_roovar)
 { 
 } 


 firstbinpdf::firstbinpdf(const firstbinpdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   roovar("roovar",this,other.roovar)
 { 
 } 



 Double_t firstbinpdf::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   return  (roovar<0.1) ? 10.0 : 0.0; 
 } 



 Int_t firstbinpdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const  
 { 
   // LIST HERE OVER WHICH VARIABLES ANALYTICAL INTEGRATION IS SUPPORTED, 
   // ASSIGN A NUMERIC CODE FOR EACH SUPPORTED (SET OF) PARAMETERS 
   // THE EXAMPLE BELOW ASSIGNS CODE 1 TO INTEGRATION OVER VARIABLE X
   // YOU CAN ALSO IMPLEMENT MORE THAN ONE ANALYTICAL INTEGRAL BY REPEATING THE matchArgs 
   // EXPRESSION MULTIPLE TIMES

   if (matchArgs(allVars,analVars,roovar)) return 1 ; 
   return 0 ; 
 } 



 Double_t firstbinpdf::analyticalIntegral(Int_t code, const char* rangeName) const  
 { 
   // RETURN ANALYTICAL INTEGRAL DEFINED BY RETURN CODE ASSIGNED BY getAnalyticalIntegral
   // THE MEMBER FUNCTION x.min(rangeName) AND x.max(rangeName) WILL RETURN THE INTEGRATION
   // BOUNDARIES FOR EACH OBSERVABLE x

    assert(code==1) ; 
    float up = std::min(10.0*roovar.max(rangeName),1.0);
    float down = std::min(10.0*roovar.min(rangeName),1.0);
    //    std::cout << "returning " << roovar.min(rangeName) << " " << roovar.max(rangeName) << " " <<  up-down << std::endl;
    return up-down;
   // return (x.max(rangeName)-x.min(rangeName)) ; 
   //return 0 ; 
 } 



