/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 
#include "RooFit.h"

#include "Riostream.h" 

#include "PdfSigMass.h" 
#include "GBRMath.h"

#include "RooDataHist.h"
#include "TMatrixDSym.h"
#include "RooFitResult.h"
ClassImp(PdfSigMass) 

using namespace RooFit;

PdfSigMass::PdfSigMass(const char *name, const char *title, 
		     RooAbsReal& _m,
         	     RooAbsReal& _mean  ,
         	     RooAbsReal& _sigma1,
         	     RooAbsReal& _alpha1,
         	     RooAbsReal& _alpha2,
         	     RooAbsReal& _n1    ,
         	     RooAbsReal& _n2    ,
		     RooAbsReal& _mFrac,
		     RooAbsReal& _rtMassTerm
		     ) :
  RooAbsPdf(name,title), 
  m("m","m",this,_m),
  mean("mean"  , "mean"  ,this,_mean  ),
  sigma1("sigma1" , "sigma1" ,this,_sigma1 ),
  alpha1("alpha1", "alpha1",this,_alpha1),
  alpha2("alpha2", "alpha2",this,_alpha2),
  n1("n1"    , "n1"    ,this,_n1    ),
  n2("n2"    , "n2"    ,this,_n2    ),
  mFrac("mFrac","mFrac",this,_mFrac),
  rtMassTerm("rtMassTerm","rtMassTerm",this,_rtMassTerm, kFALSE, kFALSE)
{
}

PdfSigMass::PdfSigMass(const char *name, const char *title, 
		     RooAbsReal& _m,
         	     RooAbsReal& _mean  ,
         	     RooAbsReal& _sigma1,
         	     RooAbsReal& _sigma2 ,
         	     RooAbsReal& _alpha1,
         	     RooAbsReal& _alpha2,
         	     RooAbsReal& _n1    ,
         	     RooAbsReal& _n2    ,
         	     RooAbsReal& _f1rt  ,
		     RooAbsReal& _mFrac,
		     RooAbsReal& _rtMassTerm
		     ) :
  RooAbsPdf(name,title), 
  m("m","m",this,_m),
  mean("mean"  , "mean"  ,this,_mean  ),
  sigma1("sigma1" , "sigma1" ,this,_sigma1 ),
  sigma2("sigma2" , "sigma2" ,this,_sigma2 ),
  alpha1("alpha1", "alpha1",this,_alpha1),
  alpha2("alpha2", "alpha2",this,_alpha2),
  n1("n1"    , "n1"    ,this,_n1    ),
  n2("n2"    , "n2"    ,this,_n2    ),
  f1rt("f1rt"  , "f1rt"  ,this,_f1rt  ),
  mFrac("mFrac","mFrac",this,_mFrac),
  rtMassTerm("rtMassTerm","rtMassTerm",this,_rtMassTerm, kFALSE, kFALSE)
{
}



PdfSigMass::PdfSigMass(const PdfSigMass& other, const char* name) :  
  RooAbsPdf(other,name), 
  m("m",this,other.m),
  mean("mean",this,other.mean),
  sigma1("sigma1",this,other.sigma1),
  sigma2("sigma2",this,other.sigma2),
  alpha1("alpha1",this,other.alpha1),
  alpha2("alpha2",this,other.alpha2),
  n1("n1",this,other.n1),
  n2("n2",this,other.n2),
  f1rt("f1rt",this,other.f1rt),
  mFrac("mFrac",this,other.mFrac),
  rtMassTerm("rtMassTerm", this, other.rtMassTerm)
  //   myrtMassTermPdf(other.myrtMassTermPdf)
//   }
{
//   myrtMassTermPdf
}



Double_t PdfSigMass::evaluate() const 
{

//   double penalty = 1;
//   std::cout << m << "  " << width << "  " << alpha1 << std::endl;
//   rtMassTerm.arg().getVariables()->Print("v");
  
  //test sara
//   RooAbsReal & marg = (RooAbsReal&)m.arg();
//   RooArgSet massSet (marg);
//   std::cout << ((RooAbsPdf&)(rtMassTerm.arg())).getVal() << std::endl;
  double ret = ((RooAbsPdf&)(rtMassTerm.arg())).getVal();
//   RooDataHist* data = ((RooAbsPdf&)(rtMassTerm.arg())).generateBinned(RooArgSet(massSet),ExpectedData());
//   RooFitResult* res = ((RooAbsPdf&)(rtMassTerm.arg())).fitTo(*data, Save(),PrintLevel(-1),Minos(kFALSE),SumW2Error(kFALSE));
//   TMatrixDSym cov = res->covarianceMatrix();
//   cov.Invert();
//   double ret =  sqrt(cov.Determinant());  
//   double ret = ((PdfCBShape&)(rtMassTerm.arg()))->getVal();
//   double ret = ((RooAbsReal*) rtMassTerm.absArg())->getVal();
//   if (isPenalised) penalty = penTermVal()->getVal();

//   double rtMassValue = rtMassTermVal()->getVal();
//   double wtMassValue = wtMassTermVal()->getVal();
  
  
//   RooAbsReal xarg = m.arg() ;
//   RooArgSet massSet (xarg);
//   RooSetProxy massSetProxy (m);
//   RooArgSet massSet ((RooAbsReal*) m.absArg() );
//   double sara = rtMassTermPdf()->createIntegral(massSet, RooFit::Range("rangename"));
//   std::cout<<"pdfSigMass: " << m <<  ": " << rtMassValue << " : " << wtMassValue  << std::endl;
  
//   RooGaussian myg = RooGaussian("myg","myg", m, mean,sigma);

//   double ret = (effCValue  * rtMassValue + mFrac * effWValue * wtMassValue) * penalty);
//   double ret = ( rtMassValue->eval() );
//   double sara = rtMassTermPdf()->eval();
//   std::cout<<"pdfSigMass: " << m <<  ": " << ret << " : " << sara  << std::endl;
//   double ret = ( rtMassValue + mFrac  * wtMassValue) * penalty;
//   return sara;


//   PdfCBShape* rt = new PdfCBShape("cbs", "Crystal Ball shape", 
//                                   *(RooAbsReal*)m.absArg(), 
//                                    *(RooAbsReal*)mean.absArg(), 
//                                    *(RooAbsReal*)width.absArg(), 
//                                    *(RooAbsReal*)alpha1.absArg(), 
//                                    *(RooAbsReal*)n1.absArg(), 
//                                    *(RooAbsReal*)alpha2.absArg(), 
//                                    *(RooAbsReal*)n2.absArg());
   
//   std::cout<<"pdfSigMass: " << m <<  ": " << myrtMassTermPdf()->eval() << " : " << sara  << std::endl;

//      
//    std::cout<<"dbc: " << myrtMassTermPdf->eval() << std::endl;
//    std::cout<<"dbc: " <<ret << std::endl;
   return ret;

}


// Int_t PdfSigMass::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
// // // Int_t PdfSigMass::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
// {
//   if ( matchArgs(allVars,analVars,m) ){
//       return 1 ;
//   }
//   return 0 ;
// }
// 
// Double_t PdfSigMass::analyticalIntegral(Int_t code, const char* rangeName) const
// {
//   assert(code>0 && code<2) ;
//   double theIntegral = ((RooAbsPdf&)(rtMassTerm.arg())).analyticalIntegral();
// 
// 
// //     //   double sara = rtMassTermPdf()->createIntegral(m, rangeName);
// //     //   double rtMass = rtMassTerm.analyticalIntegral();
// //     //   RooAbsReal& marg = m.absArg() ; 
// //     RooAbsReal & marg = (RooAbsReal&)m.arg();
// //     
// //     RooAbsReal & rtMass = (RooAbsReal&)rtMassTerm.arg();
// //     double rtMassIntegral = ((RooAbsReal* )rtMass.createIntegral(marg, RooFit::NormSet(marg)))->getVal();
// //     
// //     RooAbsReal & wtMass = (RooAbsReal&)wtMassTerm.arg();
// //     double wtMassIntegral = ((RooAbsReal* )wtMass.createIntegral(marg, RooFit::NormSet(marg)))->getVal();
// //     theIntegral = rtMassIntegral + mFrac*wtMassIntegral  ;
// //   }
// //   
// //   
// //   
// // //   else if (code ==2){
// // //     Double_t retCT =  9./(32*3.14159265) * (
// // //   					  0.75*(1-Fl)              * intCPart[0]
// // //   					  + Fl                     * intCPart[1]
// // //   					  + 0.25*(1-Fl)            * intCPart[2]
// // //   					  - Fl                     * intCPart[3]
// // //   					  + 0.5*P1*(1-Fl)          * intCPart[4]
// // //   					  + 0.5*sqrt(Fl-Fl*Fl)*P4p * intCPart[5]
// // //   					  + sqrt(Fl-Fl*Fl)*P5p     * intCPart[6]
// // //   					  - sqrt(Fl-Fl*Fl)*P6p     * intCPart[7]
// // //   					  + 0.5*sqrt(Fl-Fl*Fl)*P8p * intCPart[8]
// // //   					  + 2*(1-Fl)*P2            * intCPart[9]
// // //   					  - P3*(1-Fl)              * intCPart[10]
// // //   					  );
// // //     
// // //     Double_t retWT =  9./(32*3.14159265) * (
// // //   					  0.75*(1-Fl)              * intWPart[0]
// // //   					  + Fl                     * intWPart[1]
// // //   					  + 0.25*(1-Fl)            * intWPart[2]
// // //   					  - Fl                     * intWPart[3]
// // //   					  + 0.5*P1*(1-Fl)          * intWPart[4]
// // //   					  + 0.5*sqrt(Fl-Fl*Fl)*P4p * intWPart[5]
// // //   					  - sqrt(Fl-Fl*Fl)*P5p     * intWPart[6]
// // //   					  - sqrt(Fl-Fl*Fl)*P6p     * intWPart[7]
// // //   					  - 0.5*sqrt(Fl-Fl*Fl)*P8p * intWPart[8]
// // //   					  - 2*(1-Fl)*P2            * intWPart[9]
// // //   					  + P3*(1-Fl)              * intWPart[10]
// // //     					  );
// // //     
// // //     
// // //     if (retCT<=0) {
// // //       if (retCT<0) std::cout<<"ERROR! Negative ct pdf integral, fake value returned"<<std::endl;
// // //       else std::cout<<"ERROR! Null ct pdf integral, fake value returned"<<std::endl;
// // //       return 1e-55;
// // //     }
// // //     if (retWT<=0) {
// // //       if (retWT<0) std::cout<<"ERROR! Negative wt pdf integral, fake value returned"<<std::endl;
// // //       else std::cout<<"ERROR! Null wt pdf integral, fake value returned"<<std::endl;
// // //       return 1e-55;
// // //     }
// // //     theIntegral = retCT + mFrac*retWT ;
// // //   }
//   return theIntegral;
// }
