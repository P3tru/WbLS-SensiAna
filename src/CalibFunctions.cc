//
// Created by zsoldos on 2/24/20.
//

#include <TMath.h>

#include "CalibFunctions.hh"

static double GenEBin(){

  static int iBin=1;
  static double EBin = 0.1;
  return EBin*iBin++;

}

double ComputeECalib(const MCCalib& CalibObj, double NPE, double NHits, TGraph *grL){

  if( NPE>CalibObj.GetMaxPe() || NHits>CalibObj.GetMaxHits() ){
    return -1.;
  }

  if(!grL)
    grL = new TGraph();

  auto MaxL = DBL_MIN;
  auto MaxE = DBL_MIN;

  for(auto EBin: CalibObj.GetEBins()){

	double muPE, sigPE;
	double muHits, sigHits;

	muPE = CalibObj.GetGrMuPe()->Eval(EBin, 0 ,"S");
	sigPE = CalibObj.GetGrSigPe()->Eval(EBin, 0 ,"S");

	muHits = CalibObj.GetGrMuHits()->Eval(EBin, 0 ,"S");
	sigHits = CalibObj.GetGrSigHits()->Eval(EBin, 0 ,"S");

	double Likelihood = TMath::Gaus(NPE, muPE, sigPE)*TMath::Gaus(NHits, muHits, sigHits);
	grL->SetPoint(grL->GetN(), EBin, Likelihood);

	if(Likelihood>MaxL){
	  MaxL = Likelihood;
	  MaxE = EBin;
	}

  }

  return MaxE;

}
