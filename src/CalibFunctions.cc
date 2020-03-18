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

double ComputeECalib(const MCCalib& CalibObj, double NPE, double NHits){

  if( NPE>CalibObj.GetMaxPe() || NHits>CalibObj.GetMaxHits() ){
    return -1.;
  }

  auto *grL = new TGraph();

  for(auto EBin: CalibObj.GetEBins()){

	double muPE, sigPE;
	double muHits, sigHits;

	muPE = CalibObj.GetGrMuPe()->Eval(EBin, 0 ,"S");
	sigPE = CalibObj.GetGrSigPe()->Eval(EBin, 0 ,"S");

	muHits = CalibObj.GetGrMuHits()->Eval(EBin, 0 ,"S");
	sigHits = CalibObj.GetGrSigHits()->Eval(EBin, 0 ,"S");

	double Likelihood = TMath::Gaus(NPE, muPE, sigPE)*TMath::Gaus(NHits, muHits, sigHits);
	grL->SetPoint(grL->GetN(), EBin, Likelihood);

  }

  TF1 *fFit;
  grL->Fit("gaus", "Q0");
  fFit = grL->GetFunction("gaus");
  double mean = fFit->GetParameter(1);
  double meanErr = fFit->GetParError(1);
  double sig = fFit->GetParameter(2);
  double sigErr = fFit->GetParError(2);

  // cout << "ERec: " << mean << " +- " << sig << " MeV " << endl;

  return mean;

}
