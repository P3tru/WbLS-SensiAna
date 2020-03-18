//
// Created by zsoldos on 2/24/20.
//

#include <TMath.h>

#include "CalibFunctions.hh"

static double GenERes(){

  static int iBin=1;
  static double EBin = 0.1;
  return EBin*iBin++;

}

double ComputeLikelihood(MCCalib CalibObj, double NPE, double NHits){

  vector<double> ERess(100);
  generate(ERess.begin(), ERess.end(), GenERes);

  auto *grL = new TGraph();

  for(auto ERes: ERess){

	double muPE, sigPE;
	double muHits, sigHits;

	muPE = CalibObj.GetGrMuPe()->Eval(ERes, 0 ,"S");
	sigPE = CalibObj.GetGrSigPe()->Eval(ERes, 0 ,"S");

	muHits = CalibObj.GetGrMuHits()->Eval(ERes, 0 ,"S");
	sigHits = CalibObj.GetGrSigHits()->Eval(ERes, 0 ,"S");

	double Likelihood = TMath::Gaus(NPE, muPE, sigPE)*TMath::Gaus(NHits, muHits, sigHits);
	grL->SetPoint(grL->GetN(), ERes, Likelihood);

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
