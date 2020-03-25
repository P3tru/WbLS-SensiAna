#ifndef _FITPOSTIMECUSTOM_HH_
#define _FITPOSTIMECUSTOM_HH_

///////////////////////// STL C/C++ /////////////////////////

/////////////////////////   BOOST   /////////////////////////

/////////////////////////   ROOT   //////////////////////////
#include <TVector3.h>
#include <TH1D.h>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/SimulatedAnnealing.hh>


/////////////////////////   USER   ///////////////////////////
#include "HitClass.hh"
#include "LL.hh"

#define SoL 224.9 // mm/ns
#define SQRT5 2.2360679775
#define SQRT2 1.41421356237

using namespace std;


class FitTheia : public RAT::Minimizable {

 public:

  TH1D *hPDF;
  TH1D *hTResidsGuess;

  TVector3 POS_TRUTH;
  double T_TRUTH;

  vector<Hit> vHit;

  vector<TVector3> POS_ALL_GUESS;
  vector<double> T_ALL_GUESS;
  vector<double> NLL_ALL_GUESS;


  FitTheia() = default;

  FitTheia(TH1D *h_pdf, const TVector3 &pos_truth, double t_truth, const vector<Hit> &v_hit)
	  : POS_TRUTH(pos_truth), T_TRUTH(t_truth), vHit(v_hit) {

	hPDF = (TH1D *) h_pdf->Clone("hPDF");

	auto nBins = hPDF->GetNbinsX();
	auto minTResids = hPDF->GetXaxis()->GetXmin();
	auto maxTResids = hPDF->GetXaxis()->GetXmax();

	hTResidsGuess = new TH1D("hTResid", "",
							 nBins, minTResids, maxTResids);

  }

  virtual ~FitTheia() {

	delete hTResidsGuess;
	delete hPDF;

  };

  double operator()(double *params) { // Minimizable

	TVector3 POS_GUESS(params[0], params[1], params[2]);
	double T_GUESS(params[3]);

	if(T_GUESS < hPDF->GetXaxis()->GetXmin() || T_GUESS > hPDF->GetXaxis()->GetXmax())
	  return DBL_MAX;

	hTResidsGuess->Reset();

	for (auto hit: vHit) {

	  hTResidsGuess->Fill(hit.CalculateTResid(POS_TRUTH+POS_GUESS, SoL) - T_GUESS);

	}

	double NLL_val = 0;
	for (int iBin = 0; iBin < hPDF->GetNbinsX(); iBin++) {

	  NLL_val += EvalNLL(hTResidsGuess->GetBinContent(hTResidsGuess->GetBin(iBin)),
					   hPDF->GetBinContent(hPDF->GetBin(iBin)));
	}

	// POS_ALL_GUESS.emplace_back(POS_GUESS);
	// T_ALL_GUESS.emplace_back(T_GUESS);
	// NLL_ALL_GUESS.emplace_back(-NLL_val);

	return -NLL_val;

  }

};

#endif

