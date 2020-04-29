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
#define SQRT3 1.73205080757
#define SQRT2 1.41421356237

#define DIAMETER 10857
#define HEIGHT 10857
#define BUFFER 1000
#define MAXPOSFITSIDE BUFFER+DIAMETER/2
#define MAXPOSFITCAPS BUFFER+HEIGHT/2

using namespace std;


class FitTheia : public RAT::Minimizable {

 public:

  double MaxPosFitSide;
  double MaxPosFitCaps;

  TH1D *hPDF;
  TH1D *hTResidsGuess;

  TVector3 POS_TRUTH;
  double T_TRUTH;

  vector<Hit> vHit;

  vector<TVector3> POS_ALL_GUESS;
  vector<double> T_ALL_GUESS;
  vector<double> NLL_ALL_GUESS;

  TH1D *hTrue;

  TGraph *grTGuess;
  TGraph *grPosGuess[3];


  FitTheia() = default;

  FitTheia(TH1D *h_pdf,
		   const TVector3 &pos_truth, double t_truth,
		   const vector<Hit> &v_hit,
		   double MaxSide = MAXPOSFITSIDE, double MaxCaps = MAXPOSFITCAPS)
	  : POS_TRUTH(pos_truth), T_TRUTH(t_truth), vHit(v_hit), MaxPosFitSide(MaxSide), MaxPosFitCaps(MaxCaps) {

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

	if(abs(T_GUESS) > 9.)
	  return DBL_MAX;

	if(abs(params[0]) > MaxPosFitSide || abs(params[1]) > MaxPosFitSide || abs(params[2]) > MaxPosFitCaps)
	  return DBL_MAX;

	if(POS_GUESS.Mag() > SQRT3*MaxPosFitSide || POS_GUESS.Mag() > SQRT3*MaxPosFitCaps)
	  return DBL_MAX;

	hTResidsGuess->Reset();

	for (auto hit: vHit) {

	  hTResidsGuess->Fill(hit.CalculateTResid(POS_TRUTH+POS_GUESS, SoL) - (T_TRUTH + T_GUESS));

	}

	double NLL_val = 0;
	for (int iBin = 0; iBin < hPDF->GetNbinsX(); iBin++) {

	  NLL_val += EvalNLL(hTResidsGuess->GetBinContent(hTResidsGuess->GetBin(iBin)),
					   hPDF->GetBinContent(hPDF->GetBin(iBin)));
	}

	POS_ALL_GUESS.emplace_back(POS_GUESS);
	T_ALL_GUESS.emplace_back(T_GUESS);
	NLL_ALL_GUESS.emplace_back(-NLL_val);

	return -NLL_val;

  }

  void CreateTrueTRes(unsigned int iEvt = 0){

	auto nBins = hPDF->GetNbinsX();
	auto minTResids = hPDF->GetXaxis()->GetXmin();
	auto maxTResids = hPDF->GetXaxis()->GetXmax();

	hTrue = new TH1D(Form("hTrueEvt%d", iEvt), "",
					 nBins, minTResids, maxTResids);

	for (auto hit: vHit) {

	  hTrue->Fill(hit.CalculateTResid(POS_TRUTH, SoL) - (T_TRUTH));

	}

  }

  void CreateDebugPlots(unsigned int iEvt, TFile *foutput = NULL){

	grTGuess = new TGraph();
	grTGuess->SetName(Form("grTEvt%d", iEvt));

	vector<const char*> vcPos = {"X", "Y", "Z"};

	for(int iPos = 0; iPos<3; iPos++){

	  grPosGuess[iPos] = new TGraph();
	  grPosGuess[iPos]->SetName(Form("grPos%sEvt%d", vcPos[iPos], iEvt));

	}

	if(NLL_ALL_GUESS.size()>0){

	  for(auto i=0; i<NLL_ALL_GUESS.size(); i++){

	    grTGuess->SetPoint(grTGuess->GetN(), T_ALL_GUESS[i], NLL_ALL_GUESS[i]);

	    for(auto iPos=0; iPos<3; iPos++){

	      grPosGuess[iPos]->SetPoint(i,POS_ALL_GUESS[i](iPos), NLL_ALL_GUESS[i]);

	    } // END for iPos

	  } // END for i

	} // END if size>0

	if(foutput){
	  foutput->cd();
	  grTGuess->Write();
	  for(int iPos = 0; iPos<3; iPos++) {

	    grPosGuess[iPos]->Write();

	  }
	}

  }

  void DumpEvtInTTree(vector<TVector3> *vPos, vector<double> *vT, vector<double> *vNLL, TTree *tree = NULL){

    *vPos = POS_ALL_GUESS;
    *vT = T_ALL_GUESS;
    *vNLL = NLL_ALL_GUESS;

	if(tree)
	  tree->Fill();

  }

};

#endif

