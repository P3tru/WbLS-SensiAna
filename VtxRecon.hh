//
// Created by zsoldos on 12/5/19.
//

#ifndef _TEMPLATEANALYSIS_HH_
#define _TEMPLATEANALYSIS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <fstream>
#include <sstream>
#include <regex>
#include <string>
#include <vector>

/////////////////////////   BOOST   /////////////////////////
#include <boost/filesystem/path.hpp>

/////////////////////////   ROOT   //////////////////////////
#include <TApplication.h>
#include <TMath.h>
#include <Math/QuantFuncMathCore.h>

/////////////////////////   RAT   ///////////////////////////

/////////////////////////   USER  ///////////////////////////
#include "utils.hh"

#include <HitClass.hh>
#include <HitFunctions.hh>

#include <LL.hh>

#include "ProgressBar.hpp"


class HitShift{

 public:
  Hit hShift;

  explicit HitShift(const Hit &h_shift) : hShift(h_shift) {}

  void operator()(Hit h){
    hShift.Print();
    h - hShift;
  }

};

class HitCollection {

 protected:

  vector<Hit> hCol;

  string tag;

  TH2D *hQVST;

  TH1D *hTResid;
  TH1D *hD;

 public:

  HitCollection(){

    hQVST = nullptr;
    hTResid = nullptr;
    hD = nullptr;

  };

  HitCollection(const vector<Hit> &h_col, const string &tag) : hCol(h_col), tag(tag) {

	hQVST = nullptr;
	hTResid = nullptr;
	hD = nullptr;

  };

  virtual ~HitCollection() {
	delete hQVST;
	delete hTResid;
	delete hD;
  }

  vector<Hit> GetHCol() {
	return hCol;
  }
  void SetHCol(const vector<Hit> &h_col) {
	hCol = h_col;
  }
  string GetTag() {
	return tag;
  }
  void SetTag(const string &tag) {
	HitCollection::tag = tag;
  }

  TH2D *GethQVST() {
	return hQVST;
  }
  TH1D *GethTResid() {
	return hTResid;
  }
  TH1D *GethD() {
	return hD;
  }

  void CreateNewHQVST(const char * name,
					  const int nBinsQ, const double minQ, const double maxQ,
					  const int nBinsT, const double minT, const double maxT){

	hQVST = new TH2D(name, "Q VS T",
					 nBinsQ, minQ, maxQ,
					 nBinsT, minT, maxT);

  }

  void FillHQVST(){
    for(auto hit: hCol){
      hQVST->Fill(hit.GetQ(),hit.GetT());
    }
  }

  void CreateHTResid(const char *name,
					 const int nBins, const double min, const double max){

	hTResid = new TH1D(name, "T residuals", nBins, min, max);
  }

  void FillHTResid(const TVector3& Origin = TVector3(0,0,0)){
	for(auto hit: hCol){
	  hTResid->Fill(hit.CalculateTResid(Origin));
	}
  }


};

void ProcessVHit(vector<Hit> vHit, vector<HitCollection> &Evts, string tag,
				 double TCut=0, double TTrig=0, bool prompt = true,
				 double *TFirstHit=NULL){

  if(prompt) {

	if (vHit.size() > 0) {

	  sort(vHit.begin(), vHit.end());

	  *TFirstHit = vHit[0].GetT();

	  Hit hCut(TVector3(0, 0, 0),
			   0,
			   TCut + vHit.begin()->GetT());
	  RemoveHitsAfterCut(vHit, hCut);

	  Evts.push_back(HitCollection(vHit, tag));

	}

  } else {

	if (vHit.size() > 0) {

	  sort(vHit.begin(), vHit.end());

	  *TFirstHit = vHit[0].GetT();

 	  Hit PreTrig(TVector3(0., 0., 0.),
				  0,
				  TTrig);
	  vector<Hit> vHit_CORRECTED = CorrectDelayedHits(vHit, PreTrig);

	  Hit hCut(TVector3(0, 0, 0),
			   0,
			   TCut + vHit_CORRECTED.begin()->GetT());
	  RemoveHitsAfterCut(vHit_CORRECTED, hCut);

	  Evts.emplace_back(HitCollection(vHit_CORRECTED, tag));

	}

  }

}




void SetPoissonErrors(TH1D *h){

  const double alpha = 1 - 0.6827;

  for (int i = 1; i < h->GetNbinsX(); ++i) {
	int N = h->GetBinContent(i);
	double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
	double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
	h->SetBinError(i,std::min(L,U));
  }

}


double errBinomial(double k, double N){

  return sqrt(k*(1-k/N))/N;

}


void SetBinomialErrors(TH1D *h) {

  const double N = h->GetEntries();

  for (int i = 0; i < h->GetNbinsX(); ++i) {
	double k = h->GetBinContent(i + 1);
	h->SetBinError(i + 1, errBinomial(k, N));
  }

}

double WalkAndEvalNLL(TVector3 X_GUESS, HitCollection itEvt,
					  TH1D *hPDF){

  const int nTResidBins = hPDF->GetNbinsX();
  const double minTResid = hPDF->GetXaxis()->GetXmin();
  const double maxTResid = hPDF->GetXaxis()->GetXmax();

  // Create TResid from this guess
  TH1D *hTResidGuess = new TH1D("hTResidGuess", "",
								nTResidBins, minTResid, maxTResid);

  for (auto itHit : itEvt.GetHCol()) {
	hTResidGuess->Fill(itHit.CalculateTResid(X_GUESS),
					   1 / (double) (itEvt.GetHCol().size()));
  }

  // Evaluate LL
  double NLL_val = 0;
  for(int iBin=0; iBin<nTResidBins; iBin++) {
	NLL_val += EvalL(hTResidGuess->GetBinContent(hTResidGuess->GetBin(iBin)),
					 hPDF->GetBinContent(hPDF->GetBin(iBin)));
  }

  // Speed up code by deleting unecessary guess now
  delete hTResidGuess;

  return NLL_val;
}



void FitTResid(TH1D *hTResidPDF, vector<HitCollection> Evts,
			   TH1D *hVtxRes[3],
			   vector<TVector3> *POS_MLL){

  cout << "FIT " << endl;

  ProgressBar progressBar(Evts.size(), 70);

  for(auto itEvt: Evts){

	// record the tick
	++progressBar;

	// Set seed
	TRandom3 r(0);

	const int nWalk = 100;

	// Set NLL array
	vector<double> NLL;
	vector<TVector3> POS_GUESS;

	// Set seed between [minGuess, maxGuess]mm
	const double minGuess = -1000; // mm
	const double maxGuess = 1000; // mm

	for(int iWalk=0; iWalk<nWalk; iWalk++){

	  TVector3 X_GUESS(r.Uniform(minGuess,maxGuess),
					   r.Uniform(minGuess,maxGuess),
					   r.Uniform(minGuess,maxGuess));

	  double NLL_val = WalkAndEvalNLL(X_GUESS, itEvt, hTResidPDF);

	  NLL.emplace_back(NLL_val);
	  POS_GUESS.emplace_back(X_GUESS);

	}

	int maxNLL = distance(NLL.begin(),max_element(NLL.begin(),NLL.end()));

	for(int iPos=0; iPos<3; iPos++){
	  hVtxRes[iPos]->Fill(POS_GUESS[maxNLL](iPos));
	}

	POS_MLL->emplace_back(POS_GUESS[maxNLL]);

	// display the bar
	progressBar.display();

  }

  cout << endl;

}

void CreateDRHist(TH1D *hDR, vector<TVector3> PROMPT, vector<TVector3> DELAYED){

  for(int iEvt = 0; iEvt<PROMPT.size(); iEvt++){
	hDR->Fill(abs((PROMPT[iEvt]-DELAYED[iEvt]).Mag()));
  }

}

void FitPos(TH1D *h, TF1 *fit){

  fit->SetParameters(h->GetBinContent(h->GetMaximumBin()),
					 h->GetMean(),
					 h->GetRMS());

}

void SetFitStyle(TH1D *h, int color){

  h->GetFunction("fGaus")->SetLineColor(color);
  h->GetFunction("fGaus")->SetLineWidth(3);
  h->GetFunction("fGaus")->SetLineStyle(2);

}

void PrintFitResult(TF1 *fit){

  double mean = fit->GetParameter(1);
  double meanErr = fit->GetParError(1);
  double sigma = fit->GetParameter(2);
  double sigmaErr = fit->GetParError(2);

  cout.precision(2);

  cout << "mean: " << mean/10 << "+-" << meanErr/10 << " cm" << endl;
  cout << "sigma: " << sigma/10 << "+-" << sigmaErr/10 << " cm" << endl;

}

void ShowUsage(string name){

  cerr << "Usage: " << name << " <option(s)> SOURCES" << endl
	   << "Options:\n"

	   << "\t-h\tShow this help message\n"

	   << "\t-b\tSet batch mode\n"

	   << "\t-pc\tSet prompt cut (int)\n"
	   << "\t-nevt\tSet nEvt to process (int) (default all evts in file)\n"

	   << "\t-tbins\tSet nb bins for tresid hist (int)\n"
	   << "\t-mint\tSet min val for tresid hist (double)\n"
	   << "\t-maxt\tSet max val for tresid hist (double)\n"

	   << "\t-mc\tinput file (ROOT)\n"
	   << "\t-pdf\tpdf file (ROOT)\n"
	   << endl;

}

void ProcessArgs(TApplication *theApp, string *filename,
				 int *User_PromtCut,
				 int *User_nEvts,
				 bool *User_batch,
				 bool *User_IsDelay,
				 int *User_nTResidBins, double *User_minTResid, double *User_maxTResid,
				 string *User_fPDF) {

  // Reading user input parameters
  if (theApp->Argc() < 2) {
	ShowUsage(theApp->Argv(0));
	exit(0);
  }

  int nFiles=0;

  for (int i = 1; i < theApp->Argc(); i++) {
	string arg = theApp->Argv(i);
	if ((arg == "-h") || (arg == "--help")) {
	  ShowUsage(theApp->Argv(0));
	  exit(0);
	} else if ((arg == "-pc")) {
	  *User_PromtCut = stoi(theApp->Argv(++i));
	} else if ((arg == "-nevt")) {
	  *User_nEvts = stoi(theApp->Argv(++i));
	} else if ((arg == "-b")) {
	  *User_batch=true;
	} else if ((arg == "-del")) {
	  *User_IsDelay=true;
	} else if ((arg == "-tbins")) {
	  *User_nTResidBins = stoi(theApp->Argv(++i));
	} else if ((arg == "-mint")) {
	  *User_minTResid = stod(theApp->Argv(++i));
	} else if ((arg == "-maxt")) {
	  *User_maxTResid = stod(theApp->Argv(++i));
	} else if ((arg == "-mc")) {
	  *filename = theApp->Argv(++i);
	} else if ((arg == "-pdf")) {
	  *User_fPDF = theApp->Argv(++i);
	} else {
	  cout << "Unkown parameter" << endl;
	  continue;
	}
  }

  if(filename->empty()){
    cout << "ERROR: No MC input file provided!" << endl;
    exit(EXIT_FAILURE);
  } else if(!IsFileExist(*filename)){
	cout << "ERROR: MC file doesn't exist!!" << endl;
	exit(EXIT_FAILURE);
  }

}



#endif //_TEMPLATEANALYSIS_HH_
