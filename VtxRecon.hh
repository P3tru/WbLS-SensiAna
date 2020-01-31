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

/////////////////////////   ROOT   //////////////////////////
#include <TApplication.h>
#include <TMath.h>
#include <Math/QuantFuncMathCore.h>

/////////////////////////   RAT   ///////////////////////////

/////////////////////////   USER  ///////////////////////////
#include "utils.hh"

#define PI 3.14159
#define C 299.792458 // mm/ns

using namespace std;

class Hit {

 protected:
  TVector3 Pos;

  double Q;
  double T;

  double TResid;
  double D;

 public:
  Hit() : Pos(0.,0.,0.), Q(-1), T(-1), TResid(-1), D(-1){};
  Hit(TVector3 p, double q=0., double t=0., double tr=-1., double d=-1.) : Pos(p), Q(q), T(t), TResid(tr), D(d){};
  ~Hit(){};

  const TVector3 &GetPos() const {
	return Pos;
  }
  void SetPos(const TVector3 &pos) {
	Pos = pos;
  }

  double GetQ(){ return Q; };
  void SetQ(double x){ Q = x; };
  double GetT(){ return T; };
  void SetT(double x){ T = x; };

  double GetTResid() const {
	return TResid;
  }
  void SetTResid(double t_resid) {
	TResid = t_resid;
  }
  double GetD() const {
	return D;
  }
  void SetD(double d) {
	D = d;
  }

  double CalculateDistance(const TVector3& Origin = TVector3(0,0,0)){
    return TVector3(Pos-Origin).Mag();
  };
  double CalculateTResid(const TVector3& Origin = TVector3(0,0,0)){
    return T - (TVector3(Pos-Origin).Mag())/C;
  };

  bool operator==(const Hit &rhs) const {
	return T == rhs.T;
  }
  bool operator!=(const Hit &rhs) const {
	return !(rhs == *this);
  }
  bool operator<(const Hit &rhs) const {
	return T < rhs.T;
  }
  bool operator>(const Hit &rhs) const {
	return rhs < *this;
  }
  bool operator<=(const Hit &rhs) const {
	return !(rhs < *this);
  }
  bool operator>=(const Hit &rhs) const {
	return !(*this < rhs);
  }

  void Print(){
    cout << "X: " << Pos.x()
		 << " Y: " << Pos.y()
		 << " Z: " << Pos.z()
		 << " Q: " << Q
		 << " T: " << T << endl;
  }

};


class HitCut{

 public:
  Hit hCut;

  explicit HitCut(const Hit &h_cut) : hCut(h_cut) {}

  int operator()(Hit h){
	bool YesItIs = false;
	if(h>hCut){
	  YesItIs = true;
	}
	return YesItIs;
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

  const vector<Hit> &GetHCol() const {
	return hCol;
  }
  void SetHCol(const vector<Hit> &h_col) {
	hCol = h_col;
  }
  const string &GetTag() const {
	return tag;
  }
  void SetTag(const string &tag) {
	HitCollection::tag = tag;
  }

  TH2D *GethQVST() const {
	return hQVST;
  }
  TH1D *GethTResid() const {
	return hTResid;
  }
  TH1D *GethD() const {
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


double EvalL(double Nobs, double Npred){
  double L;
  if (Nobs>0 && Npred>0)L=Npred-Nobs+Nobs*TMath::Log(Nobs/Npred);
  else L=Npred;
  L=-L;
  return L;
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

	   << "\t-dbins\tSet nb bins for distance hist (int)\n"
	   << "\t-mind\tSet min val for distance hist (double)\n"
	   << "\t-maxd\tSet max val for distance hist (double)\n"

	   << "\t-mc\tinput file (ROOT)\n"
	   << "\t-pdf\tpdf file (ROOT)\n"
	   << endl;

}

void ProcessArgs(TApplication *theApp, string *filename,
				 int *User_PromtCut,
				 int *User_nEvts,
				 bool *User_batch,
				 int *User_nTResidBins, double *User_minTResid, double *User_maxTResid,
				 int *User_nDBins, double *User_minD, double *User_maxD,
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
	} else if ((arg == "-tbins")) {
	  *User_nTResidBins = stoi(theApp->Argv(++i));
	} else if ((arg == "-mint")) {
	  *User_minTResid = stod(theApp->Argv(++i));
	} else if ((arg == "-maxt")) {
	  *User_maxTResid = stod(theApp->Argv(++i));
	} else if ((arg == "-dbins")) {
	  *User_nDBins = stoi(theApp->Argv(++i));
	} else if ((arg == "-mind")) {
	  *User_minD = stod(theApp->Argv(++i));
	} else if ((arg == "-maxd")) {
	  *User_maxD = stod(theApp->Argv(++i));
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
