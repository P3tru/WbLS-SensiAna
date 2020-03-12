//
// Created by zsoldos on 2/24/20.
//

#ifndef _FITPOSTIME_HH_
#define _FITPOSTIME_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <vector>
#include <numeric>

/////////////////////////   BOOST   /////////////////////////

/////////////////////////   ROOT   //////////////////////////
#include <TVector3.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TMath.h>
#include <TH2D.h>

/////////////////////////   RAT   ///////////////////////////

/////////////////////////   USER   ///////////////////////////
#include "HitClass.hh"
#include "LL.hh"
// #include "../utils.hh"

using namespace std;

class FitPosTime {

 protected:

  // Truth Information
  TVector3 Pos_TRUTH;
  double T_TRUTH;

  // Set NLL array
  vector<double> NLL;
  vector<TVector3> POS_GUESS;
  vector<double> T_GUESS;

  // GRID search params
  int SizeGridPos;
  int SizeGridTime;
  int NbBinsAxis;
  int BinWidthAxis;

  // Fit results
  TVector3 Pos;
  double T;

  // Fit results (PLOTS) (DEBUG)
  TGraph *grNLL;
  TH2D *hPosVSTGuess;
  TH2D *hXVSTGuess[3];
  TH1D *hPosGuess[3];
  TH1D *hTGuess;
  TH1D *hAcptRatio;
  TH1D *hAcptRatioSuccess;

  // Collection of hit to fit
  vector<Hit> vHit;

  // Compute T offset of hit collection
  TH1D *hT;
  double TAvg;

  // PDF hit t resids
  TH1D *hPDF;

 public:

  FitPosTime(){
	NLL.clear();
	POS_GUESS.clear();
	T_GUESS.clear();

    Pos = TVector3(0.,0.,0.);
    T = -1;
    vHit.clear();

	hT = nullptr;
	TAvg = -1.;
	hPDF = nullptr;
  }

  FitPosTime(const TVector3 &pos_truth, double t_truth, const vector<Hit> &v_hit, TH1D *h_pdf)
	  : Pos_TRUTH(pos_truth), T_TRUTH(t_truth), vHit(v_hit) {

	NLL.clear();
	POS_GUESS.clear();
	T_GUESS.clear();

	hPDF = (TH1D*)h_pdf->Clone();

	Pos = TVector3(0.,0.,0.);
	T = -1;

	hT = new TH1D("hT", "hT", 1e6, 0., 1.e6);

	for(auto h:vHit){
	  hT->Fill(h.GetT());
	}

	TAvg = hT->GetMean();

	if(hPDF){
	  // FitTResid();
	  FitOnGridTResid();
	}

  }

  FitPosTime(const vector<Hit> &v_hit, TH1D *h_pdf)
	  : vHit(v_hit) {

    Pos_TRUTH = TVector3(0.,0.,0.);
    T_TRUTH = 0;

	NLL.clear();
	POS_GUESS.clear();
	T_GUESS.clear();

	hPDF = (TH1D*)h_pdf->Clone();

	Pos = TVector3(0.,0.,0.);
	T = -1;

	hT = new TH1D("hT", "hT", 1e6, 0., 1.e6);

	for(auto h:vHit){
	  hT->Fill(h.GetT());
	}

	TAvg = hT->GetMean();

	if(hPDF)
	  FitTResid();

  }

  FitPosTime(const TVector3 &pos, double t, const vector<Hit> &v_hit, TH1D *h_t, double t_avg, TH1D *h_pdf)
	  : Pos(pos), T(t), vHit(v_hit), hT(h_t), TAvg(t_avg), hPDF(h_pdf) {

  }

  virtual ~FitPosTime() {

	// delete hPDF;
	delete hT;

  }

  TH1D *CreateHTResid(){

	const int nTResidBins = hPDF->GetNbinsX();
	const double minTResid = hPDF->GetXaxis()->GetXmin();
	const double maxTResid = hPDF->GetXaxis()->GetXmax();

	// cout << "CreateHTResid()" << endl;
	// cout << "nTResidBins: " << nTResidBins << endl;
	// cout << "minTResid: " << minTResid << endl;
	// cout << "maxTResid: " <<maxTResid << endl;

	// Create TResid from this guess
	TH1D *hTResid = new TH1D("hTResid", "",
							 nTResidBins, minTResid, maxTResid);

	return hTResid;

  }

  void FillHTResid(TH1D* hTResid, const TVector3& X_SEED = TVector3(0,0,0), double T_SEED=0){

	if(hTResid->GetEntries() > 0)
	  hTResid->Reset();

	// cout << "FillHTResid()" << endl;
	for(auto hit: vHit){

	  // cout << "hit.CalculateTResid(X_SEED) - TAvg - T_SEED: "
		//    << hit.CalculateTResid(X_SEED) - TAvg - T_SEED << endl;
	  hTResid->Fill(hit.CalculateTResid(X_SEED+Pos_TRUTH) - TAvg - T_SEED,
					1 / (double) (vHit.size()));

	}

  }

  double WalkAndEvalNLL(TVector3 X_SEED, double T_SEED){

    TH1D *hTResidGuess = CreateHTResid();

	FillHTResid(hTResidGuess, X_SEED, T_SEED);

	// Evaluate LL
	double NLL_val = 0;
	for(int iBin=0; iBin<hPDF->GetNbinsX(); iBin++) {
	  NLL_val += EvalL(hTResidGuess->GetBinContent(hTResidGuess->GetBin(iBin)),
					   hPDF->GetBinContent(hPDF->GetBin(iBin)));
	}

	// Speed up code by deleting unecessary guess now
	delete hTResidGuess;

	return NLL_val;
  }

  void GenerateSeed(TVector3 *X_SEED, double *T_SEED,
					TRandom3 *r = new TRandom3(0)){

	const double sigmaPosGuess = 1000; //mm

	X_SEED->SetXYZ(r->Gaus(0, sigmaPosGuess),
				   r->Gaus(0, sigmaPosGuess),
				   r->Gaus(0, sigmaPosGuess));

	const double sigmaTGuess = 20;

	*T_SEED = r->Gaus(0, sigmaTGuess);

  }

  template <typename T>
  T GetGradientVector(vector<T> vecT){
	return vecT[vecT.size()-1] - vecT[vecT.size()-2];
  }

  template <typename T>
  T GetAcceptanceRatio(vector<T> vecT){
	return vecT[vecT.size()-1] / vecT[vecT.size()-2];
  }

  void GenerateMCMCSeed(TVector3 *X_SEED, double *T_SEED,
						TRandom3 *r = new TRandom3(0)){

	const double sigmaPosGuess = 100.; //mm
	const double sigmaTGuess = 10;

	if(NLL.size()>1){

	  // POS

	  unsigned iWalk = POS_GUESS.size()-1;
	  TVector3 X_TEST(r->Gaus(POS_GUESS[iWalk].x(), sigmaPosGuess),
					  r->Gaus(POS_GUESS[iWalk].y(), sigmaPosGuess),
					  r->Gaus(POS_GUESS[iWalk].z(), sigmaPosGuess));

	  double AcceptanceRatio = TMath::Gaus(X_TEST.x()+X_TEST.y()+X_TEST.z(),
										   Pos_TRUTH.x()+Pos_TRUTH.y()+Pos_TRUTH.z(),
										   sqrt(3)*sigmaPosGuess,
										   kTRUE);

	  double Test = r->Uniform(0,1);

	  if(AcceptanceRatio > Test){

	    *X_SEED = X_TEST;

	  } else {

	    *X_SEED = POS_GUESS[iWalk];

	  }

	  // TIME

	  iWalk = T_GUESS.size()-1;
	  double T_TEST = r->Gaus(T_GUESS[iWalk], sigmaTGuess);


	  AcceptanceRatio = TMath::Gaus(T_TEST,
									0,
									sigmaTGuess,
									kTRUE);

	  Test = r->Uniform(0,1);

	  if(AcceptanceRatio > Test){

		*T_SEED = T_TEST;

	  } else {

		*T_SEED = T_GUESS[iWalk];

	  }


	} else {

	  GenerateSeed(X_SEED, T_SEED, r);

	}

  }

  void GenerateGradDescentSeed(TVector3 *X_SEED, double *T_SEED,
							   TRandom3 *r = new TRandom3(0)){

    const double alpha = 0.5;

	if(NLL.size()>2){

	  double DNLL = GetGradientVector(NLL);
	  TVector3 DPOS = GetGradientVector(POS_GUESS);
	  double DT = GetGradientVector(T_GUESS);

	  double AcceptRatio = GetAcceptanceRatio(NLL);
	  double Test = r->Uniform(0,1);

	  hAcptRatio->Fill(AcceptRatio);

	  if(DNLL < 0){

	    // cout << "AcceptRatio: " << AcceptRatio << endl;

		*X_SEED = *X_SEED + (1-AcceptRatio)*DPOS;
		*T_SEED = *T_SEED + (1-AcceptRatio)*DT;

		hAcptRatioSuccess->Fill(AcceptRatio);

	  } else {

		GenerateSeed(X_SEED, T_SEED, r);

	  }

	} else {

	  GenerateSeed(X_SEED, T_SEED, r);

	}

  }

  void InitGridParams(){

	NbBinsAxis = 10;
	BinWidthAxis = 1000/NbBinsAxis;
	SizeGridPos = NbBinsAxis*NbBinsAxis*NbBinsAxis; // 20^3
	SizeGridTime = 50;

  }

  void InitGrid(vector<TVector3> *X_SEED, vector<double> *T_SEED){


	X_SEED->resize(SizeGridPos);
	T_SEED->resize(SizeGridTime);
	iota(T_SEED->begin(), T_SEED->end(),-SizeGridTime/2);

    for(int iX=0 ; iX<NbBinsAxis ; iX++){

	  for(int iY=0 ; iY<NbBinsAxis ; iY++) {

		for (int iZ = 0; iZ < NbBinsAxis; iZ++) {

		  unsigned index = iX*100 + iY*10 + iZ;

		  X_SEED->at(index).SetXYZ(iX*BinWidthAxis, iY*BinWidthAxis, iZ*BinWidthAxis);

		}
	  }
    }

  }

  void WalkOnGridAndEvalNLL(vector<TVector3> *X_SEED, vector<double> *T_SEED, bool isSwitch=false){

	for(auto itT: *T_SEED){

	  for(auto itX: *X_SEED){

	    if(isSwitch){
	      itX=-itX;
	    }

		double NLL_val = WalkAndEvalNLL(itX, itT);

		NLL.emplace_back(NLL_val);
		POS_GUESS.emplace_back(itX);
		T_GUESS.emplace_back(itT);

	  }

	}

  }

  void FitOnGridTResid(){

    // cout << "FIT " << endl;

    InitGridParams();

	vector<TVector3> X_SEED;
	vector<double> T_SEED;
	InitGrid(&X_SEED, &T_SEED);

	WalkOnGridAndEvalNLL(&X_SEED, &T_SEED);
	WalkOnGridAndEvalNLL(&X_SEED, &T_SEED, true);

	int minNLL = distance(NLL.begin(),min_element(NLL.begin(),NLL.end()));

	Pos = POS_GUESS[minNLL];
	T = T_GUESS[minNLL];

  }

  void FitTResid(){

	// cout << "FIT " << endl;

	// Set seed
	TRandom3 r(0);

	const int nWalk = 1e5;

	hAcptRatio = new TH1D("hAcptRatio","Acceptance Ratio",
						  1001,-0.05,10.05);
	hAcptRatioSuccess = new TH1D("hAcptRatioSuccess","Acceptance Ratio",
								 1001,-0.05,10.05);

	for(int iWalk=0; iWalk<nWalk; iWalk++){

	  TVector3 X_SEED;
	  double T_SEED;
	  // GenerateSeed(&X_SEED, &T_SEED, &r);
	  // GenerateMCMCSeed(&X_SEED, &T_SEED, &r);
	  GenerateGradDescentSeed(&X_SEED, &T_SEED, &r);

	  if(NLL.size()>1){
	    if(POS_GUESS[POS_GUESS.size()-1] == X_SEED && T_GUESS[T_GUESS.size()-1] == T_SEED)
		  continue;
	  }

	  double NLL_val = WalkAndEvalNLL(X_SEED, T_SEED);

	  NLL.emplace_back(NLL_val);
	  POS_GUESS.emplace_back(X_SEED);
	  T_GUESS.emplace_back(T_SEED);

	}

	int minNLL = distance(NLL.begin(),min_element(NLL.begin(),NLL.end()));


	Pos = POS_GUESS[minNLL];
	T = T_GUESS[minNLL];


  }

  void PlotFitResults(){

	const unsigned int nWalks = NLL.size();

    grNLL = new TGraph();

	for(int iPos = 0; iPos<3; iPos++){

	  hPosGuess[iPos] = new TH1D(Form("hPosGuess%d", iPos),
								 Form("Vtx Reconstruction Axis %d", iPos),
								 21, -1000.5, 1000.5);
	  hXVSTGuess[iPos] = new TH2D(Form("hXVST%d",iPos), "",
							  21, -1000.5, 1000.5,
							  41, -20.5, 20.5);

	}

	hTGuess = new TH1D("hTGuess", "T Reconstruction",
					   21, -100.5, 100.5);

	hPosVSTGuess = new TH2D("hPosVSTGuess", "",
							21, -0.5, 2000.5,
							21, -100.5, 100.5);

    for(unsigned int iNLL=0; iNLL<nWalks; iNLL++){
	  grNLL->SetPoint(iNLL, POS_GUESS[iNLL].Mag(), NLL[iNLL]);

	  hTGuess->Fill(T_GUESS[iNLL]);

	  hPosVSTGuess->Fill(POS_GUESS[iNLL].Mag(), T_GUESS[iNLL], -NLL[iNLL]);

	  for(int iPos = 0; iPos<3; iPos++){

		// hPosGuess[iPos]->Fill(POS_GUESS[iNLL](iPos), NLL[iNLL]);
		hPosGuess[iPos]->Fill(POS_GUESS[iNLL](iPos));
		hXVSTGuess[iPos]->Fill(POS_GUESS[iNLL](iPos), T_GUESS[iNLL], -NLL[iNLL]);
	  }

    }

    gStyle->SetPalette(kRainBow);

    auto *c1 = new TCanvas("cNLL", "cNLL", 800,600);
    c1->SetGrid();
    grNLL->Draw("AP");

    c1 = new TCanvas("cVTX", "cVTX", 800, 600);
    c1->SetGrid();
    for(auto h:hPosGuess){
      h->Draw("SAME PLC PMC");
    }

	c1 = new TCanvas("cTime", "cTime", 800, 600);
    c1->SetGrid();
    hTGuess->Draw();

	c1 = new TCanvas("cPVST", "cPVST", 800, 600);
	c1->SetGrid();
	c1->SetLogz();
	hPosVSTGuess->Draw("COLZ");

	for(auto h:hXVSTGuess){
	  c1 = new TCanvas(Form("c%s",h->GetName()),
					   Form("c%s",h->GetName()),
					   800, 600);
	  c1->SetGrid();
	  c1->SetLogz();
	  h->Draw("COLZ");
	}

	TH1D *hTResidGuess = CreateHTResid();
	FillHTResid(hTResidGuess, Pos, T);

	c1 = new TCanvas("cTResid", "cTResid", 800, 600);
	c1->SetGrid();
	hPDF->Draw();
	hTResidGuess->SetLineColor(kRed-4);
	hTResidGuess->Draw("SAME");

	c1 = new TCanvas("cTResidG", "cTResidG", 800, 600);
	c1->SetGrid();
	hTResidGuess->SetLineColor(kRed-4);
	hTResidGuess->Draw("");

	if(hAcptRatio){
	  hAcptRatio->Draw();
	  c1 = new TCanvas("cTAcc", "cAcc", 800, 600);
	  c1->SetGrid();
	  if(hAcptRatioSuccess){
		hAcptRatioSuccess->SetLineColor(kRed-4);
		hAcptRatioSuccess->Draw("SAME");
	  }
	}


  }

  const TVector3 &GetPos() const {
	return Pos;
  }
  double GetT() const {
	return T;
  }

  const vector<double> &GetNll() const {
	return NLL;
  }
  const vector<TVector3> &GetPosGuess() const {
	return POS_GUESS;
  }
  const vector<double> &GetTGuess() const {
	return T_GUESS;
  }

  TGraph *GetGrNll() const {
	return grNLL;
  }
  TH1D *const *GetHPosGuess() const {
	return hPosGuess;
  }
  TH1D *GetHtGuess() const {
	return hTGuess;
  }

  TH1D *GetHpdf() const {
	return hPDF;
  }

  double GetTAvg() const {
	return TAvg;
  }

  TH1D *GetHAcptRatio() const {
	return hAcptRatio;
  }
  TH1D *GetHAcptRatioSuccess() const {
	return hAcptRatioSuccess;
  }

  void PrintFitResults(){
	cout << " X: " << Pos.x() + Pos_TRUTH.x()
		 << " Y: " << Pos.y() + Pos_TRUTH.y()
		 << " Z: " << Pos.z() + Pos_TRUTH.z()
		 << " T: " << T
		 << " T+TAvg: " << T+TAvg << endl;
  }

};


#endif
