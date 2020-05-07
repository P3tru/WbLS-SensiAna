//
// Created by zsoldos on 2/21/20.
//

#ifndef _SELECTIBD_HH_
#define _SELECTIBD_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

/////////////////////////   BOOST   /////////////////////////

/////////////////////////   ROOT   //////////////////////////
#include <TVector3.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TFile.h>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/SimulatedAnnealing.hh>

/////////////////////////   USER   ///////////////////////////
#include "utils.hh"
#include "HitClass.hh"
#include "LL.hh"
#include "Analyzer.hh"
#include "HitClass.hh"
#include "EVFunctions.hh"
#include "HitFunctions.hh"
#include "ProgressBar.hpp"

#define SoL 224.9 // mm/ns
#define SQRT5 2.2360679775
#define SQRT2 1.41421356237


TH1D *GetHPDF(const char *filename, const char *histname) {

  TFile *fPDF = new TFile(filename);
  TH1D *hPDF = nullptr;

  if(fPDF){

	hPDF = (TH1D *) fPDF->Get(histname)->Clone();

  }

  return hPDF;

}

void ScaleSimplexCoord(vector<double> *seed){

  const double scale = 500. / (2*SQRT2);

  for(auto &s: *seed){

	s *= scale;

  }

}

double CalculateProb(TH1D *hPDF, TH1D *hExp, double *Chi2 = NULL, int *NdF = NULL){

  auto nBins = hPDF->GetNbinsX();

  auto N = hExp->Integral();
  auto W = hPDF->GetSumOfWeights();
  auto normW = 2*W*W;

  double chi2 = 0.;
  int NonNullBin = 0;

  for(int iBin = 1; iBin<=nBins; iBin++){

	double n = hExp->GetBinContent(iBin);
	double w = hPDF->GetBinContent(iBin);

	if(n == 0 || w == 0) continue;

	double s2 = pow(hPDF->GetBinError(iBin),2);
	double res = W*w - N*s2;

	double P = res + sqrt(pow(res, 2) + 4*W*W*s2*s2*n);
	if(P == 0) continue;
	P/=normW;

	chi2+=pow(n - N*P,2)/(N*P) + pow(w - W*P,2)/s2;
	NonNullBin++;

  }

  *Chi2=chi2;
  *NdF=NonNullBin-1;

  return TMath::Prob(chi2, NonNullBin-1);


}

TH1D *GetBestFitHTResids(const vector<double> vBP, TH1D *hPDF,
						 TVector3 POS_TRUTH, double T_TRUTH,
						 const vector<Hit> vHit,
						 unsigned iEvt = 0){

  TVector3 POS_GUESS(vBP[0], vBP[1], vBP[2]);
  double T_GUESS(vBP[3]);

  auto nBins = hPDF->GetNbinsX();
  auto minTResids = hPDF->GetXaxis()->GetXmin();
  auto maxTResids = hPDF->GetXaxis()->GetXmax();

  auto *hTResidsGuess = new TH1D(Form("hBFEvt#%d", iEvt), "",
								 nBins, minTResids, maxTResids);
  hTResidsGuess->Reset();

  for (auto hit: vHit) {

	hTResidsGuess->Fill(hit.CalculateTResid(POS_TRUTH+POS_GUESS, SoL) - (T_TRUTH + T_GUESS));

  }

  return hTResidsGuess;
}

void ShowUsage(string name){

  cerr << "Usage: " << name << " <option(s)> -mc FILE.root" << endl
	   << "Options:\n"

	   << "\t-h\tShow this help message\n"

	   << "\t-pc\tSet prompt cut (int)\n"
	   << "\t-DAQT\tSet DAQ Window (int)\n"

	   << "\t-nevts\tNEvts to process (int)\n"
	   << "\t-ievt\tStart at Evt #i (int)\n"

	   << "\t-minEprompt\tSet min E for prompt window (double)\n"
	   << "\t-maxEprompt\tSet max E for prompt window (double)\n"

	   << "\t-minE-nH\tSet min E for delay window n-H (double)\n"
	   << "\t-maxE-nH\tSet max E for delay window n-H (double)\n"

	   << "\t-minE-nC\tSet min E for delay window n-C (double)\n"
	   << "\t-maxE-nC\tSet max E for delay window n-C (double)\n"

	   << "\t-minE-nO\tSet min E for delay window n-O (double)\n"
	   << "\t-maxE-nO\tSet max E for delay window n-O (double)\n"

	   << "\t-xx\tSet X for true origin vector (double)\n"
	   << "\t-yy\tSet Y for true origin vector (double)\n"
	   << "\t-zz\tSet Z for true origin vector (double)\n"

	   << "\t-tt\tSet T for true time (int)\n"

	   << "\t-DR\tSet DR cutoff (double)\n"

	   << "\t-DT\tSet DT cutoff (double)\n"


	   << "\t-b\tSet batch mode\n"

	   << "\t-calib\tcalib file (ROOT)\n"
	   << "\t-mc\tinput file (ROOT)\n"
	   << "\t-pdf\tpdf file (ROOT)\n"

	   << endl;

}


void ProcessArgs(TApplication *theApp, string *filename,
				 int *User_PromptCut,
				 int *User_DAQWindow,
				 int *User_nEvts, int *User_iEvt,
				 double *User_minEPrompt, double *User_maxEPrompt,
				 double *User_MinEDelayed_nH, double *User_MaxEDelayed_nH,
				 double *User_MinEDelayed_nC, double *User_MaxEDelayed_nC,
				 double *User_MinEDelayed_nO, double *User_MaxEDelayed_nO,
				 double *User_xx, double *User_yy, double *User_zz,
				 double *User_tt,
				 double *User_DR,
				 double *User_DT,
				 double *User_E,
				 double *User_ManCalib,
				 string *User_Calib,
				 bool *User_isBatch,
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
	  *User_PromptCut = stoi(theApp->Argv(++i));

	} else if ((arg == "-DAQT")) {
	  *User_DAQWindow = stoi(theApp->Argv(++i));

	} else if ((arg == "-nevts")) {
	  *User_nEvts = stoi(theApp->Argv(++i));
	} else if ((arg == "-ievt")) {
	  *User_iEvt = stoi(theApp->Argv(++i));

	} else if ((arg == "-minEprompt")) {
	  *User_minEPrompt = stod(theApp->Argv(++i));
	} else if ((arg == "-minEprompt")) {
	  *User_minEPrompt = stod(theApp->Argv(++i));

	} else if ((arg == "-minE-nH")) {
	  *User_MinEDelayed_nH = stod(theApp->Argv(++i));
	} else if ((arg == "-maxE-nH")) {
	  *User_MaxEDelayed_nH = stod(theApp->Argv(++i));

	} else if ((arg == "-minE-nC")) {
	  *User_MinEDelayed_nC = stod(theApp->Argv(++i));
	} else if ((arg == "-maxE-nC")) {
	  *User_MaxEDelayed_nC = stod(theApp->Argv(++i));

	} else if ((arg == "-minE-nO")) {
	  *User_MinEDelayed_nO = stod(theApp->Argv(++i));
	} else if ((arg == "-maxE-nO")) {
	  *User_MaxEDelayed_nO = stod(theApp->Argv(++i));

	} else if ((arg == "-xx")) {
	  *User_xx = stod(theApp->Argv(++i));
	} else if ((arg == "-yy")) {
	  *User_yy = stod(theApp->Argv(++i));
	} else if ((arg == "-zz")) {
	  *User_zz = stod(theApp->Argv(++i));

	} else if ((arg == "-tt")) {
	  *User_tt = stod(theApp->Argv(++i));

	} else if ((arg == "-DR")) {
	  *User_DR = stod(theApp->Argv(++i));
	} else if ((arg == "-DT")) {
	  *User_DT = stod(theApp->Argv(++i));

	} else if ((arg == "-E")) {
	  *User_E = stod(theApp->Argv(++i));

	} else if ((arg == "-mcalib")) {
	  *User_ManCalib = stod(theApp->Argv(++i));

	} else if ((arg == "-b")) {
	  *User_isBatch=true;
	} else if ((arg == "-mc")) {
	  *filename = theApp->Argv(++i);
	} else if ((arg == "-calib")) {
	  *User_Calib = theApp->Argv(++i);
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

#endif //_SELECTIBD_HH_
