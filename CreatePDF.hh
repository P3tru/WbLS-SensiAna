//
// Created by zsoldos on 2/21/20.
//

#ifndef _CREATEPDF_HH_
#define _CREATEPDF_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

/////////////////////////   BOOST  ///////////////////////////
#include <boost/algorithm/string/predicate.hpp>

/////////////////////////   ROOT   ///////////////////////////
#include <TH1D.h>
#include <TH2D.h>

/////////////////////////   USER  ///////////////////////////
#include "utils.hh"

#include "Analyzer.hh"
#include "HitClass.hh"
#include "EVFunctions.hh"
#include "HitFunctions.hh"
#include "MCFunctions.hh"
#include "AnalyzerFunctions.hh"
#include "FlatParticle.hh"

#include "ProgressBar.hpp"

int GenPromptCut(){
  static int i = -20;
  return std::abs(i++);
}

static double GenCTBin(){

  static int iCTBin = 0;
  static double CTBin = 0;
  return cos(CTBin + (10 * iCTBin++));

}

// #define SOL 224.6 // WATER
#define SOL 216.5 // WbLS 1%

void LoopAndFillHistos(Analyzer *FileAnalyzer,
					   unsigned long nEvts, unsigned long iEvt,
					   int PromptWindow, int PromptCut,
					   TH1D *hNHits, TH1D *hQ,
					   TH1D *hTRes, TH1D* hCTheta, TH2D *hTResVSCT,
					   bool isScaling = true,
					   vector<TH1D*> vhTResiduals = vector<TH1D*>(), vector<int> vPromptCut = vector<int>()){

  auto nEvtToProcess = nEvts-iEvt;
  auto nEvtProcessed = 0;

  bool isVCutDefined = !vhTResiduals.empty() && !vPromptCut.empty();
  bool isVCutOK = vhTResiduals.size() == vPromptCut.size();
  auto nCuts = isVCutDefined && isVCutOK ? vhTResiduals.size() : 0;
  vector<int> NEvtsProcessedCut(nCuts, 0);

  ProgressBar progressBar(nEvtToProcess, 70);

  for(iEvt; iEvt<nEvts; iEvt++){

	if(EoF) break;

	// record the tick
	++progressBar;

	// record Prim particle info
	FlatParticle PrimPart = GetPrimaryParticleInfo(FileAnalyzer, iEvt);
	double PrimPartKE = PrimPart.GetKinE();
	const TVector3& TrueOrigin = PrimPart.GetPos();
	auto TrueDir = PrimPart.GetDir().Unit();

	// recover EV info
	double NHits, Q;
	GetNPEAndNHitsFromEV(FileAnalyzer, iEvt, 0, &NHits, &Q);
	hNHits->Fill(NHits);
	hQ->Fill(Q);

	// Recover Hit vector for 1 evt
	// vector<Hit> vHit;
	// vector<Hit> vHit = GetEVHitCollection(FileAnalyzer, iEvt);
	vector<Hit> vHit = GetMCHitCollection(FileAnalyzer, iEvt);
	SortVHits(&vHit);

	// Apply cuts if defined
	// DAQ
	if(PromptWindow>0){
	  Hit hCut(TVector3(0, 0, 0),
			   0,
			   PromptWindow + vHit.begin()->GetT());
	  RemoveHitsAfterCut(vHit, hCut);
	}

	// Prompt evt
	if(PromptCut>0){
	  Hit hCut(TVector3(0, 0, 0),
			   0,
			   PromptCut + vHit.begin()->GetT());
	  RemoveHitsAfterCut(vHit, hCut);
	}

	// Fill histograms
	if(!vHit.empty()){

	  FillAll(hTRes, hCTheta, hTResVSCT,
			  vHit,
			  TrueOrigin,
			  SOL,
			  TrueDir);

	  nEvtProcessed++;

	}

	if(isVCutDefined && isVCutOK){

	  for(auto it=vPromptCut.begin(); it!=vPromptCut.end(); it++){

		auto pos = distance(vPromptCut.begin(), it);

		SortVHits(&vHit);

		if(!vHit.empty()) {
		  RemoveHitsAfterCut(vHit, Hit(TVector3(), 0., vHit[0].GetT() + *it));
		}

		if(!vHit.empty()) {

		  FillResiduals(vhTResiduals[pos], vHit,
						TrueOrigin,
						SOL,
						false);

		  NEvtsProcessedCut[pos]++;

		}

	  }

	}


	// display the bar
	progressBar.display();

  } // END FOR iEVT

  if(isScaling){

	// #### #### #### #### #### #### #### #### #### #### #### #### //
	// ####                      RESCALE                      #### //
	// #### #### #### #### #### #### #### #### #### #### #### #### //

	hTRes->Scale(1/(double)(nEvtProcessed));
	hCTheta->Scale(1/(double)(nEvtProcessed));
	if(isVCutDefined && isVCutOK) {
	  for (auto iCut = 0; iCut < nCuts; iCut++) {
		vhTResiduals[iCut]->Scale(1 / (double) (NEvtsProcessedCut[iCut]));
	  }
	}
	hTResVSCT->Scale(1/(double)(nEvtProcessed));

  }

}

void ShowUsage(string name){

  cerr << "Usage: " << name << " <option(s)> -mc FILE.root" << endl
	   << "Options:\n"

	   << "\t-h\tShow this help message\n"

	   << "\t-DAQ\tSet a DAQ Window for the evts (int)\n"
	   << "\t-pc\tSet prompt cut from first hit (int)\n"

	   << "\t-b\tSet batch mode\n"

	   << "\t-NEvts\tNEvts to process (int)\n"
	   << "\t-iEvt\tStart at Evt #i (int)\n"

	   << "\t-TResBins\tSet nb bins for tresid hist (int)\n"
	   << "\t-TResMin\tSet min val for tresid hist (double)\n"
	   << "\t-TResMax\tSet max val for tresid hist (double)\n"

	   << "\t-CosTBins\tSet nb bins for tresid hist (int)\n"
	   << "\t-CosTMin\tSet min val for tresid hist (double)\n"
	   << "\t-CosTMax\tSet max val for tresid hist (double)\n"

	   << "\t-xx\tSet X for true origin vector (double)\n"
	   << "\t-yy\tSet Y for true origin vector (double)\n"
	   << "\t-zz\tSet Z for true origin vector (double)\n"

	   << "\t-dir-xx\tSet X for true direction vector (double)\n"
	   << "\t-dir-yy\tSet Y for true direction vector (double)\n"
	   << "\t-dir-zz\tSet Z for true direction vector (double)\n"

	   << "\t-mc\tinput file (ROOT)\n"
	   << "\t-txt\tinput file list (.txt)\n"
	   << "\t-o\toutput file (ROOT)\n"

	   << endl;

}


void ProcessArgs(TApplication *theApp,
				 string *filename, string *listname,
				 int *User_PromptWindow,
				 int *User_PromptCut,
				 int *User_nEvts, int *User_iEvt,
				 int *User_nTResidBins, double *User_minTResid, double *User_maxTResid,
				 int *User_nCThetaBins, double *User_minCTheta, double *User_maxCTheta,
				 bool *User_isBatch,
				 string *outputname) {

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

	} else if (boost::iequals(arg, "-DAQ")) {
	  *User_PromptWindow = stoi(theApp->Argv(++i));
	} else if (boost::iequals(arg, "-pc")) {
	  *User_PromptCut = stoi(theApp->Argv(++i));

	} else if (boost::iequals(arg, "-NEvts")) {
	  *User_nEvts = stoi(theApp->Argv(++i));
	} else if (boost::iequals(arg, "-iEvt")) {
	  *User_iEvt = stoi(theApp->Argv(++i));

	} else if (boost::iequals(arg,"-TResBins")) {
	  *User_nTResidBins = stoi(theApp->Argv(++i));
	} else if (boost::iequals(arg,"-TResMin")) {
	  *User_minTResid = stod(theApp->Argv(++i));
	} else if (boost::iequals(arg, "-TResMax")) {
	  *User_maxTResid = stod(theApp->Argv(++i));

	} else if (boost::iequals(arg, "-CosTBins")) {
	  *User_nCThetaBins = stoi(theApp->Argv(++i));
	} else if (boost::iequals(arg, "-CosTMin")) {
	  *User_minCTheta = stod(theApp->Argv(++i));
	} else if (boost::iequals(arg, "-CosTMax")) {
	  *User_maxCTheta = stod(theApp->Argv(++i));

	} else if (boost::iequals(arg, "-b")) {
	  *User_isBatch=true;

	} else if (boost::iequals(arg,"-mc")) {
	  *filename = theApp->Argv(++i);

	} else if (boost::iequals(arg,"-txt")) {
	  *listname = theApp->Argv(++i);

	} else if (boost::iequals(arg,"-o")) {
	  *outputname = theApp->Argv(++i);

	} else {
	  cout << "Unkown parameter" << endl;
	  continue;
	}
  }

  if(filename->empty() && listname->empty()){
	cout << "ERROR: No input file provided!" << endl;
	exit(EXIT_FAILURE);
	} else if(!IsFileExist(*filename) && !IsFileExist(*listname)){
	cout << "ERROR: file doesn't exist!!" << endl;
	exit(EXIT_FAILURE);
  }

}

#endif //_CREATEPDF_HH_
