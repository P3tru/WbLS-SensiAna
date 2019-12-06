//
// Created by zsoldos on 12/5/19.
//

#ifndef _TFILEANALYSIS_HH_
#define _TFILEANALYSIS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <cstdarg>
#include <string>

/////////////////////////   ROOT   //////////////////////////
#include <TFile.h>
#include <TTree.h>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/MC.hh>
#include <RAT/DS/Root.hh>

/////////////////////////   USER   //////////////////////////
#include "../ProgressBar.hpp"

using namespace std;

template <typename T>
class TFileAnalysis {

 protected:

  // Filename with full path + ROOT file
  // name, vtx position and E can be extracted from filename
  string filename;
  double eBin;

  // Hist
  T *Hist;

 public:

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // #### ####      CONSTRUCTORS and DESTRUCTORS       #### #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  explicit TFileAnalysis(const string &filename) : filename(filename) {}

  // #### #### #### #### #### #### #### #### #### //
  // #### ####     INIT HISTOGRAMS      #### #### //
  // #### #### #### #### #### #### #### #### #### //

  // This is a bit different than a regular setter
  // we need to clone an histogram instead of assigning pointer
  void SetHist(T *h){
	Hist = (T*)h->Clone();
  };


  // #### #### #### #### #### #### #### #### #### //
  // #### ####     DO ANALYSIS          #### #### //
  // #### #### #### #### #### #### #### #### #### //

  // MAIN FUNCTION
  // This basically load the TFile, address the TTree,
  // and loop over each events.
  // One gives in argument a function of MC parameters and
  // histogram to fill
  void DoAnalysis(void (*func)(RAT::DS::MC *mc, T *Hist)){

	auto *f = new TFile(filename.c_str());
	auto *tree = (TTree*) f->Get("T");

	auto *rds = new RAT::DS::Root();
	tree->SetBranchAddress("ds", &rds);

	int nEvents = tree->GetEntries();

	ProgressBar progressBar(nEvents, 70);

	for (int iEvt = 0; iEvt < nEvents; iEvt++) {
	  tree->GetEntry(iEvt);

	  // record the tick
	  ++progressBar;

	  // Access RAT MC info and the summary
	  RAT::DS::MC * mc = rds->GetMC();

	  // #### #### #### #### //
	  // Here we can now insert function which takes mc in input
	  // and then extract information and fill an histogram.

	  if(func){
		func(mc,Hist);
	  } else {
		cout << "NO ANALYSIS PROVIDED" << endl;
		exit(EXIT_FAILURE);
	  }

	  // display the bar
	  progressBar.display();

	} // END FOR iEvt

	// Prevent memory leaks
	// From file
	f->Close();

	// From TTree
	delete rds;

	delete f;

  };

  // #### #### #### #### #### #### #### #### #### //
  // #### ####    GETTERS AND SETTERS   #### #### //
  // #### #### #### #### #### #### #### #### #### //

  T *GetHist() const {
	return Hist;
  }

  const string &GetFilename() const {
	return filename;
  }
  void SetFilename(const string &filename) {
	TFileAnalysis::filename = filename;
  }

  double GetEBin() const {
	return eBin;
  }
  void SetEBin(double e_bin) {
	eBin = e_bin;
  }

};

#endif //_TFILEANALYSIS_HH_
