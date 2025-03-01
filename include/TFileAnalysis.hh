//
// Created by zsoldos on 12/5/19.
//

#ifndef _TFILEANALYSIS_HH_
#define _TFILEANALYSIS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <cstdarg>
#include <string>
#include <thread>

/////////////////////////   BOOST   /////////////////////////
#include <boost/filesystem/path.hpp>

/////////////////////////   ROOT   //////////////////////////
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/MC.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/Run.hh>
#include <TRandom3.h>

/////////////////////////   USER   //////////////////////////
#include "ProgressBar.hpp"

using namespace std;

template <typename T>
class TFileAnalysis {

 protected:

  // Filename with full path + ROOT file
  // name, vtx position and E can be extracted from filename
  string filename;
  double eBin;

  // Nb of Entries
  double nbEntries;

  // Unique ID
  double ID;

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

//	// First enable implicit multi-threading globally, so that the implicit parallelisation is on.
//	// The parameter of the call specifies the number of threads to use.
//	unsigned kThreadsSupported = std::thread::hardware_concurrency();
//	unsigned kNumthreads = (kThreadsSupported>0) ? kThreadsSupported : 1;
//	ROOT::EnableImplicitMT(kNumthreads);

	auto *f = new TFile(filename.c_str());

	if(f->IsOpen()){

	  ////////////////////////////////
	  // IF file is open do stuff ////
	  ////////////////////////////////

	  auto *tree = (TTree*) f->Get("T");

	  auto *rds = new RAT::DS::Root();
	  tree->SetBranchAddress("ds", &rds);

	  int nEvents = tree->GetEntries();

	  nbEntries = nEvents;

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

	} else {

	  ////////////////////////////////////////////////////
	  // Warn user if there's a problem with the file ////
	  ////////////////////////////////////////////////////

	  boost::filesystem::path p(filename);
	  cout << "FILE " << p.string() << " not open. Skip..." << endl;
	}


  };

  // #### #### #### #### #### #### #### #### #### //
  // #### ####     DO ANALYSIS          #### #### //
  // #### #### #### #### #### #### #### #### #### //

  void DoAnalysisEV(void (*func)(RAT::DS::EV *ev, T *Hist)){

	auto *f = new TFile(filename.c_str());

	if(f->IsOpen()){

	  ////////////////////////////////
	  // IF file is open do stuff ////
	  ////////////////////////////////

	  auto *tree = (TTree*) f->Get("T");

	  auto *rds = new RAT::DS::Root();
	  tree->SetBranchAddress("ds", &rds);

	  int nEvents = tree->GetEntries();

	  nbEntries = nEvents;

	  ProgressBar progressBar(nEvents, 70);

	  for (int iEvt = 0; iEvt < nEvents; iEvt++) {
		tree->GetEntry(iEvt);

		// record the tick
		++progressBar;

		// Access RAT MC info and the summary
		RAT::DS::EV *ev = rds->GetEV(0);

		// #### #### #### #### //
		// Here we can now insert function which takes mc in input
		// and then extract information and fill an histogram.

		if(func){
		  func(ev,Hist);
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

	} else {

	  ////////////////////////////////////////////////////
	  // Warn user if there's a problem with the file ////
	  ////////////////////////////////////////////////////

	  boost::filesystem::path p(filename);
	  cout << "FILE " << p.string() << " not open. Skip..." << endl;
	}


  }

  // #### #### #### #### #### #### #### #### #### //
  // #### ####        PMTs INFO         #### #### //
  // #### #### #### #### #### #### #### #### #### //

  void PrintPMTsInfo(){

	auto *f = new TFile(filename.c_str());

	if(f->IsOpen()) {

	  ////////////////////////////////
	  // IF file is open do stuff ////
	  ////////////////////////////////

	  auto *run = new RAT::DS::Run();

	  auto *runTree = (TTree*)f->Get("runT");
	  runTree->SetBranchAddress("run",&run);
	  runTree->GetEntry(0);

	  cout << "This MC was produced with " << run->GetPMTInfo()->GetPMTCount()
		   << " PMTs of model: " << run->GetPMTInfo()->GetModelName(0) << endl;

	  auto *Random = new TRandom3(0);

	  for(int i=0; i < 100; i++){
	    int id = Random->Uniform(run->GetPMTInfo()->GetPMTCount());
		cout << "Position of random PMT: "
			 << run->GetPMTInfo()->GetPosition(id).x() << "mm "
			 << run->GetPMTInfo()->GetPosition(id).y() << "mm "
			 << run->GetPMTInfo()->GetPosition(id).z() << "mm " << endl;
		cout << "Distance from center: "
			 << TMath::Sqrt(run->GetPMTInfo()->GetPosition(id).x()*run->GetPMTInfo()->GetPosition(id).x()+
			 run->GetPMTInfo()->GetPosition(id).y()*run->GetPMTInfo()->GetPosition(id).y()) << endl;
	  }


	  f->Close();

	}

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

  double GetID() const {
	return ID;
  }
  void SetID(double id) {
	ID = id;
  }

  double GetNbEntries() const {
	return nbEntries;
  }
  void SetNbEntries(double nb_entries) {
	nbEntries = nb_entries;
  }

};

#endif //_TFILEANALYSIS_HH_
