//
// Created by zsoldos on 2/21/20.
//

#ifndef _ANALYZER_HH_
#define _ANALYZER_HH_

/////////////////////////   ROOT   //////////////////////////
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/Run.hh>
#include <RAT/DS/Root.hh>


class Analyzer {

 protected:
  TFile *fMC;
  TTree *treeMC;
  RAT::DS::Root *rds;

  RAT::DS::Run *run;
  TTree *runTree;

  unsigned long int nEvts;

 public:
  explicit Analyzer(TFile *f_mc= nullptr, TTree *tree_mc= nullptr, RAT::DS::Root *rds= nullptr,
					unsigned long n_evts=0)
	  : fMC(f_mc), treeMC(tree_mc), rds(rds), nEvts(n_evts) {

  }
  explicit Analyzer(const char* filename){
	fMC = new TFile(filename);

	if (fMC->IsOpen()) {

	  ////////////////////////////////
	  // IF file is open do stuff ////
	  ////////////////////////////////

	  run = new RAT::DS::Run();

	  runTree = (TTree *) fMC->Get("runT");
	  runTree->SetBranchAddress("run", &run);
	  runTree->GetEntry(0);

	  treeMC = (TTree *) fMC->Get("T");

	  rds = new RAT::DS::Root();
	  treeMC->SetBranchAddress("ds", &rds);

	  nEvts = treeMC->GetEntries();

	}

  }

  virtual ~Analyzer() {
	fMC->Close();
	delete fMC;

	delete treeMC;
	delete rds;

	delete runTree;
	delete run;

  }

  TFile *GetFmc() const {
	return fMC;
  }
  TTree *GetTreeMc() const {
	return treeMC;
  }
  RAT::DS::Root *GetRds() const {
	return rds;
  }

  unsigned long GetNEvts() const {
	return nEvts;
  }

  RAT::DS::Run *GetRun() const {
	return run;
  }
  TTree *GetRunTree() const {
	return runTree;
  }

};

#endif //