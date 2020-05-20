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

using namespace std;

class Analyzer {

 protected:
  TFile *fMC;
  TTree *treeMC;
  string treeName;
  RAT::DS::Root *rds;

  RAT::DS::Run *run;
  TTree *runTree;

  unsigned long nEvts;

  string tag;

 public:
  Analyzer() = default;

  Analyzer(TFile *f_mc,
		   TTree *tree_mc,
		   const string &tree_name,
		   RAT::DS::Root *rds,
		   RAT::DS::Run *run,
		   TTree *run_tree,
		   unsigned long n_evts,
		   const string &tag)
	  : fMC(f_mc),
		treeMC(tree_mc),
		treeName(tree_name),
		rds(rds),
		run(run),
		runTree(run_tree),
		nEvts(n_evts),
		tag(tag) {}

  explicit Analyzer(const char* filename){

	SetAnalyzer(filename);

  }

  Analyzer(const char* filename, const char* treeName){

	SetAnalyzer(filename, treeName);

  }

  Analyzer(const char* filename, const char* treeName, string tag)
	  : tag(std::move(tag)) {

	SetAnalyzer(filename, treeName);

  }


  void SetAnalyzer(const char *filename, const char* treeName="T"){

	Analyzer::fMC = new TFile(filename);

	if (Analyzer::fMC->IsOpen()) {

	  ////////////////////////////////
	  // IF file is open do stuff ////
	  ////////////////////////////////

	  Analyzer::run = new RAT::DS::Run();

	  Analyzer::runTree = (TTree *) Analyzer::fMC->Get("runT");
	  Analyzer::runTree->SetBranchAddress("run", &run);
	  Analyzer::runTree->GetEntry(0);

	  Analyzer::treeMC = (TTree *) Analyzer::fMC->Get(treeName);

	  Analyzer::rds = new RAT::DS::Root();
	  Analyzer::treeMC->SetBranchAddress("ds", &rds);

	  Analyzer::nEvts = Analyzer::treeMC->GetEntries();

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

  TFile *GetfMC() const { return fMC; }
  TTree *GetTree() const { return treeMC; }
  RAT::DS::Root *GetDS() const { return rds; }

  unsigned long GetNEvts() const { return nEvts; }

  RAT::DS::Run *GetRun() const { return run; }
  TTree *GetRunTree() const { return runTree; }

  const string &GetTag() const { return tag; }
  void SetTag(const string &tag) { Analyzer::tag = tag; }

};

#endif //