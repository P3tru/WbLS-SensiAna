//
// Created by zsoldos on 12/5/19.
//

///////////////////////// STL C/C++ /////////////////////////
#include <string>

/////////////////////////   ROOT   //////////////////////////
#include <TApplication.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2D.h>

/////////////////////////   USER   //////////////////////////
#include "TFileAnalysis.hh"
#include "TemplateAnalysis.hh"
#include "AnalysisDefinitions.hh"
#include "utils.hh"

using namespace std;

void ShowUsage(string name){

  cerr << "Usage: " << name << " <option(s)> SOURCES" << endl
	   << "Options:\n"
	   << "\t-h\tShow this help message\n"
	   << endl
	   << "\tSOURCES\tSpecify input data file\n"
	   << endl;

}

void ProcessArgs(TApplication *theApp, string *filename) {

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
	} else {
	  if (i + 1 > theApp->Argc() && nFiles == 0) {
		cout << "NO SOURCES PROVIDED !" << endl;
		ShowUsage(theApp->Argv(0));
		exit(EXIT_FAILURE);
	  } else {
		cout << "READ " << arg << endl;
		nFiles++;
		*filename = arg;
	  }
	}
  }

}


int main(int argc, char *argv[]) {

  // Create TApp
  TApplication theApp("App", &argc, argv);

  // Set basic ROOT style
  gStyle->SetOptStat(0);
  gStyle->SetPalette(56);

  // Create ifstream for input file
  string inputName;
  ProcessArgs(&theApp,&inputName);

  TFileAnalysis<TH2D> FileAnalysis(inputName);

  TH2D *hNbPEVSHits = new TH2D("hProut","Prout",
							   101,-0.5,1000.5,
							   101,-0.5,1000.5);

  FileAnalysis.SetHist(hNbPEVSHits);

  FileAnalysis.DoAnalysis(CollectPEAndHits);

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                      DRAWING                      #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *c1 = new TCanvas("c1","c1",800,600);

  FileAnalysis.GetHist()->Draw("COLZ");


  /////////////////////////
  // ...

  cout << endl;
  cout << "Hit Ctrl+C to exit" << endl;
  theApp.Run(kTRUE);

  return 0;

}