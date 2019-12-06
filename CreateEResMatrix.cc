//
// Created by zsoldos on 11/26/19.
//

///////////////////////// STL C/C++ /////////////////////////
#include <vector>

/////////////////////////   ROOT   //////////////////////////
#include <TApplication.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>

/////////////////////////   USER   //////////////////////////
#include "CreateEResMatrix.hh"
#include "TFileAnalysis.hh"

using namespace std;

void ShowUsage(string name){

  cerr << "Usage: " << name << " <option(s)> SOURCES" << endl
			<< "Options:\n"
			<< "\t-h\tShow this help message\n"
			<< endl
			<< "\tSOURCES\tSpecify input data file (.txt)\n"
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

  // Create ifstream for input file
  string inputName;
  ProcessArgs(&theApp,&inputName);

  string line;
  ifstream file(inputName);

  vector< TFileAnalysis<TH2D> > FileAnalysis;
  int NbFileAnalysis = 0;

  while ( getline(file,line) ){
	FileAnalysis.emplace_back(line);
	FileAnalysis[NbFileAnalysis].SetEBin(ExtractEBinFromFilename(line));

	cout << "ADD " << ExtractFilenameFromPath(line) << endl;
	cout << "-> Corresponding to Ebin: " << ExtractEBinFromFilename(line) << endl;

	NbFileAnalysis++;
  }

  TH2D *hNbPEVSHits;

  vector<double> Ebins;
  vector<double> Eres;

  TCanvas *c1;

  for(auto& file : FileAnalysis){

	hNbPEVSHits = new TH2D(Form("hEbin%.1f", file.GetEBin()),"Prout",
						   1001,-0.5,1000.5,
						   1001,-0.5,1000.5);

	file.SetHist(hNbPEVSHits);

	file.DoAnalysis(CollectPEAndHits);

	///////////////////////////////
	// FIT and recover res. info //
	///////////////////////////////

	TH1D *hNbPE = (TH1D*)file.GetHist()->ProjectionX()->Clone();
	hNbPE->Fit("gaus","0LQSEM+");
	TF1 *fFit = hNbPE->GetFunction("gaus");
	double chi2 = fFit->GetChisquare();
	double c = fFit->GetParameter(0);
	double mean = fFit->GetParameter(1);
	double meanErr = fFit->GetParError(1);
	double sigma = fFit->GetParameter(2);
	double sigmaErr = fFit->GetParError(2);

	Ebins.emplace_back(file.GetEBin());
	Eres.emplace_back(sigma/mean);

	c1 = new TCanvas(Form("cEbin%f", file.GetEBin()),Form("cEbin%f", file.GetEBin()),
					 800,600);
	c1->SetGrid();

	file.GetHist()->Draw("COLZ");

  }

  TH2D *hEMatrix = new TH2D();


  /////////////////////////
  // ...

  cout << endl;
  cout << "Hit Ctrl+C to exit" << endl;
  theApp.Run(kTRUE);

  return 0;
}
