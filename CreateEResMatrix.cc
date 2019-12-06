//
// Created by zsoldos on 11/26/19.
//

///////////////////////// STL C/C++ /////////////////////////
#include <vector>
#include <algorithm>

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

  // Set basic ROOT style
  gStyle->SetOptStat(0);
  gStyle->SetPalette(56);

  // Create ifstream for input file
  string inputName;
  ProcessArgs(&theApp,&inputName);

  // Get ready to read inside inputName
  // each file name
  string line;
  ifstream file(inputName);

  // Create a FileAnalysis for each file
  // to be processed
  vector< TFileAnalysis<TH2D> > FileAnalysis;
  int NbFileAnalysis = 0;

  // Save each Ebins
  vector<double> Ebins;

  while ( getline(file,line) ){
	FileAnalysis.emplace_back(line);
	double E = ExtractEBinFromFilename(line);
	FileAnalysis[NbFileAnalysis].SetEBin(E);

	Ebins.emplace_back(E);

	cout << "ADD " << ExtractFilenameFromPath(line) << endl;
	cout << "-> Corresponding to Ebin: " << ExtractEBinFromFilename(line) << endl;

	NbFileAnalysis++;
  }

  TH2D *hNbPEVSHits;
  TH1D *hNbPE;
  TH2D *hEMatrix;

  sort(Ebins.begin(),Ebins.end());
  vector<double> corEbins = CorrectEbinsRange(Ebins);

  const int nbBinsPE = 201;
  const double minPE = -0.5;
  const double maxPE = 200.5;

  hEMatrix = new TH2D("hEMatrix", "Transition matrix describing the energy response of the detector",
					  nbBinsPE,minPE,maxPE,
					  Ebins.size(),&corEbins[0]);

  vector<double> Eres;

  TCanvas *c1;

  for(auto& file : FileAnalysis){

    cout << "PROCESSING atm " << ExtractFilenameFromPath(file.GetFilename()) << endl;

	hNbPEVSHits = new TH2D(Form("hEbin%.1f", file.GetEBin()),"Nb PE Collected VS Nb PMTs Hits",
						   nbBinsPE,minPE,maxPE,
						   1001,-0.5,1000.5);

	file.SetHist(hNbPEVSHits);

	file.DoAnalysis(CollectPEAndHits);

	///////////////////////////////
	// Recover hnPE projection   //
	///////////////////////////////

	hNbPE = (TH1D*)file.GetHist()->ProjectionX()->Clone();

	///////////////////////////////
	// FILL transition matrix    //
	///////////////////////////////

	for(int iBin=1; iBin<nbBinsPE; iBin++){

	  double BinPE = hNbPE->GetBinCenter(iBin);
	  double BinPEContent = hNbPE->GetBinContent(iBin);
	  double BinE = file.GetEBin();

	  hEMatrix->Fill(BinPE, BinE, BinPEContent);

	}

	///////////////////////////////
	// FIT and recover res. info //
	///////////////////////////////

	hNbPE->Fit("gaus","0LQSEM+");
	TF1 *fFit = hNbPE->GetFunction("gaus");
	double chi2 = fFit->GetChisquare();
	double c = fFit->GetParameter(0);
	double mean = fFit->GetParameter(1);
	double meanErr = fFit->GetParError(1);
	double sigma = fFit->GetParameter(2);
	double sigmaErr = fFit->GetParError(2);

//	c1 = new TCanvas(Form("cEbin%.1fMeV2D", file.GetEBin()),Form("cEbin%.1fMeV2D", file.GetEBin()),
//					 800,600);
//	c1->SetGrid();
//	file.GetHist()->Draw("COLZ");
//
//	c1 = new TCanvas(Form("cEbin%.1fMeV1D", file.GetEBin()),Form("cEbin%.1fMeV1D", file.GetEBin()),
//					 800,600);
//	c1->SetGrid();
//	hNbPE->Draw("COLZ");

  }


  c1 = new TCanvas("cEMatrix", "cEMatrix", 800,600);
  c1->SetGrid();
  hEMatrix->Draw("COLZ");

  /////////////////////////
  // ...

  cout << endl;
  cout << "Hit Ctrl+C to exit" << endl;
  theApp.Run(kTRUE);

  return 0;
}
