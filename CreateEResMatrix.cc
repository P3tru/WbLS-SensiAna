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
#include "TFileAnalysis.hh"
#include "AnalysisDefinitions.hh"
#include "utils.hh"

#include "CreateEResMatrix.hh"

using namespace std;

int main(int argc, char *argv[]) {

  // Create TApp
  TApplication theApp("App", &argc, argv);

  // Set basic ROOT style
  gStyle->SetOptStat(0);
  gStyle->SetPalette(56);

  // Input parameters
  // There's default value chosen by ternaries below
  string inputName;
  int User_nPEBins=-1;
  double User_minPE=-1; double User_maxPE=-1;
  int User_nPMTBins=-1;
  double User_minPMT=-1; double User_maxPMT=-1;
  string User_fOutput;

  // READ input parameters
  ProcessArgs(&theApp,&inputName,
			  &User_nPEBins, &User_minPE, &User_maxPE,
			  &User_nPMTBins, &User_minPMT, &User_maxPMT,
			  &User_fOutput);

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

    if (line.compare(0,1,"#") == 0){

	  continue;

    } else {

	  FileAnalysis.emplace_back(line);
	  double E = ExtractEBinFromFilename(line);
	  FileAnalysis[NbFileAnalysis].SetEBin(E);

	  Ebins.emplace_back(E);

	  cout << "ADD " << ExtractFilenameFromPath(line) << endl;
	  cout << "-> Corresponding to Ebin: " << ExtractEBinFromFilename(line) << endl;

	  NbFileAnalysis++;

    }

  } // END while

  TH2D *hNbPEVSHits;
  TH1D *hNbPE;
  TH2D *hEMatrix;

  sort(Ebins.begin(),Ebins.end());
  vector<double> corEbins = CorrectBinRangeArray(Ebins);

  const int nbBinsPE = (User_nPEBins>-1) ? User_nPEBins : 201;
  const double minPE = (User_minPE>-1) ? User_minPE : -0.5;
  const double maxPE = (User_maxPE>-1) ? User_maxPE : 200.5;

  const int nbBinsPMT = (User_nPMTBins>-1) ? User_nPMTBins : 201;
  const double minPMT = (User_minPMT>-1) ? User_minPMT : -0.5;
  const double maxPMT = (User_maxPMT>-1) ? User_maxPMT : 200.5;


  hEMatrix = new TH2D("hEMatrix", "Transition matrix describing the energy response of the detector",
					  Ebins.size(),&corEbins[0],
					  nbBinsPE,minPE,maxPE);

  TCanvas *c1;

  auto fOutputName = (!User_fOutput.empty()) ? User_fOutput : "output.root";
  auto *fOutput = new TFile(fOutputName.c_str(),"RECREATE");

  for(auto& file : FileAnalysis){

    cout << "PROCESSING atm " << ExtractFilenameFromPath(file.GetFilename()) << endl;

	hNbPEVSHits = new TH2D(Form("hEbin%.1f", file.GetEBin()),"Nb PE Collected VS Nb PMTs Hits",
						   nbBinsPE,minPE,maxPE,
						   nbBinsPMT,minPMT,maxPMT);

	file.SetHist(hNbPEVSHits);

	file.DoAnalysis(CollectPEAndHits);

	///////////////////////////////
	// Recover hnPE projection   //
	///////////////////////////////

	hNbPE = (TH1D*)file.GetHist()->ProjectionX()->Clone();

	///////////////////////////////
	// Save all plots in ROOT    //
	///////////////////////////////

	fOutput->cd();
	hNbPEVSHits->Write();
	hNbPE->Write();

	///////////////////////////////
	// FILL transition matrix    //
	///////////////////////////////

	for(int iBin=1; iBin<nbBinsPE; iBin++){

	  double BinPE = hNbPE->GetBinCenter(iBin);
	  double BinPEContent = hNbPE->GetBinContent(iBin);
	  double BinE = file.GetEBin();

	  hEMatrix->Fill(BinE, BinPE, BinPEContent);

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

  hEMatrix->Write();

  c1 = new TCanvas("cEMatrix", "cEMatrix", 800,600);
  c1->SetGrid();
  hEMatrix->Draw("COLZ");

  /////////////////////////
  // ...

  fOutput->Close();

  cout << endl;
  cout << "Hit Ctrl+C to exit" << endl;
  theApp.Run(kTRUE);

  return 0;
}
