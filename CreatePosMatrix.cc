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
#include "CreatePosMatrix.hh"
#include "AnalysisDefinitions.hh"
#include "utils.hh"

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

  // Save each PosBins
  vector<double> Posbins;

  while ( getline(file,line) ){

	if (line.compare(0,1,"#") == 0){

	  continue;

	} else {

	  FileAnalysis.emplace_back(line);
	  double Pos = ExtractPosBinFromFilename(line);
	  FileAnalysis[NbFileAnalysis].SetID(Pos);

	  Posbins.emplace_back(Pos);

	  cout << "ADD " << ExtractFilenameFromPath(line) << endl;
	  cout << "-> Corresponding to Position: " << Pos << "mm" << endl;

	  NbFileAnalysis++;

	}

  } // END while

  TH2D *hNbPEVSHits;
  TH1D *hNbPE;
  TH2D *hNbPEVSPos;

  sort(Posbins.begin(),Posbins.end());
  vector<double> corBins = CorrectBinRangeArray(Posbins);

  const int nbBinsPE = (User_nPEBins>-1) ? User_nPEBins : 201;
  const double minPE = (User_minPE>-1) ? User_minPE : -0.5;
  const double maxPE = (User_maxPE>-1) ? User_maxPE : 200.5;

  const int nbBinsPMT = (User_nPMTBins>-1) ? User_nPMTBins : 201;
  const double minPMT = (User_minPMT>-1) ? User_minPMT : -0.5;
  const double maxPMT = (User_maxPMT>-1) ? User_maxPMT : 200.5;


  hNbPEVSPos = new TH2D("hNbPEVSPos", "PE collection VS position of the vtx",
						nbBinsPE,minPE,maxPE,
						corBins.size()-1,&corBins[0]);

  auto fOutputName = (!User_fOutput.empty()) ? User_fOutput : "output.root";
  auto *fOutput = new TFile(fOutputName.c_str(),"RECREATE");

  TCanvas *c1;

  for(auto& itFile : FileAnalysis){

	cout << "PROCESSING atm " << ExtractFilenameFromPath(itFile.GetFilename()) << endl;

	hNbPEVSHits = new TH2D(Form("hPosBin%.0fmm", itFile.GetID()),"Nb PE Collected VS Nb PMTs Hits",
						   nbBinsPE,minPE,maxPE,
						   nbBinsPMT,minPMT,maxPMT);

	itFile.SetHist(hNbPEVSHits);

	itFile.DoAnalysis(CollectPEAndHits);

	///////////////////////////////
	// Recover hnPE projection   //
	///////////////////////////////

	hNbPE = (TH1D*)itFile.GetHist()->ProjectionX()->Clone();

	///////////////////////////////
	// Save all plots in ROOT    //
	///////////////////////////////

	fOutput->cd();
	itFile.GetHist()->Write();
	hNbPE->Write();

	///////////////////////////////
	// FILL transition matrix    //
	///////////////////////////////

	for(int iBin=1; iBin<nbBinsPE; iBin++){

	  double BinPE = hNbPE->GetBinCenter(iBin);
	  double BinPEContent = hNbPE->GetBinContent(iBin);
	  double BinPos = itFile.GetID();

	  hNbPEVSPos->Fill(BinPE, BinPos, BinPEContent);

	}

  }

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                      DRAWING                      #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  c1 = new TCanvas("c1","c1",800,600);
  hNbPEVSPos->Draw("COLZ");

  /////////////////////////
  // ...

  hNbPEVSPos->Write();
  fOutput->Close();


  cout << endl;
  cout << "Hit Ctrl+C to exit" << endl;
  theApp.Run(kTRUE);

  return 0;

}