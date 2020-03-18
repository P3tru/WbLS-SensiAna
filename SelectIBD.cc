//
// Created by zsoldos on 12/5/19.
//

///////////////////////// STL C/C++ /////////////////////////
#include <string>
#include <vector>
#include <csignal>

/////////////////////////   ROOT   //////////////////////////
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TFile.h>
#include <TF1.h>

/////////////////////////   USER   //////////////////////////
#include "utils.hh"
#include "SelectIBD.hh"

#include "Analyzer.hh"
#include "HitClass.hh"
#include "MCFunctions.hh"
#include "HitFunctions.hh"
#include "FitPosTime.hh"
#include "CalibFunctions.hh"

#include "ProgressBar.hpp"

using namespace std;

#define DIAMETER = 10857
#define HEIGHT = 10857
#define BUFFER = 500
#define SQRT2 = 1.41421
#define PRETRIG = (DIAMETER+BUFFER)*SQRT2/C

TH1D *GetHPDF(string basic_string);
int main(int argc, char *argv[]) {

  // Get Signal if user wants to interrupt loop
  EoF=0;
  signal(SIGINT,Interrupt);


  // Create TApp
  TApplication theApp("App", &argc, argv);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####             SET BASIC ROOT STYLE                  #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  SetBasicStyle();


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####            PARSE AND SET INPUT PARAMS             #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  // Input parameters
  // There's default value chosen by ternaries below
  string inputName;
  inputName.clear();

  auto User_isBatch = false;
  auto User_nEvts = -1;
  auto User_iEvt = -1;

  // Select time cut for computing residuals
  auto User_PromptCut = -1;

  // Select Prompt event in E range
  auto User_MinEPrompt = -1.; auto User_MaxEPrompt = -1.;

  // Select Delayed evts for each capture?
  // n-H 2.22 MeV
  auto User_MinEDelayed_nH = -1.; auto User_MaxEDelayed_nH = -1.;
  // n-C 4.95 MeV
  auto User_MinEDelayed_nC = -1.; auto User_MaxEDelayed_nC = -1.;
  // n-O 0.87 MeV
  auto User_MinEDelayed_nO = -1.; auto User_MaxEDelayed_nO = -1.;

  // Select DR cutoff
  auto User_DR = -1.;

  // Select DT cutoff
  auto User_DT = -1.;

  string User_Calib;
  User_Calib.clear();
  string User_fPDF;
  User_fPDF.clear();

  ProcessArgs(&theApp, &inputName,
			  &User_PromptCut,
			  &User_nEvts, &User_iEvt,
			  &User_MinEPrompt, &User_MaxEPrompt,
			  &User_MinEDelayed_nH, &User_MaxEDelayed_nH,
			  &User_MinEDelayed_nC, &User_MaxEDelayed_nC,
			  &User_MinEDelayed_nO, &User_MaxEDelayed_nO,
			  &User_DR,
			  &User_DT,
			  &User_Calib,
			  &User_isBatch,
			  &User_fPDF);

  const bool isBatch = User_isBatch;

  // Add a prompt window, because there's no separate trigger for 2 evts in rat-pac.
  // The IBD generator will not create 2 separate EV object.
  // Therefore, add a cut to physically separate prompt+delay
  const double PromptWindow = 200.; // ns

  const double MinEPrompt = SetDefValue(User_MinEPrompt, 0.);
  const double MaxEPrompt = SetDefValue(User_MaxEPrompt, 1000.);

  // TODO : Allow user to set true origin vector from input
  const TVector3 TrueOrigin(0.,0.,0.);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE HISTOGRAMS                  #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  TH1D *hPosGuess[3];
  for(int iPos = 0; iPos<3; iPos++) {

	hPosGuess[iPos] = new TH1D(Form("hPosGuess%d", iPos),
							   Form("Vtx Reconstruction Axis %d", iPos),
							   21, -1000.5, 1000.5);
  }

  TH1D *hTGuess;
  hTGuess = new TH1D("hTGuess", "T Reconstruction",
					 201, -100.5, 100.5);

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE OUTPUTS                     #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *foutput = new TFile("dummy.root", "RECREATE");


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                LOAD  PDF                          #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  TH1D *hPDF = GetHPDF(User_fPDF.c_str(), "hTResiduals");
  if(hPDF){
    cout << "PDF loaded" << endl;
  }

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                SET CALIBRATION                    #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  MCCalib Cal(User_Calib);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE ANALYZER                    #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *FileAnalyzer = new Analyzer(inputName.c_str());
  unsigned long int nEvts = User_nEvts > 0 ? User_nEvts : FileAnalyzer->GetNEvts();
  unsigned int iEvt = SetDefValue(User_iEvt, 0);
  nEvts = User_iEvt > 0 ? User_iEvt + nEvts : nEvts;
  cout << "nEvts: " << nEvts << endl;
  cout << "iEvt: " << iEvt << endl;

  unsigned int nEvtToProcess = nEvts-iEvt;

  ProgressBar progressBar(nEvtToProcess, 70);

  for(iEvt; iEvt<nEvts; iEvt++){

	if(EoF) break;

	// record the tick
	++progressBar;

	// Recover Hit vector for 1 evt
	// vector<Hit> vHit;
	vector<Hit> vHit = GetHitCollection(FileAnalyzer, iEvt);

	// Split vector hit into prompt and delay
	vector<Hit> vHitDelayed = SplitHitCollection(&vHit, 200);

	// Get E evts
	double NPE, NHits, E;
	if(!User_Calib.empty()){

	  GetNPEAndNHitsFromHits(vHit, &NPE, &NHits);
	  E = ComputeECalib(Cal, NPE, NHits);

	}

	// Fit position time both events
	if(vHit.size()>0){

	  if(hPDF){
		FitPosTime FPT_Prompt(TrueOrigin, 0, vHit, hPDF);

		for(int iPos = 0; iPos<3; iPos++) {

		  hPosGuess[iPos]->Fill(FPT_Prompt.GetPos()(iPos));

		}

		hTGuess->Fill(FPT_Prompt.GetT());

	  }

	}

	// display the bar
	progressBar.display();

  } // END FOR iEVT

  cout << " DONE" << endl;
  EoF=1;


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                      DRAWING                      #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *c1 = new TCanvas("cVTX", "cVTX", 800,600);
  c1->SetGrid();
  for(auto h:hPosGuess){
	h->Draw("SAME PLC PMC");
  }

  c1 = new TCanvas("cTime", "cTime", 800, 600);
  c1->SetGrid();
  hTGuess->Draw();

  foutput->cd(); c1->Write();
  foutput->Close();


  /////////////////////////
  // ...

  if(!isBatch&&gROOT->GetListOfCanvases()->GetEntries()>0) {

	cout << endl;
	cout << "Hit Ctrl+C to exit" << endl;
	theApp.Run(kTRUE);

  }

  return 0;

}
