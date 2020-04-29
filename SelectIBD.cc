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
#include "EVFunctions.hh"
#include "FitPosTimeSimulatedAnnealing.hh"
#include "HitClass.hh"
#include "HitFunctions.hh"
#include "CalibFunctions.hh"

#include "ProgressBar.hpp"

using namespace std;

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

  // Add DAQ window if RAT hasn't created 2evts
  auto User_DAQWindow = -1;

  auto User_xx = -1.;
  auto User_yy = -1.;
  auto User_zz = -1.;

  auto User_tt = -1.;

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

  auto User_E = -1.;

  auto User_ManCalib = -1.;
  string User_Calib;
  User_Calib.clear();

  string User_fPDF;
  User_fPDF.clear();

  ProcessArgs(&theApp, &inputName,
			  &User_PromptCut,
			  &User_DAQWindow,
			  &User_nEvts, &User_iEvt,
			  &User_MinEPrompt, &User_MaxEPrompt,
			  &User_MinEDelayed_nH, &User_MaxEDelayed_nH,
			  &User_MinEDelayed_nC, &User_MaxEDelayed_nC,
			  &User_MinEDelayed_nO, &User_MaxEDelayed_nO,
			  &User_xx, &User_yy, &User_zz,
			  &User_tt,
			  &User_DR,
			  &User_DT,
			  &User_E,
			  &User_ManCalib,
			  &User_Calib,
			  &User_isBatch,
			  &User_fPDF);

  const bool isBatch = User_isBatch;

  // Add a prompt window, because there's no separate trigger for 2 evts in rat-pac.
  // The IBD generator will not create 2 separate EV object.
  // Therefore, add a cut to physically separate prompt+delay
  const double DAQWindow = SetDefValue(User_DAQWindow, 256); // ns

  const double MinEPrompt = SetDefValue(User_MinEPrompt, 1.); // MeV
  const double MaxEPrompt = SetDefValue(User_MaxEPrompt, 10.); // MeV

  const double MinEDelayed_nH = SetDefValue(User_MinEDelayed_nH, 1.); // MeV
  const double MaxEDelayed_nH = SetDefValue(User_MaxEDelayed_nH, 10.); // MeV

  const double MinEDelayed_nC = SetDefValue(User_MinEDelayed_nC, 1.); // MeV
  const double MaxEDelayed_nC = SetDefValue(User_MaxEDelayed_nC, 10.); // MeV

  const double MinEDelayed_nO = SetDefValue(User_MinEDelayed_nO, 1.); // MeV
  const double MaxEDelayed_nO = SetDefValue(User_MaxEDelayed_nO, 10.); // MeV

  // TODO : Allow user to set true origin vector from input
  const double xx = SetDefValue(User_xx, 0.); // mm
  const double yy = SetDefValue(User_xx, 0.); // mm
  const double zz = SetDefValue(User_xx, 0.); // mm
  const TVector3 TrueOrigin(xx,yy,zz);

  const double TrueTime = SetDefValue(User_tt, 0.); //ns

  const double DR = SetDefValue(User_DR, 1000.); // mm
  const double DT = SetDefValue(User_DT, 1.e6); // ns

  const double TrueE = SetDefValue(User_E, 0.);

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                SET CALIBRATION                    #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  MCCalib Cal;
  const double ManCalib = SetDefValue(User_ManCalib, 0.);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                LOAD  PDF                          #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  TH1D *hPDF = GetHPDF(User_fPDF.c_str(), "hTResiduals");
  if(hPDF){
	cout << "PDF loaded" << endl;
  }


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


  TH1D *hEPrompt = new TH1D("hEPrompt", "E_{Rec} Prompt",
							101, -0.05, 10.05);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE OUTPUTS                     #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *foutput = new TFile(Form("%s_RECON.root",inputName.c_str()),
							"RECREATE");


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
	vector<Hit> vHit = GetEVHitCollection(FileAnalyzer, iEvt, 0);

	// Split vector hit into prompt and delay
	vector<Hit> vHitDelayed = SplitVHits(&vHit, 200);

	// Get E evts

	double NPE=-1; double NHits=-1; double E=-1;

	if(ManCalib>0){

	  GetNPEAndNHitsFromHits(vHit, &NPE, &NHits);
	  E = NHits / ManCalib;

	} else {

	  if (!User_Calib.empty()) {

		GetNPEAndNHitsFromHits(vHit, &NPE, &NHits);
		E = ComputeECalib(Cal, NPE, NHits);

	  }

	}

	if(E>0){

	  if(E>MinEPrompt && E<MaxEPrompt){

		hEPrompt->Fill(E-TrueE);

	  }

	}


	// Fit position time both events
	if(vHit.size()>0){

	  if(hPDF){
		SortVHits(&vHit);
		auto *FT = new FitTheia(hPDF, TrueOrigin, TrueTime, vHit);

		const size_t dim = 4;
		RAT::SimulatedAnnealing<4> anneal(FT);
		vector<double> point(dim), seed(dim);

		// Regular simplex in 4D
		// https://en.wikipedia.org/wiki/5-cell

		seed[0] = 1; seed[1] = 1; seed[2] = 1; seed[3] = -1/SQRT5;
		ScaleSimplexCoord(&seed);
		anneal.SetSimplexPoint(0,seed);
		seed[0] = 1; seed[1] = -1; seed[2] = -1; seed[3] = -1/SQRT5;
		ScaleSimplexCoord(&seed);
		anneal.SetSimplexPoint(1,seed);
		seed[0] = -1; seed[1] = 1; seed[2] = -1; seed[3] = -1/SQRT5;
		ScaleSimplexCoord(&seed);
		anneal.SetSimplexPoint(2,seed);
		seed[0] = -1; seed[1] = -1; seed[2] = 1; seed[3] = -1/SQRT5;
		ScaleSimplexCoord(&seed);
		anneal.SetSimplexPoint(3,seed);
		seed[0] = 0; seed[1] = 0; seed[2] = 0; seed[3] = SQRT5 - 1/SQRT5;
		ScaleSimplexCoord(&seed);
		anneal.SetSimplexPoint(4,seed);

		anneal.Anneal(10,150,50,4.0); // Minimize

		anneal.GetBestPoint(point);

		for(int iPos = 0; iPos<3; iPos++) {

		  hPosGuess[iPos]->Fill(point[iPos]);

		}

		hTGuess->Fill(point[3]);

		delete FT;

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
  c1->BuildLegend();

  c1 = new TCanvas("cTime", "cTime", 800, 600);
  c1->SetGrid();
  hTGuess->Draw();

  c1 = new TCanvas("cE", "cE", 800, 600);
  c1->SetGrid();
  hEPrompt->Draw();

  foutput->cd();
  for(auto h:hPosGuess){
	h->Write();
  }
  hTGuess->Write();
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
