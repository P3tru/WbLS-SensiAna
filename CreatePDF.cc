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
#include "CreatePDF.hh"

#include <Analyzer.hh>
#include <HitClass.hh>
#include <EVFunctions.hh>
#include <HitFunctions.hh>

#include "ProgressBar.hpp"

using namespace std;

#define DIAMETER 10857
#define HEIGHT 10857
#define BUFFER 500
#define SQRT2 1.41421
#define PRETRIG (DIAMETER+BUFFER)*SQRT2/C

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

  auto User_isBatch = false;
  auto User_nEvts = -1;
  auto User_iEvt = -1;

  int User_nTResidBins_Prompt=-1;
  double User_minTResid_Prompt=-1; double User_maxTResid_Prompt=-1;

  int User_nTResidBins_Delay=-1;
  double User_minTResid_Delay=-1; double User_maxTResid_Delay=-1;

  // Select time cut for computing residuals
  auto User_PromptCut = -1;

  double User_xx = -1;
  double User_yy = -1;
  double User_zz = -1;

  ProcessArgs(&theApp, &inputName,
			  &User_PromptCut,
			  &User_nEvts, &User_iEvt,
			  &User_nTResidBins_Prompt, &User_minTResid_Prompt, &User_maxTResid_Prompt,
			  &User_nTResidBins_Delay, &User_minTResid_Delay, &User_maxTResid_Delay,
			  &User_xx, &User_yy, &User_zz,
			  &User_isBatch);

  const bool isBatch = User_isBatch;

  // Add a prompt window, because there's no separate trigger for 2 evts in rat-pac.
  // The IBD generator will not create 2 separate EV object.
  // Therefore, add a cut to physically separate prompt+delay
  const double PromptWindow = 200.; // ns

  const int PromptCut = SetDefValue(User_PromptCut, 15);

  const int nTResidBins_Prompt = SetDefValue(User_nTResidBins_Prompt, 201);
  const double minTResid_Prompt = SetDefValue(User_minTResid_Prompt, -100.5); // ns
  const double maxTResid_Prompt = SetDefValue(User_maxTResid_Prompt, 100.5); // ns

  const int nTResidBins_Delay = SetDefValue(User_nTResidBins_Delay, 201);
  const double minTResid_Delay = SetDefValue(User_minTResid_Delay, -100.5); // ns
  const double maxTResid_Delay = SetDefValue(User_maxTResid_Delay, 100.5); // ns

  // TODO : Allow user to set true origin vector from input
  const double xx = SetDefValue(User_xx, 0.); // mm
  const double yy = SetDefValue(User_xx, 0.); // mm
  const double zz = SetDefValue(User_xx, 0.); // mm
  const TVector3 TrueOrigin(xx,yy,zz);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE HISTOGRAMS                  #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *hTResiduals = new TH1D("hTResiduals", "T Residuals",
							   nTResidBins_Prompt, minTResid_Prompt, maxTResid_Prompt);

  SetBasicTH1Style(hTResiduals, kBlue-4);
  hTResiduals->Sumw2();


  auto *hTResidualsDelayed = new TH1D("hTResidualsDelayed", "T Residuals",
									  nTResidBins_Delay, minTResid_Delay, maxTResid_Delay);

  SetBasicTH1Style(hTResidualsDelayed, kBlue-4);
  hTResidualsDelayed->Sumw2();

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE OUTPUTS                     #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *foutput = new TFile(Form("%s_PDF.root",inputName.c_str()),
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
	vector<Hit> vHit = GetEVHitCollection(FileAnalyzer, iEvt);

	// Split vector hit into prompt and delay
	vector<Hit> vHitDelayed = SplitVHits(&vHit, PromptWindow);

	SortVHits(&vHit);

	// Hit hCut(TVector3(0, 0, 0),
	// 		 0,
	// 		 PromptCut + vHit.begin()->GetT());
	// RemoveHitsAfterCut(vHit, hCut);


	FillResiduals(hTResiduals, ResetTVHits(vHit),
				  TrueOrigin,
				  224.9,
				  0,
				  false);

	if(vHitDelayed.size()>0){
	  FillResiduals(hTResidualsDelayed, ResetTVHits(vHitDelayed),
					TrueOrigin,
					224.9,
					0,
					false);
	}

	// display the bar
	progressBar.display();

  } // END FOR iEVT

  cout << " DONE" << endl;
  EoF=1;


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                      DRAWING                      #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  if(!isBatch){
	auto *c1 = new TCanvas("cTResidPrompt", "cTResidPrompt", 800,600);
	c1->SetGrid();
	hTResiduals->Scale(1/(double)(nEvtToProcess));
	hTResiduals->Draw();

	if(hTResidualsDelayed->GetEntries()>0){
	  c1 = new TCanvas("cTResidDelay", "cTResidDelay", 800,600);
	  c1->SetGrid();
	  hTResidualsDelayed->Scale(1/(double)(nEvtToProcess));
	  hTResidualsDelayed->Draw();
	}

  }

  foutput->cd();
  hTResiduals->Write();
  if(hTResidualsDelayed->GetEntries()>0) {
	hTResidualsDelayed->Write();
  }
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
