//
// Created by zsoldos on 12/5/19.
//

///////////////////////// STL C/C++ /////////////////////////
#include <string>
#include <algorithm>

/////////////////////////   ROOT   //////////////////////////
#include <TApplication.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TF1.h>

/////////////////////////   USER   //////////////////////////
#include "TFileAnalysis.hh"
#include "VtxRecon.hh"
#include "utils.hh"

using namespace std;

#define DIAMETER = 10857
#define HEIGHT = 10857
#define BUFFER = 500
#define SQRT2 = 1.41421
#define PRETRIG = (DIAMETER+BUFFER)*SQRT2/C

int main(int argc, char *argv[]) {

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

  int User_PromptCut=-1;
  int User_nEvts=-1;

  int User_nTResidBins=-1;
  double User_minTResid=-1; double User_maxTResid=-1;

  string User_fPDF;

  bool User_batch = false;

  bool User_IsDelay = false;

  // READ input parameters
  ProcessArgs(&theApp,&inputName,
			  &User_PromptCut,
			  &User_nEvts,
			  &User_batch,
			  &User_IsDelay,
			  &User_nTResidBins, &User_minTResid, &User_maxTResid,
			  &User_fPDF);

  const int PromptCut = SetDefValue(User_PromptCut, 15);

  const int nTResidBins = SetDefValue(User_nTResidBins, 101);
  const double minTResid = SetDefValue(User_minTResid, -0.5); // ns
  const double maxTResid = SetDefValue(User_maxTResid, 100.5); // ns

  // Useful to separate events if there is not time associated.
  const int PromptWindow = 256; // 64*4ns samples (typical PMT signal)

  const bool isBatch = User_batch;

  // TODO : Allow user to set true origin vector from input
  const TVector3 TrueOrigin(0.,0.,0.);

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE HISTOGRAMS                  #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *hTResiduals = new TH1D("hTResiduals", "T Residuals",
							   nTResidBins, minTResid, maxTResid);

  SetBasicTH1Style(hTResiduals, kBlue-4);
  hTResiduals->Sumw2();


  auto *hTResidualsDelayed = new TH1D("hTResidualsDelayed", "T Residuals",
							   nTResidBins, minTResid, maxTResid);

  SetBasicTH1Style(hTResidualsDelayed, kBlue-4);
  hTResidualsDelayed->Sumw2();

  auto *hDT = new TH1D("hDT", "#Delta T between prompt and delayed",
					   100, 0., 1.e6);
  SetBasicTH1Style(hDT, kBlue-4);
  hDT->Sumw2();

  const double minDRRes = -0.5; // mm
  const double maxDRRes = 10000.5; // mm
  const int NbBinsDRRes = 101; // dm resolution
  auto *hDR = new TH1D("hDR", "#Delta R between prompt and delayed",
					   NbBinsDRRes,minDRRes,maxDRRes);


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE HIT COLLECTION              #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  vector<HitCollection> Evts;
  vector<HitCollection> EvtsDelayed;


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                   CREATE PDFs                     #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *fMC = new TFile(inputName.c_str());

  if(fMC->IsOpen()){

	////////////////////////////////
	// IF file is open do stuff ////
	////////////////////////////////

	auto *run = new RAT::DS::Run();

	auto *runTree = (TTree *) fMC->Get("runT");
	runTree->SetBranchAddress("run", &run);
	runTree->GetEntry(0);

	auto *tree = (TTree *) fMC->Get("T");

	auto *rds = new RAT::DS::Root();
	tree->SetBranchAddress("ds", &rds);

	auto nEvents = (User_nEvts > 0 && User_nEvts < tree->GetEntries()) ? User_nEvts : tree->GetEntries();
	cout << "#Evts: " << nEvents << endl;
	ProgressBar progressBar(nEvents, 70);

	cout << "CREATE PDFs" << endl;

	for (int iEvt = 0; iEvt < nEvents; iEvt++) {
	  tree->GetEntry(iEvt);
	  // record the tick
	  ++progressBar;

	  // Access RAT MC info and the summary
	  // Summary useful to get nCer photons, nScint photons, etc...
	  auto *mc = rds->GetMC();

	  auto NbPMTsHits = mc->GetMCPMTCount();

	  vector<Hit> vHit;
	  vector<Hit> vHitDelayed;

	  for (int iPMT = 0; iPMT < NbPMTsHits; iPMT++) {

		auto mcPMT = mc->GetMCPMT(iPMT);
		auto PMTID = mcPMT->GetID();
		auto PMTPos = run->GetPMTInfo()->GetPosition(PMTID);

		auto NbPhotonCounts = mcPMT->GetMCPhotonCount();

		for (int iP = 0; iP < NbPhotonCounts; iP++) {

		  auto mcPhoton = mcPMT->GetMCPhoton(iP);
		  // Get Q and T
		  auto Q = mcPhoton->GetCharge();
		  auto T = mcPhoton->GetHitTime();

		  // ########## //
		  // FILL EVENT //
		  // ########## //

		  Hit hit(PMTPos, Q, T);

		  if(T<PromptWindow){

			vHit.emplace_back(hit);

			hTResiduals->Fill(hit.CalculateTResid(TrueOrigin),
							  1/(double)(NbPhotonCounts*NbPMTsHits));

		  } else {

			vHitDelayed.emplace_back(hit);

			// Fill the delayed histogram later,
			// so we can apply an offset on the hits
			// to make sense of the residuals

		  }

		}

	  } // END FOR iPMT

	  // Recover time of first hit for each collection,
	  // defining for now the trigger time
	  double TPrompt = -1;
	  double TDelayed = -1;

	  // ################ //
	  // APPLY PROMPT CUT //
	  // ################ //

	  ProcessVHit(vHit,Evts,Form("Evt%d",iEvt),
				  PromptCut, 0,true,&TPrompt);


	  // ########################################################## //
	  // Set Trigger time for delayed event                         //
	  // Take the first hit and define as new reference - pretrig ? //
	  // Take pretrig ~ 20ns for 1kT for now                        //
	  // ########################################################## //

	  ProcessVHit(vHitDelayed,EvtsDelayed,Form("DELAYED_Evt%d",iEvt),
				  PromptCut, 32, false, &TDelayed);

	  if(vHitDelayed.size() > 0) {

	    auto iE = EvtsDelayed.size()-1;

	    for(auto itHit : EvtsDelayed[iE].GetHCol()){
		  hTResidualsDelayed->Fill(itHit.CalculateTResid(TrueOrigin),
								   1/(double)(EvtsDelayed[iE].GetHCol().size()));
		}

	  }

	  hDT->Fill(TDelayed-TPrompt);

	  // display the bar
	  progressBar.display();

	} // END FOR iEvt

  } else {

	cout << "File: " << ExtractFilenameFromPath(inputName)
		 << " not open somehow!?" << endl;

  } // END IF FILE MC

  cout << endl;
  fMC->Close();

  hTResiduals->Scale(1/(double)(Evts.size()));
  hTResidualsDelayed->Scale(1/(double)(EvtsDelayed.size()));

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                        FIT                        #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  bool DoYouWannaFit = true; // pelligro ma donna

  const double minPosRes = -1000.5; // mm
  const double maxPosRes = 1000.5; // mm
  const int NbBinsPosRes = 21; // dm resolution

  TH1D *hVtxResPrompt[3];
  TH1D *hVtxResDelayed[3];

  for(int iAxis=0; iAxis<3; iAxis++){
	hVtxResPrompt[iAxis] = new TH1D(Form("hVtxResPrompt%d",iAxis),
									"Vtx reconstruction",
									NbBinsPosRes, minPosRes, maxPosRes);
	hVtxResDelayed[iAxis] = new TH1D(Form("hVtxResDelayed%d",iAxis),
									"Vtx reconstruction",
									NbBinsPosRes, minPosRes, maxPosRes);
  }

  if(DoYouWannaFit){

    vector<TVector3> POS_PROMPT;

	FitTResid(hTResiduals, Evts, hVtxResPrompt, &POS_PROMPT);

	vector<TVector3> POS_DELAYED;

    FitTResid(hTResidualsDelayed, EvtsDelayed, hVtxResDelayed, &POS_DELAYED);

    CreateDRHist(hDR, POS_PROMPT, POS_DELAYED);

  }

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                      DRAWING                      #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //


  auto *c1 = new TCanvas("cTResidsPDF","cTResidsPDF", 800,600);
  c1->SetGrid();
  hTResiduals->Draw();

  c1 = new TCanvas("cTResidsPDF_DELAYED","cTResidsPDF_DELAYED", 800,600);
  c1->SetGrid();
  hTResidualsDelayed->Draw();

  c1 = new TCanvas("cDT","cDT", 800,600);
  c1->SetGrid();
  hDT->Draw();
  TF1 *fExpo;
  fExpo = new TF1("fExpo",fitExpo,0.,1.e6,2);
  fExpo->SetParName(0,"C0");fExpo->SetParameter(0,hDT->GetMaximum());
  fExpo->SetParName(1,"tau");fExpo->SetParameter(1,200);
  hDT->Fit("fExpo", "RLE");

  c1 = new TCanvas("cDR","cDR", 800,600);
  c1->SetGrid();
  hDR->Draw();

  c1 = new TCanvas("cVTX","cVTX", 800,600);
  c1->SetGrid();
  TF1 *fGaus;
  fGaus = new TF1("fGaus",fitGaus,minPosRes,maxPosRes,3);
  fGaus->SetParName(0,"C0");
  fGaus->SetParName(1,"mu");
  fGaus->SetParName(2,"sigma");

  for(int iAxis=0; iAxis<3; iAxis++){
    if(hVtxResPrompt[iAxis]){
	  SetBasicTH1Style(hVtxResPrompt[iAxis], kBlue);
	  hVtxResPrompt[iAxis]->Draw("SAME PLC PMC");
	  FitPos(hVtxResPrompt[iAxis], fGaus);
	  TFitResultPtr r = hVtxResPrompt[iAxis]->Fit("fGaus", "RLQE");
	  if(r == 0){
		SetFitStyle(hVtxResPrompt[iAxis],kBlue-4);
		cout << iAxis << endl;
		PrintFitResult(fGaus);
	  }
    }
  }

  c1 = new TCanvas("cVTXDelayed","cVTXDelayed", 800,600);
  c1->SetGrid();

  for(int iAxis=0; iAxis<3; iAxis++){
	if(hVtxResDelayed[iAxis]){
	  SetBasicTH1Style(hVtxResDelayed[iAxis], kBlue);
	  hVtxResDelayed[iAxis]->Draw("SAME PLC PMC");
	  FitPos(hVtxResDelayed[iAxis], fGaus);
	  TFitResultPtr r = hVtxResDelayed[iAxis]->Fit("fGaus", "RLQE");
	  if(r == 0){
		SetFitStyle(hVtxResDelayed[iAxis],kBlue-4);
		cout << iAxis << endl;
		PrintFitResult(fGaus);
	  }
	}
  }

  /////////////////////////
  // ...

  if(!isBatch&&gROOT->GetListOfCanvases()->GetEntries()>0) {

	cout << endl;
	cout << "Hit Ctrl+C to exit" << endl;
	theApp.Run(kTRUE);

  }

  return 0;

}