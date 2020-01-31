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

  int User_nDBins=-1;
  double User_minD=-1; double User_maxD=-1;

  string User_fPDF;

  bool User_batch = false;

  // READ input parameters
  ProcessArgs(&theApp,&inputName,
			  &User_PromptCut,
			  &User_nEvts,
			  &User_batch,
			  &User_nTResidBins, &User_minTResid, &User_maxTResid,
			  &User_nDBins, &User_minD, &User_maxD,
			  &User_fPDF);

  const int PromptCut = SetDefValue(User_PromptCut, 15);

  const int nTResidBins = SetDefValue(User_nTResidBins, 101);
  const double minTResid = SetDefValue(User_minTResid, -0.5); // ns
  const double maxTResid = SetDefValue(User_maxTResid, 100.5); // ns

  const int nDBins = SetDefValue(User_nDBins, 1001);
  const double minD = SetDefValue(User_minD, -0.5); // mm
  const double maxD = SetDefValue(User_maxD, 16000.5); // mm

  // Useful to separate events if there is not time associated.
  const int DAQWindow = 256; // 64*4ns samples (typical PMT signal)

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


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE HIT COLLECTION              #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  vector<HitCollection> Evts;


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                   CREATE PDFs                     #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  int nEvtWithHits = 0;

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

		  if(T<DAQWindow){

			Hit hit(PMTPos, Q, T);
			vHit.emplace_back(hit);

			hTResiduals->Fill(T - (TVector3(PMTPos-TrueOrigin).Mag())/C,
							  1/(double)(NbPhotonCounts*NbPMTsHits));

		  }

		}

	  } // END FOR iPMT

	  // ################ //
	  // APPLY PROMPT CUT //
	  // ################ //

	  if(vHit.size() > 0){
		sort(vHit.begin(),vHit.end());
		auto itHit = find_if(vHit.begin(),
							 vHit.end(),
							 HitCut(Hit(TVector3(0,0,0),
										0,PromptCut+vHit.begin()->GetT())));
		vHit.erase(itHit,vHit.end());

		Evts.emplace_back(HitCollection(vHit,Form("Evt%d",iEvt)));
	  }

	  // cout << "Evt#" << iEvt << endl;
	  // for(auto h : vHit){
	  //   h.Print();
	  // }

	  // display the bar
	  progressBar.display();

	} // END FOR iEvt

  } else {

	cout << "File: " << ExtractFilenameFromPath(inputName) << " not open somehow!?" << endl;

  } // END IF FILE MC

  cout << endl;
  fMC->Close();

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                        FIT                        #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  bool DoYouWannaFit = true; // A vos risques et perils

  const double minPosRes = -1000.5;
  const double maxPosRes = 1000.5;
  const int NbBinsPosRes = 101; // cm resolution

  TH1D *hPosRes_X = nullptr;
  TH1D *hPosRes_Y = nullptr;
  TH1D *hPosRes_Z = nullptr;

  TH1D *hPosRes = nullptr;

  if(DoYouWannaFit){

	cout << "FIT " << endl;

	hPosRes_X = new TH1D("hPosRes_X","Vtx reconstruction", NbBinsPosRes, minPosRes, maxPosRes);
	hPosRes_Y = new TH1D("hPosRes_Y","Vtx reconstruction", NbBinsPosRes, minPosRes, maxPosRes);
	hPosRes_Z = new TH1D("hPosRes_Z","Vtx reconstruction", NbBinsPosRes, minPosRes, maxPosRes);

	hPosRes = new TH1D("hPosRes","VtxReconstruction",NbBinsPosRes,-0.5,maxPosRes);

	hTResiduals->Scale(1/(double)(Evts.size()));

	ProgressBar progressBar(Evts.size(), 70);

	for(auto itEvt: Evts){

	  // record the tick
	  ++progressBar;

	  // Set seed
	  TRandom3 r(0);

	  const int nWalk = 100;

	  // Set NLL array
	  vector<double> NLL;
	  vector<TVector3> POS_GUESS;

	  // Set seed between [minGuess, maxGuess]mm
	  const double minGuess = -1000; // mm
	  const double maxGuess = 1000; // mm

	  for(int iWalk=0; iWalk<nWalk; iWalk++){

		TVector3 X_GUESS(r.Uniform(minGuess,maxGuess),
						 r.Uniform(minGuess,maxGuess),
						 r.Uniform(minGuess,maxGuess));

		double NLL_val = WalkAndEvalNLL(X_GUESS, itEvt, hTResiduals);

		NLL.emplace_back(NLL_val);
		POS_GUESS.emplace_back(X_GUESS);

	  }

	  int maxNLL = distance(NLL.begin(),max_element(NLL.begin(),NLL.end()));

	  hPosRes_X->Fill(POS_GUESS[maxNLL].x());
	  hPosRes_Y->Fill(POS_GUESS[maxNLL].y());
	  hPosRes_Z->Fill(POS_GUESS[maxNLL].z());

	  hPosRes->Fill(POS_GUESS[maxNLL].Mag());

	  // display the bar
	  progressBar.display();

	}

	cout << endl;

  }

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                      DRAWING                      #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //


  auto *c1 = new TCanvas("cTResidsPDF","cTResidsPDF", 800,600);
  c1->SetGrid();
  hTResiduals->Draw();

  c1 = new TCanvas("cVTX","cVTX", 800,600);
  c1->SetGrid();
  TF1 *fGaus;
  fGaus = new TF1("fGaus",fitGaus,minPosRes,maxPosRes,3);
  fGaus->SetParName(0,"C0");
  fGaus->SetParName(1,"mu");
  fGaus->SetParName(2,"sigma");

  if(hPosRes_X){
	SetBasicTH1Style(hPosRes_X, kBlue);
	hPosRes_X->Draw();
	FitPos(hPosRes_X, fGaus);
	TFitResultPtr r = hPosRes_X->Fit("fGaus", "RLQE");
	if(r == 0){
	  SetFitStyle(hPosRes_X,kBlue-4);
	  cout << "X" << endl;
	  PrintFitResult(fGaus);
	}
  }
  if(hPosRes_Y){
	SetBasicTH1Style(hPosRes_Y, kRed);
	hPosRes_Y->Draw("SAME");
	FitPos(hPosRes_Y, fGaus);
	TFitResultPtr r = hPosRes_Y->Fit("fGaus", "RLQE");
	if(r == 0){
	  SetFitStyle(hPosRes_Y,kRed-4);
	  cout << "Y" << endl;
	  PrintFitResult(fGaus);
	}
  }
  if(hPosRes_Z){
	SetBasicTH1Style(hPosRes_Z, kGreen);
	hPosRes_Z->Draw("SAME");
	FitPos(hPosRes_Z, fGaus);
	TFitResultPtr r = hPosRes_Z->Fit("fGaus", "RLQE");
	if(r == 0){
	  SetFitStyle(hPosRes_Z,kGreen+1);
	  cout << "Z" << endl;
	  PrintFitResult(fGaus);
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