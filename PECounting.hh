//
// Created by zsoldos on 2/21/20.
//

#ifndef _PECOUNTING_HH_
#define _PECOUNTING_HH_

#include <math.h>

#include "HitClass.hh"
#include "HitFunctions.hh"

#include "ProgressBar.hpp"

// #include "utils.hh"

#define PI 3.14159

double rad2deg(double angrad){
  return angrad*180/PI;
}

double deg2rad(double angdeg){
  return angdeg*PI/180;
}

// Count hits inside the ring / outside the ring
bool IsHitInRing(TVector3 HitPos, TVector3 Dir, double cosThetaCer,
	double minRing = 15/*deg*/, double maxRing = 15/*deg*/){

  double cThetaHit = cos(Dir.Angle(HitPos));
  double up_cher = cos(acos(cosThetaCer) + deg2rad(minRing));
  double down_cher = cos(acos(cosThetaCer) - deg2rad(minRing));

  if(cThetaHit > up_cher && cThetaHit < down_cher) return true;
  else return false;

}

double GetTrueCosTCer(double n=1.33, double E=2, double m=0.511){
  return 1 / (n*sqrt(1 - pow(m,2)/pow(E,2)));
}

double GetFOM(double Cer, double Scint){
  return (Cer-Scint) / sqrt(pow(Cer,2) + pow(Scint,2));
}

double GetPID(double Cer, double Tot){
  return Cer/Tot;
}

double GetPIDErr(double Cer, double CerErr, double Tot, double TotErr){
  return (CerErr / Tot) + (TotErr * Cer / pow(Tot,2));
}

void GetCerFracInsideRing(double *CerFrac, double *CerFracErr,
						  TH1D *hCT, TH1D* hCTInside,
						  double *FOM=NULL){

  // Fit part far for the Cer hits -> scintillation
  const double minThetaScint = -1.0;
  const double maxThetaScint = - sqrt(3)/2;
  TFitResultPtr r = hCT->Fit("pol0","0","", minThetaScint, maxThetaScint);
  // Recover results of fit
  if( r == 0 ){

	double p0 = hCT->GetFunction("pol0")->GetParameter(0);
	double p0Err = hCT->GetFunction("pol0")->GetParError(0);

	// cout << "Scint Contribution: " << p0 << " +- " << p0Err << endl;

	TH1D *hCTHitsWOScint = (TH1D*)hCTInside->Clone("hWOScint");
	TF1 *fScint = hCT->GetFunction("pol0");
	hCTHitsWOScint->Add(fScint, -1.);

	double TotRingPhoton, TotRingPhotonErr;
	double CerRingPhoton, CerRingPhotonErr;
	TotRingPhoton = hCTInside->IntegralAndError(0,
												hCTInside->GetNbinsX(), TotRingPhotonErr);
	CerRingPhoton = hCTHitsWOScint->IntegralAndError(0,
													 hCTHitsWOScint->GetNbinsX(), CerRingPhotonErr);

	*CerFrac = CerRingPhoton/TotRingPhoton;
	*CerFracErr = CerRingPhotonErr/TotRingPhoton + TotRingPhotonErr*(CerRingPhoton/TotRingPhoton)/TotRingPhoton;

	cout << "Cer Frac: " << *CerFrac << " +- " << *CerFracErr << endl;

	*FOM = GetFOM(CerRingPhoton, TotRingPhoton-CerRingPhoton);

  }

}


void PECounting(const char* filename, TVector3 Dir, TGraphErrors *grErr=NULL) {

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE HISTOGRAMS                  #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *hCosThetaHits = new TH1D("hCosThetaHits",
								 "cos(#theta) between prtcle dir and hits",
								 60,-1.,1.);
  hCosThetaHits->GetXaxis()->SetTitle("cos(#theta)");

  auto *hHitTime = new TH1D("hHitTime", "Hit time", 101, -0.5, 100.5);
  hHitTime->GetXaxis()->SetTitle("Hit time (ns)");

  const int nbPromptCuts = 5;
  const int promptCuts[nbPromptCuts] = {10, 8, 6, 4, 2};

  TH1D *hCosThetaHitsVSPromptCut[nbPromptCuts];
  TH1D *hCosThetaHitsVSPromptCutInside[nbPromptCuts];

  auto *foutput = new TFile(Form("%s_plots.root",filename), "RECREATE");

  for(int iCut=0; iCut<nbPromptCuts; iCut++){

	hCosThetaHitsVSPromptCut[iCut]
		= new TH1D(Form("hCosThetaHits%dns",promptCuts[iCut]),
				   "cos(#theta) between prtcle dir and hits",
				   60,-1.,1.);
	hCosThetaHitsVSPromptCutInside[iCut]
		= new TH1D(Form("hCosThetaHitsInside%dns",promptCuts[iCut]),
				   "cos(#theta) between prtcle dir and hits",
				   60,-1.,1.);
	hCosThetaHitsVSPromptCutInside[iCut]->SetLineColor(kRed-4);

  }

  auto *fMC = new TFile(filename);

  if (fMC->IsOpen()) {

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

	auto nEvents = tree->GetEntries();
	cout << "Start Loop on nEvts: " << nEvents << endl;

	ProgressBar progressBar(nEvents, 70);

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

		hCosThetaHits->Fill(cos(Dir.Angle(PMTPos)), (double)(NbPhotonCounts)/(double)(nEvents));

		for (int iP = 0; iP < NbPhotonCounts; iP++) {

		  auto mcPhoton = mcPMT->GetMCPhoton(iP);
		  // Get Q and T
		  auto Q = mcPhoton->GetCharge();
		  auto T = mcPhoton->GetHitTime();

		  hHitTime->Fill(T);

		  // ########## //
		  // FILL EVENT //
		  // ########## //

		  vHit.emplace_back(Hit(PMTPos, Q, T));

		}

	  } // END FOR iPMT

	  // ############ //
	  // PROCESS HITS //
	  // ############ //

	  if (vHit.size() > 0) {

		sort(vHit.begin(), vHit.end());

		auto vHitCorrected = CorrectDelayedHits(vHit, Hit(TVector3(0.,0.,0.), 0., 0.));

		for(auto h: vHitCorrected){

		  for(int iCut=0; iCut<nbPromptCuts; iCut++){

		    if(h.GetT() < promptCuts[iCut]){

			  hCosThetaHitsVSPromptCut[iCut]->Fill(cos(Dir.Angle(h.GetPos())),
												   1/(double)(vHitCorrected.size()));

			  if(IsHitInRing(h.GetPos(), Dir, GetTrueCosTCer(1.33,2,0.511))){
				hCosThetaHitsVSPromptCutInside[iCut]->Fill(cos(Dir.Angle(h.GetPos())),
														   1/(double)(vHitCorrected.size()));
			  }

			}


		  }

		}

	  }

	  // display the bar
	  progressBar.display();

	} // END FOR iEVT

	cout << " DONE" << endl;

  } // END if file

  TGraphErrors *grPID = new TGraphErrors();
  grPID->SetName("grPID");
  grPID->SetTitle("Cer photon frac inside ring VS prompt Cut ; prompt cut (ns) ; Cer photon frac");
  TGraphErrors *grFOM = new TGraphErrors();
  grFOM->SetName("grFOM");
  grFOM->SetTitle("FOM: (Cer-Scint)/sqrt(Cer^2+Scint^2) ; prompt cut (ns) ; FOM");
  grFOM->SetLineColor(kRed-4);
  grFOM->SetLineStyle(2);

  for(int iCut=0; iCut<nbPromptCuts; iCut++) {
	double CerFrac, CerFracErr, FOM;
	GetCerFracInsideRing(&CerFrac, &CerFracErr,
						 hCosThetaHitsVSPromptCut[iCut], hCosThetaHitsVSPromptCutInside[iCut],
						 &FOM);

	if(grErr){

	  // grErr ->SetPoint(grErr->GetN(), promptCuts[iCut], CerFrac);
	  // grErr ->SetPointError(grErr->GetN()-1, 0., CerFracErr);

	  grErr->SetPoint(grErr->GetN(), promptCuts[iCut], FOM);
	  grErr->SetPointError(grErr->GetN()-1, 0., 0.);

	}

	grPID->SetPoint(grPID->GetN(), promptCuts[iCut], CerFrac);
	grPID->SetPointError(grPID->GetN()-1, 0., CerFracErr);

	grFOM->SetPoint(grFOM->GetN(), promptCuts[iCut], FOM);
	grFOM->SetPointError(grFOM->GetN()-1, 0., 0.);

  }


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                      DRAWING                      #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  if(!grErr){

    gStyle->SetPalette(kRainBow);

	auto *c1 = new TCanvas("cPIDVSCut","cPIDVSCut", 800,600);
	c1->SetGrid();
	grPID->Draw("APC");
	grFOM->Draw("SAME");
	grPID->GetYaxis()->SetRangeUser(0.,1.);


	c1 = new TCanvas("cCTHits","cCTHits", 800,600);
	c1->SetGrid();
	hCosThetaHits->Draw();

	c1 = new TCanvas("cHitime","cHitTime", 800,600);
	c1->SetGrid();
	hHitTime->Draw();

	for(int iCut=0; iCut<nbPromptCuts; iCut++){
	  c1 = new TCanvas(Form("cCTHitsVSPromptCut%d", promptCuts[iCut]),
					   Form("cCTHitsVSPromptCut%d", promptCuts[iCut]), 800,600);
	  c1->SetGrid();
	  hCosThetaHitsVSPromptCut[iCut]->Draw("SAME PMC PLC");
	  hCosThetaHitsVSPromptCutInside[iCut]->Draw("SAME");
	}

  }


  // foutput->cd();
  // grPID->Write();
  // grFOM->Write();
  // hCosThetaHits->Write();
  // hHitTime->Write();
  // for(auto h:hCosThetaHitsVSPromptCut)
	// h->Write();
  // for(auto h:hCosThetaHitsVSPromptCutInside)
	// h->Write();
  //
  // foutput->Close();

}

void ProcessFiles(){

  gStyle->SetPalette(kRainBow);

  vector<const char*> files;
  files.emplace_back("Theia_1kT_Cov50pct_water_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("Theia_1kT_Cov50pct_water_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("Theia_1kT_Cov50pct_wbls_1pct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("Theia_1kT_Cov50pct_wbls_1pct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("Theia_1kT_Cov50pct_wbls_5newpct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("Theia_1kT_Cov50pct_wbls_5newpct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("Theia_1kT_Cov50pct_wbls_10newpct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("Theia_1kT_Cov50pct_wbls_10newpct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  const int nFiles = files.size();

  vector<const char*> MCParams;
  MCParams.emplace_back("Water e+");
  MCParams.emplace_back("Water e-");
  MCParams.emplace_back("WbLS 1% e+");
  MCParams.emplace_back("WbLS 1% e-");
  MCParams.emplace_back("WbLS 5% e+");
  MCParams.emplace_back("WbLS 5% e-");
  MCParams.emplace_back("WbLS 10% e+");
  MCParams.emplace_back("WbLS 10% e-");
  TLegend *leg = new TLegend(0.867168,0.576389,0.99,0.99);

  TGraphErrors *grFOM[nFiles];
  TMultiGraph *mgFOM = new TMultiGraph();
  // mgFOM->SetTitle("FOM: (N_{Cer}-N_{Scint})/#sqrt{N_{Cer}^{2}+N_{Scint}^{2}}) ; Prompt cut (ns) ; FOM");
  mgFOM->SetTitle("Cer fraction inside ring ; Prompt cut (ns) ; FOM");

  auto *foutput = new TFile("e+e-_2MeV_CerFrac.root", "RECREATE");

  for(int iFile=0; iFile<nFiles; iFile++){

	grFOM[iFile] = new TGraphErrors();
	grFOM[iFile]->SetTitle("FOM: (Cer-Scint)/sqrt(Cer^2+Scint^2) ; prompt cut (ns) ; FOM");

	cout << "Processing file: " << files[iFile] << endl;
	PECounting(files[iFile], TVector3(1.,0.,0.), grFOM[iFile]);
	foutput->cd(); grFOM[iFile]->Write();
	mgFOM->Add(grFOM[iFile]);
	leg->AddEntry(grFOM[iFile], MCParams[iFile]);


  }

  auto *c1 = new TCanvas("cCerFrac","cCerFrac", 800,600);
  c1->SetGrid();
  mgFOM->Draw("A PMC PLC");
  leg->Draw();
  foutput->cd(); c1->Write();
  foutput->Close();
}

#endif //_PECOUNTING_HH_
