//
// Created by zsoldos on 11/26/19.
//

///////////////////////// STL C/C++ /////////////////////////
#include <vector>
#include <algorithm>
#include <csignal>

/////////////////////////   ROOT   //////////////////////////
#include <TApplication.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>

/////////////////////////   USER   //////////////////////////
#include "TFileAnalysis.hh"
#include "AnalysisDefinitions.hh"
#include "utils.hh"

#include "CreateEResMatrix.hh"

using namespace std;

int main(int argc, char *argv[]) {

  EoF=0;
  signal(SIGINT,Interrupt);

  // Create TApp
  TApplication theApp("App", &argc, argv);

  // Set basic ROOT style
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kDarkRainBow);

  // Input parameters
  // There's default value chosen by ternaries below
  string inputName;
  int User_nPEBins=-1;
  double User_minPE=-1; double User_maxPE=-1;
  int User_nPMTBins=-1;
  double User_minPMT=-1; double User_maxPMT=-1;
  string User_fOutput;
  int User_nThresh=-1;

  // READ input parameters
  ProcessArgs(&theApp,&inputName,
			  &User_nPEBins, &User_minPE, &User_maxPE,
			  &User_nPMTBins, &User_minPMT, &User_maxPMT,
			  &User_nThresh,
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
  TH1D *hNHits;

  TH2D *hEMatrix;
  vector< TH1D* > hNbPEVSE;
  vector< TH1D* > hNHitsVSE;

  TGraphErrors *grYieldVSE;

  TGraphErrors *grEff;

  sort(Ebins.begin(),Ebins.end());
  vector<double> corEbins = CorrectBinRangeArray(Ebins);

  const int nbBinsPE = (User_nPEBins>-1) ? User_nPEBins : 1001;
  const double minPE = (User_minPE>-1) ? User_minPE : -0.5;
  const double maxPE = (User_maxPE>-1) ? User_maxPE : 1000.5;

  const int nbBinsPMT = (User_nPMTBins>-1) ? User_nPMTBins : 1001;
  const double minPMT = (User_minPMT>-1) ? User_minPMT : -0.5;
  const double maxPMT = (User_maxPMT>-1) ? User_maxPMT : 1000.5;

  const int nThresh = (User_nThresh>-1) ? User_nThresh : 20;

  hEMatrix = new TH2D("hEMatrix", "Transition matrix describing the energy response of the detector",
					  Ebins.size(),&corEbins[0],
					  nbBinsPMT,minPMT,maxPMT);

  grYieldVSE = new TGraphErrors();
  grYieldVSE->SetName("grYieldVSE");
  grYieldVSE->SetTitle("Yield VS E ; E_{kin} positron (MeV) ; Yield (MeV^{-1})");
  grYieldVSE->SetMarkerColor(kBlue-4);
  grYieldVSE->SetMarkerStyle(kFullCrossX);
  grYieldVSE->SetMarkerSize(1);

  grEff = new TGraphErrors();
  grEff->SetName("grEff");
  grEff->SetTitle(Form("Eff VS E for nThresh=%d ; E_{kin} positron (MeV) ; Eff (%)", nThresh));
  grEff->SetMarkerColor(kRed-4);
  grEff->SetMarkerStyle(kFullCrossX);
  grEff->SetMarkerSize(1);

  TF1 *fPois = new TF1("fPois",
					   "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)",
					   0, maxPMT);
  fPois->SetParName(0,"alpha");
  fPois->SetParName(1,"lambda");
  fPois->SetParName(2,"mu");

  TCanvas *c1;

  auto fOutputName = (!User_fOutput.empty()) ? User_fOutput : "output";
  auto *fOutput = new TFile(Form("%s.root", fOutputName.c_str()),"RECREATE");

  auto *leg = new TLegend(0.83,0.38,0.99,0.99);

  for(auto& file : FileAnalysis){

    if(EoF) break;

    cout << "PROCESSING atm " << ExtractFilenameFromPath(file.GetFilename()) << endl;

	hNbPEVSHits = new TH2D(Form("hEbin%.1f", file.GetEBin()),"Nb PE Collected VS Nb PMTs Hits",
						   nbBinsPE,minPE,maxPE,
						   nbBinsPMT,minPMT,maxPMT);

	file.SetHist(hNbPEVSHits);

//	file.DoAnalysis(CollectPEAndHits);
// 	file.DoAnalysis(CollectPromptPEAndHits);
	file.DoAnalysisEV(CollectEVPEAndHits);

	//////////////////////////////////////
	// Recover hNPE/hNHits projection   //
	//////////////////////////////////////

	hNbPE = (TH1D*)file.GetHist()->ProjectionX()->Clone();
	hNbPE->SetTitle("Nb PE collected per events ; Nb PE collected ;");
	hNHits = (TH1D*)file.GetHist()->ProjectionY()->Clone();
	hNHits->SetTitle("Nb PMTHits per events ; NHit ;");

	grEff->SetPoint(grEff->GetN(),
					file.GetEBin(),
					hNHits->Integral(hNHits->FindBin(nThresh), nbBinsPMT)/hNHits->Integral());
	grEff->SetPointError(grEff->GetN()-1, 0., 0);


	hNbPEVSE.emplace_back(hNbPE);
	hNHitsVSE.emplace_back(hNHits);

	leg->AddEntry(hNHits, Form("%.1fMeV", file.GetEBin()));

	///////////////////////////////
	// Save all plots in ROOT    //
	///////////////////////////////

	fOutput->cd();
	file.GetHist()->Write();
	hNbPE->Write();
	hNHits->Write();

	///////////////////////////////
	// FILL transition matrix    //
	///////////////////////////////

	for(int iBin=1; iBin<nbBinsPMT; iBin++){

	  double BinPMT = hNHits->GetBinCenter(iBin);
	  double BinPMTContent = hNHits->GetBinContent(iBin);
	  double BinE = file.GetEBin();

	  hEMatrix->Fill(BinE, BinPMT, BinPMTContent);

	}

	///////////////////////////////
	// FIT and recover res. info //
	///////////////////////////////

	int NbEvents = file.GetNbEntries();
	fPois->SetParameters(1, 1, 1); // you MUST set non-zero initial values for parameters
	fPois->SetParLimits(0,0.,NbEvents);
	fPois->SetParLimits(1,0.,maxPE);
	fPois->SetParLimits(2,0.,NbEvents);

	if(hNHits->GetEntries() > 0){

	  TFitResultPtr r = hNHits->Fit("fPois", "R");
	  TF1 *fFit;

	  if(r == 0 && file.GetEBin() < 0.0 ) {

		hNHits->GetFunction("fPois")->SetLineColor(kBlue-4);
		hNHits->GetFunction("fPois")->SetLineWidth(1.5);
		hNHits->GetFunction("fPois")->SetLineStyle(2);
		fFit = hNHits->GetFunction("fPois");

	  } else {

		r = hNHits->Fit("gaus");

		if(r == 0 ){

		  hNHits->GetFunction("gaus")->SetLineColor(kRed-4);
		  hNHits->GetFunction("gaus")->SetLineWidth(1.5);
		  hNHits->GetFunction("gaus")->SetLineStyle(2);
		  fFit = hNHits->GetFunction("gaus");

		}

	  }

	  double chi2 = fFit->GetChisquare();
	  double alpha = fFit->GetParameter(0);
	  double mean = fFit->GetParameter(1);
	  double meanErr = fFit->GetParError(1);
	  double sigma = fFit->GetParameter(2);
	  double sigmaErr = fFit->GetParError(2);

	  double Yield = mean / file.GetEBin();
	  double YieldErr = meanErr / file.GetEBin();

	  if(r == 0){
		grYieldVSE->SetPoint(grYieldVSE->GetN(), file.GetEBin(), Yield);
		grYieldVSE->SetPointError(grYieldVSE->GetN()-1, 0., YieldErr);
	  }

	}

  }

  EoF=1;

  hEMatrix->Write();
  grYieldVSE->Write();
  grEff->Write();

  c1 = new TCanvas("cEMatrix", "cEMatrix", 800,600);
  c1->SetGrid();
  hEMatrix->Draw("COLZ");
  hEMatrix->Fit("pol1");
  TF1 *fFit = hEMatrix->GetFunction("pol1");
  double p0 = fFit->GetParameter(0);
  double p0Err = fFit->GetParError(0);
  double p1 = fFit->GetParameter(1);
  double p1Err = fFit->GetParError(1);
  TLatex *lFitResults = new TLatex();
  lFitResults->SetTextAlign(12);
  lFitResults->SetTextSize(0.04);
  lFitResults->DrawLatex(1,maxPMT*0.8,Form("%.2f #pm %.2f NHits/MeV",p1,p1Err));
  c1->Print(Form("%s_EMatrix.pdf", fOutputName.c_str()));

  c1 = new TCanvas("cYield", "cYield", 800,600);
  c1->SetGrid();
  grYieldVSE->Draw("AP");
  c1->Print(Form("%s_Yield.pdf", fOutputName.c_str()));

  c1 = new TCanvas("cEff", "cEff", 800,600);
  c1->SetGrid();
  grEff->Draw("AP");
  c1->Print(Form("%s_Eff.pdf", fOutputName.c_str()));

  c1 = new TCanvas("cNPEVSE", "cNPEVSE", 800,600);
  c1->SetGrid();
  for(auto& h : hNHitsVSE) {
    if(h->GetEntries() > 0)
	  h->Draw("SAME PLC PMC");
  }
  leg->Draw();
  c1->Write();
  c1->Print(Form("%s_NPEVSE.pdf", fOutputName.c_str()));

  /////////////////////////
  // ...

  fOutput->Close();

  cout << endl;
  cout << "Hit Ctrl+C to exit" << endl;
  theApp.Run(kTRUE);

  return 0;
}
