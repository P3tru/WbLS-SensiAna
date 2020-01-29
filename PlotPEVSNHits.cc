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
#include "PlotPEVSNHits.hh"
#include "AnalysisDefinitions.hh"
#include "utils.hh"

using namespace std;

int main(int argc, char *argv[]) {

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

  TFileAnalysis<TH2D> FileAnalysis(inputName);

  const int nbBinsPE = (User_nPEBins>-1) ? User_nPEBins : 201;
  const double minPE = (User_minPE>-1) ? User_minPE : -0.5;
  const double maxPE = (User_maxPE>-1) ? User_maxPE : 200.5;

  const int nbBinsPMT = (User_nPMTBins>-1) ? User_nPMTBins : 201;
  const double minPMT = (User_minPMT>-1) ? User_minPMT : -0.5;
  const double maxPMT = (User_maxPMT>-1) ? User_maxPMT : 200.5;

  const int nThresh = (User_nThresh>-1) ? User_nThresh : 20;

  auto fOutputName = (!User_fOutput.empty()) ? User_fOutput : "output";
  auto *fOutput = new TFile(Form("%s.root", fOutputName.c_str()),"RECREATE");


  TH2D *hNbPEVSHits = new TH2D("hNbPEVSHits","NbPE VS NbHits",
							   nbBinsPE,minPE,maxPE,
							   nbBinsPMT,minPMT,maxPMT);
  TH1D *hNbPE;
  TH1D *hNHits;

  FileAnalysis.SetHist(hNbPEVSHits);

  FileAnalysis.DoAnalysis(CollectPEAndHits);

  hNbPE = (TH1D*)FileAnalysis.GetHist()->ProjectionX()->Clone();
  hNbPE->SetTitle("Nb PE collected per events ; Nb PE collected ;");
  hNHits = (TH1D*)FileAnalysis.GetHist()->ProjectionY()->Clone();
  hNHits->SetTitle("Nb PMTHits per events ; NHit ;");

  const double eff = hNHits->Integral(hNHits->FindBin(nThresh), nbBinsPMT)/hNHits->Integral();

  cout << endl;
  cout << "Eff for nThresh: " << nThresh << " is: " << eff << endl;

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                      DRAWING                      #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *c1 = new TCanvas("c1","c1",800,600);
  hNHits->Draw();
  c1->Print(Form("%s_NHits.pdf", fOutputName.c_str()));
  fOutput->cd();
  hNHits->Write();


  /////////////////////////
  // ...

  fOutput->Close();

  cout << endl;
  cout << "Hit Ctrl+C to exit" << endl;
  theApp.Run(kTRUE);

  return 0;

}