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
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TFile.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TTree.h>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/MC.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/Run.hh>

/////////////////////////   USER   //////////////////////////
#include "PECounting.hh"
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

  bool isBatch = false;

  // TODO : Allow user to set true origin vector from input
  const TVector3 TrueOrigin(0.,0.,0.);

  vector<const char*> files;
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_water_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_water_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_1pct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_1pct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_5newpct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_5newpct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_10newpct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_10newpct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../rat-pac/Theia_1kT_Cov50pct_water_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../rat-pac/Theia_1kT_Cov50pct_water_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../rat-pac/Theia_1kT_Cov50pct_wbls_1pct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../rat-pac/Theia_1kT_Cov50pct_wbls_1pct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../rat-pac/Theia_1kT_Cov50pct_wbls_5newpct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../rat-pac/Theia_1kT_Cov50pct_wbls_5newpct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../rat-pac/Theia_1kT_Cov50pct_wbls_10newpct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../rat-pac/Theia_1kT_Cov50pct_wbls_10newpct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("../MC/Theia_1kT_Cov70pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov70pct_water_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("../MC/Theia_1kT_Cov70pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov70pct_water_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("../MC/Theia_1kT_Cov70pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov70pct_wbls_1pct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("../MC/Theia_1kT_Cov70pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov70pct_wbls_1pct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("../MC/Theia_1kT_Cov70pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov70pct_wbls_5newpct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("../MC/Theia_1kT_Cov70pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov70pct_wbls_5newpct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("../MC/Theia_1kT_Cov70pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov70pct_wbls_10newpct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  files.emplace_back("../MC/Theia_1kT_Cov70pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov70pct_wbls_10newpct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
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


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE HISTOGRAMS                  #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  TGraphErrors *grFOM[nFiles];
  TMultiGraph *mgFOM = new TMultiGraph();
  mgFOM->SetTitle("FOM: (N_{Cer}-N_{Scint})/#sqrt{N_{Cer}^{2}+N_{Scint}^{2}}) ; Prompt cut (ns) ; FOM");
  // mgFOM->SetTitle("Cer fraction inside ring ; Prompt cut (ns) ; FOM");


  auto *foutput = new TFile("e+e-_2MeV_FOM_70pct.root", "RECREATE");

  for(int iFile=0; iFile<nFiles; iFile++){

	grFOM[iFile] = new TGraphErrors();
	grFOM[iFile]->SetTitle("FOM: (Cer-Scint)/sqrt(Cer^2+Scint^2) ; prompt cut (ns) ; FOM");

	cout << "Processing file: " << files[iFile] << endl;
	PECounting(files[iFile], TVector3(1.,0.,0.), grFOM[iFile]);
	foutput->cd(); grFOM[iFile]->Write();
	mgFOM->Add(grFOM[iFile]);
	leg->AddEntry(grFOM[iFile], MCParams[iFile]);


  }


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                      DRAWING                      #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //


  auto *c1 = new TCanvas("cCerFrac","cCerFrac", 800,600);
  c1->SetGrid();
  mgFOM->Draw("A PMC PLC");
  leg->Draw();
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