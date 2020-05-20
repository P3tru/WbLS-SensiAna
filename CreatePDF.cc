//
// Created by zsoldos on 12/5/19.
//

///////////////////////// STL C/C++ /////////////////////////
#include <string>
#include <vector>
#include <csignal>
#include <numeric>
#include <climits>
#include <fstream>


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
  string listName;
  string outputName;

  auto User_isBatch = false;

  auto User_nEvts = INT_MIN;
  auto User_iEvt = INT_MIN;

  int User_nTResidBins = INT_MIN;
  auto User_minTResid = DBL_MIN; auto User_maxTResid= DBL_MIN;

  int User_nCThetaBins = INT_MIN;
  auto User_minCTheta = DBL_MIN; auto User_maxCTheta = DBL_MIN;

  // Select Prompt Window or DAQ Window
  auto User_PromptWindow = INT_MIN;

  // Select time cut for computing residuals
  auto User_PromptCut = INT_MIN;

  // Scaling?
  auto User_isScaling = true;

  ProcessArgs(&theApp,
			  &inputName,&listName,
			  &User_PromptWindow,
			  &User_PromptCut,
			  &User_nEvts, &User_iEvt,
			  &User_nTResidBins, &User_minTResid, &User_maxTResid,
			  &User_nCThetaBins, &User_minCTheta, &User_maxCTheta,
			  &User_isBatch,
			  &outputName);

  const bool isBatch = User_isBatch;

  const bool isScaling = User_isScaling;

  // Add a prompt window, because there's no separate trigger for 2 evts in rat-pac.
  // The IBD generator will not create 2 separate EV object.
  // Therefore, add a cut to physically separate prompt+delay

  // BY DEFAULT, let's not make any cuts now.
  // But one will have the option in the future

  const int PromptWindow = SetDefValue(User_PromptWindow, -1); // ns
  const int PromptCut = SetDefValue(User_PromptCut, -1); // ns

  const int nbBinsTRes = SetDefValue(User_nTResidBins, 55);
  const double minTRes = SetDefValue(User_minTResid, -5.); // ns
  const double maxTRes = SetDefValue(User_maxTResid, 50.); // ns

  const int nCThetaBins = SetDefValue(User_nCThetaBins, 24);
  const double minCTheta= SetDefValue(User_minCTheta, -1.); // costheta
  const double maxCTheta = SetDefValue(User_maxCTheta, 1.); // costheta

  const string outName = outputName.empty() ? "Output_PDF.root" : outputName;


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE OUTPUTS                     #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *foutput = new TFile(outName.c_str(),"RECREATE");


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE ANALYZER                    #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  vector<Analyzer*> vFileAnalyzer;
  AddFAnalyzers(&vFileAnalyzer, inputName, listName);
  auto nAnalyzer = vFileAnalyzer.size();

  vector<double> EBins;

  for(auto fA:vFileAnalyzer) {

	double EBin = GetPrimaryParticleInfo(fA, 0).GetKinE();
	EBins.emplace_back(EBin);
	const char *cEBin = Form("E%.1fMeV", EBin);
	fA->SetTag(cEBin);

  }

  sort(EBins.begin(),EBins.end());
  vector<double> corEbins = CorrectBinRangeArray(EBins);

  vector<double> vCTBins(12); // slices of 15deg
  generate(vCTBins.begin(), vCTBins.end(), GenCTBin);
  sort(vCTBins.begin(), vCTBins.end());


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                DEFINE HISTOGRAMS                  #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  vector<TH1D*> vHTRes;
  vector<TH1D*> vHCT;
  vector< vector<TH1D*> > vvHTresCut;
  vector<TH2D*> vHTResVSCTheta;
  vector<TH1D*> vHNHits;
  vector<TH1D*> vHQ;


  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                   ANALYSIS                        #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  for(auto fA:vFileAnalyzer){

    cout << " PROCESSING file with tag: " << fA->GetTag() << endl;
    cout << endl;

    const char *cTag = fA->GetTag().c_str();

  	unsigned long nEvts = User_nEvts > INT_MIN && User_nEvts < fA->GetNEvts() ? User_nEvts : fA->GetNEvts();
	unsigned long iEvt = SetDefValue(User_iEvt, 0);

	// #### #### #### #### #### #### #### #### #### #### #### #### //
	// ####                CREATE HISTOGRAMS                  #### //
	// #### #### #### #### #### #### #### #### #### #### #### #### //

	auto *hTResiduals = new TH1D(Form("hTResiduals%s",cTag), "T Residuals",
								 nbBinsTRes, minTRes, maxTRes);

	SetBasicTH1Style(hTResiduals, kBlue-4);
	if(isScaling)
	  hTResiduals->Sumw2();

	auto *hCTheta = new TH1D(Form("hCTheta%s",cTag), "cos#theta hits",
							 // nCThetaBins, minCTheta, maxCTheta);
							 vCTBins.size()-1,&vCTBins[0]);

	SetBasicTH1Style(hCTheta, kBlue-4);
	if(isScaling)
	  hCTheta->Sumw2();

	const std::size_t nbCuts = 20;
	vector<int> vPromptCut(nbCuts);
	generate(vPromptCut.begin(), vPromptCut.end(), GenPromptCut);
	vector<TH1D *>vhTResiduals(nbCuts);
	vector<int> NEvtsProcessedCut(nbCuts, 0);

	for(auto i=0; i<nbCuts; i++){
	  vhTResiduals[i] = new TH1D(Form("hTResiduals%sCut%d",cTag, vPromptCut[i]), "T Residuals",
								 nbBinsTRes, minTRes, maxTRes);
	  SetBasicTH1Style(vhTResiduals[i], kRed-4);
	  if(isScaling)
		vhTResiduals[i]->Sumw2();
	}

	auto *hTResVSCosTheta = new TH2D(Form("hTResVSCosTheta%s",cTag), "hTRes VS Cos #theta",
									 nbBinsTRes, minTRes, maxTRes,
									 // nCThetaBins, minCTheta, maxCTheta);
									 vCTBins.size()-1,&vCTBins[0]);

	if(isScaling)
	  hTResVSCosTheta->Sumw2();

	auto *hNHits = new TH1D(Form("hNHits%s", cTag), "NHits",
							1000, 0., 1000.);
	auto *hQ = new TH1D(Form("hQ%s", cTag), "Charge",
						1000, 0., 1000.);

	// #### #### #### #### #### #### #### #### #### #### #### #### //
	// ####           LOOP AND FILL ANALYSIS                  #### //
	// #### #### #### #### #### #### #### #### #### #### #### #### //

	LoopAndFillHistos(fA,
					  nEvts, iEvt, PromptWindow, PromptCut,
					  hNHits, hQ,
					  hTResiduals, hCTheta,
					  hTResVSCosTheta,
					  isScaling);

	vHTRes.emplace_back(hTResiduals);
	vHCT.emplace_back(hCTheta);
	vvHTresCut.emplace_back(vhTResiduals);
	vHTResVSCTheta.emplace_back(hTResVSCosTheta);
	vHNHits.emplace_back(hNHits);
	vHQ.emplace_back(hQ);

  }

  EoF = 1;

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                      DRAWING                      #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  if(!isBatch){

    if(nAnalyzer == 1){

	  PlotAHist(vHTRes[0], "");

	  PlotAHist(vHCT[0],"");

	  PlotAHist(vHTResVSCTheta[0], "COLZ");

	  PlotAHist(vHNHits[0], "");

	  PlotAHist(vHQ[0], "");

    }


  }

  foutput->cd();

  for(int iA = 0; iA<nAnalyzer; iA++){

    vHTRes[iA]->Write();
    vHCT[iA]->Write();
	// for(auto h: vvHTresCut[iA]){
	//   h->Write();
	// }
	vHTResVSCTheta[iA]->Write();
	vHNHits[iA]->Write();
	vHQ[iA]->Write();
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
