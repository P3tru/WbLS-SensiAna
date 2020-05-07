///////////////////////// STL C/C++ /////////////////////////

/////////////////////////   BOOST   /////////////////////////

/////////////////////////   ROOT   //////////////////////////
#include <TVector3.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TFile.h>

/////////////////////////   RAT   ///////////////////////////

/////////////////////////   USER   ///////////////////////////
#include "Analyzer.hh"
#include "HitClass.hh"

#include "EVFunctions.hh"
#include "HitFunctions.hh"
#include "LL.hh"
#include "CalibFunctions.hh"
#include "MCFunctions.hh"

#include "ProgressBar.hpp"
#include "utils.hh"

#define SQRT5 2.2360679775
#define SQRT2 1.41421356237

using namespace std;

R__LOAD_LIBRARY(libEVFunctions.so)
R__LOAD_LIBRARY(libHitFunctions.so)
R__LOAD_LIBRARY(libLL.so)
R__LOAD_LIBRARY(libCalibFunctions.so)
R__LOAD_LIBRARY(libMCFunctions.so)


double CalculateProb(TH2D *hPDF, TH2D *hExp, double *Chi2 = NULL, int *NdF = NULL){

  auto nBinsX = hPDF->GetNbinsX();
  auto nBinsY = hPDF->GetNbinsY();

  auto N = hExp->Integral();
  auto W = hPDF->GetSumOfWeights();
  auto normW = 2*W*W;

  double chi2 = 0.;
  int NonNullBin = 0;

  for(auto iBinX=1; iBinX<=nBinsX; iBinX++){
	for(auto iBinY=1; iBinY<=nBinsY; iBinY++) {

	  double n = hExp->GetBinContent(hExp->GetBin(iBinX, iBinY));
	  double w = hPDF->GetBinContent(hPDF->GetBin(iBinX, iBinY));

	  if(n == 0 || w == 0) continue;

	  double s2 = pow(hPDF->GetBinError(hPDF->GetBin(iBinX, iBinY)),2);
	  double res = W*w - N*s2;

	  double P = res + sqrt(pow(res, 2) + 4*W*W*s2*s2*n);
	  if(P == 0) continue;
	  P/=normW;

	  chi2+=pow(n - N*P,2)/(N*P) + pow(w - W*P,2)/s2;
	  NonNullBin++;

	}
  }

  if(Chi2)
	*Chi2=chi2;
  if(NdF)
	*NdF=NonNullBin-1;

  return TMath::Prob(chi2, NonNullBin-1);


}


void TestMacro(int User_NEvts = -1, int User_iEvt = -1){

  string sPathToFile = "../MC/";
  string sGeom = "Theia_1kT_Cov90pct";
  string sMaterial = "wbls_1pct";
  string sParticle = "e-";
  string sE = "2.6Mev";
  string sInit = "Pos_0_0_0_Dir_1_0_0";

  //../MC/Theia_1kT_Cov30pct_wbls_1pct_e-_Pos_0_0_0_Dir_1_0_0
  string sRootPath = sPathToFile + sGeom + "_" + sMaterial + "_" + sParticle + "_" + sInit + "/";
  sRootPath = sPathToFile;
  string sFilename = sRootPath + sGeom+ "_" + sMaterial + "_" + sParticle + "_" + sE + "_" + sInit + ".root";

  auto *FileAnalyzer = new Analyzer(sFilename.c_str());

  unsigned long int NEvts = User_NEvts > 0 ? User_NEvts : FileAnalyzer->GetNEvts();
  auto iEvt = SetDefValue(User_iEvt, 0);
  NEvts = User_iEvt > 0 ? User_iEvt + NEvts : NEvts;
  ProgressBar progressBar(NEvts, 70);

  const int nTResBins = 35;
  const double minTRes = -5;
  const double maxTRes = 30;

  const int nCThetaBins = 45;
  const double minCTheta = -1;
  const double maxCTheta = 1;

  auto *hTResVSCosTheta = new TH2D("hTResVSCosTheta", "hTRes VS Cos #theta",
								   nTResBins,minTRes,maxTRes,
								   nCThetaBins,minCTheta,maxCTheta);

  auto *hProb = new TH1D("hProb", "GoF", 10000, 0., 100.);

  hTResVSCosTheta->Sumw2();
  TVector3 TrueOrigin = TVector3(0.,0.,0.);
  TVector3 TrueDir = TVector3(1.,0.,0.);

  // Create 2D PDF

  for(iEvt; iEvt<NEvts; iEvt++) {

    // record the tick
    ++progressBar;

    vector<Hit> vHit = GetEVHitCollection(FileAnalyzer, iEvt, 0);

	FillTResVSCosTheta(hTResVSCosTheta, vHit, TrueOrigin, 224.9, TrueDir);

    // display the bar
    progressBar.display();

  }

  hTResVSCosTheta->Scale(1/(double)(NEvts));

  // Calculate LL

  iEvt = SetDefValue(User_iEvt, 0);
  NEvts = User_iEvt > 0 ? User_iEvt + NEvts : NEvts;
  // iEvt = 0;
  // NEvts = 1;

  ProgressBar progressBar2(NEvts, 70);
  for(iEvt; iEvt<NEvts; iEvt++) {

	// record the tick
	++progressBar2;

	vector<Hit> vHit = GetEVHitCollection(FileAnalyzer, iEvt, 0);

	auto *hTResVSCosTheta_Evt = new TH2D(Form("hTResVSCosThetaEvt%d", iEvt), "hTRes VS Cos #theta",
										 nTResBins,minTRes,maxTRes,
										 nCThetaBins,minCTheta,maxCTheta);

	FillTResVSCosTheta(hTResVSCosTheta_Evt, vHit, TrueOrigin, 224.9, TrueDir);

	// cout << CalculateProb(hTResVSCosTheta, hTResVSCosTheta_Evt) << endl;

	double Chi2 = -1;
	int NdF = -1;
	double Prob = CalculateProb(hTResVSCosTheta, hTResVSCosTheta_Evt, &Chi2, &NdF);

	hProb->Fill(Chi2/(double)(NdF));

	delete hTResVSCosTheta_Evt;

	// display the bar
	progressBar2.display();

  }


  auto *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->SetGrid();
  c1->SetLogz();
  hTResVSCosTheta->Draw("COLZ");

  c1 = new TCanvas("c2", "c2", 800, 600);
  c1->SetGrid();
  hProb->Draw();

}
