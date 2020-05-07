///////////////////////// STL C/C++ /////////////////////////
#include <climits>
#include <numeric>
#include <csignal>

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
#include <TApplication.h>
#include <TLine.h>

/////////////////////////   RAT   ///////////////////////////

/////////////////////////   USER   ///////////////////////////
// PUT UTILS.HH FIRST if you want to use ROOT CINT
// Because other headers might use function in utils.hh
#include "utils.hh"

#include "Analyzer.hh"
#include "AnalyzerFunctions.hh"
#include "HitClass.hh"

#include "EVFunctions.hh"
#include "HitFunctions.hh"
#include "LL.hh"
#include "CalibFunctions.hh"
#include "MCFunctions.hh"

#include "ProgressBar.hpp"

#define SQRT5 2.2360679775
#define SQRT2 1.41421356237

using namespace std;

static double GenEBin(){

  static int iBin=0;
  static double EBin = 10.;
  return EBin*iBin++;

}

class MCPDF{

 private:

  static double GenEBin(){

	static int iBin=1;
	static double EBin = 0.1;
	return EBin*iBin++;

  }

 protected:

  vector<double> EBins;

  vector< pair<double, TH2D*> > vpPDF;

 public:

  MCPDF() = default;

  MCPDF(string sFileName) {

    auto *fPDF = TFile::Open(sFileName.c_str());

	EBins.resize(100);
	generate(EBins.begin(), EBins.end(), GenEBin);

    if(fPDF->IsOpen()){

	  TIter next(fPDF->GetListOfKeys());
	  TKey *key;
	  while ((key = (TKey*)next())) {

		TClass *cl = gROOT->GetClass(key->GetClassName());
		if (!cl->InheritsFrom("TH2")) continue;
		// X-axis: nPE
		// Y-axis: nHits
		TH2 *h = (TH2*)key->ReadObj();

		for(auto EBin:EBins){

		  if(strcmp(Form("hTResVSCosThetaE%.1fMeV", EBin),h->GetName()) == 0){

			TH2D *h2D = (TH2D*)fPDF->Get(Form("hTResVSCosThetaE%.1fMeV", EBin));
			vpPDF.emplace_back(make_pair(EBin, h2D));

			break;

		  } // END if strcmp

		} // END for EBin

	  } // END while Key

	}

    sort(vpPDF.begin(), vpPDF.end());

  }

  TH2D *GetHPDF(double E){

    double roundedE = round(E*10)/10;

    vector< pair<double, TH2D*> >::iterator itLower;
    itLower = lower_bound(vpPDF.begin(), vpPDF.end(), make_pair(roundedE, new TH2D()));
	vector< pair<double, TH2D*> >::iterator itHigher;
	itHigher = upper_bound(vpPDF.begin(), vpPDF.end(), make_pair(roundedE, new TH2D()));

	auto lowBound = roundedE-itLower->first;
	auto highBound = itHigher->first-roundedE;

	if(lowBound<highBound){
	  return itLower->second;
	} else if (lowBound==highBound) {
	  return itLower->second;
	} else {
	  return itHigher->second;
	}

  }

};

class ZVector3 : public TVector3{

  bool operator<(const ZVector3 &rhs) const { return this->Mag2() < rhs.Mag2(); }
  bool operator>(const ZVector3 &rhs) const { return this->Mag2() > rhs.Mag2(); }

};

int GetIEvtIt(){
  static int iEvt = 0;
  return iEvt++;
}

void NextHisto(const char *command, TH1D *h)
{
   cout<<command<<endl;
   if (command[0]=='q') gApplication->Terminate(0);

}

static double GenCTBin(){

  static int iCTBin = 0;
  static double CTBin = 0;
  return cos(CTBin + (10 * iCTBin++));

}

double GetTrueCosTCer(double n=1.33, double E=2, double m=0.511){
  return 1 / (n*sqrt(1 - pow(m,2)/pow(E,2)));
}

void TestMacro(int User_NEvts = INT_MIN, int User_iEvt = INT_MIN){

  // Get Signal if user wants to interrupt loop
  EoF=0;
  signal(SIGINT,Interrupt);

  loadlibs();

  SetBasicStyle();

  string sCalFile = "Theia_1kT_Cov30pct_wbls_1pct_e-_Pos_0_0_0_Dir_0_0_1_EMatrix.root.root";
  MCCalib Cal(sCalFile);

  string sPDFFile = "Theia_1kT_Cov30pct_wbls_1pct_e-_Pos_0_0_0_Dir_0_0_1_PDF.root";
  MCPDF PDF(sPDFFile);

  string sPathToFile = "../rat-pac/";
  string sGeom = "Theia_1kT_Cov30pct";
  string sMaterial = "wbls_5newpct";
  string sParticle = "proton";
  string sE = "10.0MeV";
  string sInit = "Pos_0_0_0_Dir_0_0_1";

  //../MC/Theia_1kT_Cov30pct_wbls_1pct_e-_Pos_0_0_0_Dir_1_0_0
  string sRootPath = sPathToFile + sGeom + "_" + sMaterial + "_" + sParticle + "_" + sInit + "/";
  sRootPath = sPathToFile;
  string sFilename = sRootPath + sGeom+ "_" + sMaterial + "_" + sParticle + "_" + sE + "_" + sInit + ".root";

  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_water_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_water_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_1pct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_1pct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_5newpct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_5newpct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_10newpct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // files.emplace_back("../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_10newpct_e-_2Mev_Pos_0_0_0_Dir_1_0_0.root");
  // sFilename = "../MC/Theia_1kT_Cov30pct_2Mev_Pos_0_0_0_Dir_1_0_0/Theia_1kT_Cov30pct_wbls_5newpct_e+_2Mev_Pos_0_0_0_Dir_1_0_0.root";

  string sOutname = sGeom+ "_" + sMaterial + "_" + sParticle + "_" + sE + "_" + sInit + "_PMatrix.root";

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE OUTPUTS                     #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *foutput = new TFile(sOutname.c_str(),"RECREATE");

  auto *FileAnalyzer = new Analyzer(sFilename.c_str());

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE HISTOGRAMS                  #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  vector<double> vCTBins(12); // slices of 15deg
  generate(vCTBins.begin(), vCTBins.end(), GenCTBin);
  sort(vCTBins.begin(), vCTBins.end());
  // vector<double> corCTBins = CorrectBinRangeArray(vCTBins);
  //
  // for(auto b:corCTBins)
  //   cout << b << endl;

  auto hTRes = new TH1D("hTResids", "TResids",
						35, -5, 30);
  auto hCT = new TH1D("hCT", "Cos#theta hits",
					  vCTBins.size()-1,&vCTBins[0]);
  auto h2D = new TH2D("h2D", "TResids VS cos#theta hits",
					  35, -5, 30,
					  vCTBins.size()-1,&vCTBins[0]);

  unsigned long NEvts = User_NEvts > INT_MIN ? User_NEvts : FileAnalyzer->GetNEvts();
  unsigned long iEvt = SetDefValue(User_iEvt, 0);
  NEvts = User_iEvt > 0 ? User_iEvt + NEvts : NEvts;
  ProgressBar progressBar(NEvts, 70);

  TCanvas *c1;

  TLine *lCosTheta;

  for(iEvt; iEvt<NEvts; iEvt++) {

	// record the tick
	++progressBar;

	FlatParticle PrimPart = GetPrimaryParticleInfo(FileAnalyzer, iEvt);
	double PrimPartKE = PrimPart.GetKinE();
	TVector3 TrueOrigin = PrimPart.GetPos();
	TVector3 TrueDir = PrimPart.GetDir();

	hTRes->Reset();
	hCT->Reset();
	h2D->Reset();

	vector<Hit> vHit = GetEVHitCollection(FileAnalyzer, iEvt);
	// vector<Hit> vHit = GetVHitsFromPart(FileAnalyzer, iEvt, "proton");
	FillAll(hTRes, hCT, h2D, vHit, TrueOrigin, 224.9, TrueDir);

	lCosTheta = new TLine(GetTrueCosTCer(1.33, PrimPartKE, 0.511), 0,
						  GetTrueCosTCer(1.33, PrimPartKE, 0.511), hCT->GetMaximum());
	lCosTheta->SetLineColor(kRed-4);
	lCosTheta->SetLineWidth(2);
	lCosTheta->SetLineStyle(2);

	c1 = new TCanvas(Form("cEvt%lu", iEvt), Form("cEvt%lu", iEvt), 1200,800);
	c1->Divide(2,2);
	c1->cd(1);
	hTRes->Draw();
	c1->cd(2);
	gPad->SetGrid();
	hCT->Draw();
	// lCosTheta->Draw("SAME");
	c1->cd(3);
	c1->SetGrid();
	gPad->SetGrid();
	h2D->Draw("COLZ");

	c1->WaitPrimitive();

	delete c1;
	delete lCosTheta;

  }

  EoF = 1;

  foutput->cd();

  foutput->Close();

}
