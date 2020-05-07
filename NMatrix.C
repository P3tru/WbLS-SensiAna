///////////////////////// STL C/C++ /////////////////////////
#include <climits>
#include <numeric>

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
// PUT UTILS.HH FIRST if you want to use ROOT CINT
// Because other headers might use function in utils.hh
#include "utils.hh"

#include "Analyzer.hh"
#include "HitClass.hh"

#include "AnalyzerFunctions.hh"
#include "EVFunctions.hh"
#include "HitFunctions.hh"
#include "LL.hh"
#include "CalibFunctions.hh"
#include "MCFunctions.hh"

#include "ProgressBar.hpp"

#define SQRT5 2.2360679775
#define SQRT2 1.41421356237

using namespace std;

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

static double GenCTBin(){

  static int iCTBin = 0;
  static double CTBin = 0;
  return cos(CTBin + (10 * iCTBin++));

}

static double GenEBin(){

  static int iBin=1;
  static double EBin = 0.1;
  return EBin*iBin++;

}

double GetTrueCosTCer(double n=1.33, double E=2, double m=0.511){
  return 1 / (n*sqrt(1 - pow(m,2)/pow(E,2)));
}

class ProtonRecoilAnalysis {

 protected:
  double EBin;
  TH2D *hTResVSCosTheta;
  unsigned nEvts;

 public:
  ProtonRecoilAnalysis() {}
  ProtonRecoilAnalysis(double e_bin, TH2D *h_t_res_vs_cos_theta, unsigned int n_evts)
	  : EBin(e_bin), hTResVSCosTheta(h_t_res_vs_cos_theta), nEvts(n_evts) {}
  double GetEBin() const { return EBin; }
  void SetEBin(double e_bin) { EBin = e_bin; }
  TH2D *GetHtResVsCosTheta() const { return hTResVSCosTheta; }
  void SetHtResVsCosTheta(TH2D *h_t_res_vs_cos_theta) { hTResVSCosTheta = h_t_res_vs_cos_theta; }
  unsigned int GetNEvts() const { return nEvts; }
  void SetNEvts(unsigned int n_evts) { nEvts = n_evts; }

};

void TestMacro(int User_NEvts = INT_MIN, int User_iEvt = INT_MIN){

  loadlibs();

  SetBasicStyle();

  string sCalFile = "Theia_1kT_Cov30pct_wbls_5newpct_e+_Pos_0_0_0_Dir_0_0_1_EMatrix.root.root";
  MCCalib Cal(sCalFile);

  string sPDFFile = "Theia_1kT_Cov30pct_wbls_5newpct_e+_Pos_0_0_0_Dir_0_0_1_PDF.root";
  MCPDF PDF(sPDFFile);

  string sPathToFile = "../rat-pac/";
  string sGeom = "Theia_1kT_Cov30pct";
  string sMaterial = "wbls_5newpct";
  string sParticle = "neutron";
  string sE = "1_100MeV";
  string sInit = "Pos_0_0_0_Dir_0_0_1";

  //../MC/Theia_1kT_Cov30pct_wbls_1pct_e-_Pos_0_0_0_Dir_1_0_0
  string sRootPath = sPathToFile + sGeom + "_" + sMaterial + "_" + sParticle + "_" + sInit + "/";
  sRootPath = sPathToFile;
  string sFilename = sRootPath + sGeom+ "_" + sMaterial + "_" + sParticle + "_" + sE + "_" + sInit + ".root";

  string sOutname = sGeom+ "_" + sMaterial + "_" + sParticle + "_" + sE + "_" + sInit + "_plots.root";

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE OUTPUTS                     #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  auto *foutput = new TFile(sOutname.c_str(),"RECREATE");

  auto *FileAnalyzer = new Analyzer(sFilename.c_str());

  unsigned long NEvts = User_NEvts > INT_MIN ? User_NEvts : FileAnalyzer->GetNEvts();
  unsigned long iEvt = SetDefValue(User_iEvt, 0);
  NEvts = User_iEvt > 0 ? User_iEvt + NEvts : NEvts;
  ProgressBar progressBar(NEvts, 70);

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                CREATE HISTOGRAMS                  #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  vector<double> vCTBins(12); // slices of 15deg
  generate(vCTBins.begin(), vCTBins.end(), GenCTBin);
  sort(vCTBins.begin(), vCTBins.end());

  vector<double> EBins(100);
  generate(EBins.begin(), EBins.end(), GenEBin);
  vector<ProtonRecoilAnalysis> vPRA;

  const unsigned nBinsTRes = 35;
  const double minTRes = -5;
  const double maxTRes = 30;

  for(auto iE=0;iE<EBins.size();iE++){
	vPRA.emplace_back(ProtonRecoilAnalysis(EBins[iE],
										   new TH2D(Form("h2DEBin%.1f", EBins[iE]),
													"Cos#theta VS TResids proton hits",
													nBinsTRes, minTRes, maxTRes,
													vCTBins.size()-1,&vCTBins[0]),
										   0));
	vPRA[iE].GetHtResVsCosTheta()->Sumw2();
  }


  auto hTResidVSCTheta = new TH2D("hTResidVSCTheta",
								  "Cos#theta VS TResids proton hits ; T res (ns) ; Cos#theta",
								  nBinsTRes, minTRes, maxTRes,
								  vCTBins.size()-1,&vCTBins[0]);

  auto hProtonELike = new TH1D("hProtonELike", "Positron E-like from proton recoil ; e^{+} E-like (MeV)",
							   101,-0.05, 10.05);

  auto hEquivEMatrix = new TH2D("hEquivEMatrix", "E_{kin} equiv ; E proton recoil (MeV) ; E positron (MeV)",
								100,0.5,100.5,
								101,-0.05,10.05);

  auto hEVSNHits = new TH2D("hEVSNHits", "E VS NHits; E proton recoil (MeV) ; NHits ",
							100,0.5,100.5,
							101,-0.5,100.5);

  auto hEVSQ = new TH2D("hEVSQ", "E VS Q ; E proton recoil (MeV) ; Q",
						100,0.5,100.5,
						101,-0.5,100.5);


  auto hNHitsVSQ = new TH2D("hNHitsVSQ", "NHits VS Q from proton recoil ; Q ; NHits ",
							100,0.5,100.5,
							100,0.5,100.5);

  auto hCTProtonVSNeutron = new TH1D("hCTProtonVSNeutron",
									 "Cos#theta angle between neutron and proton daughters",
									 vCTBins.size()-1,&vCTBins[0]);

  auto hProtonEVSNeutronE = new TH2D("hProtonEVSNeutronE",
									 "proton recoil VS E neutrons ; E neutrons (MeV) ; E proton recoil (MeV)",
									 100, 0.5, 100.5,
									 100, 0.5, 100.5);

  auto hProtonMult = new TH2D("hProtonMult", "Proton Multiplicity VS E neutron ; E neutrons (MeV) ; proton multiplicity",
							  100,0.5,100.5,
							  100,0.5,1000.5);

  auto hEquivEVSNHits = new TH2D("hEquivEVSNHits", "e^{+} E-like VS NHits ; e^{+} E-like (MeV) ; NHits",
								 101,-0.05,10.05,
								 101,-0.5,100.5);

  auto hEquivEVSQ = new TH2D("hEquivEVSQ", "e^{+} E-like VS Q ; e^{+} E-like (MeV) ; Q",
							 101,-0.05,10.05,
							 101,-0.5,100.5);

  auto hProb = new TH1D("hProb", "GoF", 1000, 0., 1000.);

  vector<TGraph*> vGrL;

  for(iEvt; iEvt<NEvts; iEvt++) {

	// record the tick
	++progressBar;

	FlatParticle PrimPart = GetPrimaryParticleInfo(FileAnalyzer, iEvt);
	double PrimPartKE = PrimPart.GetKinE();
	const TVector3& TrueOrigin = PrimPart.GetPos();
	const TVector3& TrueDir = PrimPart.GetDir();
	// Proton recoil E
	vector<ComplexParticle> vProton;
	vector<Hit> vHit = GetVHitsFromPart(FileAnalyzer, iEvt, "proton", &vProton);

	// double dEdX = 0.;
	// for(auto &P:vProton){
	//   FilldEdX(P, FileAnalyzer, iEvt);
	// }

	const double nProton = vProton.size();
	hProtonMult->Fill(PrimPartKE,
					  nProton);

	double NHits, Q;
	GetNPEAndNHitsFromHits(vHit, &Q, &NHits);

	hNHitsVSQ->Fill(Q, NHits);

	auto *grL = new TGraph();
	double E = ComputeECalib(Cal, Q, NHits, grL);
	if(E>0.1){
	  grL->SetName(Form("grEvt%lu", iEvt));
	  vGrL.emplace_back(grL);
	}

	hEquivEVSNHits->Fill(E, NHits);
	hEquivEVSQ->Fill(E, Q);

	FillTResVSCosTheta(hTResidVSCTheta, vHit,
					   TVector3(0.,0.,0.), 224.9,
					   TVector3(0.,0.,1.),
					   true);

	double roundedE = round(E*10)/10;

	for(auto &p:vPRA){
	  if(round(p.GetEBin()*10) == round(roundedE*10)){
		FillTResVSCosTheta(p.GetHtResVsCosTheta(), vHit,
						   TVector3(0.,0.,0.), 224.9,
						   TVector3(0.,0.,1.),
						   true);
		p.SetNEvts(p.GetNEvts()+1);
	  }
	}


	hProtonELike->Fill(E);

	for(auto p:vProton){
	  hEquivEMatrix->Fill(p.GetKinE(), E);
	  hEVSNHits->Fill(p.GetKinE(), NHits);
	  hEVSQ->Fill(p.GetKinE(), Q);
	}



	auto *hTResVSCosTheta_Evt = new TH2D(Form("hTResVSCosThetaEvt%lu", iEvt), "hTRes VS Cos #theta",
										 nBinsTRes,minTRes,maxTRes,
										 vCTBins.size()-1,&vCTBins[0]);

	FillTResVSCosTheta(hTResVSCosTheta_Evt, vHit,
					   TVector3(0.,0.,0.), 224.9,
					   TVector3(0.,0.,1.),
					   true);

	FillProtonRecoilSpectrum(hProtonEVSNeutronE, FileAnalyzer, iEvt);

	double Chi2 = -1;
	int NdF = -1;
	double Prob = CalculateProb(PDF.GetHPDF(E), hTResVSCosTheta_Evt, &Chi2, &NdF);

	hProb->Fill(Chi2/(double)(NdF));

	delete hTResVSCosTheta_Evt;

	// display the bar
	progressBar.display();

  }

  auto *grROC = new TGraph();
  grROC->SetName(Form("gr%s", sParticle.c_str()));
  auto nBinsX = hProb->GetNbinsX();
  for(auto iBinX=1; iBinX<=nBinsX; iBinX++) {
	grROC->SetPoint(grROC->GetN(),
		hProb->GetBinCenter(hProb->GetBin(iBinX)),
		hProb->Integral(iBinX, nBinsX)/hProb->Integral());
  }

  PlotAHist(hNHitsVSQ, "COLZ");
  PlotAHist(hTResidVSCTheta, "COLZ");
  PlotAHist(hProtonELike);
  PlotAHist(hEquivEMatrix, "COLZ");
  PlotAHist(hEVSNHits, "COLZ");
  PlotAHist(hEVSQ, "COLZ");

  PlotAHist(hProb);
  PlotAHist(hProtonEVSNeutronE, "COLZ");
  PlotAHist(hProtonMult, "COLZ");
  // PlotAHist(grROC, "APC");

  PlotAHist(hEquivEVSNHits, "COLZ");
  PlotAHist(hEquivEVSQ, "COLZ");

  foutput->cd();
  hNHitsVSQ->Write();
  hTResidVSCTheta->Write();
  hProtonELike->Write();
  hEquivEMatrix->Write();
  hEVSNHits->Write();
  hEVSQ->Write();
  hProb->Write();
  hProtonEVSNeutronE->Write();
  hProtonMult->Write();

  hEquivEVSNHits->Write();
  hEquivEVSQ->Write();

  grROC->Write();

  for(auto p:vPRA){
	if(p.GetNEvts()>0){
	  p.GetHtResVsCosTheta()->Scale(1/(double)(p.GetNEvts()));
	  p.GetHtResVsCosTheta()->Write();
	}
  }

  for(auto gr:vGrL){
    gr->Write();
  }

  foutput->Close();

}
