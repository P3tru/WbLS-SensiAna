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

 protected:

  vector<double> EBins;

  vector< pair<double, TH2D*> > vpPDF;

  double minEBin;
  double maxEBin;

 public:

  MCPDF() = default;

  MCPDF(const string& sFileName) {

    minEBin=0;
    maxEBin=0;

    auto *fPDF = TFile::Open(sFileName.c_str());

	EBins.resize(100);
	int i = 0;
	generate(EBins.begin(), EBins.end(), [&](){return (double)(0.1*i++);});

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
			h2D->SetDirectory(0);
			if(h2D->GetEntries()>0)
			  vpPDF.emplace_back(make_pair(EBin, h2D));

			break;

		  } // END if strcmp

		} // END for EBin

	  } // END while Key

	}

    fPDF->Close();

    sort(vpPDF.begin(), vpPDF.end());

    minEBin = vpPDF.begin()->first;
    auto lastBin = vpPDF.size()-1;
    maxEBin = vpPDF[lastBin].first;

  }

  const vector<double> &GetEBins() const { return EBins; }
  const vector< pair<double, TH2D*> > GetVpPdf() const{ return vpPDF; }

  double GetMinEBin() const { return minEBin; }
  double GetMaxEBin() const { return maxEBin; }

  TH2D *GetHPDF(double E){

    double roundedE = round(E*10)/10;

    // if(roundedE<minEBin || roundedE>=maxEBin)
    //   return nullptr;

    if(roundedE <= minEBin)
	  return vpPDF.begin()->second;
	if( roundedE >= maxEBin)
	  return vpPDF[vpPDF.size()-1].second;

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

  SetBasicStyle();

  string sCalFile = "Theia_1kT_Cov30pct_wbls_5newpct_e+_Pos_0_0_0_Dir_0_0_1_EMatrix.root.root";
  MCCalib Cal(sCalFile);

  string sPDFFile = "Theia_1kT_Cov30pct_wbls_5newpct_e+_Pos_0_0_0_Dir_0_0_1_PDF.root";
  MCPDF PDF(sPDFFile);

  string sPDFBCKGFile = "Theia_1kT_Cov30pct_wbls_5newpct_proton_10.0_100.0MeV_Pos_0_0_0_Dir_0_0_1_PDF.root";
  MCPDF PDFBCKG(sPDFBCKGFile);

  string sPathToFile = "../rat-pac/";
  string sGeom = "Theia_1kT_Cov30pct";
  string sMaterial = "wbls_5newpct";
  string sParticle = "e-";
  string sE = "1_10MeV";
  string sInit = "Pos_0_0_0_Dir_0_0_1";

  //../MC/Theia_1kT_Cov30pct_wbls_1pct_e-_Pos_0_0_0_Dir_1_0_0
  string sRootPath = sPathToFile + sGeom + "_" + sMaterial + "_" + sParticle + "_" + sInit + "/";
  sRootPath = sPathToFile;
  string sFilename = sRootPath + sGeom+ "_" + sMaterial + "_" + sParticle + "_" + sE + "_" + sInit + ".root";

  string sOutname = sGeom+ "_" + sMaterial + "_" + sParticle + "_" + sE + "_" + sInit + "_PMatrix.root";

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
										   new TH2D(Form("hTResVSCosThetaE%.1fMeV", EBins[iE]),
													"Cos#theta VS TResids proton hits",
													nBinsTRes, minTRes, maxTRes,
													vCTBins.size()-1,&vCTBins[0]),
										   0));
	vPRA[iE].GetHtResVsCosTheta()->Sumw2();
  }


  auto hTResidVSCTheta = new TH2D("hTResidVSCTheta", "Cos#theta VS TResids proton hits",
								  nBinsTRes, minTRes, maxTRes,
								  vCTBins.size()-1,&vCTBins[0]);

  auto hProtonELike = new TH1D("hProtonELike", "Positron E-like from proton recoil",
							   101,-0.05, 10.05);

  auto hEquivEMatrix = new TH2D("hEquivEMatrix", "E_{kin} equiv ; E proton recoil (MeV) ; E positron (MeV)",
								100,0.5,100.5,
								101,-0.05,10.05);

  auto hEVSNHits = new TH2D("hEVSNHits", "E VS NHits; E proton recoil (MeV) ; NHits",
							100,0.5,100.5,
							101,-0.5,100.5);

  auto hEVSQ = new TH2D("hEVSQ", "E VS Q ; E proton recoil (MeV) ; Q",
						100,0.5,100.5,
						101,-0.5,100.5);

  auto hNHitsVSQ = new TH2D("hNHitsVSQ", "NHits VS Q from proton recoil ; Q ; NHits ",
							100,0.5,100.5,
							100,0.5,100.5);

  auto nEvtsProcessed = 0;

  vector<TGraph*> vGrL;

  vector< pair <double, TH2D*> > vhNHitsVSQ(10);
  generate(vhNHitsVSQ.begin(), vhNHitsVSQ.end(),
		   [](){
			 static int i=0;
			 return make_pair((double)(10*i), new TH2D(Form("hNHitsVSQE%d", 10*i++),
													   "NHits VS Q from proton recoil ; Q ; NHits ",
													   100,0.5,100.5,
													   100,0.5,100.5));
		   });
  sort(vhNHitsVSQ.begin(), vhNHitsVSQ.end(),
  	[](const pair<double, TH2D*> &p1, const pair<double, TH2D*> &p2){
    return ( p1.first > p2.first ||
		( p2.first <= p1.first && p1.second < p2.second ) );
	});

  auto hGOFSIG = new TH1D("hGOFSIG", "GoF SIG",
						  100, 0.5, 1000.5);
  auto hGOFBCKG = new TH1D("hGOFBCKG", "GoF BCKG",
						   100, 0.5, 1000.5);

  auto hProbSIG = new TH1D("hProbSIG", "Prob SIG",
						   100, -0.0005, 0.9995);
  auto hProbBCKG = new TH1D("hProbBCKG", "GoF BCKG",
							100, -0.0005, 0.9995);

  for(iEvt; iEvt<NEvts; iEvt++) {

	// record the tick
	++progressBar;

	FlatParticle PrimPart = GetPrimaryParticleInfo(FileAnalyzer, iEvt);
	double PrimPartKE = PrimPart.GetKinE();
	const TVector3& TrueOrigin = PrimPart.GetPos();
	const TVector3& TrueDir = PrimPart.GetDir();
	double NHits, Q;
	GetNPEAndNHitsFromEV(FileAnalyzer, iEvt, 0, &NHits, &Q);
	if(NHits == 0 && Q == 0)
	  continue;
	if(NHits>0 && Q>0)
	  nEvtsProcessed++;
	// Neutron E-Like
	auto *grL = new TGraph();
	double E = ComputeECalib(Cal, Q, NHits, grL);
	if(E>0.1){
	  grL->SetName(Form("grEvt%lu", iEvt));
	  vGrL.emplace_back(grL);
	}

	for(auto &p:vhNHitsVSQ){
	  if(PrimPartKE>p.first){
	    p.second->Fill(Q, NHits);
		break;
	  }
	}

	hNHitsVSQ->Fill(Q, NHits);

	vector<Hit> vHit = GetEVHitCollection(FileAnalyzer, iEvt);

	FillTResVSCosTheta(hTResidVSCTheta, vHit, TrueOrigin, 224.9, TrueDir);

	double roundedE = round(E*10)/10;

	for(auto &p:vPRA){
	  if(round(p.GetEBin()*10) == round(roundedE*10)){
	    FillTResVSCosTheta(p.GetHtResVsCosTheta(), vHit, TrueOrigin, 224.9, TrueDir);
		p.SetNEvts(p.GetNEvts()+1);
	  }
	}


	auto *hTResVSCosTheta_Evt = new TH2D(Form("hTResVSCosThetaEvt%lu", iEvt), "hTRes VS Cos #theta",
										 nBinsTRes,minTRes,maxTRes,
										 vCTBins.size()-1,&vCTBins[0]);

	FillTResVSCosTheta(hTResVSCosTheta_Evt, vHit,
					   TVector3(0.,0.,0.), 224.9,
					   TVector3(0.,0.,1.));

	if(E >= 1.0 && E<=3.0){

	  double Chi2 = -1;
	  int NdF = -1;
	  double Prob = CalculateProb(PDF.GetHPDF(E), hTResVSCosTheta_Evt, &Chi2, &NdF);
	  hGOFSIG->Fill(Chi2/(double)(NdF));
	  hProbSIG->Fill(Prob);

	  Chi2 = -1;
	  NdF = -1;
	  // cout << "iEvt:" << iEvt << endl;
	  Prob = CalculateProb(PDFBCKG.GetHPDF(E), hTResVSCosTheta_Evt, &Chi2, &NdF);
	  hGOFBCKG->Fill(Prob);
	  hProbBCKG->Fill(Chi2/(double)(NdF));

	}

	hProtonELike->Fill(E);

	hEquivEMatrix->Fill(PrimPartKE, E);

	hEVSNHits->Fill(PrimPartKE, NHits);

	hEVSQ->Fill(PrimPartKE, Q);

	progressBar.display();

  }

  PlotAHist(hTResidVSCTheta, "COLZ");
  PlotAHist(hProtonELike);
  PlotAHist(hEquivEMatrix, "COLZ");
  PlotAHist(hEVSNHits, "COLZ");
  PlotAHist(hEVSQ, "COLZ");
  PlotAHist(hNHitsVSQ, "COLZ");
  PlotAHist(hProbSIG);
  PlotAHist(hProbBCKG);
  PlotAHist(hGOFSIG);
  PlotAHist(hGOFBCKG);

  foutput->cd();

  hTResidVSCTheta->Write();
  hProtonELike->Write();
  hEquivEMatrix->Write();
  hEVSNHits->Write();
  hEVSQ->Write();
  hNHitsVSQ->Write();

  hProbSIG->Write();
  hProbBCKG->Write();
  hGOFSIG->Write();
  hGOFBCKG->Write();

  for(auto &p:vhNHitsVSQ) {
    p.second->Write();
  }
  for(auto p:vPRA) {
	if (p.GetNEvts() > 0) {
	  p.GetHtResVsCosTheta()->Scale(1 / (double) (p.GetNEvts()));
	  p.GetHtResVsCosTheta()->Write();
	}
  }
  for(auto gr:vGrL){
	gr->Write();
  }
  foutput->Close();

}
