//
// Created by zsoldos on 11/26/19.
//
#include "../ProgressBar.hpp"

#define PI 3.14159
#define C 299.792458 // mm/ns

class Hit {

 protected:
  double X;
  double Y;
  double Z;

  double Q;
  double T;

 public:
  Hit() : X(-1), Y(-1), Z(-1), Q(-1), T(-1){};
  Hit(double x, double y, double z, double q, double t) : X(x), Y(y), Z(z), Q(q), T(t){};
  ~Hit(){};

  double GetX(){ return X; };
  void SetX(double x){ X = x; };
  double GetY(){ return Y; };
  void SetY(double x){ Y = x; };
  double GetZ(){ return Z; };
  void SetZ(double x){ Z = x; };

  double GetQ(){ return Q; };
  void SetQ(double x){ Q = x; };
  double GetT(){ return T; };
  void SetT(double x){ T = x; };

  double GetDistance(){return sqrt(X*X + Y*Y + Z*Z);};

};

class VectorizedHit {

 protected:

  std::vector<double> T;

 public:
  VectorizedHit() : {};
  ~VectorizedHit(){};

};

double EvalL(double Nobs, double Npred){
  double L;
  if (Nobs>0 && Npred>0)L=Npred-Nobs+Nobs*TMath::Log(Nobs/Npred);
  else L=Npred;
  L=-L;
  return L;
}

double GetRMSFromShiftedMean(TH1 *h, double ShiftedMean=0){
  const unsigned int nbBins = h->GetNbinsX();

  double D2 = 0;

  for(int iBin=1; iBin<nbBins+1; iBin++){

    D2 += std::pow(h->GetBinContent(h->GetBin(iBin)) - ShiftedMean,2);

  }

  return sqrt(D2);
}

double GetTGuess(double *x, double *par){

  std::cout << sizeof(par) << std::endl;
  std::cout << sizeof(par[0]) << std::endl;

  size_t n = sizeof(*par)/sizeof(par);

  std::cout << n << std::endl;

  double T_GUESS = 0;

  for(unsigned int i=0; i<n; i++){

    std::cout << par[i] << std::endl;

	T_GUESS += (par[i] - x[0]/C) / n;

  }

  return T_GUESS;

}

void PlotIBD(const char *file=NULL){

  gROOT->SetStyle("Plain");

  gStyle->SetOptStat(1); // This determines if you want a stats box
  gStyle->SetOptFit(1); // This determines if you want a fit info box
  gStyle->GetAttDate()->SetTextColor(1);
  gStyle->SetOptTitle(1); // no title; comment out if you want a title
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetTextFont(132);
  gStyle->SetTitleFont(132,"XYZ");

  gROOT->ForceStyle();

  gStyle->SetPalette(kDarkRainBow);

  auto *f = new TFile(file);

  std::vector< std::vector<Hit> > Positrons;

  TH2D* hTVSQ = new TH2D("hTVSQ","hTVSQ",
						 201, -0.5, 200.5,
						 7, -0.5, 6.5);

  const int PromptCut = 15; //ns

  const int NbTResidBins = 101; const double minTResid = -0.5; const double maxTResid = 100.5;
  auto *hTResiduals = new TH1D("hTResiduals", "T Residuals",
							   NbTResidBins, minTResid, maxTResid);
  const int NbDBins = 101; const double minD = 5000.5; const double maxD = 8000.5;
  auto *hDistances = new TH1D("hDistances", "photon distances",
							  NbDBins, minD, maxD);
  hTResiduals->Sumw2();

  int nEvtWithHits = 0;

  std::vector<TH1D*> hTResPerEvt;
  std::vector<TH1D*> hDPerEvt;

  if(f->IsOpen()) {

	////////////////////////////////
	// IF file is open do stuff ////
	////////////////////////////////

	auto *run = new RAT::DS::Run();

	auto *runTree = (TTree *) f->Get("runT");
	runTree->SetBranchAddress("run", &run);
	runTree->GetEntry(0);

	auto *tree = (TTree *) f->Get("T");

	auto *rds = new RAT::DS::Root();
	tree->SetBranchAddress("ds", &rds);

	auto nEvents = tree->GetEntries();
//	nEvents = 10;
	std::cout << "#Evts: " << nEvents << std::endl;
	ProgressBar progressBar(nEvents, 70);

	std::cout << "CREATE PDFs" << std::endl;

	for (int iEvt = 0; iEvt < nEvents; iEvt++) {
	  tree->GetEntry(iEvt);
	  // record the tick
	  ++progressBar;

	  std::vector<Hit> Positron;

	  TH1D *hTResPerEvtPerHit = new TH1D(Form("hTResidualsPerEvt#%d", iEvt),
										 "T Residuals",
										 NbTResidBins, minTResid, maxTResid);

	  TH1D *hDPerEvtPerHit = new TH1D(Form("hDistances%d", iEvt),
									  "photon distances",
									  NbDBins, minD, maxD);


	  // Access RAT MC info and the summary
	  // Summary useful to get nCer photons, nScint photons, etc...
	  auto *mc = rds->GetMC();

	  auto NbPMTsHits = mc->GetMCPMTCount();

	  auto NbParentCounts = mc->GetMCParentCount();

	  auto IsEvtWHits = false;

	  for (int iPMT = 0; iPMT < NbPMTsHits; iPMT++) {

		auto mcPMT = mc->GetMCPMT(iPMT);
		auto PMTID = mcPMT->GetID();
		auto xHit = run->GetPMTInfo()->GetPosition(PMTID).x();
		auto yHit = run->GetPMTInfo()->GetPosition(PMTID).y();
		auto zHit = run->GetPMTInfo()->GetPosition(PMTID).z();
		Hit h;
		h.SetX(xHit);
		h.SetY(yHit);
		h.SetZ(zHit);

		auto NbPhotonCounts = mcPMT->GetMCPhotonCount();

		for (int iP = 0; iP < NbPhotonCounts; iP++) {

		  auto mcPhoton = mcPMT->GetMCPhoton(iP);
		  // Get Q and T
		  auto Q = mcPhoton->GetCharge();
		  auto T = mcPhoton->GetHitTime();
		  h.SetQ(Q);
		  h.SetT(T);

		  if (T > 0 && T < 0.5e3) {

			IsEvtWHits = true;
			auto TResid = T - (h.GetDistance() / C);
			auto w = 1 / (double) (NbPMTsHits);
			hTResiduals->Fill(TResid, w);
			hTResPerEvtPerHit->Fill(TResid, w);

			hDistances->Fill(h.GetDistance());
			hDPerEvtPerHit->Fill(h.GetDistance());

			Positron.emplace_back(h);
		  }

		  hTVSQ->Fill(T, Q);

		} // END FOR iPhoton

	  } // END FOR iPMT

	  Positrons.emplace_back(Positron);

	  hTResPerEvt.emplace_back(hTResPerEvtPerHit);
	  hDPerEvt.emplace_back(hDPerEvtPerHit);

	  if (IsEvtWHits)
		nEvtWithHits++;

	  // display the bar
	  progressBar.display();

	} // END for iEvt


	hTResiduals->Scale(1/(double)(nEvtWithHits));

	// NOW FIT

	std::cout << std::endl;
	std::cout << "### ##### ###" << std::endl;
	std::cout << "###  FIT  ###" << std::endl;
	std::cout << "### ##### ###" << std::endl;

	const int NbDBinsMin = 101;
	const double minDMin = 6000;
	const double maxDMin = 7000;
	const double StepDMin = (double)(maxDMin - minDMin)/(double)(NbDBinsMin);

	TH1D* hPosRes = new TH1D("hPosRes","PosRes",
							 2*NbDBinsMin, -NbDBinsMin*StepDMin, NbDBinsMin*StepDMin);

	auto NbEvtsToFit = Positrons.size();
//	NbEvtsToFit = 1000;

	ProgressBar progressBarFIT(NbEvtsToFit, 70);


	for(unsigned int iEvt=0; iEvt<NbEvtsToFit; iEvt++){

//	  // record the tick
//	  ++progressBarFIT;

	  if(iEvt % 100 == 0){
		std::cout << "Evt#" << iEvt << std::endl;
	  }

	  const unsigned int NHits = Positrons[iEvt].size();

	  std::vector<double> NLL;
	  for(int iD = 0; iD<NbDBinsMin; iD++){

	    double D_GUESS = minDMin + iD*StepDMin;

	    TH1D* hTResPerEvt_GUESS = new TH1D(Form("hTResidualsPerEvt_GUESS_%d_%d",iEvt,iD),"",
										   NbTResidBins, minTResid, maxTResid);

		for(int iHit=0; iHit<NHits; iHit++){
		  hTResPerEvt_GUESS->Fill(Positrons[iEvt][iHit].GetT() - (D_GUESS / C),
								  1/(double)(NHits));
		}

//		std::cout << "D_GUESS: " << D_GUESS << std::endl;
		double NLL_val = 0;
		for(int iBin=0; iBin<NbDBinsMin; iBin++){
		  NLL_val += EvalL(hTResPerEvt_GUESS->GetBinContent(hTResPerEvt_GUESS->GetBin(iBin)),
								 hTResiduals->GetBinContent(hTResiduals->GetBin(iBin)));
//		  std::cout << "#Bin" << iBin
//					<< " TRes_OBS: " << hTResPerEvt_GUESS->GetBinContent(hTResPerEvt_GUESS->GetBin(iBin))
//					<< " TRes_PRED: " << hTResiduals->GetBinContent(hTResiduals->GetBin(iBin))
//					<< std::endl;
		}
//		std::cout << "NLL: " << NLL_val << std::endl;
		NLL.emplace_back(NLL_val);

		delete hTResPerEvt_GUESS;

	  }

	  int minElementIndex = std::min_element(NLL.begin(),NLL.end()) - NLL.begin();

	  hPosRes->Fill(minDMin + StepDMin*minElementIndex - hDistances->GetMean());

//	  std::cout << "FIT d: "
//				<< minDMin + StepDMin*minElementIndex << std::endl;
//	  std::cout << "FIT d-d_mean: "
//				<< minDMin + StepDMin*minElementIndex - hDistances->GetMean() << std::endl;

//	  // display the bar
//	  progressBarFIT.display();

	}

	//////////////////////////////
	// #### #### DRAW #### #### //
	//////////////////////////////

	auto *c1 = new TCanvas("c1","c1",800,600);
	c1->SetGrid();
	hTVSQ->Draw("COLZ");

	auto *c2 = new TCanvas("c2","c2",800,600);
	c2->SetGrid();
	hTResiduals->Draw("HIST E");

	auto *c3 = new TCanvas("c3","c3",800,600);
	c3->SetGrid();
	hDistances->Draw();

	auto *c4 = new TCanvas("c4","c4",800,600);
	c4->SetGrid();
	for(int iH=0; iH<10; iH++){
	  hTResPerEvt[iH]->Draw("SAME PLC PMC");
	}

	auto *c5 = new TCanvas("c5","c5",800,600);
	c5->SetGrid();
	hPosRes->Draw();
	hPosRes->Fit("gaus");

  } // END if f

} // END PlotIBD()
