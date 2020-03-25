//
// Created by zsoldos on 2/21/20.
//

///////////////////////// STL C/C++ /////////////////////////

/////////////////////////   RAT   ///////////////////////////

/////////////////////////   ROOT  ///////////////////////////
#include <TGraph.h>

/////////////////////////   USER  ///////////////////////////
#include <MCFunctions.hh>

using namespace std;

RAT::DS::MC * GetRATMCOnEvt(Analyzer *fAnalyzer, unsigned int iEvt){

  fAnalyzer->GetTreeMc()->GetEntry(iEvt);

  // Access RAT MC info and the summary
  // Summary useful to get nCer photons, nScint photons, etc...
  return fAnalyzer->GetRds()->GetMC();

}

vector<Hit> GetMCHitCollection(Analyzer *fAnalyzer, unsigned int iEvt){


  RAT::DS::MC * mc = GetRATMCOnEvt(fAnalyzer, iEvt);

  vector<Hit> vHit;

  for (int iPMT = 0; iPMT < mc->GetMCPMTCount(); iPMT++) {

	auto mcPMT = mc->GetMCPMT(iPMT);
	auto PMTID = mcPMT->GetID();
	auto PMTPos = fAnalyzer->GetRun()->GetPMTInfo()->GetPosition(PMTID);

	auto NbPhotonCounts = mcPMT->GetMCPhotonCount();

	for (int iP = 0; iP < NbPhotonCounts; iP++) {

	  auto mcPhoton = mcPMT->GetMCPhoton(iP);
	  // Get Q and T
	  auto Q = mcPhoton->GetCharge();
	  auto T = mcPhoton->GetHitTime();

	  // ########## //
	  // FILL EVENT //
	  // ########## //

	  vHit.emplace_back(Hit(PMTPos, Q, T));

	}

  } // END FOR iPMT

  return vHit;

}

static bool IsParticle(RAT::DS::MCTrack *mctrack, const char *name){

  if(mctrack->GetParticleName() == name){
	return true;
  } else {
	return false;
  }

}

double GetBirks(RAT::DS::MC *mc){

  auto * mcs = mc->GetMCSummary();

  return mcs->GetTotalScintEdepQuenched()/mcs->GetNumScintPhoton();

}

void PlotStuff(Analyzer *fAnalyzer){

  TH1D *hPlot = new TH1D("hPlot","",1001,-0.5,1000.5);

  for(int iEvt=0; iEvt<fAnalyzer->GetNEvts(); iEvt++){

	hPlot->Fill(GetBirks(GetRATMCOnEvt(fAnalyzer, iEvt)));

  }

  hPlot->Draw();

}

// Loop for each proton
// Get all optical photon child


static double MeV2lambda(double MeV){

  const double hc = 1.23984193e-3; //MeV.nm
  return hc/MeV;

}


vector<TProton> GetProtonLY(Analyzer *fAnalyzer, unsigned int iEvt){

  auto * mc = GetRATMCOnEvt(fAnalyzer, iEvt);
  auto NbTracks = mc->GetMCTrackCount();

  vector<TProton> protons;

  for(auto iTrack=0; iTrack<NbTracks; iTrack++){

	auto * mctrack = mc->GetMCTrack(iTrack);

	if(IsParticle(mctrack, "proton")){

	  TProton proton;
	  proton.SetId(mctrack->GetID());

	  auto *mctrackstep = mctrack->GetMCTrackStep(0);
	  auto dE = mctrackstep->GetKE();
	  auto dX = mctrackstep->GetLength();


	  auto NbTrackSteps = mctrack->GetMCTrackStepCount();

	  for(auto iStep=1; iStep<NbTrackSteps; iStep++){

		auto *mctrackstep = mctrack->GetMCTrackStep(iStep);

		dE -= mctrackstep->GetKE();
		dX += abs(mctrackstep->GetLength());

	  }

	  proton.SetDe(dE);
	  proton.SetDx(dX);
	  protons.emplace_back(proton);

	}

  }

  for(auto iTrack=0; iTrack<NbTracks; iTrack++) {

	auto * mctrack = mc->GetMCTrack(iTrack);

	if(IsParticle(mctrack, "opticalphoton")) {

	  for(auto &p:protons){

	    if (mctrack->GetParentID() == p.GetId()){

		  p.AddPhoton(TPhoton(mctrack->GetID(),
							  mctrack->GetParentID(),
							  MeV2lambda(mctrack->GetMCTrackStep(0)->GetKE())));

		  break;

	    }


	  }


	}

  }

  return protons;

}

void PrintTrackInfo(RAT::DS::MC *mc, unsigned int iTrack){

  auto * mctrack = mc->GetMCTrack(iTrack);
  auto NbTrackSteps = mctrack->GetMCTrackStepCount();

  cout << mctrack->GetParticleName() << endl;
  cout << " With ID: " << mctrack->GetID() << " and parent ID: " << mctrack->GetParentID() << endl;
  cout << "TOT Nb Steps: " << NbTrackSteps << endl;

  for(auto iStep=0; iStep<NbTrackSteps; iStep++){

    auto *mctrackstep = mctrack->GetMCTrackStep(iStep);

    cout << mctrackstep->GetLength() <<"mm " << mctrackstep->GetKE() << "MeV " << endl;

  }


}

void FillBirksLaw(Analyzer *fAnalyzer, unsigned int iEvt, TH2D *h2D){

  vector<TProton> protons = GetProtonLY(fAnalyzer, iEvt);

  for(auto p:protons){

    double dL = p.GetPhotonChildren().size();
    double dX = abs(p.GetDx());
    double dE = abs(p.GetDe());

	if(h2D)
	  h2D->Fill(dL/dX, (dE/dX)/(1+dE/dX));

  }

}

void PlotBirksLaw(Analyzer *fAnalyzer){

  TH2D *hBirks = new TH2D("hBirks", "Birks Law according to rat-pac",
						  100,-0.05,10.05,
						  100,0.55,1.05);

  hBirks->SetTitle("10MeV protons; dLY/dX; (dE/dX)/(1+dE/dX))");

  for(unsigned int iEvt=0; iEvt<fAnalyzer->GetNEvts(); iEvt++){

    FillBirksLaw(fAnalyzer, iEvt, hBirks);

  }

  hBirks->Draw("COLZ");

}

void GetMCNumPEAndHits(double *NHits, double *NPE, Analyzer *fAnalyzer, unsigned int iEvt){


  RAT::DS::MC * mc = GetRATMCOnEvt(fAnalyzer, iEvt);

  *NHits=mc->GetMCPMTCount();
  *NPE=mc->GetNumPE();

}

