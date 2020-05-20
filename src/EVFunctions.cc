//
// Created by zsoldos on 2/21/20.
//

/////////////////////////   USER  ///////////////////////////
#include "EVFunctions.hh"

RAT::DS::EV * GetRATEVOnEvt(Analyzer *fAnalyzer, unsigned int iEvt, unsigned int iEV){

  fAnalyzer->GetTree()->GetEntry(iEvt);

  // Access RAT MC info and the summary
  // Summary useful to get nCer photons, nScint photons, etc...
  return fAnalyzer->GetDS()->GetEV(iEV);

}

vector<Hit> GetEVHitCollection(Analyzer *fAnalyzer, unsigned int iEvt, unsigned int iEV){

  auto *EV = GetRATEVOnEvt(fAnalyzer, iEvt, iEV);

  vector<Hit> vHit;

  for (int iPMT = 0; iPMT < EV->GetPMTCount(); iPMT++) {

	auto PMT = EV->GetPMT(iPMT);
	auto PMTID = PMT->GetID();
	auto PMTPos = fAnalyzer->GetRun()->GetPMTInfo()->GetPosition(PMTID);

	// Get Q and T
	auto Q = PMT->GetCharge();
	auto T = PMT->GetTime();

	// ########## //
	// FILL EVENT //
	// ########## //

	vHit.emplace_back(Hit(PMTPos, Q, T));


  } // END FOR iPMT

  return vHit;

}

void GetNPEAndNHitsFromEV(Analyzer *fAnalyzer, unsigned int iEvt, unsigned int iEV,
						  double *NHits, double *Q){

  auto *EV = GetRATEVOnEvt(fAnalyzer, iEvt, iEV);

  *Q = EV->GetTotalCharge();
  *NHits = EV->Nhits();

}