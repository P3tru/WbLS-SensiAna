/////////////////////////   ROOT   //////////////////////////
#include <TVector3.h>

/////////////////////////   USER  ///////////////////////////
#include <AnalyzerFunctions.hh>

void GetHitTime(Analyzer *Ana, unsigned int iEvt, TH1D* hHitTime){

  Ana->GetTreeMc()->GetEntry(iEvt);
  auto *mc = Ana->GetRds()->GetMC();

  auto NbPMTsHits = mc->GetMCPMTCount();

  for (int iPMT = 0; iPMT < NbPMTsHits; iPMT++) {

	auto mcPMT = mc->GetMCPMT(iPMT);
	auto PMTID = mcPMT->GetID();
	auto PMTPos = Ana->GetRun()->GetPMTInfo()->GetPosition(PMTID);

	auto NbPhotonCounts = mcPMT->GetMCPhotonCount();

	for (int iP = 0; iP < NbPhotonCounts; iP++) {

	  auto mcPhoton = mcPMT->GetMCPhoton(iP);
	  // Get Q and T
	  auto Q = mcPhoton->GetCharge();
	  auto T = mcPhoton->GetHitTime();

	  // ########## //
	  // FILL EVENT //
	  // ########## //

	  if(hHitTime)
		hHitTime->Fill(T);

	}

  } // END FOR iPMT

}
