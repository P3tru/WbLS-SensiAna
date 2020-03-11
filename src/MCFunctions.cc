#include <vector>

/////////////////////////   USER  ///////////////////////////
#include <MCFunctions.hh>

using namespace std;

RAT::DS::MC * GetRATMCOnEvt(Analyzer *fAnalyzer, unsigned int iEvt){

  fAnalyzer->GetTreeMc()->GetEntry(iEvt);

  // Access RAT MC info and the summary
  // Summary useful to get nCer photons, nScint photons, etc...
  return fAnalyzer->GetRds()->GetMC();

}


vector<Hit> GetHitCollection(Analyzer *fAnalyzer, unsigned int iEvt){


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

vector<Hit> SplitHitCollection(vector<Hit> *vHit, double PromptWindow){

  vector<Hit> Empty;

  unsigned int SizeInit = vHit->size();

  if(vHit->size() == 0){

	return Empty;

  } else if ( vHit->size() > 0 ){

	sort(vHit->begin(), vHit->end());
	for(unsigned int iHit=0; iHit<vHit->size()-1; iHit++){

	  const double dT = vHit->at(iHit+1).GetT() - vHit->at(iHit).GetT();

	  if(dT>PromptWindow){

		if(vHit->begin() + iHit + 1 < vHit->end()){

		  try {

			vector<Hit> vHitDelayed(vHit->begin() + iHit + 1, vHit->end());

			vHit->erase(vHit->begin() + iHit + 1, vHit->end());

			return vHitDelayed;

		  } catch (vector<Hit> const *vHit){

			for(auto h : *vHit)
			  h.Print();

			return Empty;

		  }


		} else if (vHit->begin() + iHit + 1 >= vHit->end()){

		  return Empty;

		}
	  }

	}

  }

  if(vHit->size() == SizeInit){
	return Empty;
  }

}
