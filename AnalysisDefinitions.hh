//
// Created by zsoldos on 11/26/19.
//

#ifndef _ANALYSISDEFINITIONS_HH_
#define _ANALYSISDEFINITIONS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <fstream>
#include <sstream>
#include <regex>
#include <numeric>
#include <vector>
#include <array>
#include <iostream>
#include <functional>
#include <iterator>

/////////////////////////   BOOST   /////////////////////////
#include <boost/filesystem/path.hpp>

/////////////////////////   ROOT   //////////////////////////
#include <TH2D.h>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/Root.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCPhoton.hh>

using namespace std;

void CollectPEAndHits(RAT::DS::MC *mc, TH2D *Hist){

  // Get Nb of PMTs which at least 1hit
  const int nbPE = mc->GetNumPE();
  const int nbPMTsHits = mc->GetMCPMTCount();

  if(Hist)
	Hist->Fill(nbPE, nbPMTsHits);

}

void CollectEVPEAndHits(RAT::DS::EV *ev, TH2D *Hist){

  const double nbPE = ev->GetTotalCharge();
  const double nHits = ev->Nhits();

  if(Hist)
	Hist->Fill(nbPE, nHits);

}

void CollectPromptPEAndHits(RAT::DS::MC *mc, TH2D *Hist){

  // Define container for PMT Hit Times
  vector<float> HitTime;
  vector<float>::iterator itHitTime;  // declare an iterator to a vector of strings

  // Define prompt hits
  const float threshPromptHits = 9.; // ns

  const int nbPMTsHits = mc->GetMCPMTCount();
  int PromptHits = 0;

  if(nbPMTsHits>0){

	// Loop over PMT hits
	for(int iPMT = 0; iPMT<nbPMTsHits; iPMT++) {
	  RAT::DS::MCPMT *PMT = mc->GetMCPMT(iPMT);

	  if(PMT){

		const int nbPhotons = PMT->GetMCPhotonCount();
		for(int iPhoton = 0; iPhoton<nbPhotons; iPhoton++){
		  RAT::DS::MCPhoton *photon = PMT->GetMCPhoton(iPhoton);

		  if(photon) {
			HitTime.emplace_back(photon->GetFrontEndTime());
		  } // END if photon

		} // END for iPhoton

	  } // END if PMT

	} // END for iPMT

	// Search for prompts hits
	sort(HitTime.begin(),HitTime.end());
	adjacent_difference(HitTime.begin(),HitTime.end(),HitTime.begin());
	HitTime[0] = 0.;

	for(itHitTime = HitTime.begin(); itHitTime != HitTime.end(); itHitTime++){
	  if(accumulate(HitTime.begin(),itHitTime,0) > threshPromptHits){
		break;
	  }
	  PromptHits++;
	}

  }

  const int nbPE = mc->GetNumPE();

  if(Hist)
	Hist->Fill(nbPE, PromptHits);

}


#endif // _ANALYSISDEFINITIONS_HH_

