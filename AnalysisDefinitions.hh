//
// Created by zsoldos on 11/26/19.
//

#ifndef _ANALYSISDEFINITIONS_HH_
#define _ANALYSISDEFINITIONS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <fstream>
#include <sstream>
#include <regex>

/////////////////////////   BOOST   /////////////////////////
#include <boost/filesystem/path.hpp>

/////////////////////////   ROOT   //////////////////////////
#include <TH2D.h>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/MC.hh>
#include <RAT/DS/Root.hh>

void CollectPEAndHits(RAT::DS::MC *mc, TH2D *Hist){

  // Get Nb of PMTs which at least 1hit
  const int nbPE = mc->GetNumPE();
  const int nbPMTsHits = mc->GetMCPMTCount();

  if(Hist)
	Hist->Fill(nbPE, nbPMTsHits);

}


#endif // _ANALYSISDEFINITIONS_HH_

