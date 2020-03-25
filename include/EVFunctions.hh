//
// Created by zsoldos on 2/21/20.
//

#ifndef _EVFUNCTIONS_HH_
#define _EVFUNCTIONS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <vector>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/EV.hh>

/////////////////////////   ROOT  ///////////////////////////

/////////////////////////   USER  ///////////////////////////
#include "Analyzer.hh"
#include "HitClass.hh"

using namespace std;

// Access the EV object from RAT
// at iEvt (from MCTree) and for iEV (from EVCount)
RAT::DS::EV * GetRATEVOnEvt(Analyzer *fAnalyzer, unsigned int iEvt=0, unsigned int iEV=0);

// Fill vector<Hit> with EV Info
vector<Hit> GetEVHitCollection(Analyzer *fAnalyzer, unsigned int iEvt=0, unsigned int iEV=0);

#endif