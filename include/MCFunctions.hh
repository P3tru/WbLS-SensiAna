//
// Created by zsoldos on 2/24/20.
//

#ifndef _MCFUNCTIONS_HH_
#define _MCFUNCTIONS_HH_


///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/Run.hh>
#include <RAT/DS/Root.hh>

#include "Analyzer.hh"
#include "HitClass.hh"

using namespace std;

RAT::DS::MC * GetRATMCOnEvt(Analyzer *fAnalyzer, unsigned int iEvt=0);
vector<Hit> GetHitCollection(Analyzer *fAnalyzer, unsigned int iEvt=0);
vector<Hit> SplitHitCollection(vector<Hit> *vHit, double PromptWindow);

#endif