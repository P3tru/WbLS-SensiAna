//
// Created by zsoldos on 2/21/20.
//

#ifndef _EVFUNCTIONS_HH_
#define _EVFUNCTIONS_HH_

#include <RAT/DS/EV.hh>

#include <Analyzer.hh>

RAT::DS::EV * GetRATEVOnEvt(Analyzer *fAnalyzer, unsigned int iEvt=0, unsigned int iEV=0);

#endif