//
// Created by zsoldos on 2/21/20.
//

#ifndef _ANALYZERFUNCTIONS_HH_
#define _ANALYZERFUNCTIONS_HH_

#include <Analyzer.hh>

#include <TH1D.h>

// PROTOTYPES FOR ANALYZER FUNCTIONS
void GetHitTime(Analyzer *Ana, unsigned int iEvt = 0, TH1D* hHitTime = nullptr);

#endif