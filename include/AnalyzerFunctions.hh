//
// Created by zsoldos on 2/21/20.
//

#ifndef _ANALYZERFUNCTIONS_HH_
#define _ANALYZERFUNCTIONS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <string>
#include <vector>

/////////////////////////   ROOT   //////////////////////////

/////////////////////////   USER   //////////////////////////
#include "Analyzer.hh"
#include "FlatParticle.hh"

using namespace std;

// Add list of analyzer to vector
void AddFAnalyzers(vector<Analyzer*> *vFAnalyzer, const string& inputName, const string& listName = "");

// Dump vector<Hit> inside a npz file
void GetVHitAndDumpFlatNPZ(Analyzer *fAnalyzer, unsigned iEvt, const string& NPZName, const string& mode="a");

// Get primary particle info from
FlatParticle GetPrimaryParticleInfo(Analyzer *fAnalyzer, unsigned int iEvt);

// Plot quenching spectrum
void PlotQuenchingSpectrum(Analyzer *fAnalyzer);

#endif