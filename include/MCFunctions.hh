//
// Created by zsoldos on 2/24/20.
//

#ifndef _MCFUNCTIONS_HH_
#define _MCFUNCTIONS_HH_


///////////////////////// STL C/C++ /////////////////////////
#include <vector>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/MC.hh>

/////////////////////////   ROOT  ///////////////////////////
#include <TH2D.h>

/////////////////////////   USER  ///////////////////////////
#include "Analyzer.hh"
#include "HitClass.hh"
#include "FlatParticle.hh"

using namespace std;

// Print basic track info
// Useful for debugging
void PrintTrackInfo(RAT::DS::MC *mc, unsigned int iTrack);

// Access the MC object from RAT
// at iEvt (from MCTree)
RAT::DS::MC * GetRATMCOnEvt(Analyzer *fAnalyzer, unsigned int iEvt=0);

// Get vector<FlatPhoton> of photons in evt, storing their wl and create proc
vector<FlatPhoton> GetPhotonsFromEvt(RAT::DS::MC * mc);

// Fill vector<Hit> with MC Info
vector<Hit> GetMCHitCollection(Analyzer *fAnalyzer, unsigned int iEvt=0, bool isSource=false);

// Fill vector<Hit> with MC Info from one mother particle
vector<Hit> GetVHitsFromPart(Analyzer *fAnalyzer, unsigned int iEvt=0,
							 string sPartName="e+",
							 vector<ComplexParticle> *vCPart = NULL);

// Get info from MC particle in a flat way
void GetPartInfoFromTrackStep(FlatParticle *fp, RAT::DS::MCTrack *mctrack);

// Get ID from MC particle in a flat way
void GetPartIDFromTrackStep(G4Particle *p, RAT::DS::MCTrack *mctrack);

// Get vector<ComplexParticle> of MC object
vector<ComplexParticle> GetVPart(RAT::DS::MC * mc, const string& sPartName);

// Fill vector<ComplexParticle> with daughters photons
vector<FlatPhoton> GetAndFillPhotonsToVPart(RAT::DS::MC * mc, vector<ComplexParticle> *vPart);

bool IsParticle(RAT::DS::MCTrack *mctrack, const string& name);

// Fill dEdX for a ComplexParticle
void FilldEdX(ComplexParticle &CP, Analyzer *fAnalyzer, unsigned int iEvt=0);

// Sum total dEdX for a vector<ComplexParticle> from by an event
double SumdEdX(vector<ComplexParticle> vCP);

// Get vector protons with photons
vector<ComplexParticle> GetProtonLY(Analyzer *fAnalyzer, unsigned int iEvt);

// Create spectrum E incoming VS proton recoil
void FillProtonRecoilSpectrum(TH2D* h2D, Analyzer *fAnalyzer, unsigned iEvt);

// Fill Quenching spectrum (Check birks law for example)
void FillQuenchingSpectrum(TH2D *h2D, Analyzer *fAnalyzer, unsigned int iEvt);

#endif