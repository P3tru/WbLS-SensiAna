//
// Created by zsoldos on 2/24/20.
//

#ifndef _HITFUNCTIONS_HH_
#define _HITFUNCTIONS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

/////////////////////////   USER  ///////////////////////////
#include "HitClass.hh"

void PrintVHits(vector<Hit> const Hits);
vector<Hit> CorrectDelayedHits(vector<Hit> rawHits,
							   const Hit& PreTrig = Hit(TVector3(0.,0.,0.), 0., 0.));
void RemoveHitsAfterCut(vector<Hit> &Hits,
						const Hit& hCut = Hit(TVector3(0.,0.,0.), 0., 0.));
void GetNPEAndNHitsFromHits(vector<Hit> Hits, double *NPE, double *NHits);

#endif //_HITFUNCTIONS_HH_