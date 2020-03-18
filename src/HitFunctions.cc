//
// Created by zsoldos on 2/24/20.
//

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <TH2D.h>

/////////////////////////   USER  ///////////////////////////
#include "HitClass.hh"


Hit operator-(Hit h1, Hit h2){
  h1.SetT(h1.GetT()-h2.GetT());
  return h1;
}

Hit operator+(Hit h1, Hit h2){
  h1.SetT(h1.GetT()+h2.GetT());
  return h1;
}

void PrintVHits(vector<Hit> const Hits){
  for(auto itHit : Hits){
	itHit.Print();
  }
}

void RemoveHitsAfterCut(vector<Hit> &Hits, const Hit& hCut){

  auto itHit = find_if(Hits.begin(),
					   Hits.end(),
					   HitCut(hCut));
  Hits.erase(itHit,Hits.end());

}

vector<Hit> CorrectDelayedHits(vector<Hit> rawHits, const Hit& PreTrig){

  vector<Hit> vHitDelayed_CORRECTED;
  for(auto itHit : rawHits){

	itHit = itHit - rawHits[0] + PreTrig;
	vHitDelayed_CORRECTED.emplace_back(itHit);

  }

  return vHitDelayed_CORRECTED;

}

TH2D *GetHQVST(vector<Hit> vHit){

  auto *hQVST = new TH2D("hQVST", "Q VS T",
						 100, -0.05, 10.05,
						 100, -0.05, 10.05);

  hQVST->GetXaxis()->SetTitle("Q (PE)");
  hQVST->GetYaxis()->SetTitle("T (ns)");

  sort(vHit.begin(), vHit.end());

  vector<Hit> vHitCor = CorrectDelayedHits(vHit, Hit(TVector3(0.,0.,0.), 0., 0.));

  for(auto hit: vHitCor){
    hQVST->Fill(hit.GetQ(), hit.GetT());
  }

  return hQVST;


}

void GetNPEAndNHitsFromHits(vector<Hit> Hits, double *NPE, double *NHits){

  double mNPE=0.;
  int mNHits=0;

  for(auto h:Hits){
	mNPE+=h.GetQ();
	mNHits++;
  }

  *NPE=mNPE;
  *NHits=mNHits;

}