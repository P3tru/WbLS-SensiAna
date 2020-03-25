//
// Created by zsoldos on 2/24/20.
//

/////////////////////////   USER  ///////////////////////////
#include "HitFunctions.hh"

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

void SortVHits(vector<Hit> *vHit){

  sort(vHit->begin(), vHit->end());

}

vector<Hit> SplitVHits(vector<Hit> *vHit, double DAQWindow) {

  vector<Hit> Empty;

  unsigned int SizeInit = vHit->size();

  if (vHit->size() == 0) {

	return Empty;

  } else if (vHit->size() > 0) {

	SortVHits(vHit);
	for (unsigned int iHit = 1; iHit < vHit->size(); iHit++) {

	  const double dT = vHit->at(iHit).GetT() - vHit->at(0).GetT();

	  if (dT > DAQWindow) {

		if (vHit->begin() + iHit < vHit->end()) {

		  try {

			vector<Hit> vHitDelayed(vHit->begin() + iHit, vHit->end());

			vHit->erase(vHit->begin() + iHit, vHit->end());

			return vHitDelayed;

		  } catch (vector<Hit> const *vHit) {

			for (auto h : *vHit)
			  h.Print();

			return Empty;

		  }

		} else if (vHit->begin() + iHit >= vHit->end()) {

		  return Empty;

		}
	  }

	}

  }

  if (vHit->size() == SizeInit) {
	return Empty;
  }

}

void RemoveHitsAfterCut(vector<Hit> &Hits, const Hit& hCut){

  auto itHit = find_if(Hits.begin(),
					   Hits.end(),
					   HitCut(hCut));
  Hits.erase(itHit,Hits.end());

}

vector<Hit> ResetTVHits(vector<Hit> rawHits, const Hit& PreTrig){

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

  vector<Hit> vHitCor = ResetTVHits(vHit, Hit(TVector3(0.,0.,0.), 0., 0.));

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