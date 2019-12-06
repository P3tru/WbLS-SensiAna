//
// Created by zsoldos on 12/5/19.
//

#ifndef _TEMPLATEANALYSIS_HH_
#define _TEMPLATEANALYSIS_HH_

void CollectPEAndHits(RAT::DS::MC *mc, TH2D *Hist){

  // Get Nb of PMTs which at least 1hit
  const int nbPMTsHits = mc->GetMCPMTCount();
  const int nbPE = mc->GetNumPE();

  if(Hist)
	Hist->Fill(nbPE, nbPMTsHits);

}


#endif //_TEMPLATEANALYSIS_HH_
