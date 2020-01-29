//
// Created by zsoldos on 12/5/19.
//

#ifndef _TEMPLATEANALYSIS_HH_
#define _TEMPLATEANALYSIS_HH_

///////////////////////// STL C/C++ /////////////////////////

/////////////////////////   BOOST   /////////////////////////

/////////////////////////   ROOT   //////////////////////////
#include <TH2D.h>
#include <TRandom3.h>

/////////////////////////   RAT   ///////////////////////////

/////////////////////////   USER  ///////////////////////////
#include "TFileAnalysis.hh"
#include "AnalysisDefinitions.hh"

void ThreadAnalysis(TFileAnalysis<TH2D> file, TH2D *h2D){

  file.SetHist(h2D);

  file.DoAnalysis(CollectPromptPEAndHits);

}

void TestThreadAnalysis(TH2D *h2D){

  TRandom3 *r3;

  TH2D *h = dynamic_cast<TH2D *>(h2D->Clone());

  h->Fill(r3->Gaus(100,10), r3->Gaus(100,10));

}

#endif //_TEMPLATEANALYSIS_HH_
