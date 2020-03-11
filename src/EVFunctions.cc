#include <EVFunctions.hh>

RAT::DS::EV * GetRATEVOnEvt(Analyzer *fAnalyzer, unsigned int iEvt, unsigned int iEV){

  fAnalyzer->GetTreeMc()->GetEntry(iEvt);

  // Access RAT MC info and the summary
  // Summary useful to get nCer photons, nScint photons, etc...
  return fAnalyzer->GetRds()->GetEV(iEV);

}
