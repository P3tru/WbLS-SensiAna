//
// Created by zsoldos on 11/26/19.
//

#include <TFile.h>
#include <TTree.h>
#include <RAT/DS/Root.hh>
#include "PlotCollectedPE.hh"

int main(int argc, const char **argv) {

  TFile *f = new TFile("");
  TTree *tree = (TTree*) f->Get("T");

  RAT::DS::Root *rds = new RAT::DS::Root();
  tree->SetBranchAddress("ds", &rds);

  int nEvents = tree->GetEntries();

  for (int i = 0; i < nEvents; i++) {
	tree->GetEntry(i);
	RAT::DS::MC *mc = rds->GetMC();
	RAT::DS::MCParticle *prim = mc->GetMCParticle(0);
  }

  return 0;
}
