///////////////////////// STL C/C++ /////////////////////////
#include <string>
#include <vector>
#include <fstream>

/////////////////////////   ROOT   //////////////////////////

/////////////////////////   USER  ///////////////////////////
#include "AnalyzerFunctions.hh"
#include "HitClass.hh"
#include "HitFunctions.hh"
#include "EVFunctions.hh"
#include "MCFunctions.hh"

#include "cnpy.h"

void AddFAnalyzers(vector<Analyzer*> *vFAnalyzer, const string& inputName, const string& listName){

  if(!inputName.empty()){
	vFAnalyzer->push_back(new Analyzer(inputName.c_str()));
  }

  if(!listName.empty()) {

	// Get ready to read inside inputName
	// each file name
	string line;
	ifstream file(listName);

	while (getline(file, line)) {

	  if (line.compare(0, 1, "#") == 0) {

		continue;

	  } else {

		vFAnalyzer->push_back(new Analyzer(line.c_str()));

	  }

	}

  }


}

void GetVHitAndDumpFlatNPZ(Analyzer *fAnalyzer, unsigned iEvt, const string& NPZName, const string& mode){

  vector<Hit> vHit = GetMCHitCollection(fAnalyzer, iEvt, true);
  auto MCID = GetRATMCOnEvt(fAnalyzer, iEvt)->GetID();

  const auto NHits = vHit.size();
  vector<double> vNPY = FlatenVHit(vHit, MCID);

  cnpy::npz_save(NPZName,
				 Form("Evt%d", iEvt),
				 &vNPY[0],
				 {NHits, 6}, // NHits vector of 6 dimension {X, Y, Z, T, Source, MCID}
				 mode);

}

FlatParticle GetPrimaryParticleInfo(Analyzer *fAnalyzer, unsigned int iEvt){

  RAT::DS::MC * mc = GetRATMCOnEvt(fAnalyzer, iEvt);
  RAT::DS::MCParticle *prim = mc->GetMCParticle(0);

  return FlatParticle(prim->GetParticleName(), prim->GetPosition(), prim->GetMomentum().Unit(), prim->GetKE());

}

void PlotQuenchingSpectrum(Analyzer *fAnalyzer){

  TH2D *hBirks = new TH2D("hBirks", "Birks Law according to rat-pac",
						  100,-0.05,10.05,
						  100,0.55,1.05);

  hBirks->SetTitle("10MeV protons; dLY/dX; (dE/dX)/(1+dE/dX))");

  for(unsigned int iEvt=0; iEvt<fAnalyzer->GetNEvts(); iEvt++){

	FillQuenchingSpectrum(hBirks, fAnalyzer, iEvt);

  }

  hBirks->Draw("COLZ");

}
