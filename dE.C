#include <utils.hh>

#include <Analyzer.hh>
#include <FlatParticle.hh>
#include <MCFunctions.hh>

void CalculatedE(const char *filename, unsigned long nEvts = 0){

  auto fA = new Analyzer(filename);

  auto hdEdX = new TH1D("hdEdX", "dE/dX ; dE/dX (MeV/mm) ; Counts", 1000, 0., 1000.);

  auto hdEdXPerStep = new TH1D("hdEdXPerStep", "dE/dX per Step ; dE/dX (MeV/mm) ; Counts",
							   100, 0., 50.);

  nEvts = nEvts > 0? nEvts : fA->GetNEvts();

  for(auto iEvt=0; iEvt<nEvts; iEvt++){

	// auto * mc = GetRATMCOnEvt(fA, iEvt);
	// auto NbTracks = mc->GetMCTrackCount();
	//
	// vector<ComplexParticle> vPart;
	//
	// for(auto iTrack=0; iTrack<NbTracks; iTrack++) {
	//
	//   auto *mctrack = mc->GetMCTrack(iTrack);
	//
	//   if (IsParticle(mctrack, "e-")) {
	//
	// 	auto NbTrackSteps = mctrack->GetMCTrackStepCount();
	//
	// 	auto Ekin = mctrack->GetMCTrackStep(0)->GetKE();
	// 	auto X =  mctrack->GetMCTrackStep(0)->GetLength();
	//
	// 	for (auto iStep = 1; iStep < NbTrackSteps; iStep++) {
	//
	// 	  auto *mctrackstep = mctrack->GetMCTrackStep(iStep);
	//
	// 	  hdEdXPerStep->Fill(abs(mctrackstep->GetKE() - Ekin)/abs(mctrackstep->GetLength()) - X);
	// 	}
	//   }
	// }


	vector<ComplexParticle> vPart = GetdEdX(fA, iEvt, "e-");

	for(auto part:vPart){
	  hdEdX->Fill(part.GetDe()/part.GetDx());
	  for(auto S:part.GetVStep()){
		hdEdXPerStep->Fill(S.GetDe()/S.GetDx());
	  }
	}

  }

  PlotAHist(hdEdX);

  PlotAHist(hdEdXPerStep);

}