//
// Created by zsoldos on 11/26/19.
//

void PlotTracks(const char *file=NULL){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(56);

  TFile *f = new TFile(file);
  TTree *tree = (TTree*) f->Get("T");

  RAT::DS::Root *rds = new RAT::DS::Root();
  tree->SetBranchAddress("ds", &rds);

  int nEvents = tree->GetEntries();

  TCanvas *c1;

  for (int iEvt = 0; iEvt < nEvents; iEvt++) {
    tree->GetEntry(iEvt);

    // Access RAT MC info and the summary
    // Summary useful to get nCer photons, nScint photons, etc...
    RAT::DS::MC * mc = rds->GetMC();

    const int nbTracks = mc->GetMCTrackCount();

    for(int iTrack = 0; iTrack<nbTracks; iTrack++){

      RAT::DS::MCTrack * mctrack = mc->GetMCTrack(iTrack);

      int originalID = -1;

      if(mctrack->GetParentID() == 0){
	std::cout << "Evt#" << iEvt
		  << " Primary particle: " << mctrack->GetParticleName()
		  << " With trackID: " << mctrack->GetID() << std::endl;

	originalID=mctrack->GetID();

	mctrack->PruneIntermediateMCTrackSteps();
	std::cout<< "Process at end of track: " << mctrack->GetMCTrackStep(1)->GetProcess() << std::endl;	
      } else if(mctrack->GetParentID() == 1){
	std::cout << "Evt#" << iEvt
		  << " Primary particle: " << mctrack->GetParticleName()
		  << " With trackID: " << mctrack->GetID() << std::endl;

	mctrack->PruneIntermediateMCTrackSteps();
	std::cout<< "Process at end of track: " << mctrack->GetMCTrackStep(1)->GetProcess() << std::endl;  
      }


      
      // const int nbTrackSteps = mctrack->GetMCTrackStepCount();
      // for(int iTrackStep = 0; iTrackStep<nbTrackSteps; iTrackStep++){	
      // }
      
    } // END FOR iTrack
    
  } // END FOR iEvt

  //////////////////////////////
  // #### #### DRAW #### #### //
  //////////////////////////////

  c1 = new TCanvas("c1","c1",1600,1200);
  c1->SetGrid();
  
  
}
