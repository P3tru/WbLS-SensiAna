//
// Created by zsoldos on 11/26/19.
//

using namespace std;

void PrintParticleInfo(RAT::DS::MCParticle *mcp);
void PrintTrackStepInfo(RAT::DS::MCTrackStep *InitStep);

#define PI 3.14159

void PlotTruePos(const char *file=NULL){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(56);

  TFile *f = new TFile(file);
  TTree *tree = (TTree*) f->Get("T");

  RAT::DS::Root *rds = new RAT::DS::Root();
  tree->SetBranchAddress("ds", &rds);

  auto nEvents = tree->GetEntries();
  // nEvents = 1;

  // #### #### #### #### #### #### #### #### #### #### #### #### //
  // ####                 SETUP HISTOGRAMS                  #### //
  // #### #### #### #### #### #### #### #### #### #### #### #### //

  const int nPrimParticle = 2;
  const int nAxis = 3;

  TH1D *hDiffPos[nAxis];

  for(auto iAx=0; iAx<nAxis; iAx++){

	hDiffPos[iAx] = new TH1D(Form("hDiffPos%d",iAx),Form("Diff Pos Ax %d",iAx),
							 201,-100.5,100.5);

  }


  TH2D *hTVSPhiDiff = new TH2D("hTVSPhiDiff", "#theta VS #phi",
							   18, 0, PI,
							   36, 0, 2*PI);

  TH1D *hPos[nPrimParticle][nAxis];
  TH2D *hTVSPhi[nPrimParticle];

  for(auto iPart = 0; iPart<nPrimParticle; iPart++){

    hTVSPhi[iPart] = new TH2D(Form("hTVSPhi%d",iPart),"#theta VS #phi",
							  18, 0, PI,
							  36, -PI, PI);

	for(auto iAx=0; iAx<nAxis; iAx++){

	  hPos[iPart][iAx] = new TH1D(Form("hPos%d%d",iPart,iAx),
								  Form("Prtcle %d Ax %d",iPart,iAx),
								  201,-100.5,100.5);

	}

  }

  TCanvas *c1;

  for (int iEvt = 0; iEvt < nEvents; iEvt++) {
    tree->GetEntry(iEvt);

    if (iEvt % 1000 == 0)
      cout << "iEvt#" << iEvt << endl;

    // Access RAT MC info and the summary
    // Summary useful to get nCer photons, nScint photons, etc...
    RAT::DS::MC * mc = rds->GetMC();

    auto nbTracks = mc->GetMCTrackCount();

	for(int iTrack = 0; iTrack<nbTracks; iTrack++){

	  RAT::DS::MCTrack * mctrack = mc->GetMCTrack(iTrack);

	  if(mctrack->GetParentID() == 0){

		mctrack->PruneIntermediateMCTrackSteps();
		RAT::DS::MCTrackStep *InitStep = mctrack->GetMCTrackStep(0);
		RAT::DS::MCTrackStep *EndStep  = mctrack->GetMCTrackStep(1);

		if(mctrack->GetParticleName() == "e+"){

		  if(EndStep->GetProcess() == "annihil"){

			for(auto iAx=0; iAx<nAxis; iAx++) {
			  hPos[0][iAx]->Fill(EndStep->GetEndpoint()(iAx));
			}

			hTVSPhi[0]->Fill(EndStep->GetEndpoint().Theta(),
							 EndStep->GetEndpoint().Phi());

		  }

		} else if(mctrack->GetParticleName() == "neutron"){

		  if(EndStep->GetProcess() == "nCapture"){

			for(auto iAx=0; iAx<nAxis; iAx++) {
			  hPos[1][iAx]->Fill(EndStep->GetEndpoint()(iAx));
			}

			hTVSPhi[1]->Fill(EndStep->GetEndpoint().Theta(),
							 EndStep->GetEndpoint().Phi());

		  }

		}

	  }

    } // END FOR iTrack
    
  } // END FOR iEvt

  //////////////////////////////
  // #### #### DRAW #### #### //
  //////////////////////////////

  c1 = new TCanvas("cPositron","cPositron",1600,1200);
  c1->Divide(2,2);
  c1->cd(1);
  c1->SetGrid();
  hTVSPhi[0]->Draw("COLZ");
  for(int iAx = 0 ; iAx < nAxis ; iAx++){
    c1->cd(2+iAx);
	c1->SetGrid();
    hPos[0][iAx]->Draw();
  }

  c1 = new TCanvas("cNeutron","cNeutron",1600,1200);
  c1->Divide(2,2);
  c1->cd(1);
  c1->SetGrid();
  hTVSPhi[1]->Draw("COLZ");
  for(int iAx = 0 ; iAx < nAxis ; iAx++){
	c1->cd(2+iAx);
	c1->SetGrid();
	hPos[1][iAx]->Draw();
  }

  c1 = new TCanvas("cDiff","cDiff",1600,200);
  c1->Divide(2,2);
  c1->cd(1);
  c1->SetGrid();

}

void PrintParticleInfo(RAT::DS::MCParticle *mcp){

  cout << "GEN Particle name: " << mcp->GetParticleName() << endl;
  cout << " With KinE:" << mcp->GetKE() << endl;
  cout << " At pos:"
	   << mcp->GetPosition().X() << " "
	   << mcp->GetPosition().Y() << " "
	   << mcp->GetPosition().Z() << " " << endl;
  cout << " At t:" << mcp->GetTime();
  cout << endl;

}

void PrintTrackStepInfo(RAT::DS::MCTrackStep *Step){

  cout << "INIT STEP" << endl;
  cout << " With length: " << Step->GetLength() << endl;
  cout << " With endpoint: "
	   << " X:" << Step->GetEndpoint().X()
	   << " Y:" << Step->GetEndpoint().Y()
	   << " Z:" << Step->GetEndpoint().Z() << endl;
  cout << " At time:" << Step->GetLocalTime() << endl;
  cout << " At proper time:" << Step->GetProperTime() << endl;
  cout << " With KinE:" << Step->GetKE() << endl;

}

void FillTruePos(TH2D *hTVSPhi, RAT::DS::MCTrackStep *Step){

  hTVSPhi->Fill(Step->GetEndpoint().Theta(), Step->GetEndpoint().Phi());

}