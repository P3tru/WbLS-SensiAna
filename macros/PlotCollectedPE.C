//
// Created by zsoldos on 11/26/19.
//

void PlotCollectedPE(const char *file=NULL){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(56);

  TFile *f = new TFile(file);
  TTree *tree = (TTree*) f->Get("T");

  RAT::DS::Root *rds = new RAT::DS::Root();
  tree->SetBranchAddress("ds", &rds);

  int nEvents = tree->GetEntries();

  // Histogram definition
  const int NbPEMin=0;
  const int NbPEMax=2e3;
  const int NbPEBins=1e2;
  const int NbHitsMin=0;
  const int NbHitsMax=2e3;
  const int NbHitsBins=1e2;
  TH2D *hNbPEVSNbHits = new TH2D("hNbPEVSNbHits","Nb PE Collected VS Nb PMT hits",
				 NbPEBins, NbPEMin, NbPEMax,
				 NbHitsBins, NbHitsMin, NbHitsMax);
  hNbPEVSNbHits->GetXaxis()->SetTitle("Nb PE");
  hNbPEVSNbHits->GetYaxis()->SetTitle("Nb PMTs hits");

  TCanvas *c1;

  for (int iEvt = 0; iEvt < nEvents; iEvt++) {
    tree->GetEntry(iEvt);

    // Access RAT MC info and the summary
    // Summary useful to get nCer photons, nScint photons, etc...
    RAT::DS::MC * mc = rds->GetMC();
    RAT::DS::MCSummary * mcs = mc->GetMCSummary();

    // Get Nb of PMTs which at least 1hit
    const int nbPMTsHits = mc->GetMCPMTCount();
    
    for(int iPMT = 0; iPMT<nbPMTsHits; iPMT++){
      RAT::DS::MCPMT * pmt = mc->GetMCPMT(iPMT);

      /* // Do stuff to PMT */
      /* int nPEPMT = pmt->GetMCPhotonCount(); */
      /* double QPMT = pmt->GetTotalCharge(); */

    } // END FOR iPMT

    const int nbParticle = mc->GetMCParticleCount();
    for(int iParticle = 0; iParticle<nbParticle; iParticle++){

      RAT::DS::MCParticle * particle = mc->GetMCParticle(iParticle);

      // Do stuff to particle
	  
    } // END FOR iParticle

    const int nbPE = mc->GetNumPE();

    hNbPEVSNbHits->Fill(nbPE, nbPMTsHits);

  } // END FOR iEvt

  //////////////////////////////
  // #### #### DRAW #### #### //
  //////////////////////////////

  c1 = new TCanvas("c1","c1",1600,1200);
  c1->SetGrid();
  hNbPEVSNbHits->Draw("COLZ");

  c1 = new TCanvas("c2","c2",1600,1200);
  TH1D *hNbPE = hNbPEVSNbHits->ProjectionX();
  hNbPE->SetTitle("Nb PE collected per event");
  hNbPE->Draw();
  hNbPE->Fit("gaus");
  
  
}
