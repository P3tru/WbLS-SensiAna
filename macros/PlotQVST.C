//
// Created by zsoldos on 11/26/19.
//

void PlotQVST(const char *file=NULL){

  // gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  TFile *f = new TFile(file);
  TTree *tree = (TTree*) f->Get("T");

  RAT::DS::Root *rds = new RAT::DS::Root();
  tree->SetBranchAddress("ds", &rds);

  int nEvents = tree->GetEntries();

  // auto Calib = 13.28; // Hits/MeV for water
  auto Calib = 18.21; // Hits/MeV for wbls 1%
  // auto Calib = 50.02; // Hits/MeV for wbls 5%
  // auto Calib = 85.79; // Hits/MeV for wbls 10%
  

  // Histogram definition
  auto NbPEMin=-0.5;
  auto NbPEMax=100.5;
  auto NbPEBins=101;
  auto NbHitsMin=-0.5;
  auto NbHitsMax=100.5;
  auto NbHitsBins=101;
  auto NbTMin=-0.5;
  auto NbTMax=200.5;
  auto NbTBins=201;
  TH2D *hPEVST = new TH2D("hPEVST","Nb PE Collected VS T",
			  NbTBins, NbTMin, NbTMax,
			  NbPEBins, NbPEMin, NbPEMax);
  hPEVST->GetXaxis()->SetTitle("T");
  hPEVST->GetYaxis()->SetTitle("Nb PE");

  vector<TH1D*> vecTPromptLikeEvt;

  TH2D *hHitsVST = new TH2D("hHitsVST","Nb Hits Collected VS T",
			  NbTBins, NbTMin, NbTMax,
			  NbHitsBins, NbHitsMin, NbHitsMax);
  hHitsVST->GetXaxis()->SetTitle("T");
  hHitsVST->GetYaxis()->SetTitle("Nb Hits");

  TH1D *hNeutronPromptE = new TH1D("hNeutronPromptE","E equivalent prompt evt",
				   101,-0.05,10.05);

  TCanvas *c1;

  for (int iEvt = 0; iEvt < nEvents; iEvt++) {
    tree->GetEntry(iEvt);

    // Access RAT MC info and the summary
    // Summary useful to get nCer photons, nScint photons, etc...
    RAT::DS::MC * mc = rds->GetMC();
    RAT::DS::MCSummary * mcs = mc->GetMCSummary();

    // Get Nb of PMTs which at least 1hit
    auto nbPMTsHits = mc->GetMCPMTCount();
    auto nbPromptE = 0.;

    for(int iPMT = 0; iPMT<nbPMTsHits; iPMT++){

      RAT::DS::MCPMT * pmt = mc->GetMCPMT(iPMT);

      auto nbPhotonCounts = pmt->GetMCPhotonCount();
      for(int iPhoton=0;iPhoton<nbPhotonCounts;iPhoton++){

	RAT::DS::MCPhoton * photon = pmt->GetMCPhoton(iPhoton);

	hPEVST->Fill(photon->GetFrontEndTime(), photon->GetCharge());

	if(photon->GetFrontEndTime() < 200)
	  nbPromptE += (double)(photon->GetCharge())/(double)(Calib);
	
      } // END for iPHOTON
      
      

    } // END FOR iPMT
    
    hNeutronPromptE->Fill(nbPromptE);

    if(nbPromptE > 1.){
      
      TH1D *hT = new TH1D(Form(""),"",
			  NbTBins, NbTMin, NbTMax); 
      
      for(int iPMT = 0; iPMT<nbPMTsHits; iPMT++){
	
	RAT::DS::MCPMT * pmt = mc->GetMCPMT(iPMT);

	auto nbPhotonCounts = pmt->GetMCPhotonCount();
	for(int iPhoton=0;iPhoton<nbPhotonCounts;iPhoton++){

	  RAT::DS::MCPhoton * photon = pmt->GetMCPhoton(iPhoton);

	  hT->Fill(photon->GetFrontEndTime());
	
	} // END for iPHOTON
      }

      vecTPromptLikeEvt.emplace_back(hT);

    }

  } // END FOR iEvt

  //////////////////////////////
  // #### #### DRAW #### #### //
  //////////////////////////////

  c1 = new TCanvas("c1","c1",800,600);
  c1->SetGrid();
  hPEVST->Draw("COLZ");

  c1 = new TCanvas("c2","c2",800,600);
  c1->SetGrid();
  hNeutronPromptE->Draw("");

  int count=0;
  for(auto &h: vecTPromptLikeEvt){
    c1 = new TCanvas(Form("c3%d",count++),"cEvt",800,600);
    c1->SetGrid();
    h->Draw("SAME PLC PMC");
  }
  cout << vecTPromptLikeEvt.size() << endl;
  
}
