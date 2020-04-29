//
// Created by zsoldos on 2/21/20.
//

///////////////////////// STL C/C++ /////////////////////////

/////////////////////////   RAT   ///////////////////////////

/////////////////////////   ROOT  ///////////////////////////

/////////////////////////   USER  ///////////////////////////
#include "MCFunctions.hh"

using namespace std;

// Convert photon MeV to nm
static double MeV2lambda(double MeV){

  const double hc = 1.23984193e-3; //MeV.nm
  return hc/MeV;

}

bool operator==(G4Particle const& a, G4Particle const& b){
  return a.GetTrackId() == b.GetTrackId() && a.GetParentId() == b.GetParentId();
}

RAT::DS::MC * GetRATMCOnEvt(Analyzer *fAnalyzer, unsigned int iEvt){

  fAnalyzer->GetTreeMc()->GetEntry(iEvt);

  // Access RAT MC info and the summary
  // Summary useful to get nCer photons, nScint photons, etc...
  return fAnalyzer->GetRds()->GetMC();

}

static int ConvertProcess(string proc){
  if(proc == "Cerenkov"){
	return 0;
  } else if (proc == "Scintillation") {
	return 1;
  } else if(proc.rfind("Reemission", 0) == 0) {
	return 2;
  } else {
	return -1;
  }
}

vector<FlatPhoton> GetPhotonsFromEvt(RAT::DS::MC * mc){

  vector<FlatPhoton> vFP;

  for(auto iTrack=0; iTrack<mc->GetMCTrackCount(); iTrack++){

    auto track = mc->GetMCTrack(iTrack);

    if(IsParticle(track, "opticalphoton")){

	  auto first = track->GetMCTrackStep(0);
	  FlatPhoton p;
	  GetPartIDFromTrackStep(&p, track);
	  p.SetWl(MeV2lambda(first->GetKE()));
	  p.SetProcess(ConvertProcess(first->GetProcess()));

	  vFP.emplace_back(p);

    }

  }

  return vFP;

}

vector<Hit> GetMCHitCollection(Analyzer *fAnalyzer, unsigned int iEvt, bool isSource){

  RAT::DS::MC * mc = GetRATMCOnEvt(fAnalyzer, iEvt);

  vector<Hit> vHit;

  vector<FlatPhoton> vFP;
  if(isSource){
	vFP = GetPhotonsFromEvt(mc);
  }

  for (int iPMT = 0; iPMT < mc->GetMCPMTCount(); iPMT++) {

	auto mcPMT = mc->GetMCPMT(iPMT);
	auto PMTID = mcPMT->GetID();
	auto PMTPos = fAnalyzer->GetRun()->GetPMTInfo()->GetPosition(PMTID);

	auto NbPhotonCounts = mcPMT->GetMCPhotonCount();

	for (int iP = 0; iP < NbPhotonCounts; iP++) {

	  auto mcPhoton = mcPMT->GetMCPhoton(iP);
	  // Get Q and T
	  auto Q = mcPhoton->GetCharge();
	  auto T = mcPhoton->GetHitTime();

	  auto isDarkHit = mcPhoton->IsDarkHit();
	  if(isDarkHit)
		cout << "DARK HIT" << endl;

	  if(isSource){

	    auto trackID = mcPhoton->GetTrackID();

	    for(auto p:vFP){

	      if(trackID == p.GetTrackId()){

			// ########## //
			// FILL EVENT //
			// ########## //

			auto CreaProc = p.GetProcess();
			auto WL = p.GetWl();

			vHit.emplace_back(Hit(PMTPos, Q, T,
								  TVector3(0.,0.,0.),
								  TVector3(0.,0.,0.),
								  CreaProc, WL));

			break;

	      }

	    }

	  } else {

		// ########## //
		// FILL EVENT //
		// ########## //

		vHit.emplace_back(Hit(PMTPos, Q, T));

	  }

	}

  } // END FOR iPMT

  return vHit;

}

void GetPartInfoFromTrackStep(FlatParticle *fp, RAT::DS::MCTrack *mctrack){

  auto *mctrackstep = mctrack->GetMCTrackStep(0);

  fp->SetPos(mctrackstep->GetEndpoint());
  fp->SetKinE(mctrackstep->GetKE());
  fp->SetDir(mctrackstep->GetMomentum().Unit());

}

void GetPartIDFromTrackStep(G4Particle *p, RAT::DS::MCTrack *mctrack){

  p->SetParentId(mctrack->GetParentID());
  p->SetTrackId(mctrack->GetID());

}

vector<ComplexParticle> GetVPart(RAT::DS::MC * mc, const string& sPartName){

  auto NbTracks = mc->GetMCTrackCount();

  vector<ComplexParticle> vPart;

  // First, identify how many particle with name sPartName are created in the evt.
  // Loop over tracks and records the one called sPartName.
  // Save only information of interest.

  for(auto iTrack=0; iTrack<NbTracks; iTrack++) {

	auto *mctrack = mc->GetMCTrack(iTrack);

	if (IsParticle(mctrack, sPartName)) {

	  ComplexParticle cPart;
	  cPart.SetName(mctrack->GetName());
	  GetPartIDFromTrackStep(&cPart, mctrack);
	  GetPartInfoFromTrackStep(&cPart, mctrack);

	  vPart.emplace_back(cPart);

	}

  }

  return vPart;

}

vector<FlatPhoton> GetAndFillPhotonsToVPart(RAT::DS::MC * mc, vector<ComplexParticle> *vPart){

  vector<FlatPhoton> vFP;

  auto NbTracks = mc->GetMCTrackCount();

  for(auto iTrack=0; iTrack<NbTracks; iTrack++) {

	auto * mctrack = mc->GetMCTrack(iTrack);

	if(IsParticle(mctrack, "opticalphoton")) {

	  FlatPhoton fp;
	  GetPartIDFromTrackStep(&fp, mctrack);
	  fp.SetWl(MeV2lambda(mctrack->GetMCTrackStep(0)->GetKE()));

	  for(auto &p:*vPart){

		if(fp.GetParentId() == p.GetTrackId()){

		  p.AddPhoton(fp);
		  vFP.emplace_back(fp);

		} // END if daughter

	  } // ENF for vPart

	} // END if photon

  } // END for iTrack

  return vFP;

}

vector<Hit> GetVHitsFromPart(Analyzer *fAnalyzer, unsigned int iEvt,
							 string sPartName,
							 vector<ComplexParticle> *vCPart){

  RAT::DS::MC * mc = GetRATMCOnEvt(fAnalyzer, iEvt);

  vector<ComplexParticle> vPart = GetVPart(mc, sPartName);
  vector<FlatPhoton> vFP = GetAndFillPhotonsToVPart(mc, &vPart);

  // cout << "#DEBUGGING " << vPart.size() << " protons." << endl;
  //
  // for (auto proton = vPart.begin(); proton != vPart.end(); proton++) {
  //
	// proton->GetPos().Print();
  //   proton->GetDir().Print();
  //   cout << proton->GetKinE() << "MeV" << endl;
	// cout << proton->GetVp().size() << endl;
  //
  // }

	vector<Hit> vHit;

  for (int iPMT = 0; iPMT < mc->GetMCPMTCount(); iPMT++) {

	auto mcPMT = mc->GetMCPMT(iPMT);
	auto PMTID = mcPMT->GetID();
	auto PMTPos = fAnalyzer->GetRun()->GetPMTInfo()->GetPosition(PMTID);

	auto NbPhotonCounts = mcPMT->GetMCPhotonCount();

	for (int iP = 0; iP < NbPhotonCounts; iP++) {

	  auto mcPhoton = mcPMT->GetMCPhoton(iP);
	  // Get Q and T
	  auto Q = mcPhoton->GetCharge();
	  auto T = mcPhoton->GetHitTime();

	  // cout << mcPhoton->GetTrackID() << endl;

	  // ########## //
	  // FILL EVENT //
	  // ########## //

	  bool isPhotonProcessed = false;

	  for (auto proton = vPart.begin(); proton != vPart.end(); proton++) {

		for (auto photon = proton->GetVp().begin(); photon != proton->GetVp().end(); photon++) {

		  if (mcPhoton->GetTrackID() == photon->GetTrackId()) {

			// proton->GetPos().Print();
			// proton->GetDir().Print();
			// cout << proton->GetKinE() << "MeV" << endl;
			// cout << proton->GetVp().size() << endl;

			vHit.emplace_back(Hit(PMTPos, Q, T, proton->GetPos(), proton->GetDir()));

			isPhotonProcessed = true;
			proton->RemovePhoton(photon);
			break;

		  }

		  if (isPhotonProcessed)
			break;

		} // END for photon

	  } // END for proton

	} // END FOR iPhoton

  } // END for iPMT

  if(vCPart)
    *vCPart = vPart;

  return vHit;

}


bool IsParticle(RAT::DS::MCTrack *mctrack, const string& name){

  return mctrack->GetParticleName() == name;

}

void FilldEdX(ComplexParticle &CP, Analyzer *fAnalyzer, unsigned int iEvt) {

  auto *mc = GetRATMCOnEvt(fAnalyzer, iEvt);
  auto NbTracks = mc->GetMCTrackCount();

  for (auto iTrack = 0; iTrack < NbTracks; iTrack++) {
	auto *mctrack = mc->GetMCTrack(iTrack);

	if (CP.GetTrackId() == mctrack->GetID()) {

	  auto NbTrackSteps = mctrack->GetMCTrackStepCount();

	  for(auto iStep=1; iStep<NbTrackSteps; iStep++){

		auto *mctrackstep = mctrack->GetMCTrackStep(iStep);

		CP.ELoss(mctrackstep->GetKE());
		CP.AddTrackLength(abs(mctrackstep->GetLength()));

	  }

	}

  }

}

double SumdEdX(vector<ComplexParticle> vCP){

  double dEdX = 0.;
  for(auto CP:vCP){
	dEdX+=CP.GetDe()/CP.GetDx();
  }
  return dEdX;
}

vector<ComplexParticle> GetProtonLY(Analyzer *fAnalyzer, unsigned int iEvt){

  auto * mc = GetRATMCOnEvt(fAnalyzer, iEvt);
  auto NbTracks = mc->GetMCTrackCount();

  vector<ComplexParticle> vProtons;

  for(auto iTrack=0; iTrack<NbTracks; iTrack++){

	auto * mctrack = mc->GetMCTrack(iTrack);

	if(IsParticle(mctrack, "proton")){

	  ComplexParticle proton;
	  proton.SetName("proton");
	  GetPartInfoFromTrackStep(&proton, mctrack);
	  GetPartIDFromTrackStep(&proton, mctrack);

	  auto NbTrackSteps = mctrack->GetMCTrackStepCount();

	  for(auto iStep=1; iStep<NbTrackSteps; iStep++){

		auto *mctrackstep = mctrack->GetMCTrackStep(iStep);

		proton.ELoss(mctrackstep->GetKE());
		proton.AddTrackLength(abs(mctrackstep->GetLength()));

	  }

	  vProtons.emplace_back(proton);

	}

  }

  vector<FlatPhoton> vFP = GetAndFillPhotonsToVPart(mc, &vProtons);

  return vProtons;

}

void PrintTrackInfo(RAT::DS::MC *mc, unsigned int iTrack){

  auto * mctrack = mc->GetMCTrack(iTrack);
  auto NbTrackSteps = mctrack->GetMCTrackStepCount();

  cout << mctrack->GetParticleName() << endl;
  cout << " With ID: " << mctrack->GetID() << " and parent ID: " << mctrack->GetParentID() << endl;
  cout << "TOT Nb Steps: " << NbTrackSteps << endl;

  for(auto iStep=0; iStep<NbTrackSteps; iStep++){

	auto *mctrackstep = mctrack->GetMCTrackStep(iStep);

	cout << mctrackstep->GetLength() <<"mm " << mctrackstep->GetKE() << "MeV " << endl;

  }


}

void FillProtonRecoilSpectrum(TH2D* h2D, Analyzer *fAnalyzer, unsigned iEvt){

  auto * mc = GetRATMCOnEvt(fAnalyzer, iEvt);
  auto nTracks = mc->GetMCTrackCount();

  for(auto iTrack=0; iTrack<nTracks; iTrack++){

	auto * mctrack = mc->GetMCTrack(iTrack);

	if(IsParticle(mctrack, "proton")){

	  auto *mctrackstep = mctrack->GetMCTrackStep(0);

	  h2D->Fill(mc->GetMCParticle(0)->GetKE(),
				mctrackstep->GetKE());

	}

  }

}


void FillQuenchingSpectrum(TH2D *h2D, Analyzer *fAnalyzer, unsigned int iEvt){

  vector<ComplexParticle> vProtons = GetProtonLY(fAnalyzer, iEvt);

  for(auto p:vProtons){

	double dL = p.GetVp().size();
	double dX = abs(p.GetDx());
	double dE = abs(p.GetDe());

	if(h2D)
	  h2D->Fill(dL/dX, (dE/dX)/(1+dE/dX));

  }

}

void GetMCNumPEAndHits(double *NHits, double *NPE, Analyzer *fAnalyzer, unsigned int iEvt){


  RAT::DS::MC * mc = GetRATMCOnEvt(fAnalyzer, iEvt);

  *NHits=mc->GetMCPMTCount();
  *NPE=mc->GetNumPE();

}

