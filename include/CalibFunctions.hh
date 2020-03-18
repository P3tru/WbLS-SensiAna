//
// Created by zsoldos on 2/24/20.
//

#ifndef _MCCALIB_HH_
#define _MCCALIB_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

/////////////////////////   ROOT  ///////////////////////////
#include <TGraph.h>
#include <TFile.h>
#include <TROOT.h>
#include <TKey.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>

using namespace std;

void FitAndRecoverGaussParams(TH1D *h, double *mu, double *sigma);
void FillPair(TH1D *h, double EBin,
			  vector< pair<double, double> > *vPMu, vector< pair<double, double> > *vPSig);
TGraph *CreateGraph(vector< pair<double, double> > *vP);
TGraphErrors *CreateGraphErrors(vector< pair<double, double> > *vP,
								vector< pair<double, double> > *vPErr);
vector<double> GetArray(vector< pair<double, double> > *vP);

class MCCalib{

 protected:

  char *filename;

  TGraph *grMuPE;
  TGraph *grSigPE;
  TGraph *grMuHits;
  TGraph *grSigHits;

 public:

  static double GenEBin(){

	static int iBin=1;
	static double EBin = 0.1;
	return EBin*iBin++;

  }


  explicit MCCalib(char *filename) : filename(filename) {

	auto *FileCalib = TFile::Open(filename);

	vector<double> EBins(100);
	generate(EBins.begin(), EBins.end(), GenEBin);

	vector< pair<double, double> > vMuPE;
	vector< pair<double, double> > vSigPE;

	vector< pair<double, double> > vMuHits;
	vector< pair<double, double> > vSigHits;

	TIter next(FileCalib->GetListOfKeys());
	TKey *key;
	while ((key = (TKey*)next())) {

	  TClass *cl = gROOT->GetClass(key->GetClassName());
	  if (!cl->InheritsFrom("TH1")) continue;
	  // X-axis: nPE
	  // Y-axis: nHits
	  TH1 *h = (TH1*)key->ReadObj();

	  for(auto EBin:EBins){

		if(strcmp(Form("hEbin%.1f", EBin),h->GetName()) == 0){

		  TH1D *hNPE = (TH1D*)FileCalib->Get(Form("hEbin%.1f_px", EBin));
		  FillPair(hNPE, EBin, &vMuPE, &vSigPE);

		  TH1D *hNHits = (TH1D*)FileCalib->Get(Form("hEbin%.1f_py", EBin));
		  FillPair(hNHits, EBin, &vMuHits, &vSigHits);

		  break;

		} // END if strcmp

	  } // END for EBin

	} // END while Key

	grMuPE=new TGraph();
	grSigPE=new TGraph();
	grMuHits=new TGraph();
	grSigHits=new TGraph();

	if(grMuPE){
	  grMuPE = CreateGraph(&vMuPE);
	}
	if(grSigPE){
	  grSigPE = CreateGraph(&vSigPE);
	}

	if(grMuHits){
	  grMuHits = CreateGraph(&vMuHits);
	}
	if(grSigHits){
	  grSigHits = CreateGraph(&vSigHits);
	}

	grMuPE->SetBit(TGraph::kIsSortedX);
	grSigPE->SetBit(TGraph::kIsSortedX);

	grMuHits->SetBit(TGraph::kIsSortedX);
	grSigHits->SetBit(TGraph::kIsSortedX);

  }

  TGraph *GetGrMuPe() const {
	return grMuPE;
  }
  TGraph *GetGrSigPe() const {
	return grSigPE;
  }
  TGraph *GetGrMuHits() const {
	return grMuHits;
  }
  TGraph *GetGrSigHits() const {
	return grSigHits;
  }

  char *GetFilename() const {
	return filename;
  }

};

void FitAndRecoverGaussParams(TH1D *h, double *mu, double *sigma){

  TF1 *fFit;
  h->Fit("gaus", "Q0");
  fFit = h->GetFunction("gaus");
  double mean = fFit->GetParameter(1);
  double meanErr = fFit->GetParError(1);
  double sig = fFit->GetParameter(2);
  double sigErr = fFit->GetParError(2);

  *mu=mean;
  *sigma=sig;


}

void FillPair(TH1D *h, double EBin,
			  vector< pair<double, double> > *vPMu, vector< pair<double, double> > *vPSig){

  double mu, sig;
  FitAndRecoverGaussParams(h, &mu, &sig);
  vPMu->push_back( make_pair(EBin, mu) );
  vPSig->push_back( make_pair(EBin, sig) );

}

TGraph *CreateGraph(vector< pair<double, double> > *vP){

  sort(vP->begin(), vP->end());

  auto *gr = new TGraph();

  for(auto p: *vP){

	gr->SetPoint(gr->GetN(), p.first, p.second);

  }

  return gr;

}

TGraphErrors *CreateGraphErrors(vector< pair<double, double> > *vP,
								vector< pair<double, double> > *vPErr){

  sort(vP->begin(), vP->end());
  sort(vPErr->begin(), vPErr->end());

  auto *gr = new TGraphErrors();

  for(auto p: *vP){

	gr->SetPoint(gr->GetN(), p.first, p.second);

  }
  vector<double> sig = GetArray(vPErr);

  for(auto i=0; i<sig.size(); i++){

	gr->SetPointError(i, 0., sig[i]);

  }

  return gr;

}

vector<double> GetArray(vector< pair<double, double> > *vP){

  sort(vP->begin(), vP->end());

  vector<double> vD;

  for(auto p: *vP){

	vD.push_back(p.second);

  }

  return vD;

}


double ComputeLikelihood(MCCalib CalibObj, double NPE, double NHits);

#endif