//
// Created by zsoldos on 2/21/20.
//

#ifndef _CREATEPDF_HH_
#define _CREATEPDF_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

/////////////////////////   USER  ///////////////////////////
#include "HitClass.hh"

double GetTAvg(vector<Hit> vHit){
  TH1D *hdummy = new TH1D("hdummy", "dummy", 1e6, 0., 1e6);

  for(auto h:vHit){
    hdummy->Fill(h.GetT());
  }

  double TAvg = hdummy->GetMean();
  delete hdummy;
  return TAvg;
}

void FillResiduals(TH1D *hResid, vector<Hit> vHit,
				   TVector3 Origin = TVector3(0.,0.,0.),
				   double SoL = 224.9,
				   double TAvg = 0.,
				   bool isWeight=true){

  for(auto h: vHit){
    if(isWeight){
	  hResid->Fill(h.CalculateTResid(Origin, SoL) - TAvg,
				   1/(double)(vHit.size()));
    } else{
	  hResid->Fill(h.CalculateTResid(Origin, SoL) - TAvg);
    }
  }

}




void ShowUsage(string name){

  cerr << "Usage: " << name << " <option(s)> -mc FILE.root" << endl
	   << "Options:\n"

	   << "\t-h\tShow this help message\n"

	   << "\t-b\tSet batch mode\n"

	   << "\t-nevts\tNEvts to process (int)\n"
	   << "\t-ievt\tStart at Evt #i (int)\n"

	   << "\t-tbins-prompt\tSet nb bins for tresid hist (int)\n"
	   << "\t-mint-prompt\tSet min val for tresid hist (double)\n"
	   << "\t-maxt-prompt\tSet max val for tresid hist (double)\n"

	   << "\t-tbins-delay\tSet nb bins for tresid hist (int)\n"
	   << "\t-mint-delay\tSet min val for tresid hist (double)\n"
	   << "\t-maxt-delay\tSet max val for tresid hist (double)\n"

	   << "\t-xx\tSet X for true origin vector (double)\n"
	   << "\t-yy\tSet X for true origin vector (double)\n"
	   << "\t-zz\tSet X for true origin vector (double)\n"

	   << "\t-mc\tinput file (ROOT)\n"

	   << endl;

}


void ProcessArgs(TApplication *theApp, string *filename,
				 int *User_PromptCut,
				 int *User_nEvts, int *User_iEvt,
				 int *User_nTResidBins_Prompt, double *User_minTResid_Prompt, double *User_maxTResid_Prompt,
				 int *User_nTResidBins_Delay, double *User_minTResid_Delay, double *User_maxTResid_Delay,
				 double *User_xx, double *User_yy, double *User_zz,
				 bool *User_isBatch) {

  // Reading user input parameters
  if (theApp->Argc() < 2) {
	ShowUsage(theApp->Argv(0));
	exit(0);
  }

  int nFiles=0;

  for (int i = 1; i < theApp->Argc(); i++) {
	string arg = theApp->Argv(i);
	if ((arg == "-h") || (arg == "--help")) {
	  ShowUsage(theApp->Argv(0));
	  exit(0);

	} else if ((arg == "-pc")) {
	  *User_PromptCut = stoi(theApp->Argv(++i));

	} else if ((arg == "-nevts")) {
	  *User_nEvts = stoi(theApp->Argv(++i));
	} else if ((arg == "-ievt")) {
	  *User_iEvt = stoi(theApp->Argv(++i));

	} else if ((arg == "-tbins-prompt")) {
	  *User_nTResidBins_Prompt = stoi(theApp->Argv(++i));
	} else if ((arg == "-mint-prompt")) {
	  *User_minTResid_Prompt = stod(theApp->Argv(++i));
	} else if ((arg == "-maxt-prompt")) {
	  *User_maxTResid_Prompt = stod(theApp->Argv(++i));

	} else if ((arg == "-tbins-delay")) {
	  *User_nTResidBins_Delay = stoi(theApp->Argv(++i));
	} else if ((arg == "-mint-delay")) {
	  *User_minTResid_Delay = stod(theApp->Argv(++i));
	} else if ((arg == "-maxt-delay")) {
	  *User_maxTResid_Delay = stod(theApp->Argv(++i));

	} else if ((arg == "-xx")) {
	  *User_xx = stod(theApp->Argv(++i));
	} else if ((arg == "-yy")) {
	  *User_yy = stod(theApp->Argv(++i));
	} else if ((arg == "-zz")) {
	  *User_zz = stod(theApp->Argv(++i));

	} else if ((arg == "-b")) {
	  *User_isBatch=true;

	} else if ((arg == "-mc")) {
	  *filename = theApp->Argv(++i);

	} else {
	  cout << "Unkown parameter" << endl;
	  continue;
	}
  }

  if(filename->empty()){
	cout << "ERROR: No MC input file provided!" << endl;
	exit(EXIT_FAILURE);
  // } else if(!IsFileExist(*filename)){
	// cout << "ERROR: MC file doesn't exist!!" << endl;
	// exit(EXIT_FAILURE);
  }

}

#endif //_CREATEPDF_HH_
