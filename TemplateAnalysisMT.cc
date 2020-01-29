//
// Created by zsoldos on 12/5/19.
//

///////////////////////// STL C/C++ /////////////////////////
#include <thread>
#include <csignal>
#include <fstream>

/////////////////////////   ROOT   //////////////////////////
#include <TApplication.h>
#include <TH2D.h>

/////////////////////////   USER   //////////////////////////
#include "TFileAnalysis.hh"
#include "CreateEResMatrix.hh"
#include "TemplateAnalysisMT.hh"
#include "utils.hh"

using namespace std;

int main(int argc, char *argv[]) {

  EoF=0;
  signal(SIGINT,Interrupt);

  const unsigned int CORES = std::thread::hardware_concurrency();

  // Create TApp
  TApplication theApp("App", &argc, argv);

  // Input parameters
  // There's default value chosen by ternaries below
  string inputName;
  int User_nPEBins=-1;
  double User_minPE=-1; double User_maxPE=-1;
  int User_nPMTBins=-1;
  double User_minPMT=-1; double User_maxPMT=-1;
  string User_fOutput;
  int User_nThresh=-1;

  // READ input parameters
  ProcessArgs(&theApp,&inputName,
			  &User_nPEBins, &User_minPE, &User_maxPE,
			  &User_nPMTBins, &User_minPMT, &User_maxPMT,
			  &User_nThresh,
			  &User_fOutput);

  // Get ready to read inside inputName
  // each file name
  string line;
  ifstream file(inputName);

  // Create a FileAnalysis for each file
  // to be processed
  vector< TFileAnalysis<TH2D> > FileAnalysis;
  unsigned int NbFileToAnalyze = 0;

  while ( getline(file,line) ){

	if (line.compare(0,1,"#") == 0){

	  continue;

	} else {

	  if(IsFileExist(line)){

		FileAnalysis.emplace_back(line);
		double E = ExtractEBinFromFilename(line);
		cout << E << endl;
		FileAnalysis[FileAnalysis.size()-1].SetEBin(E);

		cout << "ADD " << ExtractFilenameFromPath(line) << endl;

	  } else {

		cerr << "FILE: " << ExtractFilenameFromPath(line) << " don't exist. Skip..." <<  endl;

	  }


	}

  } // END while

  NbFileToAnalyze = FileAnalysis.size();
  cout << "Nb of Files to process: " << NbFileToAnalyze << endl;

  TH2D *hNbPEVSHits[NbFileToAnalyze];
  const int nbBinsPE = (User_nPEBins>-1) ? User_nPEBins : 201;
  const double minPE = (User_minPE>-1) ? User_minPE : -0.5;
  const double maxPE = (User_maxPE>-1) ? User_maxPE : 200.5;

  const int nbBinsPMT = (User_nPMTBins>-1) ? User_nPMTBins : 201;
  const double minPMT = (User_minPMT>-1) ? User_minPMT : -0.5;
  const double maxPMT = (User_maxPMT>-1) ? User_maxPMT : 200.5;

  const int nThresh = (User_nThresh>-1) ? User_nThresh : 20;

  vector<thread> threads;

  for(unsigned int i=0; i < NbFileToAnalyze; i++) {

	cout << FileAnalysis[i].GetEBin() << endl;

	hNbPEVSHits[i] = new TH2D(Form("hEbin%.1f", FileAnalysis[i].GetEBin()),
							  "Nb PE Collected VS Nb PMTs Hits",
							  nbBinsPE, minPE, maxPE,
							  nbBinsPMT, minPMT, maxPMT);

//	threads.emplace_back(thread(ThreadAnalysis, FileAnalysis[i], hNbPEVSHits));
	threads.emplace_back(thread(TestThreadAnalysis, hNbPEVSHits[i]));

  }

  if(threads.size() > 0){
	for_each(threads.begin(), threads.end(), mem_fn(&std::thread::join));
  }

  /////////////////////////
  // ...

  cout << endl;
  cout << "Hit Ctrl+C to exit" << endl;
  theApp.Run(kTRUE);

  return 0;

}