//
// Created by zsoldos on 11/26/19.
//

#ifndef _CREATEERSMATRIX_HH_
#define _CREATEERSMATRIX_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <fstream>
#include <sstream>
#include <regex>
#include <string>
#include <vector>

/////////////////////////   BOOST   /////////////////////////

/////////////////////////   ROOT   //////////////////////////
#include <TApplication.h>

/////////////////////////   RAT   ///////////////////////////


void ShowUsage(string name){

  cerr << "Usage: " << name << " <option(s)> SOURCES" << endl
	   << "Options:\n"
	   << "\t-h\tShow this help message\n"
	   << "\t-nPEBins\tSet nPE bins (default 201)\n"
	   << "\t-minPE\tSet low edge bin PE (default -0.5)\n"
	   << "\t-maxPE\tSet max edge bin PE (default 200.5)\n"
	   << "\t-nPMTBins\tSet nPMT bins (default 201)\n"
	   << "\t-minPMT\tSet low edge bin PMT (default -0.5)\n"
	   << "\t-maxPMT\tSet max edge bin PMT (default 200.5)\n"
	   << "\t-nThresh\tSet threshold for Eff (default nHits=20)\n"
	   << "\t-foutputname\tSet name output file (default \"output.root\")\n"
	   << endl
	   << "\tSOURCES\tSpecify input data file (.txt)\n"
	   << endl;

}

void ProcessArgs(TApplication *theApp, string *filename,
				 int *User_nPEBins, double *User_minPE, double *User_maxPE,
				 int *User_nPMTBins, double *User_minPMT, double *User_maxPMT,
				 int *User_nThresh,
				 string *User_fOutput) {

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
	} else if ((arg == "-nPEBins")) {
	  *User_nPEBins = stoi(theApp->Argv(++i));
	} else if ((arg == "-minPE")) {
	  *User_minPE = stod(theApp->Argv(++i));
	} else if ((arg == "-maxPE")) {
	  *User_maxPE = stod(theApp->Argv(++i));
	} else if ((arg == "-nPMTBins")) {
	  *User_nPMTBins = stoi(theApp->Argv(++i));
	} else if ((arg == "-minPMT")) {
	  *User_minPMT = stod(theApp->Argv(++i));
	} else if ((arg == "-maxPMT")) {
	  *User_maxPMT = stod(theApp->Argv(++i));
	} else if ((arg == "-nThresh")) {
	  *User_nThresh = stod(theApp->Argv(++i));
	} else if ((arg == "-foutputname")) {
	  *User_fOutput = theApp->Argv(++i);
	} else {
	  if (i + 1 > theApp->Argc() && nFiles == 0) {
		cout << "NO SOURCES PROVIDED !" << endl;
		ShowUsage(theApp->Argv(0));
		exit(EXIT_FAILURE);
	  } else {
		cout << "READ " << arg << endl;
		nFiles++;
		*filename = arg;
	  }
	}
  }

}

double ExtractEBinFromFilename(const std::string& filename){

  std::vector <std::string> tokens;
  std::stringstream sname(filename);
  std::string intermediate;

  // Tokenizing w.r.t. space ' '
  while(getline(sname, intermediate, '_')){
	tokens.push_back(intermediate);
  }

  // Extract float from token using regex
  std::smatch m;
  std::regex e (R"(\d+\.\d+)");   // matches floats inside filename
  std::regex_search (tokens[tokens.size()-3],m,e);
  return std::stod(m.str());

}



#endif // _CREATEERSMATRIX_HH_

