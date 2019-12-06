//
// Created by zsoldos on 11/26/19.
//

#ifndef _CREATEERSMATRIX_HH_
#define _CREATEERSMATRIX_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <fstream>
#include <sstream>
#include <regex>

/////////////////////////   BOOST   /////////////////////////
#include <boost/filesystem/path.hpp>

/////////////////////////   ROOT   //////////////////////////
#include <TH2D.h>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/MC.hh>
#include <RAT/DS/Root.hh>

void CollectPEAndHits(RAT::DS::MC *mc, TH2D *Hist){

  // Get Nb of PMTs which at least 1hit
  const int nbPE = mc->GetNumPE();
  const int nbPMTsHits = mc->GetMCPMTCount();

  if(Hist)
	Hist->Fill(nbPE, nbPMTsHits);

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

std::string ExtractFilenameFromPath(std::string pathname){

  boost::filesystem::path p(pathname);
  return p.filename().string();

}

std::vector<double> CorrectEbinsRange(std::vector<double> inputArray){

  std::vector<double> output;
  output.resize(inputArray.size()+1);

  double DeltaBin0 = inputArray[1] - inputArray[0];

  output[0] = inputArray[0] - DeltaBin0/2;

  for(int i = 1; i<inputArray.size(); i++){

    output[i] = inputArray[i-1] + (inputArray[i] - inputArray[i-1])/2;

  }

  output[inputArray.size()] =
	  inputArray[inputArray.size()-1] +
		  (inputArray[inputArray.size()-1] - inputArray[inputArray.size()-2])/2;

  return output;

}


#endif // _CREATEERSMATRIX_HH_

