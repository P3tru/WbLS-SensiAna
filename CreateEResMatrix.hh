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
  const int nbPMTsHits = mc->GetMCPMTCount();
  const int nbPE = mc->GetNumPE();

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

#endif // _CREATEERSMATRIX_HH_

