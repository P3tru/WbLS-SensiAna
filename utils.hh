//
// Created by zsoldos on 11/26/19.
//

#ifndef _UTILS_HH_
#define _UTILS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <string>

/////////////////////////   BOOST   /////////////////////////
#include <boost/filesystem/path.hpp>

/////////////////////////   ROOT   //////////////////////////

/////////////////////////   RAT   ///////////////////////////

std::string ExtractFilenameFromPath(std::string pathname){

  boost::filesystem::path p(pathname);
  return p.filename().string();

}

std::vector<double> CorrectBinRangeArray(std::vector<double> inputArray){

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


#endif // _UTILS_HH_
