///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <vector>

/////////////////////////   BOOST   /////////////////////////

/////////////////////////   ROOT   //////////////////////////

/////////////////////////   RAT   ///////////////////////////

/////////////////////////   USER   ///////////////////////////
#include "utils.hh"

#include "Analyzer.hh"
#include "AnalyzerFunctions.hh"
#include "HitClass.hh"
#include "HitFunctions.hh"
#include "MCFunctions.hh"
#include "EVFunctions.hh"
#include "LL.hh"
#include "FlatParticle.hh"

#include "ProgressBar.hpp"

#include "cnpy.h"

using namespace std;

void TestMacro(int User_NEvts = INT_MIN, int User_iEvt = INT_MIN){

  string sPathToFile = "/home/zsoldos/theia/MC/Theia_1kT_Cov30pct_e-_2.6MeV_Pos_0_0_0_Dir_0_0_1/";
  string sFilename = sPathToFile + "Theia_1kT_Cov30pct_wbls_5pct_WM_0420_e-_2.6MeV_Pos_0_0_0_Dir_0_0_0.root";

  auto *FileAnalyzer = new Analyzer(sFilename.c_str());

  unsigned long NEvts = User_NEvts > 0 ? User_NEvts : FileAnalyzer->GetNEvts();
  auto iEvt = SetDefValue(User_iEvt, 0);
  NEvts = User_iEvt > 0 ? User_iEvt + NEvts : NEvts;
  ProgressBar progressBar(NEvts, 70);

  string NPZName = sPathToFile + "Theia_1kT_Cov30pct_wbls_5pct_WM_0420_e-_2.6MeV_Pos_0_0_0_Dir_0_0_0.npz";
  GetVHitAndDumpFlatNPZ(FileAnalyzer, iEvt, NPZName, "w");
  iEvt++;

  for(iEvt; iEvt<NEvts; iEvt++) {

    // record the tick
    ++progressBar;

    GetVHitAndDumpFlatNPZ(FileAnalyzer, iEvt, NPZName);

    // display the bar
    progressBar.display();

  }

}
