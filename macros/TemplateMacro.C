///////////////////////// STL C/C++ /////////////////////////

/////////////////////////   BOOST   /////////////////////////

/////////////////////////   ROOT   //////////////////////////
#include <TVector3.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TFile.h>

/////////////////////////   RAT   ///////////////////////////

/////////////////////////   USER   ///////////////////////////
// PUT UTILS.HH FIRST if you want to use ROOT CINT
// Because other headers might use function in utils.hh
#include "utils.hh"

#include "Analyzer.hh"
#include "HitClass.hh"

#include "EVFunctions.hh"
#include "HitFunctions.hh"
#include "LL.hh"
#include "CalibFunctions.hh"
#include "MCFunctions.hh"

#include "ProgressBar.hpp"

#define SQRT5 2.2360679775
#define SQRT2 1.41421356237

using namespace std;

void TestMacro(){

  loadlibs();

}
