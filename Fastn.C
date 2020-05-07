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
#include "HitClass.hh"
#include "LL.hh"
#include "Analyzer.hh"
#include "HitClass.hh"
#include "EVFunctions.hh"
#include "HitFunctions.hh"
#include "ProgressBar.hpp"
#include "utils.hh"

#define SQRT5 2.2360679775
#define SQRT2 1.41421356237

using namespace std;

R__LOAD_LIBRARY(libEVFunctions.so)
R__LOAD_LIBRARY(libHitFunctions.so)
R__LOAD_LIBRARY(libLL.so)

void TestMacro(int User_NEvts = -1, int User_iEvt = -1){

  string sPathToFile = "../rat-pac/RAW_OUTPUTS/";
  string sGeom = "Theia_1kT_Cov50perCent";
  string sMaterial = "wbls_10newpct";
  string sParticle = "neutrons";
  string sE = "10_100MeV";
  string sInit = "Hor0mm_Ver0mm";

  //../MC/Theia_1kT_Cov30pct_wbls_1pct_e-_Pos_0_0_0_Dir_1_0_0
  string sRootPath = sPathToFile + sGeom + "_" + sMaterial + "_" + sParticle + "_" + sInit + "/";
  sRootPath = sPathToFile;
  string sFilename = sRootPath + sGeom+ "_" + sMaterial + "_" + sParticle + "_" + sE + "_" + sInit + ".root";

  auto *FileAnalyzer = new Analyzer(sFilename.c_str());

  unsigned long int NEvts = User_NEvts > 0 ? User_NEvts : FileAnalyzer->GetNEvts();
  auto iEvt = SetDefValue(User_iEvt, 0);
  NEvts = User_iEvt > 0 ? User_iEvt + NEvts : NEvts;
  ProgressBar progressBar(NEvts, 70);

  string NPZName = sRootPath + sGeom+ "_" + sMaterial + "_" + sParticle + "_" + sE + "_" + sInit + ".npz";
  GetVHitAndDumpFlatNPZ(FileAnalyzer, iEvt, NPZName, "w");
  iEvt++;

  for(iEvt; iEvt<NEvts; iEvt++) {

    // record the tick
    ++progressBar;

    // display the bar
    progressBar.display();

  }

}
