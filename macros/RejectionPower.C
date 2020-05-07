///////////////////////// STL C/C++ /////////////////////////
#include <vector>
#include <climits>
#include <numeric>
#include <algorithm>

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

using namespace std;

pair<double, double> GetMaxRejection(TGraph *gr1, TGraph *gr2, vector<double> vPtX){

  pair<double, double> pMax = make_pair(DBL_MIN, DBL_MIN);
  
  for(auto X:vPtX){

    pMax = abs(gr1->Eval(X) - gr2->Eval(X)) > pMax.second ? make_pair(X, abs(gr1->Eval(X) - gr2->Eval(X))) : pMax;
    
  }

  return pMax;

}

void TestMacro(){

  TFile *_file0 = TFile::Open("Theia_1kT_Cov30pct_wbls_1pct_e+_1MeV_Pos_0_0_0_Dir_0_0_1_plots.root");
  TGraph *grPositron = (TGraph*)_file0->Get("gre+");
  grPositron->SetLineColor(kBlue-4);
  grPositron->SetLineWidth(2);
  TFile *_file1 = TFile::Open("Theia_1kT_Cov30pct_wbls_1pct_neutron_1MeV_Pos_0_0_0_Dir_0_0_1_plots.root");
  TGraph *grNeutron = (TGraph*)_file1->Get("grneutron");
  grNeutron->SetLineColor(kRed-4);
  grNeutron->SetLineWidth(2);
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("PID e+/neutron ; #chi^{2}/NdF ; Ratio of Integral");
  mg->Add(grPositron, "lp");
  mg->Add(grNeutron, "lp");
  
  vector<double> vBinsX(1000);
  iota(vBinsX.begin(), vBinsX.end(), 1);

  
  pair<double, double> pMax = GetMaxRejection(grPositron, grNeutron, vBinsX);

  cout << "MAX Positron: " << 1-grPositron->Eval(pMax.first) << endl;
  cout << "MAX Neutron: " << 1-grNeutron->Eval(pMax.first) << endl;

  TLine *lMax = new TLine(pMax.first, grPositron->Eval(pMax.first),
			  pMax.first, grNeutron->Eval(pMax.first));
    
  lMax->SetLineStyle(2);
  lMax->SetLineWidth(2);
  
  auto *c1 = new TCanvas("c1","c1", 800, 600);
  mg->Draw("A");
  c1->BuildLegend();
  lMax->Draw("SAME");
  
}
