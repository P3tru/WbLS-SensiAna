//
// Created by zsoldos on 11/26/19.
//

#ifndef _UTILS_HH_
#define _UTILS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h>

/////////////////////////   BOOST   /////////////////////////
#include <boost/filesystem/path.hpp>

/////////////////////////   ROOT   //////////////////////////
#include <TStyle.h>
#include <TVector3.h>
#include <TROOT.h>
#include <TMath.h>

/////////////////////////   RAT   ///////////////////////////

using namespace std;

static string ExtractFilenameFromPath(string pathname){

  boost::filesystem::path p(pathname);
  return p.filename().string();

}

static vector<double> CorrectBinRangeArray(vector<double> inputArray){

  vector<double> output;
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

int EoF = 0;

static void Interrupt(int arg){

  if(EoF==0) { printf("got a control-C, stop\n"); EoF=1; return; }
  else { printf("got a control-C, exiting\n"); exit(0); }

}

static void SetStyleVariables(TStyle *t2kStyle){

  t2kStyle->SetFrameBorderMode(0);
  t2kStyle->SetCanvasBorderMode(0);
  t2kStyle->SetPadBorderMode(0);
  t2kStyle->SetPadColor(0);
  t2kStyle->SetCanvasColor(0);
  t2kStyle->SetStatColor(0);
  t2kStyle->SetFillColor(0);
  t2kStyle->SetLegendBorderSize(1);

  t2kStyle->SetPaperSize(20,26);
  t2kStyle->SetPadTopMargin(0.15);
  t2kStyle->SetPadRightMargin(0.15); //0.05
  t2kStyle->SetPadBottomMargin(0.16);
  t2kStyle->SetPadLeftMargin(0.13);

  t2kStyle->SetTextFont(132);
  t2kStyle->SetTextSize(0.08);
  t2kStyle->SetLabelFont(132,"x");
  t2kStyle->SetLabelFont(132,"y");
  t2kStyle->SetLabelFont(132,"z");
  t2kStyle->SetLabelSize(0.05,"x");
  t2kStyle->SetTitleSize(0.06,"x");
  t2kStyle->SetLabelSize(0.05,"y");
  t2kStyle->SetTitleSize(0.06,"y");
  t2kStyle->SetLabelSize(0.05,"z");
  t2kStyle->SetTitleSize(0.06,"z");
  t2kStyle->SetLabelFont(132,"t");
  t2kStyle->SetTitleFont(132,"x");
  t2kStyle->SetTitleFont(132,"y");
  t2kStyle->SetTitleFont(132,"z");
  t2kStyle->SetTitleFont(132,"t");
  t2kStyle->SetTitleFillColor(0);
  t2kStyle->SetTitleX(0.05);
  t2kStyle->SetTitleFontSize(0.08);
  t2kStyle->SetTitleFont(132,"pad");

  // t2kStyle->SetPadGridX(true);
  // t2kStyle->SetPadGridY(true);


  t2kStyle->SetHistLineWidth(1.85);
  t2kStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes


  // t2kStyle->SetOptTitle(0);
  t2kStyle->SetOptStat(0);
  t2kStyle->SetOptFit(0);

  t2kStyle->SetPadTickX(1);
  t2kStyle->SetPadTickY(1);

  t2kStyle->SetPalette(1,0);  // use the nice red->blue palette
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
								   NCont);
  t2kStyle->SetNumberContours(NCont);

}

static void SetBasicStyle(){
  gROOT->SetStyle("Plain");

  gStyle->SetOptStat(0); // This determines if you want a stats box
  gStyle->SetOptFit(0); // This determines if you want a fit info box
  gStyle->GetAttDate()->SetTextColor(1);
  gStyle->SetOptTitle(1); // no title; comment out if you want a title
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetTextFont(132);
  gStyle->SetTitleFont(132,"XYZ");

  gROOT->ForceStyle();

  gStyle->SetPalette(kDarkRainBow);
}

static void SetBasicTH1Style(TH1 *h,
					  int Color=kBlue-4,
					  int LineWidth=1, int LineStyle=1,
					  int MarkerSize=1, int MarkerStyle=kPlus){

  h->SetLineColor(Color);
  h->SetLineWidth(LineWidth);
  h->SetLineStyle(LineStyle);

  h->SetMarkerColor(Color);
  h->SetMarkerSize(MarkerSize);
  h->SetMarkerStyle(MarkerStyle);

}

static bool IsFileExist(const string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

template <typename T>
T SetDefValue(const T& User, const T& def){
  return (User > 0) ? User : def;
}

// define a function with 3 parameters
static Double_t fitGaus(Double_t *x,Double_t *par) {

  Double_t arg = 0;
  if (par[2]!=0) arg = (x[0] - par[1])/par[2];
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;

}

// define a function with 3 parameters
static Double_t fitExpo(Double_t *x,Double_t *par) {

  return par[0]*TMath::Exp(-x[0]/par[1]);

}

void PrintPos(TVector3 Pos){
  cout << " X: " << Pos.x()
	   << " Y: " << Pos.y()
	   << " Z: " << Pos.z() << endl;

}

#endif // _UTILS_HH_
