#include <algorithm>
#include <iostream>
#include <vector>

#include <boost/filesystem/path.hpp>

using namespace std;

void SetStyleVariables(TStyle *t2kStyle);

double gen(){
  
  static double i = 1;
  return i++/10;
  
}

void LoadLibs(){

  gSystem->Load("/usr/lib/libboost_filesystem.so");

}

void EffVSCov(const char *f0pct=NULL,
	      const char *f1pct=NULL,
	      const char *f5pct=NULL,
	      const char *f10pct=NULL){

  TStyle *t2kstyle = new TStyle("T2K","T2K approved plots style");
  SetStyleVariables(t2kstyle);
  gROOT->SetStyle("T2K");

  vector<string> sFiles;
  sFiles.emplace_back(f0pct);
  sFiles.emplace_back(f1pct);
  sFiles.emplace_back(f5pct);
  sFiles.emplace_back(f10pct);

  const unsigned int nbFiles = sFiles.size();

  vector<string> Frac = {"0%", "1%", "5%", "10%"};
  const char *IDName = "Cov30pct";

  TFile *f = NULL;

  const int nEBins = 100;
  vector<double> v(nEBins);
  generate(v.begin(), v.end(), gen);

  TH1D *h = NULL;
  std::vector<TH1D*> hNHits;

  const int nThresh = 3;
  const int Thresh[nThresh] = {10, 15, 20};

  TGraphErrors *grEff[nbFiles][nThresh];

  unsigned int iFile = 0;
  for(auto& s : sFiles) {

    f = new TFile(s.c_str(),"READ");

    if(f){

      for(int iThresh = 0; iThresh<nThresh; iThresh++){
	grEff[iFile][iThresh] = new TGraphErrors();
	grEff[iFile][iThresh]->SetName(Form("%s_grEff%dFile%d",IDName, Thresh[iThresh], iFile));
	grEff[iFile][iThresh]->SetTitle(Form("WbLS %s; E_{kin} positron (MeV) ; Eff (%)",
					     Frac[iFile].c_str()));
	if(iFile == 0){
	  grEff[iFile][iThresh]->SetTitle("H2O; E_{kin} positron (MeV) ; Eff (%)");
	}
	
	grEff[iFile][iThresh]->SetLineWidth(2.);
	
      }
    
      for (auto iv: v) {

	h = (TH1D*)f->Get(Form("hEbin%.1f_py",iv));
	
	if(h){
	  cout << "Open hist for EBin: " << iv << endl;
	  
	  if(h->GetEntries()>0){
	    
	    h->SetTitle(Form("%.1f MeV",iv));
	    hNHits.emplace_back((TH1D*)h->Clone());
	    
	    for(int iThresh = 0; iThresh<nThresh; iThresh++){
	      
	      grEff[iFile][iThresh]->SetPoint(grEff[iFile][iThresh]->GetN(),
				       iv,
				       h->Integral(h->FindBin(Thresh[iThresh]),
						   h->GetNbinsX())/h->Integral());
	      grEff[iFile][iThresh]->SetPointError(grEff[iFile][iThresh]->GetN()-1, 0., 0);
	      
	      
	    } // END FOR iThresh

	    if(iv == 2.5) h->Fit("gaus","0");
	    
	  } // END if h->GetEntries()>0
      
	} // END if h
	
      } // END for iv

      iFile++;
      
    } // END if f

  } // END FOR f
  

  TCanvas *c1;
  TMultiGraph *mg[nThresh];
  for(int iThresh = 0; iThresh<nThresh; iThresh++){

    mg[iThresh] = new TMultiGraph();
    mg[iThresh]->SetTitle(Form("Positron detection efficiency for Thresh=%d NHits ; E_{kin} e^{+} (MeV); Eff (%)",Thresh[iThresh]));

    for(int iFile = 0; iFile<nbFiles; iFile++){
      mg[iThresh]->Add(grEff[iFile][iThresh],"pl");
    }

    c1 = new TCanvas(Form("c%s_Eff_%d",IDName,iThresh),Form("c%s_Eff_%d",IDName,iThresh),800,600);
    mg[iThresh]->Draw("A PMC PLC");
    mg[iThresh]->GetYaxis()->SetRangeUser(0,1.1);
    c1->BuildLegend();
    c1->Print(Form("%s_Eff_%dHits.pdf",IDName,Thresh[iThresh]));
    
  }

    // c1 = new TCanvas(Form("c%s_NHits",IDName),Form("c%s_NHits",IDName), 800,600);
  // for(auto& h : hNHits) {
  //   if(h->GetEntries() > 0){
  //     h->Draw("SAME PLC PMC");
  //   }
  // }
  // // gStyle->SetOptTitle(0);
  // // c1->BuildLegend();

}

void SetStyleVariables(TStyle *t2kStyle){
    
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
