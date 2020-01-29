#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
void SetStyleVariables(TStyle *t2kStyle);

void Neutrons(const char *f0pct,
	      const char *f1pct,
	      const char *f5pct,
	      const char *f10pct,
	      const char *fdoped,
	      const char *fdoped1pct=NULL){


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

  const unsigned int nFrac = 4;
  const double wblsFrac[nFrac] = {0., 0.01, 0.05, 0.10};
  const double dwblsFrac[nFrac] = {0., 0., 0., 0.};

  const int nThresh = 3;
  const int Thresh[nThresh] = {10, 15, 20};

  TFile *f = NULL;

  TH1D *h = NULL;
  std::vector<TH1D*> hNHits;

  TGraphErrors *grEff[nThresh];
  TGraphErrors *grDopedEff[nThresh];

  for(int iThresh = 0; iThresh<nThresh; iThresh++){
    
    grEff[iThresh] = new TGraphErrors();
    grEff[iThresh]->SetName(Form("grThresh%d",
				 Thresh[iThresh]));
    grEff[iThresh]->SetTitle(Form("nThresh = %d; WbLS (%) ; Eff (%)",
				  Thresh[iThresh]));

    grEff[iThresh]->SetMarkerStyle(20);

    grDopedEff[iThresh] = new TGraphErrors();
    grDopedEff[iThresh]->SetName(Form("grDopedThresh%d",
				 Thresh[iThresh]));
    grDopedEff[iThresh]->SetTitle(Form("H2O+Gd nThresh = %d; WbLS (%) ; Eff (%)",
				  Thresh[iThresh]));

    grDopedEff[iThresh]->SetMarkerStyle(23);
    
  }

  f = new TFile(fdoped,"READ");
  h = (TH1D*)f->Get("hNbPEVSHits_py");
  h->SetTitle("H2O+Gd");
  h->Scale(1/h->GetMaximum());
  
  hNHits.emplace_back((TH1D*)h->Clone());
  for(int iThresh = 0; iThresh<nThresh; iThresh++){

    grDopedEff[iThresh]->SetPoint(grDopedEff[iThresh]->GetN(),
				  0,
				  h->Integral(h->FindBin(Thresh[iThresh]),
					      h->GetNbinsX())/h->Integral());
    grDopedEff[iThresh]->SetPointError(grDopedEff[iThresh]->GetN()-1, 0., 0.);

  } // END for iThresh
	      
  
  if(fdoped1pct){

    f = new TFile(fdoped1pct,"READ");
    h = (TH1D*)f->Get("hNbPEVSHits_py");
    h->SetTitle("WbLS 1% + Gd");
    h->Scale(1/h->GetMaximum());
  
    hNHits.emplace_back((TH1D*)h->Clone());
    for(int iThresh = 0; iThresh<nThresh; iThresh++){

      grDopedEff[iThresh]->SetPoint(grDopedEff[iThresh]->GetN(),
				    0.01,
				    h->Integral(h->FindBin(Thresh[iThresh]),
						h->GetNbinsX())/h->Integral());
      grDopedEff[iThresh]->SetPointError(grDopedEff[iThresh]->GetN()-1, 0., 0.);

    } // END for iThresh

  }
  
  TLegend *leg = new TLegend(0.7,0.7,0.99,0.99);

  unsigned int iFile = 0;

  for(auto& s : sFiles) {

    f = new TFile(s.c_str(),"READ");

    if(f){

      cout << "File " << s.c_str() << " OPEN " << endl;

      h = (TH1D*)f->Get("hNbPEVSHits_py");
      h->SetTitle(Form("WbLS %s",Frac[iFile].c_str()));
	if(iFile == 0){
	  h->SetTitle("H2O");
	}
      h->Scale(1/h->GetMaximum());
      hNHits.emplace_back((TH1D*)h->Clone());

      for(int iThresh = 0; iThresh<nThresh; iThresh++){

	cout << "#WbLSFrac: " << wblsFrac[iFile]
	     << " Eff: " << h->Integral(h->FindBin(Thresh[iThresh]),h->GetNbinsX())/h->Integral()
	     << endl;

      	grEff[iThresh]->SetPoint(grEff[iThresh]->GetN(),
      				 wblsFrac[iFile],
      				 h->Integral(h->FindBin(Thresh[iThresh]),
      					     h->GetNbinsX())/h->Integral());
      	grEff[iThresh]->SetPointError(grEff[iThresh]->GetN()-1, 0., 0.);

      } // END for iThresh

      iFile++;

    } // END if f

  } // END iFile

  auto mg = new TMultiGraph();
  mg->SetTitle("Neutrons detection efficiency for 1kT 30% Coverage; WbLS (%) ; Eff (%)");
  for(int iThresh = 0; iThresh<nThresh; iThresh++){
    mg->Add(grEff[iThresh],"pl");
  }
  for(int iThresh = 0; iThresh<nThresh; iThresh++){
    mg->Add(grDopedEff[iThresh],"pl");
  }
  
  const char *IDName = "Cov30pct";
  auto c1 = new TCanvas("cEff","cEff", 800,600);
  mg->Draw("A PMC PLC");
  mg->GetYaxis()->SetRangeUser(0.,1.01);
  c1->BuildLegend();
  c1->Print(Form("%s_Eff.pdf",IDName));

  // for(unsigned int i=0; i<5; i++){
  //   leg->AddEntry(hNHits[i],Form("WbLS %s",Frac[i].c_str()));
  // }

  c1 = new TCanvas("c1","c1", 800,600);
  for(auto& h : hNHits) {
    if(h->GetEntries() > 0){
      h->Draw("SAME HIST PLC PMC");
    }
  }
  gStyle->SetOptTitle(0);
  c1->BuildLegend();
  
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

