#include "utils.hh"

TH1D *GetHPos(const char *filename, const char *histname){

  TFile *f = new TFile(filename);

  return (TH1D*)f->Get(histname);

}

// define a function with 3 parameters
Double_t fitGaus(Double_t *x,Double_t *par) {

  Double_t arg = 0;
  if (par[2]!=0) arg = (x[0] - par[1])/par[2];
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;

}

// define a function with 3 parameters
Double_t fitExpo(Double_t *x,Double_t *par) {
  
  if(x[0]>0)
    return par[0]*TMath::Exp(-x[0]/par[1]);
  else
    return 0;

}

// define a function with 3 parameters
Double_t fitNegExpo(Double_t *x,Double_t *par) {
  
  if(x[0]<0)
    return par[0]*TMath::Exp(x[0]/par[1]);
  else
    return 0;
  
}

Double_t fitFunc(Double_t *x, Double_t *par){

  return fitGaus(x, par) + fitExpo(x, &par[3]);
  
}


void GetFitResults(TH1D *h,
		   double *mean, double *sigma,
		   double *meanErr=NULL, double *sigmaErr=NULL,
		   bool isLargeDrift=false){

  TF1 *fitFcn;
  // create a TF1 with the range from -1000 to 1000 and 5 parameters
  if(isLargeDrift){
    fitFcn = new TF1("fitFcn",fitFunc,-1000,1000,5);
    fitFcn->SetParameter(0, h->GetMaximum());
    fitFcn->SetParameter(1, h->GetMean());
    fitFcn->SetParameter(2, h->GetRMS());
    fitFcn->SetParameter(3, h->GetMaximum());
    fitFcn->SetParameter(4, 200.);
      
    fitFcn->SetParLimits(0, 0.,1e5);
    fitFcn->SetParLimits(1, 0.,1000);
    fitFcn->SetParLimits(2, 0.,1000);
    fitFcn->SetParLimits(3, 0.,1e5);
    fitFcn->SetParLimits(4, 1.e-6,1000);
  } else {
    fitFcn = new TF1("fitFcn",fitGaus,-1000,1000,3);
    fitFcn->SetParameter(0, h->GetMaximum());
    fitFcn->SetParameter(1, h->GetMean());
    fitFcn->SetParameter(2, h->GetRMS());
    fitFcn->SetParLimits(0, 0.,1e5);
    fitFcn->SetParLimits(1, -1000.,1000);
    fitFcn->SetParLimits(2, -1000.,1000);
  }
  
  TFitResultPtr r = h->Fit("fitFcn", "LSEMR0+");
  TF1 *fFit;
  fFit = h->GetFunction("fitFcn");

  *mean = fFit->GetParameter(1);
  if(meanErr)
    *meanErr = fFit->GetParError(1);
  *sigma = fFit->GetParameter(2);
  if(sigmaErr)
    *sigmaErr = fFit->GetParError(2);

  delete fFit;

}

void SetGrPoint(TH1D *h, TGraph *gr, double WbLSFrac, bool isLargeDrift=false){
  double mean, sigma;
  GetFitResults(h, &mean, &sigma);
  gr->SetPoint(gr->GetN(), WbLSFrac, mean);

}

void SetGrErrPoint(TH1D *h, TGraphErrors *gr, double WbLSFrac, bool isLargeDrift=false){
  double mean, sigma;
  double meanErr, sigmaErr;
  GetFitResults(h, &mean, &sigma, &meanErr, &sigmaErr, isLargeDrift);
  gr->SetPoint(gr->GetN(), WbLSFrac, abs(mean));
  gr->SetPointError(gr->GetN()-1, 0., abs(meanErr));

}

void SetCondensedGrErrPoint(TH1D *h, TGraphErrors *gr, double WbLSFrac, bool isLargeDrift=false){

  double mean, sigma;
  double meanErr, sigmaErr;
  GetFitResults(h, &mean, &sigma, &meanErr, &sigmaErr, isLargeDrift);
  gr->SetPoint(gr->GetN(), WbLSFrac, mean);
  gr->SetPointError(gr->GetN()-1, 0., sigma);

}

void RebinHist(TH1D *h, const double n = 3){

  h->Rebin(n);
  
}

int GetPosHistMax(vector<TH1D*> vH){

  double min = DBL_MIN;
  int posmax=-1;

  for(auto ith = vH.begin(); ith != vH.end(); ith++){

    auto pos = distance(vH.begin(), ith);

    if((*ith)->GetMaximum()>min){

      posmax = pos;
      min = (*ith)->GetMaximum();
      
    }
    
  }

  return posmax;
  
}

void SetHistColor(TH1D *h, int color){

  h->SetLineColor(color);
  h->SetMarkerColor(color);
  
}

void PlotHPos(TCanvas *c, vector<TH1D*> vH){

  c->SetGrid();
  auto pos = GetPosHistMax(vH);
  vH[pos]->Draw();
  for(auto ith = vH.begin(); ith != vH.end(); ith++){

    if(pos != distance(vH.begin(), ith)){

      (*ith)->Draw("SAME");
      
    }

  }

  // c->BuildLegend();

}

void PlotFitResults(const char *f0pct,
		    const char *f1pct = NULL,
		    const char *f5pct = NULL,
		    const char *f10pct= NULL){

  vector<TH1D*> vhPos_0pct(3);
  vector<TH1D*> vhPos_1pct(3);
  vector<TH1D*> vhPos_5pct(3);
  vector<TH1D*> vhPos_10pct(3);
  // TGraph *gr[3];
  TGraphErrors *gr[3];
  for(unsigned iPos=0; iPos<3;iPos++){

    // gr[iPos] = new TGraph();
    gr[iPos] = new TGraphErrors();

  }
  TCanvas *c1;

  auto *mg = new TMultiGraph();

  vector<double> WbLSFrac = {0., 0.01, 0.05, 0.10};
  vector<int> Color = {kBlue-4, kRed-4, kGreen+1};
  vector<const char*> Leg = {"X", "Y", "Z"};
  TLegend *leg = new TLegend(0.92,0.7,0.99,0.99);

  vector<bool> vIsLargeDrift = {true, false, false};
  // vIsLargeDrift[0] = true;

  for(unsigned iPos=0; iPos<3;iPos++){

    gr[iPos]->SetLineColor(Color[iPos]);
    gr[iPos]->SetLineWidth(2);
    gr[iPos]->SetMarkerColor(Color[iPos]);
    gr[iPos]->SetMarkerStyle(kFullCrossX);
    gr[iPos]->SetMarkerSize(1);

    if(f0pct){
      vhPos_0pct[iPos] = GetHPos(f0pct, Form("hPosGuess%d",iPos));
      RebinHist(vhPos_0pct[iPos]);
      SetHistColor(vhPos_0pct[iPos], Color[iPos]);
      vhPos_0pct[iPos]->SetTitle(Leg[iPos]);
      SetGrErrPoint(vhPos_0pct[iPos], gr[iPos], WbLSFrac[0], vIsLargeDrift[iPos]);
      // SetCondensedGrErrPoint(vhPos_0pct[iPos], gr[iPos], WbLSFrac[0], vIsLargeDrift[iPos]);
    }

    if(f1pct){
      vhPos_1pct[iPos] = GetHPos(f1pct, Form("hPosGuess%d",iPos));
      RebinHist(vhPos_1pct[iPos]);
      SetHistColor(vhPos_1pct[iPos], Color[iPos]);
      vhPos_1pct[iPos]->SetTitle(Leg[iPos]);
      SetGrErrPoint(vhPos_1pct[iPos], gr[iPos], WbLSFrac[1], vIsLargeDrift[iPos]);
      // SetCondensedGrErrPoint(vhPos_1pct[iPos], gr[iPos], WbLSFrac[1], vIsLargeDrift[iPos]);
    }

    if(f5pct){
      vhPos_5pct[iPos] = GetHPos(f5pct, Form("hPosGuess%d",iPos));
      RebinHist(vhPos_5pct[iPos]);
      SetHistColor(vhPos_5pct[iPos], Color[iPos]);
      vhPos_5pct[iPos]->SetTitle(Leg[iPos]);
      SetGrErrPoint(vhPos_5pct[iPos], gr[iPos], WbLSFrac[2], vIsLargeDrift[iPos]);
      // SetCondensedGrErrPoint(vhPos_5pct[iPos], gr[iPos], WbLSFrac[2], vIsLargeDrift[iPos]);
    }

    if(f10pct){
      vhPos_10pct[iPos] = GetHPos(f10pct, Form("hPosGuess%d",iPos));
      RebinHist(vhPos_10pct[iPos]);
      SetHistColor(vhPos_10pct[iPos], Color[iPos]);
      vhPos_10pct[iPos]->SetTitle(Leg[iPos]);
      SetGrErrPoint(vhPos_10pct[iPos], gr[iPos], WbLSFrac[3], vIsLargeDrift[iPos]);
      // SetCondensedGrErrPoint(vhPos_10pct[iPos], gr[iPos], WbLSFrac[3], vIsLargeDrift[iPos]);
    }

    // c1 = new TCanvas(Form("c%d",iPos),Form("c%d",iPos),800,600);
    // c1->SetGrid();
    // gr[iPos]->Draw("APL");
    mg->Add(gr[iPos], "lp");
    leg->AddEntry(gr[iPos], Leg[iPos] ,"lp");

  }

  c1 = new TCanvas("cmg","cmg", 800,600);
  c1->SetGrid();
  mg->SetTitle("5 MeV e- Pos (0,0,0) Dir (1,0,0) ; WbLS Frac ; Drift (mm)");
  mg->Draw("A");
  mg->GetYaxis()->SetRangeUser(-50,100);
  leg->Draw();
  c1->Print("mg.png","png");

  c1 = new TCanvas("c0pct", "c0pct", 800, 600);
  PlotHPos(c1, vhPos_0pct);
  leg->Draw();
  c1->Print("c0pct.png","png");
  
  c1 = new TCanvas("c1pct", "c1pct", 800, 600);
  PlotHPos(c1, vhPos_1pct);
  leg->Draw();
  c1->Print("c1pct.png","png");
  
  c1 = new TCanvas("c5pct", "c5pct", 800, 600);
  PlotHPos(c1, vhPos_5pct);
  leg->Draw();
  c1->Print("c5pct.png","png");

  c1 = new TCanvas("c10pct", "c10pct", 800, 600);
  PlotHPos(c1, vhPos_10pct);
  leg->Draw();
  c1->Print("c10pct.png","png");

}
