
void PlotAll(){

  gStyle->SetPalette(kDarkRainBow);

  const unsigned int nFrac = 4;
  const double wblsFrac[nFrac] = {0., 0.01, 0.05, 0.10};
  const double dwblsFrac[nFrac] = {0., 0., 0., 0.};
  const unsigned int nCov = 3;
  const double PECoverage[nCov] = {0.3, 0.5, 0.7};
  const double dPECoverage[nCov] = {0., 0., 0.};

  // Cov30percent
  // const double NHitsMeV30percent[nFrac] = {13.28, 18.21, 50.02, 85.79};
  // const double dNHitsMeV30percent[nFrac] = {0.11, 0.14, 0.19, 0.26};
  const double NHitsMeV30percent[nFrac] = {2.75162e+01,4.87803e+01, 1.61785e+02, 2.89749e+02 };
  const double dNHitsMeV30percent[nFrac] = {6.23532e+00, 7.75317e+00, 1.33405e+01, 1.74006e+01};
  TGraphErrors *gr30percent = new TGraphErrors(nFrac, wblsFrac, NHitsMeV30percent,
					       dwblsFrac, dNHitsMeV30percent);
  gr30percent->SetMarkerStyle(20);

  // Cov50percent
  // const double NHitsMeV50percent[nFrac] = {22.67, 31.27, 86.47, 86.68};
  // const double dNHitsMeV50percent[nFrac] = {0.12, 0.15, 0.33, 0.28};
  const double NHitsMeV50percent[nFrac] = {4.50312e+01, 8.23783e+01, 2.73371e+02};
  const double dNHitsMeV50percent[nFrac] = {8.61573e+00, 1.02169e+01, 1.74615e+01};
  TGraphErrors *gr50percent = new TGraphErrors(nFrac-1, wblsFrac, NHitsMeV50percent,
					       dwblsFrac, dNHitsMeV50percent);
  gr50percent->SetMarkerStyle(21);

  // Cov70percent
  // const double NHitsMeV70percent[nFrac] = {33.85, 47.23, 131.71, 230.50};
  // const double dNHitsMeV70percent[nFrac] = {0.18, 0.21, 0.27, 0.42};
  const double NHitsMeV70percent[nFrac] = {6.97052e+01, 1.27646e+02, 4.20850e+02, 7.58632e+02};
  const double dNHitsMeV70percent[nFrac] = {1.18162e+01,1.34051e+01, 2.26087e+01, 3.06383e+01};
  TGraphErrors *gr70percent = new TGraphErrors(nFrac, wblsFrac, NHitsMeV70percent,
					       dwblsFrac, dNHitsMeV70percent);
  gr70percent->SetMarkerStyle(23);

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("NHits/MeV for a 2.5 MeV positron in 1kT Theia; WbLS Fraction (%); NHits");
  mg->Add(gr30percent,"lp");
  mg->Add(gr50percent,"lp");
  mg->Add(gr70percent,"lp");

  TLegend *leg = new TLegend(0.13,0.56,0.42,0.87);
  leg->AddEntry(gr30percent,"30% Cov");
  leg->AddEntry(gr50percent,"50% Cov");
  leg->AddEntry(gr70percent,"70% Cov");

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  // c1->SetGrid();
  mg->Draw("A PMC PLC");
  leg->Draw();

  // ############################ //
  // -------- NEUTRONS ---------- //
  // ############################ //

  const unsigned int DopedWaternFrac = 2;
  const double DopedWaterLSFrac[DopedWaternFrac] = {0., 0.01};

  // Cov30percent
  const double NNHitsMeV30percent[nFrac] = {0.14, 0.94, 0.999, 0.999};
  TGraph *grN30percent = new TGraph(nFrac, wblsFrac, NNHitsMeV30percent);
  grN30percent->SetMarkerStyle(20);
  const double NNHitsMeVDopedWater30percent[DopedWaternFrac] = {0.87, 0.99};
  TGraph *grNDopedWater30percent =
    new TGraph(DopedWaternFrac, DopedWaterLSFrac, NNHitsMeVDopedWater30percent);
  grNDopedWater30percent->SetMarkerStyle(20);

  // Cov50percent
  const double NNHitsMeV50percent[nFrac] = {0.61, 0.999, 0.999, 0.999};
  TGraph *grN50percent = new TGraph(nFrac, wblsFrac, NNHitsMeV50percent);
  grN50percent->SetMarkerStyle(21);
  const double NNHitsMeVDopedWater50percent[DopedWaternFrac] = {0.95, 0.999};
  TGraph *grNDopedWater50percent =
    new TGraph(DopedWaternFrac, DopedWaterLSFrac, NNHitsMeVDopedWater50percent);
  grNDopedWater50percent->SetMarkerStyle(21);

  // Cov70percent
  const double NNHitsMeV70percent[nFrac] = {0.86, 0.999, 0.999, 0.999};
  TGraph *grN70percent = new TGraph(nFrac, wblsFrac, NNHitsMeV70percent);
  grN70percent->SetMarkerStyle(23);
  const double NNHitsMeVDopedWater70percent[DopedWaternFrac] = {0.98, 0.999};
  TGraph *grNDopedWater70percent =
    new TGraph(DopedWaternFrac, DopedWaterLSFrac, NNHitsMeVDopedWater70percent);
  grNDopedWater70percent->SetMarkerStyle(23);

  TMultiGraph *mgN = new TMultiGraph();
  mgN->SetTitle("Detection eff n-capture for 1kT Theia; WbLS Fraction (%); Eff (%)");
  mgN->Add(grN30percent,"lp");
  mgN->Add(grN50percent,"lp");
  mgN->Add(grN70percent,"lp");

  mgN->Add(grNDopedWater30percent,"lp");
  mgN->Add(grNDopedWater50percent,"lp");
  mgN->Add(grNDopedWater70percent,"lp");

  TLegend *legN = new TLegend(0.13,0.56,0.42,0.87);
  legN->AddEntry(grN30percent,"30% Cov");
  legN->AddEntry(grN50percent,"50% Cov");
  legN->AddEntry(grN70percent,"70% Cov");

  legN->AddEntry(grNDopedWater30percent,"30% Cov + Gd");
  legN->AddEntry(grNDopedWater50percent,"50% Cov + Gd");
  legN->AddEntry(grNDopedWater70percent,"70% Cov + Gd");

  TCanvas *c1n = new TCanvas("c1n","c1n",800,600);
  // c1n->SetGrid();
  mgN->Draw("A PMC PLC");
  legN->Draw();


}
