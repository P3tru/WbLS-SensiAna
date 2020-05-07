using namespace std;

void Plot(){

  vector<double> WbLSFrac = {0.004, 0.01, 1.};
  vector<double> Quenching = {0.70, 0.44, 0.07};
  vector<double> Err = {TMath::Sqrt(0.12*0.12+0.07*0.07), TMath::Sqrt(0.01*0.01+0.04*0.04), TMath::Sqrt(0.01*0.01+0.01*0.01)};

  TGraphErrors *gr = new TGraphErrors();
  for(int iPt=0; iPt<3; iPt++){

    gr->SetPoint(gr->GetN(), WbLSFrac[iPt], Quenching[iPt]);
    gr->SetPointError(gr->GetN()-1, 0., Err[iPt]);
    
  }

  auto *c1 = new TCanvas("c1", "c1", 800,600);
  c1->SetGrid();
  c1->SetLogx();
  gr->SetTitle("Birks' quenching parameter VS WbLS Frac; WbLS Frac; Quenching Parameter");
  gr->Draw("APL");
  gr->SetLineWidth(2);
  gr->GetYaxis()->SetRangeUser(0.,1.);

  cout << "Birks constant for WbLS1% : " << gr->Eval(0.01) << " mm/MeV"  << endl;
  cout << "Birks constant for WbLS3% : " << gr->Eval(0.03) << " mm/MeV" << endl;
  cout << "Birks constant for WbLS5% : " << gr->Eval(0.05) << " mm/MeV" << endl;
  cout << "Birks constant for WbLS10% : " << gr->Eval(0.10) << " mm/MeV" << endl;
  
}
