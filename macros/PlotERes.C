#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

double errBinomial(double k, double N){
  return sqrt(k*(1-k/N))/N;
}

void PlotERes(std::string filename=""){

  TFile *file = new TFile(filename.c_str(), "READ");

  std::string histprefix = "hEbin";
  std::string histsuffix = "_px";

  std::vector<double> eBin(100);
  std::iota(eBin.begin(), eBin.end(), 1);

  TH1D *h1D = NULL;

  TGraphErrors *grERes = new TGraphErrors();
  grERes->SetTitle("Energy Resolution ; Positron Kin E (MeV) ; E Res (MeV)");

  unsigned int nbPtsProcessed = 0;

  TF1 *fitERes = new TF1("fitSqrt","sqrt([0]*x)",0,10);
  fitERes->SetParName(0,"Res");

  TH1D *hChi2 = new TH1D("hChi2", "Chi2 reduced", 50, 0., 1.);
  
  for(int iBin = 0; iBin<eBin.size(); iBin++){

    std::string histEBin = Form("%.1f",eBin[iBin]/10);
    std::string histname = histprefix + histEBin + histsuffix;

    h1D = (TH1D*)file->Get(histname.c_str());

    if(h1D){
      std::cout << histname << " OPEN" << std::endl;

      h1D->Fit("gaus","0LQSEM+");
      TF1 *fFit = h1D->GetFunction("gaus");
      
      if(fFit){
      
	double chi2 = fFit->GetChisquare()/(h1D->GetXaxis()->GetNbins()-1);
	double c = fFit->GetParameter(0);
	double mean = fFit->GetParameter(1);
	double meanErr = fFit->GetParError(1);
	double sigma = fFit->GetParameter(2);
	double sigmaErr = fFit->GetParError(2);

	std::cout << eBin[iBin]/10 << " " << eBin[iBin]/10 * (sigma/mean) << std::endl;

	std::cout << mean << " " << meanErr << std::endl;
	std::cout << sigma << " " << sigmaErr << std::endl;

	grERes->SetPoint(nbPtsProcessed,
			 eBin[iBin]/10,
			 eBin[iBin]/10 * (sigma/mean));
	grERes->SetPointError(nbPtsProcessed,
			      0.,
			      errBinomial(eBin[iBin]/10 * sigmaErr,mean));
			      // eBin[iBin]/10 * (sigmaErr/mean));

	hChi2->Fill(chi2);
	nbPtsProcessed++;
	
      }
      
    }
    
  }

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->SetGrid();
  grERes->Draw("APC");
  grERes->Fit("fitSqrt");

  c1 = new TCanvas("cChi2","cChi2",800,600);
  c1->SetGrid();
  hChi2->Draw("HIST");
  
}
