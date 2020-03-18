using namespace std;

vector<double> GetArray(vector< pair<double, double> > *vP);

double GenEBin(){

  static int iBin=1;
  static double EBin = 0.1;
  return EBin*iBin++;
  
}

void FitAndRecoverGaussParams(TH1D *h, double *mu, double *sigma){

  TF1 *fFit;  
  h->Fit("gaus", "Q0");
  fFit = h->GetFunction("gaus");
  double mean = fFit->GetParameter(1);
  double meanErr = fFit->GetParError(1);
  double sig = fFit->GetParameter(2);
  double sigErr = fFit->GetParError(2);

  *mu=mean;
  *sigma=sig;

  
}

void FillPair(TH1D *h, double EBin,
	      vector< pair<double, double> > *vPMu, vector< pair<double, double> > *vPSig){

  double mu, sig;
  FitAndRecoverGaussParams(h, &mu, &sig);
  vPMu->push_back( make_pair(EBin, mu) );
  vPSig->push_back( make_pair(EBin, sig) );
  
}

TGraph *CreateGraph(vector< pair<double, double> > *vP){

  sort(vP->begin(), vP->end());

  auto *gr = new TGraph();

  for(auto p: *vP){

    gr->SetPoint(gr->GetN(), p.first, p.second);
    
  }

  return gr;
  
}

TGraphErrors *CreateGraphErrors(vector< pair<double, double> > *vP,
				vector< pair<double, double> > *vPErr){

  sort(vP->begin(), vP->end());
  sort(vPErr->begin(), vPErr->end());

  auto *gr = new TGraphErrors();

  for(auto p: *vP){

    gr->SetPoint(gr->GetN(), p.first, p.second);
    
  }
  vector<double> sig = GetArray(vPErr);

  for(auto i=0; i<sig.size(); i++){

    gr->SetPointError(i, 0., sig[i]);
    
  }

  return gr;
  
}

vector<double> GetArray(vector< pair<double, double> > *vP){

  sort(vP->begin(), vP->end());

  vector<double> vD;

  for(auto p: *vP){

    vD.push_back(p.second);
    
  }

  return vD;

}

void CreateCalibSplines(const char *filename,
			TGraph *grMuPE=NULL, TGraph *grSigPE=NULL,
			TGraph *grMuHits=NULL, TGraph *grSigHits=NULL){

  auto *FileCalib = TFile::Open(filename);

  vector<double> EBins(100);
  generate(EBins.begin(), EBins.end(), GenEBin);

  vector< pair<double, double> > vMuPE;
  vector< pair<double, double> > vSigPE;
  
  vector< pair<double, double> > vMuHits;
  vector< pair<double, double> > vSigHits;

  TIter next(FileCalib->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {

    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    // X-axis: nPE
    // Y-axis: nHits
    TH1 *h = (TH1*)key->ReadObj();

    for(auto EBin:EBins){

      if(strcmp(Form("hEbin%.1f", EBin),h->GetName()) == 0){

	TH1D *hNPE = (TH1D*)FileCalib->Get(Form("hEbin%.1f_px", EBin));
	FillPair(hNPE, EBin, &vMuPE, &vSigPE);
	
	TH1D *hNHits = (TH1D*)FileCalib->Get(Form("hEbin%.1f_py", EBin));
	FillPair(hNHits, EBin, &vMuHits, &vSigHits);

	break;
	
      } // END if strcmp      
      
    } // END for EBin
    
  } // END while Key

  if(grMuPE){
    grMuPE = (TGraph*)CreateGraph(&vMuPE)->Clone();
  }
  if(grSigPE){
    grSigPE = (TGraph*)CreateGraph(&vSigPE)->Clone();
  }

  if(grMuHits){
    grMuHits = (TGraph*)CreateGraph(&vMuHits)->Clone();
  }
  if(grSigHits){
    grSigHits = (TGraph*)CreateGraph(&vSigHits)->Clone();
  }

  // auto *grMuPE = CreateGraph(&vMuPE);
  // grMuPE->SetLineColor(kBlue-4);
  // grMuPE->SetLineWidth(2);

  // auto *grMuPE = CreateGraphErrors(&vMuPE, &vSigPE);
  // grMuPE->SetLineColor(kBlue-4);
  // grMuPE->SetLineWidth(2);

  // auto *grMuHits = CreateGraph(&vMuHits);
  // grMuHits->SetLineColor(kRed-4);
  // grMuHits->SetLineWidth(2);

  // auto *grMuHits = CreateGraphErrors(&vMuHits, &vSigHits);
  // grMuHits->SetLineColor(kRed-4);
  // grMuHits->SetLineWidth(2);

  // auto *mg = new TMultiGraph();
  // mg->Add(grMuPE,"cp");
  // mg->Add(grMuHits,"cp");
  // mg->Draw("A");

}

double GenERes(){

  static int iBin=1;
  static double EBin = 0.1;
  return EBin*iBin++;

}

void ComputeLikelihood(const char *filename, double NPE, double NHits){
  
  vector<double> ERess(100);
  generate(ERess.begin(), ERess.end(), GenERes);

  auto *FileCalib = TFile::Open(filename);

  vector<double> EBins(100);
  generate(EBins.begin(), EBins.end(), GenEBin);

  vector< pair<double, double> > vMuPE;
  vector< pair<double, double> > vSigPE;
  
  vector< pair<double, double> > vMuHits;
  vector< pair<double, double> > vSigHits;

  TIter next(FileCalib->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {

    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    // X-axis: nPE
    // Y-axis: nHits
    TH1 *h = (TH1*)key->ReadObj();

    for(auto EBin:EBins){

      if(strcmp(Form("hEbin%.1f", EBin),h->GetName()) == 0){

	TH1D *hNPE = (TH1D*)FileCalib->Get(Form("hEbin%.1f_px", EBin));
	FillPair(hNPE, EBin, &vMuPE, &vSigPE);
	
	TH1D *hNHits = (TH1D*)FileCalib->Get(Form("hEbin%.1f_py", EBin));
	FillPair(hNHits, EBin, &vMuHits, &vSigHits);

	break;
	
      } // END if strcmp      
      
    } // END for EBin
    
  } // END while Key

  auto *grMuPE=new TGraph();
  auto *grSigPE=new TGraph();
  auto *grMuHits=new TGraph();
  auto *grSigHits=new TGraph();

  // if(grMuPE){
  //   grMuPE = CreateGraph(&vMuPE);
  // }
  // if(grSigPE){
  //   grSigPE = CreateGraph(&vSigPE);
  // }

  // if(grMuHits){
  //   grMuHits = CreateGraph(&vMuHits);
  // }
  // if(grSigHits){
  //   grSigHits = CreateGraph(&vSigHits);
  // }

  CreateCalibSplines(filename,
  		     grMuPE, grSigPE,
  		     grMuHits, grSigHits);

  auto *c1 = new TCanvas("c1", "c1", 800, 600);

  auto *mg = new TMultiGraph();
  mg->Add(grMuPE,"cp");
  mg->Add(grMuHits,"cp");
  mg->Draw("A");

  grMuPE->SetBit(TGraph::kIsSortedX);
  grSigPE->SetBit(TGraph::kIsSortedX);

  auto *grL = new TGraph();

  for(auto ERes: ERess){
									 
    double muPE, sigPE;
    double muHits, sigHits;

    muPE = grMuPE->Eval(ERes, 0 ,"S");
    sigPE = grSigPE->Eval(ERes, 0 ,"S");

    muHits = grMuHits->Eval(ERes, 0 ,"S");
    sigHits = grSigHits->Eval(ERes, 0 ,"S");

    double Likelihood = TMath::Gaus(NPE, muPE, sigPE)*TMath::Gaus(NHits, muHits, sigHits);
    grL->SetPoint(grL->GetN(), ERes, Likelihood);
    
  }

  c1 = new TCanvas("c2", "c2", 800, 600);
  grL->Draw();
  TF1 *fFit;  
  grL->Fit("gaus", "Q0");
  fFit = grL->GetFunction("gaus");
  double mean = fFit->GetParameter(1);
  double meanErr = fFit->GetParError(1);
  double sig = fFit->GetParameter(2);
  double sigErr = fFit->GetParError(2);

  cout << "ERec: " << mean << " +- " << sig << " MeV " << endl;
  
  
}
