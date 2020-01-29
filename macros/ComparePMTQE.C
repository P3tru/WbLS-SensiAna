#include "ComparePMTQE.hh"
#include "data.hh"

void ComparePMTQE(){

  gStyle->SetPalette(kDarkRainBow);

  //Colors
  const int SNO_Color = kBlue;
  const int Theia_Color = kRed;

  // Dimensions
  const double SNO_Radius = 8500; //distance in mm from center to PMT

  // const double Theia_Buffer = 543; // 
  const double Theia_Buffer = 3000; // 
  const double Theia_Height = 10857+2*Theia_Buffer; // FV+buffer
  const double Theia_Diameter = 10857+2*Theia_Buffer; // FV+buffer
  const double Theia_MeanLength = ComputeMeanLengthTheia(Theia_Height/2); //mean distance travel by photons in mm

  std::cout << "Theia_MeanLength: " << Theia_MeanLength << std::endl;

  const double SNO_PE_Coverage = 0.31;
  const double SNO_PE_Coverage_W_Collector = 0.59;
  const double SNO_PE_Coverage_W_Collector_Ref = 0.54;

  const double Theia_PE_Coverage = 0.30;

  const double AreaRatio = ComputeAreaRatio(SNO_Radius, SNO_PE_Coverage_W_Collector_Ref,
					    Theia_Height, Theia_Diameter/2, Theia_PE_Coverage);
  const double SolidAngleRatio = ComputeSolidAngleRatio(Theia_Height, Theia_Diameter/2);

  std::cout << "AreaRatio: " << AreaRatio << std::endl;
  std::cout << "SolidAngleRatio: " << SolidAngleRatio << std::endl;
    
  //  Photocathode optics for the Hamamatsu R1408
  //  SNO PMTs
  std::vector<double> R1408_WL = GetR1408_WL();
  std::vector<double> R1408_QE = GetR1408_QE();
  const unsigned int R1408_size = R1408_WL.size();
  TGraph *grR1408 = new TGraph(R1408_size, &R1408_WL[0], &R1408_QE[0]);
  grR1408->SetLineColor(SNO_Color-4); // 
  grR1408->SetLineStyle(1);
  grR1408->SetLineWidth(1);

  //  Photocathode optics for the Hamamatsu R11780
  //  Theia PMTs
  std::vector<double> R11780_WL = GetR11780_WL();
  std::vector<double> R11780_QE = GetR11780_QE();
  const unsigned int R11780_size = R11780_WL.size();
  TGraph *grR11780 = new TGraph(R11780_size, &R11780_WL[0], &R11780_QE[0]);
  grR11780->SetLineColor(Theia_Color-4);
  grR11780->SetLineStyle(1);
  grR11780->SetLineWidth(1);
  
  //  Water optics from rat-pac
  //  Theia lightwater
  std::vector<double> H2O_WL = GetH2O_WL();
  std::vector<double> H2O_ATT = GetH2O_ATT();
  const unsigned int H2O_size = H2O_WL.size();
  std::vector<double> H2O_RINDEX_WL = GetH2O_RINDEX_WL();
  std::vector<double> H2O_RINDEX = GetH2O_RINDEX();
  const unsigned int H2O_RINDEX_size = H2O_RINDEX_WL.size();
  TGraph *grH2O_RINDEX = new TGraph(H2O_RINDEX_size, &H2O_RINDEX_WL[0], &H2O_RINDEX[0]);

  
  double H2O_ABS[H2O_size];
  double H2O_NCER[H2O_size];
  double Theia_CE[H2O_size];
  for(unsigned int ipt=0;ipt<H2O_size;ipt++){
    H2O_ABS[ipt]=ATT2ABS(Theia_MeanLength/10, H2O_ATT[ipt]); // convert to cm
    H2O_NCER[ipt]=nCerenkov(5,1e-2, grH2O_RINDEX->Eval(H2O_WL[ipt]), H2O_WL[ipt]);
    Theia_CE[ipt]=grR11780->Eval(H2O_WL[ipt])*H2O_ABS[ipt]*H2O_NCER[ipt];
    // Theia_CE[ipt]=grR11780->Eval(H2O_WL[ipt])*H2O_ABS[ipt];
  }

  //  Water optics from SNO rat
  //  SNO heavywater
  std::vector<double> D2O_WL = GetD2O_WL();  
  std::vector<double> D2O_ATT = GetD2O_ATT();
  const unsigned int D2O_size = D2O_WL.size();
  std::vector<double> D2O_RINDEX = GetD2O_RINDEX();
  //  SNO lightwater
  std::vector<double> H2O_SNO_WL = GetH2O_SNO_WL();
  std::vector<double> H2O_SNO_ATT = GetH2O_SNO_ATT();
  const unsigned int H2O_SNO_size = H2O_SNO_WL.size();

  double D2O_ABS[D2O_size];
  double D2O_NCER[D2O_size];
  double SNO_CE[D2O_size];
  for(int ipt=0;ipt<D2O_size;ipt++){
    D2O_ABS[ipt]=ATT2ABS(SNO_Radius, D2O_ATT[ipt]);
    D2O_NCER[ipt]=nCerenkov(5,1e-2, D2O_RINDEX[ipt], D2O_WL[ipt]);
    SNO_CE[ipt]=grR1408->Eval(D2O_WL[ipt])*D2O_ABS[ipt]*D2O_NCER[ipt];
    // SNO_CE[ipt]=grR1408->Eval(D2O_WL[ipt])*D2O_ABS[ipt];
  }

  TGraph *grD2O = new TGraph(D2O_size, &D2O_WL[0], &D2O_ABS[0]);
  grD2O->SetLineColor(SNO_Color);
  grD2O->SetLineStyle(2);
  grD2O->SetLineWidth(2);
  TGraph *grD2OnCer = new TGraph(D2O_size, &D2O_WL[0], &D2O_NCER[0]);
  grD2OnCer->SetLineColor(SNO_Color+4);
  grD2OnCer->SetLineStyle(3);
  grD2OnCer->SetLineWidth(2);
  TGraph *grSNO = new TGraph(D2O_size, &D2O_WL[0], &SNO_CE[0]);
  grSNO->SetLineColor(SNO_Color-4);
  grSNO->SetLineStyle(3);
  grSNO->SetLineWidth(2);

  TGraph *grH2O = new TGraph(H2O_size, &H2O_WL[0], &H2O_ABS[0]);
  grH2O->SetLineColor(Theia_Color);
  grH2O->SetLineStyle(2);
  grH2O->SetLineWidth(2);
  TGraph *grH2OnCer = new TGraph(H2O_size, &H2O_WL[0], &H2O_NCER[0]);
  grH2OnCer->SetLineColor(Theia_Color+4);
  grH2OnCer->SetLineStyle(3);
  grH2OnCer->SetLineWidth(2);
  TGraph *grTheia = new TGraph(H2O_size, &H2O_WL[0], &Theia_CE[0]);
  grTheia->SetLineColor(Theia_Color-4);
  grTheia->SetLineStyle(3);
  grTheia->SetLineWidth(2);
  
  TLegend *leg = new TLegend(0.13,0.6,0.43,0.8);
  leg->AddEntry(grR11780,"R11780 QE","PL");
  leg->AddEntry(grH2O,"Beer-Lambert H2O","PL");
  leg->AddEntry(grR1408,"R1408 QE","PL");
  leg->AddEntry(grD2O,"Beer-Lambert D2O","PL");
  leg->AddEntry(grH2OnCer,"Nb Cer Theia","PL");
  leg->AddEntry(grD2OnCer,"Nb Cer SNO","PL");
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("QE for Hamamatsu R1408 VS R11780 ; Wavelength (nm) ; QE (%)");
  mg->Add(grR1408,"PL");
  mg->Add(grR11780,"PL");
  mg->Add(grH2O,"PL");
  mg->Add(grD2O,"PL");
  mg->Add(grH2OnCer,"PL");
  mg->Add(grD2OnCer,"PL");

  TLegend *legPE = new TLegend(0.6,0.75,0.9,0.9);
  legPE->AddEntry(grH2OnCer,"Nb Cer Theia","PL");
  legPE->AddEntry(grD2OnCer,"Nb Cer SNO","PL");
  TMultiGraph *mgPE = new TMultiGraph();
  mgPE->SetTitle("nPE Cerenkov created in Heavywater SNO VS Lightwater Theia ; Wavelength (nm) ; QE (%)");
  mgPE->Add(grH2OnCer,"PL");
  mgPE->Add(grD2OnCer,"PL");

  TLegend *legCE = new TLegend(0.7,0.75,0.9,0.9);
  legCE->AddEntry(grSNO,"CE SNO","PL");
  legCE->AddEntry(grTheia,"CE Theia","PL");
  TMultiGraph *mgCE = new TMultiGraph();
  mgCE->SetTitle("CE for Hamamatsu R1408 (SNO) VS R11780 (Theia) ; Wavelength (nm) ; A.U");
  mgCE->Add(grSNO,"PL");
  mgCE->Add(grTheia,"PL");
    
  TCanvas *c1 = new TCanvas("cQE","cQE",800,600);
  c1->SetGrid();
  mg->Draw("A");
  leg->Draw();

  c1 = new TCanvas("cCE","cCE",800,600);
  c1->SetGrid();
  mgCE->Draw("A");
  legCE->Draw();

  c1 = new TCanvas("cPE","cPE",800,600);
  c1->SetGrid();
  mgPE->Draw("A");
  legPE->Draw();

  double inf=300.0;
  double sup=500.0;

  TF1 * fSNO = new TF1("fSNO",[&](double*x, double *p){ return grSNO->Eval(x[0]); }, 60.0, 800.0, 1);
  fSNO->FixParameter(0,1.);
  double integralSNO = fSNO->Integral(inf,sup);
  TF1 * fTheia = new TF1("fTheia",[&](double*x, double *p){ return grTheia->Eval(x[0]); }, 60.0, 800.0, 1);
  fTheia->FixParameter(0,1.);
  double integralTheia = fTheia->Integral(inf,sup);

  std::cout << "Ratio of integral: " << integralTheia/integralSNO << std::endl;
  std::cout << "Overall increase factor from SNO: " << AreaRatio * SolidAngleRatio * integralTheia/integralSNO << std::endl;
}
