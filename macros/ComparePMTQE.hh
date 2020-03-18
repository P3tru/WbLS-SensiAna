double nCerenkov(double E, double l, double n, double lambda){

  const double mass_e=0.511;//MeV
  const double beta = sqrt(1 - (mass_e*mass_e)/(E*E));

  const double thetaC = acos(1/(beta*n));

  const double pi = 3.14159;
  const double twopi = 2*pi;

  const double epsilon0=55.26349406e3; //e^2/MeV/nm
  const double eC=1.; //J
  const double hbarc=1.23984193e-3/twopi; //MeV.nm;
  const double alpha = 1/(2*twopi*epsilon0) * (eC*eC)/(hbarc);

  return twopi*alpha*l*1e6*std::pow(sin(thetaC),2)/lambda;
  
}

double MeV2lambda(double MeV){

    const double hc = 1.23984193e-3; //MeV.nm
  return hc/lambda;

}

double lambda2MeV(double lambda){

  const double hc = 1.23984193e-3; //MeV.nm
  return hc/lambda;

}

double ATT2ABS(double x, double att){

  return TMath::Exp(-x/att);

}

double bethe(double E){

  return 0;  

}

double ComputeMeanLengthTheia(double Theia_R){

  return 0.5 * Theia_R * (TMath::Sqrt(2)+TMath::ASinH(1));
  
}

double ComputeAreaRatio(double SNO_R, double SNO_PE, double Theia_H, double Theia_R, double Theia_PE){

  const double pi = 3.14159;  
  return (2*pi*Theia_R*(Theia_R+Theia_H)*Theia_PE) / (4*pi*SNO_R*SNO_R*SNO_PE) ;

}

double ComputeSolidAngleRatio(double Theia_H, double Theia_R){

  const double pi = 3.14159;

  const double Omega_wall = 4*pi*Theia_H / (TMath::Sqrt(4*Theia_R*Theia_R + Theia_H*Theia_H));

  const double Theta = TMath::ATan(Theia_R/Theia_H);
  const double Omega_caps = 4*pi*TMath::Power(TMath::Sin(Theta/2),2);
  
  return (Omega_wall+2*Omega_caps) / (4*pi);

}
