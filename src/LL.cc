#include <TMath.h>

#include <LL.hh>

double EvalL(double Nobs, double Npred){
  double L;
  if(Nobs>0 && Npred>0)
	L=Npred-Nobs+Nobs*TMath::Log(Nobs/Npred);
  else
	L=Npred;
  L = -L;
  return L;
}
