//
// Created by zsoldos on 2/24/20.
//

#ifndef _HITCLASS_HH_
#define _HITCLASS_HH_

///////////////////////// STL C/C++ /////////////////////////
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include <TVector3.h>


#define PI 3.14159
#define C 299.792458 // mm/ns

using namespace std;

class Hit {

 protected:
  TVector3 Pos;

  double Q;
  double T;

  double TResid;
  double D;

 public:
  Hit() : Pos(0.,0.,0.), Q(-1), T(-1), TResid(-1), D(-1){};
  Hit(TVector3 p, double q=0., double t=0., double tr=-1., double d=-1.)
	  : Pos(p), Q(q), T(t), TResid(tr), D(d){};
  ~Hit(){};

  const TVector3 &GetPos() const {
	return Pos;
  }
  void SetPos(const TVector3 &pos) {
	Pos = pos;
  }

  double GetQ(){ return Q; };
  void SetQ(double x){ Q = x; };
  double GetT(){ return T; };
  void SetT(double x){ T = x; };

  double GetTResid() const {
	return TResid;
  }
  void SetTResid(double t_resid) {
	TResid = t_resid;
  }
  double GetD() const {
	return D;
  }
  void SetD(double d) {
	D = d;
  }

  double CalculateDistance(const TVector3& Origin = TVector3(0,0,0)){
	return TVector3(Pos-Origin).Mag();
  };
  double CalculateTResid(const TVector3& Origin = TVector3(0,0,0), double SoL = C){
	return T - (TVector3(Pos-Origin).Mag())/SoL;
  };

  bool operator==(const Hit &rhs) const {
	return T == rhs.T;
  }
  bool operator!=(const Hit &rhs) const {
	return !(rhs == *this);
  }
  bool operator<(const Hit &rhs) const {
	return T < rhs.T;
  }
  bool operator>(const Hit &rhs) const {
	return rhs < *this;
  }
  bool operator<=(const Hit &rhs) const {
	return !(rhs < *this);
  }
  bool operator>=(const Hit &rhs) const {
	return !(*this < rhs);
  }

  void operator+=(Hit rhs) {
	SetT(T+rhs.GetT());
  }
  void operator-=(Hit rhs) {
	SetT(T-rhs.GetT());
  }

  void Print() const {
	cout << "X: " << Pos.x()
		 << " Y: " << Pos.y()
		 << " Z: " << Pos.z()
		 << " Q: " << Q
		 << " T: " << T << endl;
  }

};

Hit operator-(Hit h1, Hit h2);

Hit operator+(Hit h1, Hit h2);

class HitCut{

 public:
  Hit hCut;

  explicit HitCut(const Hit &h_cut) : hCut(h_cut) {}

  int operator()(Hit h){
	bool YesItIs = false;
	if(h>hCut){
	  YesItIs = true;
	}
	return YesItIs;
  }

};


#endif //_HITCLASS_HH_
