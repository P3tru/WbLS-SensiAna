//
// Created by zsoldos on 2/24/20.
//

#ifndef _MCFUNCTIONS_HH_
#define _MCFUNCTIONS_HH_


///////////////////////// STL C/C++ /////////////////////////
#include <vector>

/////////////////////////   RAT   ///////////////////////////
#include <RAT/DS/MC.hh>

/////////////////////////   ROOT  ///////////////////////////
#include <TH2D.h>

/////////////////////////   USER  ///////////////////////////
#include "Analyzer.hh"
#include "HitClass.hh"

using namespace std;

// Access the MC object from RAT
// at iEvt (from MCTree)
RAT::DS::MC * GetRATMCOnEvt(Analyzer *fAnalyzer, unsigned int iEvt=0);

// Fill vector<Hit> with MC Info
vector<Hit> GetMCHitCollection(Analyzer *fAnalyzer, unsigned int iEvt=0);


void FillBirksLaw(Analyzer *fAnalyzer, unsigned int iEvt, TH2D *h2D= nullptr);

class TPhoton{

 protected:
  int ID;
  int ParentID;
  double WL;

 public:
  TPhoton(int id, int parent_id, double wl) : ID(id), ParentID(parent_id), WL(wl) {

  }
  virtual ~TPhoton() {

  }
  int GetId() const {
	return ID;
  }
  void SetId(int id) {
	ID = id;
  }
  int GetParentId() const {
	return ParentID;
  }
  void SetParentId(int parent_id) {
	ParentID = parent_id;
  }
  double GetWl() const {
	return WL;
  }
  void SetWl(double wl) {
	WL = wl;
  }

};

class TProton{

 protected:
  int ID;
  vector<TPhoton> photonChildren;
  double dE;
  double dX;

 public:
  TProton() {

	dE=-1;
	dX=-1;
	// photonChildren = vector<TPhoton>();

  }
  virtual ~TProton() {

  }

  vector<TPhoton> GetPhotonChildren() const {
	return photonChildren;
  }
  void AddPhoton(TPhoton photon){
	photonChildren.push_back(photon);
  }
  double GetDe() const {
	return dE;
  }
  void SetDe(double d_e) {
	dE = d_e;
  }
  double GetDx() const {
	return dX;
  }
  void SetDx(double d_x) {
	dX = d_x;
  }
  int GetId() const {
	return ID;
  }
  void SetId(int id) {
	ID = id;
  }
};


#endif