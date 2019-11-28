//
// Created by zsoldos on 11/26/19.
//

#ifndef EVENTWRAPPER__TPARTICLE_HH_
#define EVENTWRAPPER__TPARTICLE_HH_

#include <iostream>
#include <iomanip>

#include <TVector3.h>

using namespace std;

class TParticle {

 public:
  TParticle();
  virtual ~TParticle();

  void setProperties(string name, string gen_process, string end_process, TVector3 *steps, double *times, double *energies, int nSteps);
  void addChild(TParticle *child);

  int numRemaining();
  int numChildren();

  TParticle *getNext();
  TParticle *getChildren();

  void dumpTree(ostream &out);
  void dumpList(ostream &out);

  string name;
  string gen_process;
  string end_process;
  TVector3 *steps; //take ownership of this
  double *times;
  double *energies;
  int nSteps;

 protected:
  void dumpTree(int depth, ostream &out);

  //don't take ownership of any of these pointers
  TParticle *next;
  TParticle *children;
  TParticle *child_tail;

};

#endif //EVENTWRAPPER__TPARTICLE_HH_
