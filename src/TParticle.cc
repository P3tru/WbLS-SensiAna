//
// Created by zsoldos on 11/26/19.
//

#include "TParticle.hh"

TParticle::TParticle() : name("UNSET"), gen_process("UNSET"), end_process("UNSET"), next(NULL), children(NULL), child_tail(NULL), steps(NULL), nSteps(0) {
}

TParticle::~TParticle() {
  if (steps) delete steps;
  if (times) delete times;
  if (energies) delete energies;
}

void TParticle::setProperties(string name, string gen_process, string end_process, TVector3 *steps, double *times, double *energies, int nSteps) {
  this->name = name;
  this->gen_process = gen_process;
  this->end_process = end_process;
  this->steps = steps;
  this->times = times;
  this->energies = energies;
  this->nSteps = nSteps;
}

void TParticle::addChild(TParticle *child) {
  if (children) {
	child_tail->next = child;
  } else {
	children = child;
  }
  child_tail = child;
}

int TParticle::numChildren() {
  if (children) return 1+children->numRemaining();
  return 0;
}

int TParticle::numRemaining() {
  int remaining = 0;
  TParticle *cur = this;
  while (cur->next) {
	remaining++;
	cur = cur->next;
  }
  return remaining;
}

TParticle *TParticle::getNext() {
  return next;
}

TParticle *TParticle::getChildren() {
  return children;
}

void TParticle::dumpList(ostream &out) {
  out << scientific << setprecision(5);
  out << '[';
  out << '\"' << name << "\",";
  out << '\"' << gen_process << "\",";
  out << '\"' << end_process << "\",";
  out << '[';
  for (int i = 0; i < nSteps-1; i++) {
	out << '[' << steps[i].X() << ',' << steps[i].Y() << ',' << steps[i].Z() << ',' << times[i] << ',' << energies[i] << "],";
  }
  out << '[' << steps[nSteps-1].X() << ',' << steps[nSteps-1].Y() << ',' << steps[nSteps-1].Z() << ',' << times[nSteps-1] << ',' << energies[nSteps-1] << "]";
  out << "],";
  out << '[';
  for (TParticle *child = children; child; child = child->next) {
	child->dumpList(out);
  }
  out << ']';
  out << ']';
  if (next) out << ",\n";
}

void TParticle::dumpTree(ostream &out) {
  dumpTree(0,out);
}

void TParticle::dumpTree(int depth, ostream &out) {
  for (int i = 0; i < depth; i++) {
	out << '|';
  }
  out << name << ' ';
  out << gen_process << ' ';
  out << end_process << ' ';
  out << '\n';
  for (TParticle *child = children; child; child = child->next) {
	child->dumpTree(depth+1,out);
  }
}
