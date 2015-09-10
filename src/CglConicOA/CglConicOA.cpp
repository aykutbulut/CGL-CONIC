// Name:     CglConicOA.cpp
// Author:   Aykut Bulut
//           Lehigh University
//           email: aykut@lehigh.edu, aykutblt@gmail.com
//-----------------------------------------------------------------------------
// Copyright (C) 2015, Lehigh University.  All Rights Reserved.

#include "CglConicOA.hpp"
#include <string>
#include <sstream>
#include <iomanip>

#define EPS 1e-5

CglConicOA::CglConicOA()
  : param_(0) {
  param_ = new CglConicOAParam();
}

// copy constructor
CglConicOA::CglConicOA(const CglConicOA & other) {
  // copy param_
  param_ = new CglConicOAParam(*(other.getParam()));
}

// copy assignment operator
CglConicOA & CglConicOA::operator=(const CglConicOA & rhs) {
  // copy param_
  param_ = new CglConicOAParam(*(rhs.getParam()));
  return *this;
}

CglConicOA::~CglConicOA() {
  if (param_)
    delete param_;
}

void CglConicOA::setParam(CglConicOAParam const & param) {
  param_ = new CglConicOAParam(param);
}

// generate outer approximating hyperplanes for conic constraints
// todo(aykut): approximates Lorentz cones only for now.
void CglConicOA::generateCuts(OsiConicSolverInterface const & si,
			       OsiCuts & cuts,
			       const CglTreeInfo info) {
  // get solution
  double const * sol = si.getColSolution();
  int num_cones = si.getNumCones();
  //OsiLorentzConeType * type = new OsiLorentzConeType[num_cones];
  OsiLorentzConeType type;
  int size;
  int * members;
  for (int i=0; i<num_cones; ++i) {
    si.getConicConstraint(i, type, size, members);
    // generate support for the cone
    OsiRowCut * rc = new OsiRowCut();
    int feas = generate_support(size, type, members, sol, rc);
    if (!feas) {
      cuts.insert(rc);
    }
    delete[] members;
  }
}

void CglConicOA::generateCuts(OsiConicSolverInterface const & si,
			       OsiConicCuts & cuts,
			       const CglTreeInfo info) {
}

/// Clone
CglConicCutGenerator * CglConicOA::clone() const {
  CglConicOA * new_cutg = new CglConicOA(*this);
  return new_cutg;
}

/// Create C++ lines to get to current state
std::string CglConicOA::generateCpp( FILE * fp) {
  std::cerr << "GenerateCpp is not implemented yet!" << std::endl;
  throw std::exception();
  return std::string();
}

int CglConicOA::generate_support(int size,
                                 OsiLorentzConeType type,
                                 int const * members,
                                 double const * sol,
                                 OsiRowCut * rc) const {
  int feas;
  if (type==OSI_QUAD) {
    feas = generate_support_lorentz(size, members, sol, rc);
  }
  else {
    feas = generate_support_rotated_lorentz(size,
                                            members, sol, rc);
  }
  return feas;
}

int CglConicOA::generate_support_lorentz(int size,
                                 int const * members,
                                 double const * sol,
                                 OsiRowCut * rc) const {
  int feas;
  double * par_point = new double[size];
  for(int j=0; j<size; ++j) {
    par_point[j] = sol[members[j]];
  }
  double * p = par_point;
  double term1;
  double term2;
  term1 = p[0];
  term2 = std::inner_product(p+1, p+size, p+1, 0.0);
  term2 = sqrt(term2);
  double activity = term1-term2;
  delete[] par_point;
  if (activity<-1e-5) {
    // current solution is infeasible to conic constraint i.
    double * coef = new double[size];
    double sum_rest;
    double x1 = term2;
    // cone is in canonical form
    for (int j=1; j<size; ++j) {
      coef[j] = 2.0*p[j];
    }
    coef[0] = -2.0*x1;
    // insert constraint (coef,0) to constraint pool
    rc->setRow(size, members, coef);
    // todo(aykut): fix setting infinity
    rc->setLb(-1e80);
    rc->setUb(0.0);
    delete[] coef;
    feas = 0;
  }
  else {
    feas = 1;
  }
  return feas;
}

int CglConicOA::generate_support_rotated_lorentz(int size,
                                                 int const * members,
                                                 double const * sol,
                                                 OsiRowCut * rc) const {
  int feas;
  double activity;
  double * par_point = new double[size];
  for(int j=0; j<size; ++j) {
    par_point[j] = sol[members[j]];
  }
  double * p = par_point;
  // check feasibility of the given solution
  double sum_rest = 0.0;
  sum_rest = std::inner_product(p+2, p+size, p+2, 0.0);
  activity = 2.0*p[0]*p[1]-sum_rest;
  if (activity<-1e-5) {
    //  at the end, set coef and lhs
    // map point from RLORENTZ space to LORENTZ space, find the projection on LORENTZ,
    // project this point to RLORENTZ and generate cut
    //sum_rest = std::inner_product(p+2, p+size, p+2, 0.0);
    double * coef = new double[size];
    double x1 = 0.0;
    double x2 = 0.0;
    // cone is a rotated cone
    // from point we move along [2point_2 2point_1 0 ... 0] until we hit
    // boundary. Then from this point in boundry we generate coef.
    // first compute u, step length
    double p1 = p[0];
    double p2 = p[1];
    p1 = p[0];
    p2 = p[1];
    x1 = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) - (-p1+p2)) / 2.0;
    x2 = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) + (-p1+p2)) / 2.0;
    // generate cut from xbar
    coef[0] = -2.0*x2;
    coef[1] = -2.0*x1;
    for (int i=2; i<size; ++i) {
      coef[i] = 2.0*p[i];
    }
    // insert constraint (coef,0) to constraint pool
    rc->setRow(size, members, coef);
    // todo(aykut): fix setting infinity
    rc->setLb(-1e80);
    rc->setUb(0.0);
    delete[] coef;
    feas = 0;
  }
  else {
    feas = 1;
  }
  delete[] p;
  return feas;
}
