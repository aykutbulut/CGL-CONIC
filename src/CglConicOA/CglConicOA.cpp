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

CglConicOA::CglConicOA(OsiConicSolverInterface * solver)
  : param_(0), solver_(solver) {
  param_ = new CglConicOAParam();
}

// copy constructor
CglConicOA::CglConicOA(const CglConicOA & other) {
  // copy param_
  param_ = new CglConicOAParam(*(other.getParam()));
  solver_ = other.solver();
}

// copy assignment operator
CglConicOA & CglConicOA::operator=(const CglConicOA & rhs) {
  // copy param_
  param_ = new CglConicOAParam(*(rhs.getParam()));
  solver_ = rhs.solver();
  return *this;
}

CglConicOA::~CglConicOA() {
  if (param_)
    delete param_;
}

void CglConicOA::setParam(CglConicOAParam const & param) {
  param_ = new CglConicOAParam(param);
}

void CglConicOA::generateCuts(const OsiConicSolverInterface & si,
			       OsiCuts & cuts,
			       const CglTreeInfo info) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
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

OsiConicSolverInterface * CglConicOA::solver() const {
  return solver_;
}

