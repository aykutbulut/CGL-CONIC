#include "CglConicOAParam.hpp"
#include <iostream>
#include <exception>

// Constructor.
CglConicOAParam::CglConicOAParam(double coneTol) {
  coneTol_ = coneTol;
}

// destructor
CglConicOAParam::~CglConicOAParam() {
}

/// Clone
CglParam * CglConicOAParam::clone() const {
  CglParam * par = new CglConicOAParam(*this);
  return par;
}
