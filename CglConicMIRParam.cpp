#include "CglConicMIRParam.hpp"

// constructor
CglConicMIRParam::CglConicMIRParam() {
}

// copy constructor
CglConicMIRParam::CglConicMIRParam(const CglConicMIRParam & other) {
}

// copy assignment operator
CglConicMIRParam & CglConicMIRParam::operator=(
		       const CglConicMIRParam & other) {
  return *this;
}

// destructor
CglConicMIRParam::~CglConicMIRParam() {
}

/// Clone
CglParam * CglConicMIRParam::clone() const {
  CglParam * par = new CglConicMIRParam(*this);
  return par;
}
