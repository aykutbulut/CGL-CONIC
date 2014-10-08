#include "CglConicMIR.hpp"

CglConicMIR::CglConicMIR(): param_(0) {
  param_ = new CglConicMIRParam();
}

// copy constructor
CglConicMIR::CglConicMIR(const CglConicMIR & other) {
  // copy param_
  param_ = new CglConicMIRParam(*(other.getParam()));
}

// copy assignment operator
CglConicMIR & CglConicMIR::operator=(const CglConicMIR & rhs) {
  // copy param_
  param_ = new CglConicMIRParam(*(rhs.getParam()));
  return *this;
}

CglConicMIR::~CglConicMIR() {
  if (param_)
    delete param_;
}

void CglConicMIR::setParam(const CglConicMIRParam & param) {
  param_ = new CglConicMIRParam(param);
}

// generates MIR cut.
void CglConicMIR::generateCuts(const OsiSolverInterface & si,
			       OsiCuts & cs,
			       const CglTreeInfo info) {
}

/// Clone
CglCutGenerator * CglConicMIR::clone() const {
  //
  CglConicMIR * new_cut = new CglConicMIR(*this);
  return new_cut;
}

/// Create C++ lines to get to current state
std::string CglConicMIR::generateCpp( FILE * fp) {
  return std::string();
}
