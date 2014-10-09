#include "CglConicGD2.hpp"

CglConicGD2::CglConicGD2(): param_(0) {
  param_ = new CglConicGD2Param();
}

// copy constructor
CglConicGD2::CglConicGD2(const CglConicGD2 & other) {
  // copy param_
  param_ = new CglConicGD2Param(*(other.getParam()));
}

// copy assignment operator
CglConicGD2 & CglConicGD2::operator=(const CglConicGD2 & rhs) {
  // copy param_
  param_ = new CglConicGD2Param(*(rhs.getParam()));
  return *this;
}

CglConicGD2::~CglConicGD2() {
  if (param_)
    delete param_;
}

void CglConicGD2::setParam(const CglConicGD2Param & param) {
  param_ = new CglConicGD2Param(param);
}

// generates MIR cut.
void CglConicGD2::generateCuts(const OsiSolverInterface & si,
			       OsiCuts & cs,
			       const CglTreeInfo info) {
}

/// Clone
CglCutGenerator * CglConicGD2::clone() const {
  //
  CglConicGD2 * new_cut = new CglConicGD2(*this);
  return new_cut;
}

/// Create C++ lines to get to current state
std::string CglConicGD2::generateCpp( FILE * fp) {
  return std::string();
}
