#include "CglConicMIR.hpp"

CglConicMIR::CglConicMIR(): param_(0) {
}

CglConicMIR::~CglConicMIR() {
  if (param_)
    delete param_;
}

CglConicMIR::setParam(const CglConicMIRParam & param) {
  param_ = new CglConicMIRParam(param);
}

// generates MIR cut.
void CglConicMIR::generateCuts(const OsiSolverInterface & si,
			       OsiCuts & cs,
			       const CglTreeInfo info = CglTreeInfo());

/// Clone
CglCutGenerator * CglConicMIR::clone() const {
  return this;
}

/// Create C++ lines to get to current state
std::string CglConicMIR::generateCpp( FILE * fp) {
  return std::string();
}
