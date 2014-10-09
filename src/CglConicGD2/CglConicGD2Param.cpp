#include "CglConicGD2Param.hpp"

// constructor
CglConicGD2Param::CglConicGD2Param() {
}

// copy constructor
CglConicGD2Param::CglConicGD2Param(const CglConicGD2Param & other) {
}

// copy assignment operator
CglConicGD2Param & CglConicGD2Param::operator=(
		       const CglConicGD2Param & other) {
  return *this;
}

// destructor
CglConicGD2Param::~CglConicGD2Param() {
}

/// Clone
CglParam * CglConicGD2Param::clone() const {
  CglParam * par = new CglConicGD2Param(*this);
  return par;
}
