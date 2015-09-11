#include "CglConicOAParam.hpp"
#include <iostream>
#include <exception>

// constructor
CglConicOAParam::CglConicOAParam() {
}

// copy constructor
CglConicOAParam::CglConicOAParam(const CglConicOAParam & other) {
}

// copy assignment operator
CglConicOAParam & CglConicOAParam::operator=(
		       const CglConicOAParam & other) {
  std::cerr << "operator= is not implemented yet!" << std::endl;
  throw std::exception();
  return *this;
}

// destructor
CglConicOAParam::~CglConicOAParam() {
}

/// Clone
CglParam * CglConicOAParam::clone() const {
  CglParam * par = new CglConicOAParam(*this);
  return par;
}
