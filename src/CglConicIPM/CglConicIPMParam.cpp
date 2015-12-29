#include "CglConicIPMParam.hpp"
#include <iostream>
#include <exception>

// constructor
CglConicIPMParam::CglConicIPMParam() {
}

// copy constructor
CglConicIPMParam::CglConicIPMParam(const CglConicIPMParam & other) {
}

// copy assignment operator
CglConicIPMParam & CglConicIPMParam::operator=(
		       const CglConicIPMParam & other) {
  std::cerr << "operator= is not implemented yet!" << std::endl;
  throw std::exception();
  return *this;
}

// destructor
CglConicIPMParam::~CglConicIPMParam() {
}

/// Clone
CglParam * CglConicIPMParam::clone() const {
  CglParam * par = new CglConicIPMParam(*this);
  return par;
}
