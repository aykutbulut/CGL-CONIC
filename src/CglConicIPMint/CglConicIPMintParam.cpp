#include "CglConicIPMintParam.hpp"
#include <iostream>
#include <exception>

// constructor
CglConicIPMintParam::CglConicIPMintParam() {
}

// copy constructor
CglConicIPMintParam::CglConicIPMintParam(const CglConicIPMintParam & other) {
}

// copy assignment operator
CglConicIPMintParam & CglConicIPMintParam::operator=(
		       const CglConicIPMintParam & other) {
  std::cerr << "operator= is not implemented yet!" << std::endl;
  throw std::exception();
  return *this;
}

// destructor
CglConicIPMintParam::~CglConicIPMintParam() {
}

/// Clone
CglParam * CglConicIPMintParam::clone() const {
  CglParam * par = new CglConicIPMintParam(*this);
  return par;
}
