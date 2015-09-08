#include "CglConicGD1Param.hpp"
#include <iostream>
#include <exception>

// constructor
CglConicGD1Param::CglConicGD1Param() {
}

// copy constructor
CglConicGD1Param::CglConicGD1Param(const CglConicGD1Param & other) {
}

// copy assignment operator
CglConicGD1Param & CglConicGD1Param::operator=(
		       const CglConicGD1Param & other) {
  std::cerr << "operator= is not implemented yet!" << std::endl;
  throw std::exception();
  return *this;
}

// destructor
CglConicGD1Param::~CglConicGD1Param() {
}

/// Clone
CglParam * CglConicGD1Param::clone() const {
  CglParam * par = new CglConicGD1Param(*this);
  return par;
}
