#include "CglConicCutGenerator.hpp"

// default constructor
CglConicCutGenerator::CglConicCutGenerator():
  aggressive_(0),
  canDoGlobalCuts_(0) {
}

// copy constructor
CglConicCutGenerator::CglConicCutGenerator(CglConicCutGenerator const & other) {
  setAggressiveness(other.aggressiveness());
  setGlobalCuts(other.canDoGlobalCuts());
}

// copy assignment operator
CglConicCutGenerator & CglConicCutGenerator::operator=(CglConicCutGenerator const & rhs) {
  setAggressiveness(rhs.aggressiveness());
  setGlobalCuts(rhs.canDoGlobalCuts());
  return *this;
}

CglConicCutGenerator::~CglConicCutGenerator() {
}

int CglConicCutGenerator::aggressiveness() const {
  return aggressive_;
}

void CglConicCutGenerator::setAggressiveness(int value) {
  aggressive_=value;
}

void CglConicCutGenerator::setGlobalCuts(bool ability) {
  canDoGlobalCuts_ = ability;
}

bool CglConicCutGenerator::canDoGlobalCuts() const {
  return canDoGlobalCuts_;
}

bool CglConicCutGenerator::mayGenerateRowCutsInTree() const {
  return true;
}

// Return true if needs optimal basis to do cuts
bool CglConicCutGenerator::needsOptimalBasis() const {
  return false;
}
