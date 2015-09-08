// Name:     CglConicOA.hpp
// Author:   Aykut Bulut
//           Lehigh University
//           email: aykut@lehigh.edu, aykutblt@gmail.com
//-----------------------------------------------------------------------------
// Copyright (C) 2015, Lehigh University.  All Rights Reserved.

//
// creates and adds cuts given in Belotti et al. to OsiConicSolverInterface.
// This is more practical. This operation is irreversable.
//
// Implements CglConicCutGenerator abstract base class.

#ifndef CglConicOA_H
#define CglConicOA_H
// CCGL headers
#include "CglConicOAParam.hpp"
#include "CglConicCutGenerator.hpp"
// STDLIB headers
#include <string>

// generate cuts using the current solution stored in the solver interface.
class CglConicGD1: public CglConicCutGenerator {
  CglConicOAParam * param_;
  OsiConicSolverInterface * solver_;
public:
  // default constructor
  CglConicOA(OsiConicSolverInterface * solver);
  // copy constructor
  CglConicOA(const CglConicOA & other);
  // copy assignment operator
  CglConicOA & operator=(const CglConicOA & rhs);
  /// Destructor
  virtual ~CglConicOA();
  // Set the parameters
  void setParam(const CglConicOAParam & param);
  // Return parameter object
  CglConicOAParam * getParam() const {return param_;}
  // Virtual functions
  virtual void generateCuts(const OsiConicSolverInterface & si,
			    OsiCuts & cuts,
			    const CglTreeInfo info = CglTreeInfo());
  /// Return true if needs optimal basis to do cuts
  virtual bool needsOptimalBasis() const { return false; }
  /// Clone
  virtual CglConicCutGenerator * clone() const;
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  // generate linear/ordinary cuts.
  virtual void generateCuts(const OsiConicSolverInterface & si,
			    OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());
  // return pointer to solver
  OsiConicSolverInterface * solver() const;
};

#endif
