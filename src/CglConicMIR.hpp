// Last edit: 10/05/2014
//
// Name:     CglConicMIR.hpp
// Author:   Aykut Bulut
//           Lehigh University
//           email: aykut@lehigh.edu, aykutblt@gmail.com
// Date:     10/05/2014
//-----------------------------------------------------------------------------
// Copyright (C) 2009, Aykut Bulut.  All Rights Reserved.

#ifndef CglConicMIR_H
#define CglConicMIR_H
// CGL headers
#include "CglCutGenerator.hpp"
// STDLIB headers
#include <string>

class CglConicMIR: public CglCutGenerator {
  CglConicMIRParam * param_;
public:
  // default constructor
  CglConicMIR();
  /// Destructor
  virtual ~CglConicMIR();
  // Set the parameters
  void setParam(const CglConicMIRParam & param);
  // Return parameter object
  CglConicMIRParam getParam() const {return param_;}
  // Virtual functions
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo());
  /// Return true if needs optimal basis to do cuts (will return true)
  virtual bool needsOptimalBasis() const { return false; }
  /// Clone
  virtual CglCutGenerator * clone() const;
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
};

#endif
