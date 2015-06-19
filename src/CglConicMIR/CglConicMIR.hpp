// Last edit: 10/05/2014
//
// Name:     CglConicMIR.hpp
// Author:   Aykut Bulut
//           Lehigh University
//           email: aykut@lehigh.edu, aykutblt@gmail.com
// Date:     10/05/2014
//-----------------------------------------------------------------------------
// Copyright (C) 2009, Aykut Bulut.  All Rights Reserved.

//
// creates and adds cuts to OsiConicSolverInterface.
// This is more practical. This operation is irreversable.
//

#ifndef CglConicMIR_H
#define CglConicMIR_H
// CCGL headers
#include "CglConicMIRParam.hpp"
#include "CglConicCutGenerator.hpp"
// STDLIB headers
#include <string>

///
/// A container for MIR cuts can be generated as follows,
/// keep remodeling rows, (CoinPackedMatrix)
/// keep cut, (CoinPackedMatrix)
/// keep which member in which cone will be replaced with what variable
//
//

class CglConicMIR: public CglConicCutGenerator {
  CglConicMIRParam * param_;
  std::set<int> choose_cut_row(OsiConicSolverInterface const & si);
  double phi(double a, double f);
  int add_cut(OsiConicSolverInterface & si, int cut_cone, int cut_var,
	      int cut_row);
  OsiConicSolverInterface const * solver_;
  // keeps track of cone members that are already remodeled (lifted) using
  // auxilary t variables.
  int * lifted_;
  // binding_[i] is 1 if cone i is binding, 0 otherwise
  int * binding_;
  int num_cuts_;
  void update_binding_cones();
  void add_slacks(OsiConicSolverInterface & si);
public:
  // default constructor
  CglConicMIR();
  // copy constructor
  CglConicMIR(const CglConicMIR & other);
  // copy assignment operator
  CglConicMIR & operator=(const CglConicMIR & rhs);
  /// Destructor
  virtual ~CglConicMIR();
  // Set the parameters
  void setParam(const CglConicMIRParam & param);
  // Return parameter object
  CglConicMIRParam * getParam() const {return param_;}
  // Virtual functions
  virtual void generateCuts(const OsiConicSolverInterface & si,
			    OsiConicCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());
  /// Return true if needs optimal basis to do cuts (will return true)
  virtual bool needsOptimalBasis() const { return false; }
  /// Clone
  virtual CglConicCutGenerator * clone() const;
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  // generate and add cuts
  void generateAndAddCuts(OsiConicSolverInterface & si,
			  const CglTreeInfo info = CglTreeInfo());
  int getNumCutsAdded() const;

};

#endif
