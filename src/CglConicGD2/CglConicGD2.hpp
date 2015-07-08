// Last edit: 10/05/2014
// Implements cutting strategy defined by Karzan and Yildiz.
//
// Name:     CglConicGD2.hpp
// Author:   Aykut Bulut
//           Lehigh University
//           email: aykut@lehigh.edu, aykutblt@gmail.com
// Date:     10/05/2014
//-----------------------------------------------------------------------------
// Copyright (C) 2009, Aykut Bulut.  All Rights Reserved.

#ifndef CglConicGD2_H
#define CglConicGD2_H
// CCGL headers
#include "CglConicGD2Param.hpp"
#include "CglConicCutGenerator.hpp"
// STDLIB headers
#include <string>


class CglConicGD2: public CglConicCutGenerator {
  CglConicGD2Param * param_;
  int cut_generating_cone_;
  Disjunction * disjunction_;
  double beta_;
  // holds generated cuts
  OsiConicCuts * cuts_;
  // compute cut generating cone using param_
  void compute_cut_generating_cone();
  // compute disjunction using param_
  void compute_disjunction();
  // compute cut parameter beta using param_
  void compute_beta();
public:
  // default constructor
  CglConicGD2();
  // copy constructor
  CglConicGD2(const CglConicGD2 & other);
  // copy assignment operator
  CglConicGD2 & operator=(const CglConicGD2 & rhs);
  /// Destructor
  virtual ~CglConicGD2();
  // Set the parameters
  void setParam(const CglConicGD2Param & param);
  // Return parameter object
  CglConicGD2Param * getParam() const {return param_;}
  // Virtual functions
  // this one is inherited from CglConicCutGenerator
  virtual void generateCuts(const OsiConicSolverInterface & si,
			    OsiConicCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());
  /// Return true if needs optimal basis to do cuts (will return true)
  virtual bool needsOptimalBasis() const { return false; }
  /// Clone
  virtual CglConicCutGenerator * clone() const;
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
};

#endif
