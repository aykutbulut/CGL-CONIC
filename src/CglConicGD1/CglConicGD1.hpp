// Name:     CglConicGD1.hpp
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

#ifndef CglConicGD1_H
#define CglConicGD1_H

// CCGL headers
#include "CglConicGD1Param.hpp"
#include "CglConicCutGenerator.hpp"
#include "CglConicGD1Cut.hpp"
// STDLIB headers
#include <string>

/// generates and adds cuts to solver_
class CglConicGD1: public CglConicCutGenerator {
  CglConicGD1Param * param_;
  CglConicGD1Cut * ccut_;
  int num_cuts_;
  std::vector<CglConicGD1Cut*> cuts_;
  std::vector<int> cuts_cone_ind_;
  // get equality constraints the cut cone members are in
  // get Ax=b
  void get_rows(OsiConicSolverInterface const & si,
                int cut_cone, int & num_eq_rows, int *& rows) const;
  void add_cut(OsiConicSolverInterface * solver, CglConicGD1Cut const * cut);
  // adds generated cuts to the model.
  void add_cuts(OsiConicSolverInterface * solver);
  // frees memory
  void clear_cuts();
  // compute disjunction var and the cone var is in.
  std::vector<std::pair<int, int> > compute_dis_var_cone(
                                    OsiConicSolverInterface const & si) const;
  void add_cone_from_cut(OsiConicSolverInterface & solver,
                         CglConicGD1Cut const * cut,
                         int cut_cone_index);
  void print_cut(CglConicGD1Cut const * cut) const;
  void get_input_set(OsiConicSolverInterface const * solver,
                     int dis_var,
                     int cut_cone,
                     int num_eq_rows,
                     int const * rows,
                     CoinPackedMatrix *& A,
                     double *& b,
                     double *& sol,
                     int & rel_dis_var) const;

public:
  // default constructor
  CglConicGD1(OsiConicSolverInterface * solver);
  // copy constructor
  CglConicGD1(const CglConicGD1 & other);
  // copy assignment operator
  CglConicGD1 & operator=(const CglConicGD1 & rhs);
  /// Destructor
  virtual ~CglConicGD1();
  // Set the parameters
  void setParam(const CglConicGD1Param & param);
  // Return parameter object
  CglConicGD1Param * getParam() const {return param_;}
  // Virtual functions
  virtual void generateCuts(const OsiConicSolverInterface & si,
                            OsiConicCuts & cs,
                            const CglTreeInfo info = CglTreeInfo());
  /// Return true if needs optimal basis to do cuts
  virtual bool needsOptimalBasis() const { return false; }
  /// Clone
  virtual CglConicCutGenerator * clone() const;
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  // generate and add cuts, return conic interface, primal form
  OsiConicSolverInterface * generateAndAddCuts(OsiConicSolverInterface const & si,
                                    const CglTreeInfo info = CglTreeInfo());
  OsiConicSolverInterface * generateAndAddBestCut(
                                    OsiConicSolverInterface const & si,
                                    const CglTreeInfo info = CglTreeInfo());
  // generate and add cuts, return conic interface, dual form
  OsiConicSolverInterface * generateAndAddCuts(
                              OsiConicSolverInterface const & si,
                              std::vector<CoinPackedMatrix const *> AA,
                              std::vector<double const *> bb);
  OsiConicSolverInterface * generateAndAddBestCut(
                              OsiConicSolverInterface const & si,
                              std::vector<CoinPackedMatrix const *> AA,
                              std::vector<double const *> bb);
  int getNumCutsAdded() const;
};

#endif
