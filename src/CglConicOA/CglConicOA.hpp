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
class CglConicOA: public CglConicCutGenerator {
  CglConicOAParam * param_;
  // private functions
  // generate support given a point on the cone.
  // assumes the point is on the cone boundry.
  void generate_support(int size, OsiLorentzConeType type,
		       int const * members,
		       double const * point,
		       OsiRowCut * rc) const;
  void generate_support_lorentz(int size,
			       int const * members,
			       double const * sol,
			       OsiRowCut * rc) const;
  void generate_support_rotated_lorentz(int size,
				       int const * members,
				       double const * sol,
				       OsiRowCut * rc) const;
  // project given infeasible solution to an array of points on the cones
  // of the problem
  void project(int n, int num_cones, int const * cone_size,
	       OsiLorentzConeType const * cone_type,
	       int const * const * members, double const * sol,
	       double * point, int * const feasible) const;
  // project given infeasible solution to the cone
  void project_one(int n, int num_cones, int const * cone_size,
	       OsiLorentzConeType const * cone_type,
	       int const * const * members, double const * sol,
	       double * point, int * const feasible) const;
  // project given infeasible solution to an array of points the cones
  // of the problem.
  void project_random(int n, int num_cones, int const * cone_size,
	       OsiLorentzConeType const * cone_type,
	       int const * const * members, double const * sol,
	       double * const * point, int * const feasible,
	       int num_points) const;
  // project given infeasible solution to an array of points the cones
  // of the problem using polar coordinates.
  void project_trig(int n, int num_cones, int const * cone_size,
	       OsiLorentzConeType const * cone_type,
	       int const * const * members, double const * sol,
	       double * const * point, int * const feasible,
	       int num_points) const;
public:
  // default constructor
  CglConicOA();
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
  // generate linear/ordinary cuts.
  virtual void generateCuts(OsiConicSolverInterface const & si,
			    OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());
  // generate cuts for a linear solver interface
  virtual void generateCuts(OsiSolverInterface const & si, OsiCuts & cuts,
			    int num_cones, OsiLorentzConeType const * cone_type,
			    int const * cone_size, int const * const * members);
  virtual void generateCuts(OsiSolverInterface const & si, OsiCuts & cuts,
			    int num_cones,
			    OsiLorentzConeType const * cone_type,
			    int const * cone_size,
			    int const * const * members, int num_points);
  /// Return true if needs optimal basis to do cuts
  virtual bool needsOptimalBasis() const { return false; }
  /// Clone
  virtual CglConicCutGenerator * clone() const;
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
};

#endif
