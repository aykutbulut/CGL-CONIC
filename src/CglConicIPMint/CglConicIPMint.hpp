// Name:     CglConicIPMint.hpp
// Author:   Aykut Bulut
//           Lehigh University
//           email: aykut@lehigh.edu, aykutblt@gmail.com
//-----------------------------------------------------------------------------
// Copyright (C) 2015, Lehigh University.  All Rights Reserved.

//
// Cuts current solution from feasible set. This is done by solving a
// SOCO problem. Moreover conic feasible set is approximated around
// optimal solution of this problem.
//

#ifndef CglConicIPMint_H
#define CglConicIPMint_H
// CCGL headers
#include "CglConicIPMintParam.hpp"
#include "CglConicCutGenerator.hpp"
// STDLIB headers
#include <string>

// generate cuts using the current solution stored in the solver interface.
class CglConicIPMint: public CglConicCutGenerator {
  CglConicIPMintParam * param_;
  // fix integer variables store the model at solver_
  void fixIntegerVars(OsiSolverInterface const & si);
  // solver is preferable an IPM implementation for efficiency
  OsiConicSolverInterface * solver_;
  // generate cuts for a linear solver interface
  void method1(OsiSolverInterface const & si, OsiCuts & cuts,
	       int num_cones, OsiLorentzConeType const * cone_type,
	       int const * cone_size, int const * const * members,
	       int num_points);
  void method2(OsiSolverInterface const & si, OsiCuts & cuts,
	       int num_cones, OsiLorentzConeType const * cone_type,
	       int const * cone_size, int const * const * members,
	       int num_points);
  void method3(OsiSolverInterface const & si, OsiCuts & cuts,
	       int num_cones, OsiLorentzConeType const * cone_type,
	       int const * cone_size, int const * const * members,
	       int num_points);
  void add_supports(int size, double const * sol, int num_cones,
		OsiLorentzConeType const * cone_type,
		int const * cone_size, int const * const * members,
		OsiCuts & cuts) const;
  void add_cuts(int size, double const * sol, int num_cones,
		OsiLorentzConeType const * cone_type,
		int const * cone_size, int const * const * members,
		OsiCuts & cuts) const;
  // generate cuts where  point is on the cone
  void add_cuts2(int num_cols, double const * point, int num_cones,
		 OsiLorentzConeType const * cone_type,
		 int const * cone_size,
		 int const * const * members, OsiCuts & cuts);
  // private functions
  // return 0 if sol is infeasible for the cone, nonzero otherwise
  int generate_support(int size, OsiLorentzConeType type,
		       int const * members,
		       double const * sol,
		       OsiRowCut * rc) const;
  int generate_support_lorentz(int size,
			       int const * members,
			       double const * sol,
			       OsiRowCut * rc) const;
  int generate_support_rotated_lorentz(int size,
				       int const * members,
				       double const * sol,
				       OsiRowCut * rc) const;
  void create_rand_points(int num_cols, double const * sol,
			  int num_cones,
			  OsiLorentzConeType const * cone_type,
			  int const * cone_size,
			  int const * const * members,
			  double ** points, int num_points) const;

  void create_rand_point2(int num_cols, double const * sol,
			  int num_cones,
			  OsiLorentzConeType const * cone_type,
			  int const * cone_size,
			  int const * const * members,
			  double * point) const;
  void create_rand_point3(int cone_size, double const * par_sol,
			  OsiLorentzConeType cone_type,
			  int const * members,
			  double * par_point) const;
public:
  // default constructor
  CglConicIPMint();
  // copy constructor
  CglConicIPMint(const CglConicIPMint & other);
  // copy assignment operator
  CglConicIPMint & operator=(const CglConicIPMint & rhs);
  /// Destructor
  virtual ~CglConicIPMint();
  // Set the parameters
  void setParam(const CglConicIPMintParam & param);
  // Return parameter object
  CglConicIPMintParam * getParam() const {return param_;}
  // return solver
  OsiConicSolverInterface * getSolver() const;
  // generate linear/ordinary cuts.
  virtual void generateCuts(OsiConicSolverInterface const & si,
			    OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());
  // generate cuts for a linear solver interface
  virtual void generateCuts(OsiSolverInterface const & si, OsiCuts & cuts,
			    int num_cones, OsiLorentzConeType const * cone_type,
			    int const * cone_size, int const * const * members,
			    int num_points);
  /// Return true if needs optimal basis to do cuts
  virtual bool needsOptimalBasis() const { return false; }
  /// Clone
  virtual CglConicCutGenerator * clone() const;
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
};

#endif
