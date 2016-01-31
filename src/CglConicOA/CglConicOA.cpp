// Name:     CglConicOA.cpp
// Author:   Aykut Bulut
//           Lehigh University
//           email: aykut@lehigh.edu, aykutblt@gmail.com
//-----------------------------------------------------------------------------
// Copyright (C) 2015, Lehigh University.  All Rights Reserved.

#include "CglConicOA.hpp"
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <numeric>

#define CONE_EPS 1e-6
#define COEF_EPS 1e-5

CglConicOA::CglConicOA()
  : param_(0) {
  param_ = new CglConicOAParam();
}

// copy constructor
CglConicOA::CglConicOA(const CglConicOA & other) {
  // copy param_
  param_ = new CglConicOAParam(*(other.getParam()));
}

// copy assignment operator
CglConicOA & CglConicOA::operator=(const CglConicOA & rhs) {
  // copy param_
  param_ = new CglConicOAParam(*(rhs.getParam()));
  return *this;
}

CglConicOA::~CglConicOA() {
  if (param_)
    delete param_;
}

void CglConicOA::setParam(CglConicOAParam const & param) {
  param_ = new CglConicOAParam(param);
}

// generate outer approximating hyperplanes for conic constraints
void CglConicOA::generateCuts(OsiConicSolverInterface const & si,
			       OsiCuts & cuts,
			       const CglTreeInfo info) {
  int num_cones = si.getNumCones();
  OsiLorentzConeType * cone_type = new OsiLorentzConeType[num_cones];
  int ** members = new int*[num_cones];
  int * cone_size = new int[num_cones];
  for (int i=0; i<num_cones; ++i) {
    si.getConicConstraint(i, cone_type[i], cone_size[i], members[i]);
  }
  generateCuts(si, cuts, num_cones, cone_type, cone_size, members);
  delete[] cone_type;
  for (int i=0; i<num_cones; ++i) {
    delete[] members[i];
  }
  delete[] members;
  delete[] cone_size;
}

// generate cuts for a linear solver interface
// this function ask to project solution to num_points different points
// on the conic constraints.
void CglConicOA::generateCuts(OsiSolverInterface const & si, OsiCuts & cuts,
		  int num_cones, OsiLorentzConeType const * cone_type,
		  int const * cone_size, int const * const * members,
		  int num_points) {
  int n = si.getNumCols();
  double const * sol = si.getColSolution();
  double ** point = new double*[num_points];
  for (int k=0; k<num_points; ++k) {
    point[k] = new double[n];
  }
  int * feasible = new int[num_cones];
  // project sol to the conic constraints
  //project(n, num_cones, cone_size, cone_type, members, sol, point, feasible);
  //int num_points = 5;
  project_random(n, num_cones, cone_size, cone_type, members, sol, point, feasible, num_points);
  // iterate over points on the cone
  for (int k=0; k<num_points; ++k) {
    // iterate over cones and generate support
    for (int i=0; i<num_cones; ++i) {
      if (feasible[i]) {
	continue;
      }
      // generate support for cone i
      OsiRowCut * rc = new OsiRowCut();
      double * par_point = new double[cone_size[i]];
      for (int j=0; j<cone_size[i]; ++j) {
	par_point[j] = point[k][members[i][j]];
      }
      generate_support(cone_size[i], cone_type[i], members[i], par_point, rc);
      cuts.insert(rc);
      delete[] par_point;
      delete rc;
    }
  }
  for (int k=0; k<num_points; ++k) {
    delete[] point[k];
  }
  delete[] point;
  delete[] feasible;
}

// generate cuts for a given linear solver interface
void CglConicOA::generateCuts(OsiSolverInterface const & si, OsiCuts & cuts,
		  int num_cones, OsiLorentzConeType const * cone_type,
		  int const * cone_size, int const * const * members) {
  generateCuts(si, cuts, num_cones, cone_type, cone_size, members, 1);
}


/// Clone
CglConicCutGenerator * CglConicOA::clone() const {
  CglConicOA * new_cutg = new CglConicOA(*this);
  return new_cutg;
}

/// Create C++ lines to get to current state
std::string CglConicOA::generateCpp( FILE * fp) {
  std::cerr << "GenerateCpp is not implemented yet!" << std::endl;
  throw std::exception();
  return std::string();
}

void CglConicOA::project(int n, int num_cones, int const * cone_size,
			 OsiLorentzConeType const * cone_type,
			 int const * const * members, double const * sol,
			 double * point, int * const feasible) const {
  project_one(n, num_cones, cone_size, cone_type, members, sol, point,
	      feasible);
}

// project given infeasible solution to the cone
// set feasible[i] to 0 if solution violates cone i, 1 otherwise.
// n is size of solution
// store projected point in point array.
void CglConicOA::project_one(int n, int num_cones, int const * cone_size,
			 OsiLorentzConeType const * cone_type,
			 int const * const * members, double const * sol,
			 double * point, int * const feasible) const {
  // copy solution to the point
  std::copy(sol, sol+n, point);
  // iterate over cones and project them
  for (int i=0; i<num_cones; ++i) {
    // get data for cone i
    int size = cone_size[i];
    int const * mem = members[i];
    OsiLorentzConeType type = cone_type[i];
    // get relevent part of the solution
    double * par_sol = new double[size];
    for (int j=0; j<size; ++j) {
      par_sol[j] = sol[mem[j]];
    }
    // check cone feasibility
    double activity = 0.0;
    double sum_rest = 0.0;
    int start;
    if (type==OSI_QUAD) {
      start = 1;
    }
    else if (type==OSI_RQUAD) {
      start = 2;
    }
    else {
      std::cerr << "Unknown cone type!"  << std::endl;
      throw std::exception();
    }
    sum_rest = std::inner_product(par_sol+start, par_sol+size,
				  par_sol+start, 0.0);
    if (type==OSI_QUAD) {
      activity = par_sol[0] - sqrt(sum_rest);
    }
    else if (type==OSI_RQUAD) {
      activity = 2.0*par_sol[0]*par_sol[1] - sum_rest;
    }
    else {
      std::cerr << "Unknown cone type!"  << std::endl;
      throw std::exception();
    }
    if (activity<-CONE_EPS) {
      // current solution is infeasible to conic constraint i.
      feasible[i] = 0;
    }
    else {
      // solution is feasible for this cone, keep going with the next one
      feasible[i] = 1;
      continue;
    }
    // project sol to the cone
    if (type==OSI_QUAD) {
      point[mem[0]] = sqrt(sum_rest);
    }
    else if (type==OSI_RQUAD) {
      double p1 = par_sol[0];
      double p2 = par_sol[1];
      point[mem[0]] = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) - (-p1+p2)) / 2.0;
      point[mem[1]] = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) + (-p1+p2)) / 2.0;
    }
    else {
      std::cerr << "Unknown cone type!"  << std::endl;
      throw std::exception();
    }
    delete[] par_sol;
  }
}

// project given infeasible solution to an array of points on the cones
// of the problem. This function first projects the point using the
// method in the paper then creates random points around the projection.
// set feasible[i] to 0 if solution violates cone i, 1 otherwise.
// n is size of solution
// store projected points in point array.
void CglConicOA::project_random(int n, int num_cones, int const * cone_size,
			 OsiLorentzConeType const * cone_type,
			 int const * const * members, double const * sol,
			 double * const * point, int * const feasible,
			 int num_points) const {
  // copy solution to the point
  for (int i=0; i<num_points; ++i) {
    std::copy(sol, sol+n, point[i]);
  }
  // store projected point in point[0]
  double * p = point[0];
  // iterate over cones and project relevant part
  for (int i=0; i<num_cones; ++i) {
    // get data for cone i
    int size = cone_size[i];
    int const * mem = members[i];
    OsiLorentzConeType type = cone_type[i];
    // get relevent part of the solution
    double * par_sol = new double[size];
    for (int j=0; j<size; ++j) {
      par_sol[j] = sol[mem[j]];
    }
    // check cone feasibility
    double activity = 0.0;
    double sum_rest = 0.0;
    int start;
    if (type==OSI_QUAD) {
      start = 1;
    }
    else if (type==OSI_RQUAD) {
      start = 2;
    }
    else {
      std::cerr << "Unknown cone type!"  << std::endl;
      throw std::exception();
    }
    sum_rest = std::inner_product(par_sol+start, par_sol+size,
				  par_sol+start, 0.0);
    if (type==OSI_QUAD) {
      activity = par_sol[0] - sqrt(sum_rest);
    }
    else if (type==OSI_RQUAD) {
      activity = 2.0*par_sol[0]*par_sol[1] - sum_rest;
    }
    else {
      std::cerr << "Unknown cone type!"  << std::endl;
      throw std::exception();
    }
    if (activity<-CONE_EPS) {
      // current solution is infeasible to conic constraint i.
      feasible[i] = 0;
    }
    else {
      // solution is feasible for this cone, keep going with the next one
      feasible[i] = 1;
      continue;
    }
    // project sol to the cone
    if (type==OSI_QUAD) {
      p[mem[0]] = sqrt(sum_rest);
    }
    else if (type==OSI_RQUAD) {
      double p1 = par_sol[0];
      double p2 = par_sol[1];
      p[mem[0]] = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) - (-p1+p2)) / 2.0;
      p[mem[1]] = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) + (-p1+p2)) / 2.0;
    }
    else {
      std::cerr << "Unknown cone type!"  << std::endl;
      throw std::exception();
    }
    delete[] par_sol;
  }
  // create random points around point[0] for parts that are
  // infeasible to a cone.
  double eps = 0.0001;
  // for each cone
  for (int j=0; j<num_cones; ++j) {
    // check if cone j is feasible
    if (feasible[j]) {
      continue;
    }
    int start;
    if (cone_type[j]==OSI_QUAD) {
      start = 1;
    }
    else if (cone_type[j]==OSI_RQUAD) {
      start = 2;
    }
    else {
      std::cerr << "Unknown cone type!" << std::endl;
      throw std::exception();
    }
    // for each point
    for (int i=1; i<num_points; ++i) {
      for (int k=start; k<cone_size[j]; ++k) {
	int flip = rand()%2;
	int perc = rand()%100;
	if (flip==0) {
	  point[i][members[j][k]] = point[0][members[j][k]] + eps*perc*point[0][members[j][k]];
	}
	else {
	  point[i][members[j][k]] = point[0][members[j][k]] - eps*perc*point[0][members[j][k]];
	}
      }
      // compute leading points
      // todo(aykut) copying to a contiguous memory might speed things up here
      double rest = 0.0;
      for (int k=start; k<cone_size[j]; ++k) {
	rest += point[i][members[j][k]]*point[i][members[j][k]];
      }
      if (cone_type[j]==OSI_QUAD) {
	point[i][members[j][0]] = sqrt(rest);
      }
      else if (cone_type[j]==OSI_RQUAD) {
	double val = sqrt(rest/2.0);
	point[i][members[j][0]] = val;
	point[i][members[j][1]] = val;
      }
      else {
	std::cerr << "Unknown cone type!" << std::endl;
	throw std::exception();
      }
    }
  }
}

// project given infeasible solution to points on the cones
// of the problem. This function first projects the point using the
// method in the paper then creates more points closer to the projection.
// set feasible[i] to 0 if solution violates cone i, 1 otherwise.
// n is size of solution
// store projected points in point array.
//
// use trigonometric functions to generate points around sol on the boundry
// this function works for 3 dimensional Lorentz cones only for now (no
// rotated cones).
// we generate two points around projected point. these points are
// (x1, x1 cos(theta+epsilon), x1 sin(theta+epsilon))
// (x1, x1 cos(theta-epsilon), x1 sin(theta-epsilon))
//
// where theta is acos(x2/x1)
// at the end points will look like following
// point[0] -> (x1, x2, x3)
// point[1] -> (x1, x1 cos(theta+epsilon), x1 sin(theta+epsilon))
// point[2] -> (x1, x1 cos(theta-epsilon), x1 sin(theta-epsilon))
// point[3] -> (x1, x1 cos(theta+2epsilon), x1 sin(theta+2epsilon))
// point[4] -> (x1, x1 cos(theta-2epsilon), x1 sin(theta-2epsilon))
// point[5] -> (x1, x1 cos(theta+3epsilon), x1 sin(theta+3epsilon))
// point[6] -> (x1, x1 cos(theta-3epsilon), x1 sin(theta-3epsilon))
// .
// .
// .
//
void CglConicOA::project_trig(int n, int num_cones, int const * cone_size,
			 OsiLorentzConeType const * cone_type,
			 int const * const * members, double const * sol,
			 double * const * point, int * const feasible,
			 int num_points) const {
  // check cone sizes
  if (num_points>1) {
    for (int i=0; i<num_cones; ++i) {
      if (cone_size[i]!=3) {
	std::cerr << "This is implemented for cones of size 3 only." << std::endl;
	throw std::exception();
      }
    }
  }
  // copy solution to the point
  for (int i=0; i<num_points; ++i) {
    std::copy(sol, sol+n, point[i]);
  }
  // store projected point in point[0]
  double * p = point[0];
  // iterate over cones and project relevant part
  for (int i=0; i<num_cones; ++i) {
    // get data for cone i
    int size = cone_size[i];
    int const * mem = members[i];
    OsiLorentzConeType type = cone_type[i];
    // get relevent part of the solution
    double * par_sol = new double[size];
    for (int j=0; j<size; ++j) {
      par_sol[j] = sol[mem[j]];
    }
    // check cone feasibility
    double activity = 0.0;
    double sum_rest = 0.0;
    int start;
    if (type==OSI_QUAD) {
      start = 1;
    }
    else if (type==OSI_RQUAD) {
      start = 2;
    }
    else {
      std::cerr << "Unknown cone type!"  << std::endl;
      throw std::exception();
    }
    sum_rest = std::inner_product(par_sol+start, par_sol+size,
				  par_sol+start, 0.0);
    if (type==OSI_QUAD) {
      activity = par_sol[0] - sqrt(sum_rest);
    }
    else if (type==OSI_RQUAD) {
      activity = 2.0*par_sol[0]*par_sol[1] - sum_rest;
    }
    else {
      std::cerr << "Unknown cone type!"  << std::endl;
      throw std::exception();
    }
    if (activity<-CONE_EPS) {
      // current solution is infeasible to conic constraint i.
      feasible[i] = 0;
    }
    else {
      // solution is feasible for this cone, keep going with the next one
      feasible[i] = 1;
      continue;
    }
    // project sol to the cone
    if (type==OSI_QUAD) {
      p[mem[0]] = sqrt(sum_rest);
    }
    else if (type==OSI_RQUAD) {
      double p1 = par_sol[0];
      double p2 = par_sol[1];
      p[mem[0]] = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) - (-p1+p2)) / 2.0;
      p[mem[1]] = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) + (-p1+p2)) / 2.0;
    }
    else {
      std::cerr << "Unknown cone type!"  << std::endl;
      throw std::exception();
    }
    delete[] par_sol;
  }
  // create points around point[0] for parts that are infeasible to a cone.
  double eps = 0.001;
  // for each cone
  for (int j=0; j<num_cones; ++j) {
    // check if cone j is feasible
    if (feasible[j]) {
      continue;
    }
    if (cone_type[j]!=OSI_QUAD) {
      std::cerr << "Not implemented for rotated cones!" << std::endl;
      throw std::exception();
    }
    // print point[0]
    //std::cout << "point 0" << std::endl;
    //std::cout << point[0][members[j][0]] << " " << point[0][members[j][1]] << " " << point[0][members[j][2]] << std::endl;
    double x1 = point[0][members[j][0]];
    double x2 = point[0][members[j][1]];
    double x3 = point[0][members[j][2]];
    double theta = acos(x2/x1);
    //std::cout << "x2 " << x2 << " x3 " << x3 << " theta " << theta << std::endl;
    // for each point
    for (int i=1; i<num_points; ++i) {
      // point[0] -> (x1, x2, x3)
      // point[1] -> (x1, x1 cos(theta+epsilon), x1 sin(theta+epsilon))
      // point[2] -> (x1, x1 cos(theta-epsilon), x1 sin(theta-epsilon))
      // point[3] -> (x1, x1 cos(theta+2epsilon), x1 sin(theta+2epsilon))
      // point[4] -> (x1, x1 cos(theta-2epsilon), x1 sin(theta-2epsilon))
      int coef = (i-1)/2 +1;
      int sign = i%2==0 ? -1: +1;
      point[i][members[j][0]] = x1;
      point[i][members[j][1]] = x1*cos(theta+sign*coef*eps);
      if (x3<0.0) {
	point[i][members[j][2]] = -1.0*x1*sin(theta+sign*coef*eps);
      }
      else {
	point[i][members[j][2]] = x1*sin(theta+sign*coef*eps);
      }
      //std::cout << "point " << i << std::endl;
      //std::cout << point[i][members[j][0]] << " " << point[i][members[j][1]] << " " << point[i][members[j][2]] << std::endl;
    }
  }
}

void CglConicOA::generate_support(int size,
				  OsiLorentzConeType type,
				  int const * members,
				  double const * point,
				  OsiRowCut * rc) const {
  if (type==OSI_QUAD) {
    generate_support_lorentz(size, members, point, rc);
  }
  else {
    generate_support_rotated_lorentz(size, members, point, rc);
  }
}

void CglConicOA::generate_support_lorentz(int size,
					  int const * members,
					  double const * point,
					  OsiRowCut * rc) const {
  // todo(aykut) check whether p is 0,
  // the method will return to 0 coef if p is 0.
  double * coef = new double[size];
  // cone is in canonical form
  for (int j=1; j<size; ++j) {
    if ((point[j]<CONE_EPS) && (point[j]>-CONE_EPS)) {
      coef[j] = 0.0;
    }
    else {
      coef[j] = 2.0*point[j];
    }
  }
  coef[0] = -2.0*point[0];
  // insert constraint (coef,0) to constraint pool
  rc->setRow(size, members, coef);
  // todo(aykut): fix setting infinity
  rc->setLb(-1e80);
  rc->setUb(0.0);
  delete[] coef;
}

void CglConicOA::generate_support_rotated_lorentz(int size,
						 int const * members,
						 double const * point,
						 OsiRowCut * rc) const {
  // map point from RLORENTZ space to LORENTZ space, find the projection on LORENTZ,
  // project this point to RLORENTZ and generate cut
  double * coef = new double[size];
  // cone is a rotated cone
  // from point we move along [2point_2 2point_1 0 ... 0] until we hit
  // boundary. Then from this point in boundry we generate coef.
  // first compute u, step length
  // generate cut from xbar
  coef[0] = -2.0*point[1];
  coef[1] = -2.0*point[0];
  for (int i=2; i<size; ++i) {
    coef[i] = 2.0*point[i];
  }
  // insert constraint (coef,0) to constraint pool
  rc->setRow(size, members, coef);
  // todo(aykut): fix setting infinity
  rc->setLb(-1e80);
  rc->setUb(0.0);
  delete[] coef;
}
