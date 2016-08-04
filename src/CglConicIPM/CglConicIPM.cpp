// Name:     CglConicIPM.cpp
// Author:   Aykut Bulut
//           Lehigh University
//           email: aykut@lehigh.edu, aykutblt@gmail.com
//-----------------------------------------------------------------------------
// Copyright (C) 2015, Lehigh University.  All Rights Reserved.

#include "CglConicIPM.hpp"
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
//#include <random>

#if defined(__OSI_MOSEK__)
  // use mosek as IPM solver
  #include <OsiMosekSolverInterface.hpp>
  //#define IPM_SOLVER OsiMosekSolverInterface
  typedef OsiMosekSolverInterface IPM_SOLVER;
#elif defined(__OSI_CPLEX__)
  // use cplex as IPM solver
  #include <OsiCplexSolverInterface.hpp>
  typedef OsiCplexSolverInterface IPM_SOLVER;
#elif defined(__OSI_IPOPT__)
  // use ipopt
  #include <OsiIpoptSolverInterface.hpp>
  typedef OsiIpoptSolverInterface IPM_SOLVER;
#endif

#define CONE_EPS 1e-6
#define COEF_EPS 1e-5

CglConicIPM::CglConicIPM()
  : param_(0) {
  param_ = new CglConicIPMParam();
  solver_ = 0;
  srand(0);
}

// copy constructor
CglConicIPM::CglConicIPM(const CglConicIPM & other) {
  // copy param_
  param_ = new CglConicIPMParam(*(other.getParam()));
  OsiConicSolverInterface * os = other.getSolver();
  if (os!=0) {
    /// notes(aykut) solver interface should provide a constructor
    /// that copies from a OsiConicSolverInterface object.
    /// todo(aykut) mosek does not have this copy constructor, add
    /// it.
    solver_ = new IPM_SOLVER(os);
    solver_->setHintParam(OsiDoReducePrint,true,OsiHintDo, 0);
  }
  else {
    solver_ = 0;
  }
  srand(0);
}

// copy assignment operator
CglConicIPM & CglConicIPM::operator=(const CglConicIPM & rhs) {
  // copy param_
  param_ = new CglConicIPMParam(*(rhs.getParam()));
  // copy solver
  OsiConicSolverInterface * os = rhs.getSolver();
  if (os!=0) {
    solver_ = new IPM_SOLVER(os);
    solver_->setHintParam(OsiDoReducePrint,true,OsiHintDo, 0);
  }
  else {
    solver_ = 0;
  }
  return *this;
}

CglConicIPM::~CglConicIPM() {
  if (param_)
    delete param_;
  if (solver_)
    delete solver_;
}

void CglConicIPM::setParam(CglConicIPMParam const & param) {
  param_ = new CglConicIPMParam(param);
}

// generate outer approximating hyperplanes for conic constraints
// todo(aykut): approximates Lorentz cones only for now.
void CglConicIPM::generateCuts(OsiConicSolverInterface const & si,
                               OsiCuts & cuts,
                               const CglTreeInfo info) {
  // this function does nothing.
}

// generate cuts for a linear solver interface
void CglConicIPM::generateCuts(OsiSolverInterface const & si, OsiCuts & cuts,
                  int num_cones, OsiLorentzConeType const * cone_type,
                               int const * cone_size, int const * const * members,
                               int num_points) {
  method3(si, cuts, num_cones, cone_type, cone_size, members, num_points);
}

// optimize given SOCO. Find points around optimal solution and add cuts
// on these points.
void CglConicIPM::method1(OsiSolverInterface const & si, OsiCuts & cuts,
                          int num_cones, OsiLorentzConeType const * cone_type,
                          int const * cone_size, int const * const * members,
                          int num_points) {
  // si solution
  double const * sol = si.getColSolution();
  // unboundedness direction, if si is unbounded.
  double * u_dir = 0;
  // if linear problem is unbounded restrict unboundedness direction
  if (si.isProvenDualInfeasible()) {
    // check if primal is infeasible, then the problem is infeasible
    if (si.isProvenPrimalInfeasible()) {
      // Both LP primal and dual is infeasible, conic problem is infeasible
      std::cerr << "CglConic: Conic problem is infeasible."
                << std::endl;
    }
    // get one ray
    // todo(aykut) for now we get only one ray
    std::vector<double*> rays = si.getPrimalRays(1);
    double const * vec = 0;
    if (!rays.empty() and rays[0]!=0) {
      vec = rays[0];
    }
    else {
      std::cerr << "CglConic: Warning! "
                << "LP relaxation is unbounded but solver did not return a "
        "direction of unboundedness." << std::endl
                << "CglConic: Trying to generate supports using objective "
        "function coefficients..." << std::endl;
      vec = si.getObjCoefficients();
    }
    int num_cols = si.getNumCols();
    u_dir = new double[num_cols];
    std::copy(vec, vec+num_cols, u_dir);
    // delete all rays not just first one.
    if (!rays.empty()) {
      for(int i=0; i<rays.size(); ++i)
        delete[] rays[i];
      rays.clear();
    }
  }
  // if there is an unboundedness direction replace sol with it
  if (u_dir) {
    sol = u_dir;
  }
  int feasible = 1;
  for (int i=0; i<num_cones; ++i) {
    double activity;
    double * par_sol = new double[cone_size[i]];
    for (int j=0; j<cone_size[i]; ++j) {
      par_sol[j] = sol[members[i][j]];
    }
    if (cone_type[i]==OSI_QUAD) {
      activity = par_sol[0] - sqrt(std::inner_product(par_sol+1,
                                par_sol+cone_size[i], par_sol+1, 0.0));
    }
    else if (cone_type[i]==OSI_RQUAD) {
      activity = 2*par_sol[0]*par_sol[1] - std::inner_product(par_sol+2,
                                par_sol+cone_size[i], par_sol+2, 0.0);
    }
    else {
      std::cerr << "Unknown cone." << std::endl;
      throw std::exception();
    }
    if (activity>-1e-5) {
      // solution is already feasible
      delete[] par_sol;
      continue;
    }
    else {
      feasible = 0;
      delete[] par_sol;
      break;
    }
  }
  if (feasible) {
    return;
  }
  // free solver
  if (solver_) {
    delete solver_;
  }
  //solver_ = new OsiMosekSolverInterface();
  solver_ = new IPM_SOLVER();
  solver_->setHintParam(OsiDoReducePrint,true,OsiHintDo, 0);

  // load data to solver
  CoinPackedMatrix const * matrix = si.getMatrixByCol();
  double const * rowlb = si.getRowLower();
  double const * rowub = si.getRowUpper();
  double const * collb = si.getColLower();
  double const * colub = si.getColUpper();
  double const * obj = si.getObjCoefficients();
  solver_->loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
  // add conic constraints
  for (int i=0; i<num_cones; ++i) {
    solver_->addConicConstraint(cone_type[i], cone_size[i], members[i]);
  }
  // solve problem
  solver_->initialSolve();
  if (solver_->isProvenPrimalInfeasible() || solver_->isProvenDualInfeasible()) {
    // problem is infeasible, this means we can fathom the node in BB
    //std::cerr << "Cut generation problem is infeasible!" << std::endl;
    // add infinity as a lower bound for objective
    double infinity = solver_->getInfinity();
    int num_cols = solver_->getNumCols();
    int * ind = new int[num_cols];
    double * val = new double[num_cols];
    int s = 0;
    for (int i=0; i<num_cols; ++i) {
      if (obj[i]!=0.0) {
        ind[s] = i;
        val[s] = obj[i];
        s++;
      }
    }
    OsiRowCut * rc = new OsiRowCut();
    // insert constraint (coef,0) to constraint pool
    rc->setRow(s, ind, val);
    // todo(aykut): fix setting infinity
    rc->setLb(infinity);
    cuts.insert(rc);
    delete[] ind;
    delete[] val;
    delete rc;
    return;
  }
  else if (!(solver_->isProvenOptimal())) {
    std::cerr << "Cut problem could not be solved!" << std::endl;
    std::cerr << "No cuts generated!" << std::endl;
    return;
  }
  // get solution
  double const * ipm_sol = solver_->getColSolution();
  // create random points
  int num_cols = solver_->getNumCols();
  // create random points around solution, on the cone
  double ** points = new double*[num_points];
  for (int i=0; i<num_points; ++i) {
    // allocate memory for each point
    points[i] = new double[num_cols]();
  }
  // create random points around sol, on the cone
  create_rand_points(num_cols, ipm_sol, num_cones, cone_type, cone_size,
                     members, points, num_points);
  // add cuts from the points that are on the cones.
  for (int j=0; j<num_points; ++j) {
    add_cuts2(num_cols, points[j], num_cones, cone_type, cone_size,
             members, cuts);
  }
  // add cuts that actually cuts the point.
  // add cut if it actually cuts the lp solution
  int num_cuts = cuts.sizeRowCuts();
  std::vector<int> cuts_to_del(num_cuts, -1);
  for (int i=0; i<num_cuts; ++i) {
    OsiRowCut const * rc = cuts.rowCutPtr(i);
    double infeasibility = rc->violated(sol);
    if (infeasibility<1e-5) {
      // rc does not cut sol
      cuts_to_del.push_back(i);
    }
  }
  // remove cuts
  std::vector<int>::reverse_iterator rit;
  for (rit=cuts_to_del.rbegin(); rit!=cuts_to_del.rend(); ++rit) {
    if (*rit!=-1) {
      cuts.eraseRowCut(*rit);
    }
  }
  // free memory allocated to points
  for (int i=0; i<num_points; ++i) {
    delete[] points[i];
  }
  delete[] points;
  if (u_dir) {
    delete[] u_dir;
  }
}

// find maximally violated cut by solving the following
// min ||x-xbar||
// s.t. x in X
// This can be reformulated as
// min z
// s.t. y -x = -xbar
//      x in X
//      (z,y) in L^(n+1)
// a point is in (x,y,z) form.
// where xbar is approximating LP solution and X is feasible region of SOCO.
// Then add cut on the optimal solution of this problem.
void CglConicIPM::method2(OsiSolverInterface const & si, OsiCuts & cuts,
                          int num_cones, OsiLorentzConeType const * cone_type,
                          int const * cone_size, int const * const * members,
                          int num_points) {
  // free solver
  if (solver_) {
    delete solver_;
  }
  solver_ = new IPM_SOLVER();
  solver_->setHintParam(OsiDoReducePrint,true,OsiHintDo, 0);

  double infinity = solver_->getInfinity();
  // add linear constraints and bounds
  int num_cols = si.getNumCols();
  int num_rows = si.getNumRows();
  double const * sol = si.getColSolution();
  // if solution is conic feasible return
  int feasible = 1;
  for (int i=0; i<num_cones; ++i) {
    double activity;
    double * par_sol = new double[cone_size[i]];
    for (int j=0; j<cone_size[i]; ++j) {
      par_sol[j] = sol[members[i][j]];
    }
    if (cone_type[i]==OSI_QUAD) {
      activity = par_sol[0] - sqrt(std::inner_product(par_sol+1,
                                par_sol+cone_size[i], par_sol+1, 0.0));
    }
    else if (cone_type[i]==OSI_RQUAD) {
      activity = 2*par_sol[0]*par_sol[1] - std::inner_product(par_sol+2,
                                par_sol+cone_size[i], par_sol+2, 0.0);
    }
    else {
      std::cerr << "Unknown cone." << std::endl;
      throw std::exception();
    }
    if (activity>-1e-5) {
      // solution is already feasible
      continue;
    }
    else {
      feasible = 0;
      break;
    }
    delete[] par_sol;
  }
  if (feasible) {
    return;
  }
  double * neg_sol = new double[num_cols];
  for (int i=0; i<num_cols; ++i) {
    neg_sol[i] = -sol[i];
  }
  // == compute problem matrix
  // ==== get matrix of input problem
  CoinPackedMatrix const * cm = si.getMatrixByRow();
  // ==== add additional rows of projection problem
  CoinPackedMatrix * pm = new CoinPackedMatrix(*cm);
  CoinPackedVectorBase ** add_rows = new CoinPackedVectorBase*[num_cols+1];
  for (int i=0; i<num_cols; ++i) {
    int ind[2] = {i, num_cols+i};
    double val[2] = {-1.0, 1.0};
    add_rows[i] = new CoinPackedVector(2, ind, val);
  }
  // ==== add additional rows
  pm->appendRows(num_cols, add_rows);
  // == compute row lower bound
  double const * si_rowlb = si.getRowLower();
  double * rowlb = new double[num_rows+num_cols];
  std::copy(si_rowlb, si_rowlb+num_rows, rowlb);
  std::copy(neg_sol, neg_sol+num_cols, rowlb+num_rows);
  // == compute row upper bound
  double const * si_rowub = si.getRowUpper();
  double * rowub = new double[num_rows+num_cols];
  std::copy(si_rowub, si_rowub+num_rows, rowub);
  std::copy(neg_sol, neg_sol+num_cols, rowub+num_rows);
  // == compute col lower bound
  double const * si_collb = si.getColLower();
  double * collb = new double[2*num_cols+1];
  std::copy(si_collb, si_collb+num_cols, collb);
  std::fill(collb+num_cols, collb+2*num_cols+1, -infinity);
  // == compute col upper bound
  double const * si_colub = si.getColUpper();
  double * colub = new double[2*num_cols+1];
  std::copy(si_colub, si_colub+num_cols, colub);
  std::fill(colub+num_cols, colub+2*num_cols+1, infinity);
  // == compute objective coef
  double * obj = new double[2*num_cols+1]();
  obj[2*num_cols] = 1.0;
  // == load data to solver
  solver_->loadProblem(*pm, collb, colub, obj, rowlb, rowub);
  // add dummy constraint z - y_1 >= 0;
  // int dum_ind[2] = {2*num_cols, num_cols};
  // double dum_val[2] = {1.0, -1.0};
  // solver_->addRow(2, dum_ind, dum_val, 0.0, infinity);
  // add column z
  solver_->addCol(0, 0, 0, 0.0, infinity, 1.0);
  // supress messages
  solver_->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  // end of dummy row add
  delete pm;
  for (int i=0; i<num_cols; ++i) {
    delete add_rows[i];
  }
  delete[] add_rows;
  delete[] rowlb;
  delete[] rowub;
  delete[] collb;
  delete[] colub;
  delete[] obj;
  delete[] neg_sol;
  // add conic constraints
  for (int i=0; i<num_cones; ++i) {
    solver_->addConicConstraint(cone_type[i], cone_size[i], members[i]);
  }
  // add (z,y) in L^{n+1}
  int * new_cone = new int[num_cols+1];
  for (int i=0; i<num_cols; ++i) {
    new_cone[i+1] = num_cols+i;
  }
  new_cone[0] = 2*num_cols;
  solver_->addConicConstraint(OSI_QUAD, num_cols+1, new_cone);
  delete[] new_cone;
  // solve problem
  solver_->initialSolve();
  // get solution and add cuts
  double const * proj_sol = solver_->getColSolution();
  // create random points around solution, on the cone
  double ** points = new double*[num_points];
  for (int i=0; i<num_points; ++i) {
    // allocate memory for each point
    points[i] = new double[num_cols]();
  }
  // create random points around sol, on the cone
  create_rand_points(num_cols, proj_sol, num_cones, cone_type, cone_size,
                     members, points, num_points);
  // add cuts from the points that are on the cones.
  for (int j=0; j<num_points; ++j) {
    add_cuts2(num_cols, points[j], num_cones, cone_type, cone_size,
             members, cuts);
  }
  // free memory allocated to points
  for (int i=0; i<num_points; ++i) {
    delete[] points[i];
  }
  delete[] points;
}


void CglConicIPM::method3(OsiSolverInterface const & si, OsiCuts & cuts,
             int num_cones, OsiLorentzConeType const * cone_type,
             int const * cone_size, int const * const * members,
             int num_points) {
  // si solution
  double const * sol = si.getColSolution();
  // unboundedness direction, if si is unbounded.
  double * u_dir = 0;
  // if linear problem is unbounded restrict unboundedness direction
  if (si.isProvenDualInfeasible()) {
    // check if primal is infeasible, then the problem is infeasible
    if (si.isProvenPrimalInfeasible()) {
      // Both LP primal and dual is infeasible, conic problem is infeasible
      std::cerr << "CglConic: Conic problem is infeasible."
                << std::endl;
    }
    // get one ray
    // todo(aykut) for now we get only one ray
    std::vector<double*> rays = si.getPrimalRays(1);
    double const * vec = 0;
    if (!rays.empty() and rays[0]!=0) {
      vec = rays[0];
    }
    else {
      std::cerr << "CglConic: Warning! "
                << "LP relaxation is unbounded but solver did not return a "
        "direction of unboundedness." << std::endl
                << "CglConic: Trying to generate supports using objective "
        "function coefficients..." << std::endl;
      vec = si.getObjCoefficients();
    }
    int num_cols = si.getNumCols();
    u_dir = new double[num_cols];
    std::copy(vec, vec+num_cols, u_dir);
    // delete all rays not just first one.
    if (!rays.empty()) {
      for(int i=0; i<rays.size(); ++i)
        delete[] rays[i];
      rays.clear();
    }
  }
  // if there is an unboundedness direction replace sol with it
  if (u_dir) {
    sol = u_dir;
  }
  int feasible = 1;
  for (int i=0; i<num_cones; ++i) {
    double activity;
    double * par_sol = new double[cone_size[i]];
    for (int j=0; j<cone_size[i]; ++j) {
      par_sol[j] = sol[members[i][j]];
    }
    if (cone_type[i]==OSI_QUAD) {
      activity = par_sol[0] - sqrt(std::inner_product(par_sol+1,
                                par_sol+cone_size[i], par_sol+1, 0.0));
    }
    else if (cone_type[i]==OSI_RQUAD) {
      activity = 2*par_sol[0]*par_sol[1] - std::inner_product(par_sol+2,
                                par_sol+cone_size[i], par_sol+2, 0.0);
    }
    else {
      std::cerr << "Unknown cone." << std::endl;
      throw std::exception();
    }
    if (activity>-1e-5) {
      // solution is already feasible
      delete[] par_sol;
      continue;
    }
    else {
      feasible = 0;
      delete[] par_sol;
      break;
    }
  }
  if (feasible) {
    return;
  }
  // free solver
  if (solver_) {
    delete solver_;
  }
  solver_ = new IPM_SOLVER();
  solver_->setHintParam(OsiDoReducePrint,true,OsiHintDo, 0);

  // load data to solver
  CoinPackedMatrix const * matrix = si.getMatrixByCol();
  double const * rowlb = si.getRowLower();
  double const * rowub = si.getRowUpper();
  double const * collb = si.getColLower();
  double const * colub = si.getColUpper();
  double const * obj = si.getObjCoefficients();
  double * new_collb = new double[si.getNumCols()];
  double * new_colub = new double[si.getNumCols()];
  char const * col_type = si.getColType();
  std::copy(collb, collb+si.getNumCols(), new_collb);
  std::copy(colub, colub+si.getNumCols(), new_colub);
  for (int i=0; i<si.getNumCols(); ++i) {
    if (col_type[i]=='\000') {
      continue;
    }
    else {
      new_collb[i] = sol[i];
      new_colub[i] = sol[i];
    }
  }
  solver_->loadProblem(*matrix, new_collb, new_colub, obj, rowlb, rowub);
  // add conic constraints
  for (int i=0; i<num_cones; ++i) {
    solver_->addConicConstraint(cone_type[i], cone_size[i], members[i]);
  }
  // solve problem
  solver_->initialSolve();
  if (solver_->isProvenPrimalInfeasible() || solver_->isProvenDualInfeasible()) {
    // problem is infeasible, do not generate any cuts.
    return;
  }
  else if (!(solver_->isProvenOptimal())) {
    std::cerr << "Cut problem could not be solved!" << std::endl;
    std::cerr << "No cuts generated!" << std::endl;
    return;
  }
  // get solution
  double const * ipm_sol = solver_->getColSolution();
  // create random points
  int num_cols = solver_->getNumCols();
  // create random points around solution, on the cone
  double ** points = new double*[num_points];
  for (int i=0; i<num_points; ++i) {
    // allocate memory for each point
    points[i] = new double[num_cols]();
  }
  // create random points around sol, on the cone
  create_rand_points(num_cols, ipm_sol, num_cones, cone_type, cone_size,
                     members, points, num_points);
  // add cuts from the points that are on the cones.
  for (int j=0; j<num_points; ++j) {
    add_cuts2(num_cols, points[j], num_cones, cone_type, cone_size,
             members, cuts);
  }
  // add cuts that actually cuts the point.
  // add cut if it actually cuts the lp solution
  int num_cuts = cuts.sizeRowCuts();
  std::vector<int> cuts_to_del(num_cuts, -1);
  for (int i=0; i<num_cuts; ++i) {
    OsiRowCut const * rc = cuts.rowCutPtr(i);
    double infeasibility = rc->violated(sol);
    if (infeasibility<1e-5) {
      // rc does not cut sol
      cuts_to_del.push_back(i);
    }
  }
  // remove cuts
  std::vector<int>::reverse_iterator rit;
  for (rit=cuts_to_del.rbegin(); rit!=cuts_to_del.rend(); ++rit) {
    if (*rit!=-1) {
      cuts.eraseRowCut(*rit);
    }
  }
  // free memory allocated to points
  for (int i=0; i<num_points; ++i) {
    delete[] points[i];
  }
  delete[] points;
  if (u_dir) {
    delete[] u_dir;
  }
}


// add cuts from the point on the cone
// we know that points are on the cone.
void CglConicIPM::add_cuts2(int num_cols, double const * point, int num_cones,
                            OsiLorentzConeType const * cone_type,
                            int const * cone_size,
                            int const * const * members, OsiCuts & cuts) {
  double * par_point;
  int start;
  for (int i=0; i<num_cones; ++i) {
    par_point = new double[cone_size[i]];
    for(int j=0; j<cone_size[i]; ++j) {
      par_point[j] = point[members[i][j]];
    }
    double * p = par_point;
    // skip the point if it is close to 0.
    double sq = std::inner_product(p, p+cone_size[i], p, 0.0);
    if (sq<1e-5) {
      delete[] par_point;
      continue;
    }
    if (cone_type[i]==OSI_QUAD) {
      start = 1;
    }
    else if (cone_type[i]==OSI_RQUAD) {
      start = 2;
    }
    else {
      std::cerr << "Unknown cone type." << std::endl;
      throw std::exception();
    }
    // check whether p is 0
    // double norm;
    // norm = std::inner_product(p, p+cone_size[i], p, 0.0);
    // if (norm<1e-5) {
    //   continue;
    // }
    // the method will return to 0 coef if p is o.
    double sum_sq = std::inner_product(p+start, p+cone_size[i], p+start, 0.0);
    double * coef = new double[cone_size[i]];
    // cone is in canonical form
    for (int j=start; j<cone_size[i]; ++j) {
      coef[j] = 2.0*p[j];
    }
    if (cone_type[i]==OSI_QUAD) {
      coef[0] = -2.0*p[0];
    }
    else if (cone_type[i]==OSI_RQUAD) {
      //std::cout << "Not implemented yet." << std::endl;
      //throw std::exception();
      double p1 = p[0];
      double p2 = p[1];
      double x1 = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_sq) - (-p1+p2)) / 2.0;
      double x2 = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_sq) + (-p1+p2)) / 2.0;
      // generate cut from xbar
      coef[0] = -2.0*x2;
      coef[1] = -2.0*x1;
    }
    else {
      std::cerr << "Unknown cone type." << std::endl;
      throw std::exception();
    }
    OsiRowCut * rc = new OsiRowCut();
    // insert constraint (coef,0) to constraint pool
    rc->setRow(cone_size[i], members[i], coef);
    // todo(aykut): fix setting infinity
    rc->setLb(-1e80);
    rc->setUb(0.0);
    cuts.insert(rc);
    delete[] coef;
    delete[] par_point;
    delete rc;
  }
}


void CglConicIPM::create_rand_points(int num_cols, double const * sol,
                                     int num_cones,
                                     OsiLorentzConeType const * cone_type,
                                     int const * cone_size,
                                     int const * const * members,
                                     double ** points,
                                     int num_points) const {
  std::copy(sol, sol+num_cols, points[0]);
  for (int i=1; i<num_points; ++i) {
    // generate a random point around sol
    create_rand_point2(num_cols, sol, num_cones, cone_type, cone_size,
                       members, points[i]);
  }
}


// create a random point around par_sol that is on the tight cones.
void CglConicIPM::create_rand_point2(int num_cols, double const * sol,
                                     int num_cones,
                                     OsiLorentzConeType const * cone_type,
                                     int const * cone_size,
                                     int const * const * members,
                                     double * point) const {
  double * par_sol = 0;
  double * par_point = 0;
  // check whether cone is tight.
  for (int i=0; i<num_cones; ++i) {
    par_sol = new double[cone_size[i]];
    // store partial solution
    for (int j=0; j<cone_size[i]; ++j) {
      par_sol[j] = sol[members[i][j]];
    }
    // generate point on a tight cone
    par_point = new double[cone_size[i]];
    create_rand_point3(cone_size[i], par_sol, cone_type[i], members[i],
                       par_point);
    // store par_point to points[i]
    for (int j=0; j<cone_size[i]; ++j) {
      point[members[i][j]] = par_point[j];
    }
    delete[] par_point;
    delete[] par_sol;
  }
}

void CglConicIPM::create_rand_point3(int cone_size, double const * par_sol,
                                     OsiLorentzConeType cone_type,
                                     int const * members,
                                     double * par_point) const {
  double eps = 1e-1;
  int rand_sign;
  double rand_number;
  int start;
  double sum_sq;
  if (cone_type==OSI_QUAD) {
    start = 1;
  }
  else if (cone_type==OSI_RQUAD) {
    start = 2;
  }
  else {
    std::cerr << "Unknown cone type." << std::endl;
    throw std::exception();
  }
  for (int i=0; i<cone_size; ++i) {
    rand_sign = rand()%2;
    rand_number = eps*(rand()/double(RAND_MAX));
    if (rand_sign==0) {
      par_point[i] = par_sol[i] + rand_number;
    }
    else {
      par_point[i] = par_sol[i] - rand_number;
    }
  }
  sum_sq = std::inner_product(par_point+start, par_point+cone_size,
                              par_point+start, 0.0);
  if (cone_type==OSI_QUAD) {
    par_point[0] = sqrt(sum_sq);
  }
  else if (cone_type==OSI_RQUAD) {
    double lead = sqrt(sum_sq)/sqrt(2.0);
    par_point[0] = lead;
    par_point[1] = lead;
  }
  else {
    std::cerr << "Unknown cone type." << std::endl;
    throw std::exception();
  }
  // scale point so that it is on the unit circle.
  if (cone_type==OSI_QUAD) {
    // scale such that par_point[0] is par_sol[0]
    double scale = par_sol[0]/par_point[0];
    for (int i=0; i<cone_size; ++i) {
      par_point[i] = scale*par_point[i];
    }
  }
  else if (cone_type==OSI_RQUAD) {
    double scale = (par_sol[0]*par_sol[1])/(par_point[0]*par_point[1]);
    for (int i=0; i<cone_size; ++i) {
      par_point[i] = scale*par_point[i];
    }
  }
  else {
    std::cerr << "Unknown cone type." << std::endl;
    throw std::exception();
  }
}

void CglConicIPM::add_cuts(int size, double const * sol, int num_cones,
                           OsiLorentzConeType const * cone_type,
                           int const * cone_size, int const * const * members,
                           OsiCuts & cuts) const {
  for (int i=0; i<num_cones; ++i) {
    OsiRowCut * rc = new OsiRowCut();
    double feas = generate_support(cone_size[i], cone_type[i], members[i], sol, rc);
    if (!feas) {
      cuts.insert(rc);
    }
    else {
      delete rc;
    }
  }
}

int CglConicIPM::generate_support(int size,
                                 OsiLorentzConeType type,
                                 int const * members,
                                 double const * sol,
                                 OsiRowCut * rc) const {
  int feas;
  if (type==OSI_QUAD) {
    feas = generate_support_lorentz(size, members, sol, rc);
  }
  else {
    feas = generate_support_rotated_lorentz(size,
                                            members, sol, rc);
  }
  return feas;
}

int CglConicIPM::generate_support_lorentz(int size,
                                 int const * members,
                                 double const * sol,
                                 OsiRowCut * rc) const {
  int feas;
  double * par_point = new double[size];
  for(int j=0; j<size; ++j) {
    par_point[j] = sol[members[j]];
  }
  double * p = par_point;
  // check the point
  // set variables almost 0 to 0.
  for (int i=0; i<size; ++i) {
    if ((p[i]<CONE_EPS) && (p[i]>-CONE_EPS)) {
      p[i] = 0.0;
    }
  }
  // todo(aykut) check whether p is 0,
  // the method will return to 0 coef if p is o.
  double term1;
  double term2;
  term1 = p[0];
  term2 = std::inner_product(p+1, p+size, p+1, 0.0);
  term2 = sqrt(term2);
  double activity = term1-term2;
  if (activity<-CONE_EPS) {
    // current solution is infeasible to conic constraint i.
    double * coef = new double[size];
    double x1 = term2;
    // cone is in canonical form
    for (int j=1; j<size; ++j) {
      coef[j] = 2.0*p[j];
    }
    coef[0] = -2.0*x1;
    // insert constraint (coef,0) to constraint pool
    rc->setRow(size, members, coef);
    // todo(aykut): fix setting infinity
    rc->setLb(-1e80);
    rc->setUb(0.0);
    delete[] coef;
    feas = 0;
  }
  else {
    feas = 1;
  }
  delete[] par_point;
  return feas;
}

int CglConicIPM::generate_support_rotated_lorentz(int size,
                                                 int const * members,
                                                 double const * sol,
                                                 OsiRowCut * rc) const {
  int feas;
  double activity;
  double * par_point = new double[size];
  for(int j=0; j<size; ++j) {
    par_point[j] = sol[members[j]];
  }
  double * p = par_point;
  // check feasibility of the given solution
  double sum_rest = 0.0;
  sum_rest = std::inner_product(p+2, p+size, p+2, 0.0);
  activity = 2.0*p[0]*p[1]-sum_rest;
  if (activity<-CONE_EPS) {
    //  at the end, set coef and lhs
    // map point from RLORENTZ space to LORENTZ space, find the projection on LORENTZ,
    // project this point to RLORENTZ and generate cut
    //sum_rest = std::inner_product(p+2, p+size, p+2, 0.0);
    double * coef = new double[size];
    double x1 = 0.0;
    double x2 = 0.0;
    // cone is a rotated cone
    // from point we move along [2point_2 2point_1 0 ... 0] until we hit
    // boundary. Then from this point in boundry we generate coef.
    // first compute u, step length
    double p1 = p[0];
    double p2 = p[1];
    p1 = p[0];
    p2 = p[1];
    x1 = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) - (-p1+p2)) / 2.0;
    x2 = (sqrt((-p1+p2)*(-p1+p2)+2.0*sum_rest) + (-p1+p2)) / 2.0;
    // generate cut from xbar
    coef[0] = -2.0*x2;
    coef[1] = -2.0*x1;
    for (int i=2; i<size; ++i) {
      coef[i] = 2.0*p[i];
    }
    // insert constraint (coef,0) to constraint pool
    rc->setRow(size, members, coef);
    // todo(aykut): fix setting infinity
    rc->setLb(-1e80);
    rc->setUb(0.0);
    delete[] coef;
    feas = 0;
  }
  else {
    feas = 1;
  }
  delete[] p;
  return feas;
}


/// Clone
CglConicCutGenerator * CglConicIPM::clone() const {
  CglConicIPM * new_cutg = new CglConicIPM(*this);
  return new_cutg;
}

/// Create C++ lines to get to current state
std::string CglConicIPM::generateCpp( FILE * fp) {
  std::cerr << "GenerateCpp is not implemented yet!" << std::endl;
  throw std::exception();
  return std::string();
}

OsiConicSolverInterface * CglConicIPM::getSolver() const {
  return solver_;
}
