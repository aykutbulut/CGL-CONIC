// Name:     CglConicGD1.cpp
// Author:   Aykut Bulut
//           Lehigh University
//           email: aykut@lehigh.edu, aykutblt@gmail.com
//-----------------------------------------------------------------------------
// Copyright (C) 2015, Lehigh University.  All Rights Reserved.

#include "CglConicGD1.hpp"
#include <string>
#include <sstream>
#include <iomanip>

#define EPS 1e-5

CglConicGD1::CglConicGD1(OsiConicSolverInterface * solver)
  : param_(0) {
  param_ = new CglConicGD1Param();
  num_cuts_ = 0;
  cuts_.clear();
}

// copy constructor
CglConicGD1::CglConicGD1(const CglConicGD1 & other) {
  // copy param_
  param_ = new CglConicGD1Param(*(other.getParam()));
  num_cuts_ = 0;
  cuts_.clear();
}

// copy assignment operator
CglConicGD1 & CglConicGD1::operator=(const CglConicGD1 & rhs) {
  // copy param_
  param_ = new CglConicGD1Param(*(rhs.getParam()));
  num_cuts_ = 0;
  cuts_.clear();
  return *this;
}

CglConicGD1::~CglConicGD1() {
  if (param_)
    delete param_;
  int size = cuts_.size();
  for (int i=0; i<size; ++i) {
    delete cuts_[i];
  }
  cuts_.clear();
}

void CglConicGD1::setParam(const CglConicGD1Param & param) {
  param_ = new CglConicGD1Param(param);
}

void CglConicGD1::generateCuts(const OsiConicSolverInterface & si,
                               OsiConicCuts & cs,
                               const CglTreeInfo info) {
  std::cerr << "Not implemented yet!" << std::endl;
  throw std::exception();
}

// return true if all constraints are equality constraints.
bool CglConicGD1::constraint_check(
          OsiConicSolverInterface const & solver) const {
  double const * rowlb = solver.getRowLower();
  double const * rowub = solver.getRowUpper();
  char const * row_sense = solver.getRowSense();
  int num_rows = solver.getNumRows();
  for (int i=0; i<num_rows; ++i) {
    if (row_sense[i]!='E' and rowlb[i]!=rowub[i]) {
      return false;
    }
  }
  return true;
}

OsiConicSolverInterface * CglConicGD1::generateAndAddCuts(
                   OsiConicSolverInterface const & si,
                   const CglTreeInfo info) {
  // create a copy of the input
  OsiConicSolverInterface * solver = si.clone();
  // decide disjunction var and cut cone
  std::vector<std::pair<int,int> > candidates;
  candidates = compute_dis_var_cone(si);
  std::vector<std::pair<int,int> >::const_iterator it;
  // print cut candidates
  for (it=candidates.begin(); it!=candidates.end(); it++) {
    int dis_var = it->first;
    int cut_cone = it->second;
    std::cout << "Cut var " << dis_var
              << " cone " << cut_cone << std::endl;
  }
  // add cuts for all (var,cone) pair
  for (it=candidates.begin(); it!=candidates.end(); it++) {
    int dis_var = it->first;
    int cut_cone = it->second;
    std::cout << "Adding cut using var " << dis_var
              << " cone " << cut_cone << std::endl;
    // get equality constraints the cut cone members are in
    // get Ax=b
    int num_eq_rows;
    int * rows;
    get_rows(si, cut_cone, num_eq_rows, rows);
    if (num_eq_rows==0) {
      delete[] rows;
      continue;
    }
    // model.addConicCutGenerator(rows1, 2, 0);
    CglConicGD1Cut * cut = new CglConicGD1Cut(&si, num_eq_rows, rows,
                                              cut_cone, dis_var);
    num_cuts_++;
    if(!cut->valid()) {
      std::cerr << "Generated cut is not valid." << std::endl;
      delete cut;
    }
    else {
      cuts_.push_back(cut);
    }
    //add_cut(cut);
    delete[] rows;
  }
  add_cuts(solver);
  // clear cuts after adding them
  clear_cuts();
  //std::copy(cut_row.begin(), cut_row.end(), cut_row_in);
  // once we determined the cut row, we can generate cut using any cone and
  // any cone member we want.
  return solver;
}

// generates cuts and adds them to si.
OsiConicSolverInterface * CglConicGD1::generateAndAddCuts2(
                                     OsiConicSolverInterface const & si,
                                     const CglTreeInfo info) {
  OsiConicSolverInterface * solver = si.clone();

  {
    char const * row_sense = solver->getRowSense();
    double const * lb = solver->getRowLower();
    double const * ub = solver->getRowUpper();
    // print row sense
    std::cout << "============================== lb ub =============================="
              << std::endl;
    for (int i=0; i<solver->getNumRows(); ++i) {
      std::cout << i << " " << row_sense[i] << " " << lb[i] << " " << ub[i] << std::endl;
    }
  }
  preprocess(*solver, si);

  {
    char const * row_sense = solver->getRowSense();
    double const * lb = solver->getRowLower();
    double const * ub = solver->getRowUpper();
    double const * rhs = solver->getRightHandSide();
    // print row sense
    std::cout << "============================== lb ub =============================="
              << std::endl;
    for (int i=0; i<solver->getNumRows(); ++i) {
      std::cout << i << " " << row_sense[i] << " " << lb[i] << " " << ub[i] << std::endl;
    }
  }

  // check whether all constraints are equality
  if (!constraint_check(*solver)) {
    std::cout << "Only equality constraints are allowed." << std::endl;
    throw std::exception();
  }
  //OsiConicSolverInterface * solver = si.clone();


  // Summary of this function
  // 1. compute null space for the whole constraint matrix
  // Ax = b, x = x^0 + Hw
  // then (x^0 + Hw)^i in L^{n_i} for each x^i in L^{n_i}
  // 2. we are interested in cones where (x^0 + Hw)^i _j
  // is an integer variable with a fractional value
  // 3. Create disjunction based on (x^0 + Hw)^i _j.
  // 4. Compute quadric for the cone we are interested.
  // Q^i <- H^iT J H^i, q^i <- H^iT Jx^0i , \rho^i <- x^0i^T J X^0i.
  //
  //

  // == 1. compute null space for coefficient matrix
  int num_cols = solver->getNumCols();
  int num_rows = solver->getNumRows();
  // col ordered dense A
  double * denseA = new double[num_cols*num_rows]();
  // copy A to dense A
  CoinPackedMatrix const * matA = solver->getMatrixByCol();
  int const * starts = matA->getVectorStarts();
  int const * lengths = matA->getVectorLengths();

  for (int i=0; i<num_cols; ++i) {
    int const * indices = matA->getIndices()+starts[i];
    double const * values = matA->getElements()+starts[i];
    double * col = denseA+i*num_rows;
    for (int j=0; j<lengths[i]; ++j) {
      col[indices[j]] = values[j];
    }
  }

  // print dense A
  {
    std::cout << "==================== dense A ===================="
              << std::endl;
    for (int i=0; i<num_rows; ++i) {
      for (int j=0; j<num_cols; ++j) {
        std::cout << denseA[j*num_rows+i] << " ";
      }
      std::cout << std::endl;
    }
  }
  //

  //Right hand side singular vectors of A
  double * VT = new double[num_cols*num_cols];
  svDecompICL(num_rows, num_cols, denseA, VT);
  double * matH = new double[(num_cols-num_rows)*num_cols]();
  // Take the last n-m columns of V, lapack returns V^T
  for(int i=0; i<(num_cols-num_rows); ++i) {
    cblas_dcopy(num_cols, (VT+num_rows+i), num_cols, (matH+i*num_cols), 1);
  }
  delete[] denseA;
  delete[] VT;


  // print H
  {
    std::cout << "==================== H ===================="
              << std::endl;
    int m = num_cols;
    int n = num_cols-num_rows;
    for (int i=0; i<m; ++i) {
      for (int j=0; j<n; ++j) {
        std::cout << matH[j*m+i] << " ";
          }
      std::cout << std::endl;
    }
  }
  //
  double const * x0 = si.getColSolution();
  //print x0
  {
    std::cout << "==================== x0 ===================="
              << std::endl;
    int m = num_cols;
    for (int i=0; i<m; ++i)
      std::cout << x0[i] << " ";
    std::cout << std::endl;
  }
  // == x2. we are interested in cones where (x^0 + Hw)^i _j
  // is an integer variable with a fractional value
  // ==== compute list of fractional variables and their cone
  std::vector<std::pair<int,int> > candidates;
  candidates = compute_dis_var_cone(si);
  std::vector<std::pair<int,int> >::const_iterator it;
  // add cuts for all (var,cone) pair
  for (it=candidates.begin(); it!=candidates.end(); it++) {
    int dis_var = it->first;
    int cut_cone = it->second;
    // print cut candidates
    std::cout << "Cut var " << dis_var
              << " value " << x0[dis_var]
              << " cone " << cut_cone << std::endl;
    // 3. Create disjunction based on (x^0 + Hw)^i _j.
    // Lorentz cone type of cut generating cone, rotated or not.
    OsiLorentzConeType lctype;
    // Size of cut generating cone.
    int csize;
    // Members of cut generating cone.
    int * cmembers;
    solver->getConicConstraint(cut_cone, lctype, csize, cmembers);
    {
      std::cout << "==================== cone members ===================="
                << std::endl;
      for (int i=0; i<csize; ++i)
        std::cout << cmembers[i] << " ";
      std::cout << std::endl;
    }
    // get relted portion of H, this means rows corresponding to the
    // variables in the cut cone.
    // partH is column ordered.
    int m = num_cols;
    int n = num_cols-num_rows;
    int parm = csize;
    double * partH = new double[csize*n];
    for (int i=0; i<csize; ++i) {
      // get row cmembers[j] from matH.
      // row i of partH is row cmembers[i] of matH.
      for (int j=0; j<n; ++j) {
        partH[j*csize+i] = matH[j*m+cmembers[i]];
      }
    }
    {
      std::cout << "==================== partial H ===================="
                << std::endl;
      for (int i=0; i<parm; ++i) {
        for (int j=0; j<n; ++j) {
          std::cout << partH[j*parm+i] << " ";
        }
        std::cout << std::endl;
      }
    }
    // end of getting partial H
    double * a = new double[csize]();
    int new_dis_index = 0;
    for (int i=0; i<csize; ++i) {
      if (cmembers[i]!=dis_var)
        new_dis_index++;
      else
        break;
    }
    a[new_dis_index] = 1.0;
    // compute u <- H^Ta
    double * u = new double[n]();
    cblas_dgemv(CblasColMajor, CblasTrans, parm, n, 1.0, partH, parm, a, 1, 0.0, u, 1);
    double norm_of_u = cblas_dnrm2(n, u, 1);
    if (norm_of_u<EPS) {
      std::cout << "Numerical problems. Go to next candidate." << std::endl;
      delete[] partH;
      delete[] a;
      delete[] u;
      continue;
    }
    // todo(aykut) do this using blas
    for (int i=0; i<n; ++i) {
      u[i] = u[i]/norm_of_u;
    }

    // compute partial x0 and x0J
    double * vecx0 = new double[csize];
    for (int i=0; i<csize; ++i) {
      vecx0[i] = x0[cmembers[i]];
    }
    double * vecx0J = new double[csize];
    std::copy(vecx0, vecx0+csize, vecx0J);
    vecx0J[0] = -1.0*vecx0J[0];
    // end of partial x0

    // define alpha and beta
    double aT_x0 = cblas_ddot(parm, a, 1, vecx0, 1);
    double alpha = floor(x0[dis_var])+1.0 - aT_x0;
    double beta =  floor(x0[dis_var]) - aT_x0;
    alpha = alpha / norm_of_u;
    beta = beta / norm_of_u;
    // create coefficient matrix
    Disjunction * disjunction = new Disjunction(csize, u, alpha, u, beta);
    delete[] u;
    delete[] a;

    // 4. Compute quadric for the cone we are interested.
    // Q^i <- H^iT J H^i, q^i <- H^iT Jx^0i , \rho^i <- x^0i^T J X^0i.
    // == matrix Q
    double * matQ = new double[n*n]();
    // Temporary array, stores the first row of H
    double * d = new double[n];
    cblas_dcopy(n, partH, parm, d, 1);
    // Temporary array, exlcudes the first row of H
    double * tempH = new double[(parm-1)*n];
    for(int i=0; i<n; i++)
      cblas_dcopy(parm-1, partH+i*parm+1, 1, tempH+i*(parm-1), 1);
    //computes Q = tempH^T tempH - dd^T =  H^T J H
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, n, parm-1,
                1.0, tempH, parm-1, 0.0, matQ, n);
    cblas_dsyr(CblasColMajor, CblasUpper, n, -1.0, d, 1, matQ, n);
    delete[] d;
    delete[] tempH;
    // print Q
    {
      std::cout << "==================== Q ===================="
                << std::endl;
      for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
          if (i<=j)
            std::cout << matQ[j*n+i] << " ";
          else
            std::cout << matQ[i*n+j] << " ";
        }
        std::cout << std::endl;
      }
    }
    // vector q
    double * vecq = new double[n]();
    cblas_dgemv (CblasColMajor, CblasTrans, parm, n, 1.0, partH, parm, vecx0J, 1,
                 0.0, vecq, 1);
    // print vector q
    {
      std::cout << "==================== q ===================="
                << std::endl;
      for (int i=0; i<n; ++i)
        std::cout << vecq[i] << " ";
      std::cout << std::endl;
    }
    // scalar rho
    double rho = - (vecx0J[0]*vecx0J[0]);
    rho += cblas_ddot(parm-1, vecx0J+1, 1, vecx0J+1, 1);
    {
      std::cout << "==================== rho ===================="
                << std::endl;
      std::cout << rho << std::endl;
    }

    // OsiLorentzConeType lctype;
    // int csize;
    // int * cmembers;
    // solver->getConicConstraint(cut_cone, lctype, csize, cmembers);
    CglConicGD1Cut * cut = new CglConicGD1Cut(matQ, vecq, rho,
                                              parm, n, partH,
                                              vecx0,
                                              cmembers,
                                              dis_var,
                                              *disjunction);
    delete[] cmembers;
    num_cuts_++;
    if(!cut->valid()) {
      std::cerr << "Generated cut is not valid." << std::endl;
      delete cut;
    }
    else {
      cuts_.push_back(cut);
    }
    delete disjunction;
    delete[] partH;
    delete[] matQ;
    delete[] vecx0;
    delete[] vecx0J;
    delete[] vecq;
    break;
  }

  delete[] matH;

  add_cuts(solver);
  // clear cuts after adding them
  clear_cuts();
  //std::copy(cut_row.begin(), cut_row.end(), cut_row_in);
  // once we determined the cut row, we can generate cut using any cone and
  // any cone member we want.
  return solver;
}

void CglConicGD1::preprocess(OsiConicSolverInterface & si, OsiConicSolverInterface const & solver) const {
  // check number of equality constraints
  int num_rows = solver.getNumRows();
  char const * row_sense = solver.getRowSense();
  double const * rhs = solver.getRightHandSide();
  int num_eq_rows = 0;
  for (int i=0; i<num_rows; ++i) {
    if (row_sense[i]=='E')
      num_eq_rows++;
  }
  // modify rows by adding slacks.
  // CoinPackedMatrix const * mat = solver.getMatrixByRow();
  // double const * elements = mat->getElements();
  // int const * cols = mat->getIndices();
  int rows[1];
  double vals[1];
  for (int i=0; i<num_rows; ++i) {
    if (row_sense[i]=='L') {
      // add slack column
      rows[0] = i;
      vals[0] = 1.0;
      si.addCol(1, rows, vals, 0.0, si.getInfinity(),
                0.0, std::string("s"));
      si.setRowLower(i, rhs[i]);
      si.setRowUpper(i, rhs[i]);
      std::cout << "setting row "<< i << " to " << rhs[i] << std::endl;
    }
    else if (row_sense[i]=='G') {
      // add slack column
      rows[0] = i;
      vals[0] = -1.0;
      si.addCol(1, rows, vals, 0.0, si.getInfinity(),
                0.0, std::string("s"));
      si.setRowLower(i, rhs[i]);
      si.setRowUpper(i, rhs[i]);
      std::cout << "setting row "<< i << " to " << rhs[i] << std::endl;
    }
    else if (row_sense[i]=='R') {
      std::cout << "Range constraint, do not know what to do!" << std::endl;
      throw std::exception();
    }
  }
}

void CglConicGD1::clear_cuts() {
  int size = cuts_.size();
  for (int i=0; i<size; ++i) {
    delete cuts_[i];
  }
  cuts_.clear();
}

void CglConicGD1::add_cuts(OsiConicSolverInterface * solver) {
  std::vector<CglConicGD1Cut*>::const_iterator it;
  for (it=cuts_.begin(); it!=cuts_.end(); ++it) {
    add_cut(solver, *it);
  }
}

// when is cut invalid? (1) when no hyperplane intersects quadric.
// we have 3 options, cut_type_ 1, 2, 3.
void CglConicGD1::add_cut(OsiConicSolverInterface * solver,
                          CglConicGD1Cut * cut) {
  if(!cut->valid()) {
    std::cerr << "Generated cut is not valid." << std::endl;
    return;
  }
  int cut_type = cut->cutType();
  if (cut_type==1)  {
    add_cone_form_cut(solver, cut);
  }
  else if (cut_type==2) {
    int const * ind = cut->linear_cut_ind();
    double const * coef = cut->linear_cut_coef();
    double rhs = cut->linear_cut_rhs();
    int size = cut->linear_cut_size();
    // solver->addRow(size, ind, coef, rhs, solver->getInfinity());
    // for now assume disjunction is variable disjunction and
    // update column bound.
    solver->setColLower(ind[0], rhs);
  }
  else if (cut_type==3) {
    int const * ind = cut->linear_cut_ind();
    double const * coef = cut->linear_cut_coef();
    double rhs = cut->linear_cut_rhs();
    int size = cut->linear_cut_size();
    // solver->addRow(size, ind, coef, -solver->getInfinity(), rhs);
    solver->setColUpper(ind[0], rhs);
  }
  else {
    std::cerr << "Unknown conic cut type!" << std::endl;
    throw std::exception();
  }
}

void CglConicGD1::add_cone_form_cut(OsiConicSolverInterface * solver, CglConicGD1Cut * cut) {
  int * members = cut->getMembers();
  int num_cut_rows = cut->getNumCols();
  // num_cut_rows is same as cut cone size
  int num_cut_cols = cut->getNumRows();
  double const * newMatA = cut->getNewMatA();
  double const * newRhs =  cut->getNewRhs();
  int num_orig_rows = solver->getNumRows();
  int num_orig_cols = solver->getNumCols();
  // add Ix + Ay = b and y \in L to the model. x is old variable.
  // y is new variable.
  // add rows Ix = b first.
  int * ind = new int[1];
  double * val = new double[1];
  val[0] = 1.0;
  for (int i=0; i<num_cut_rows; ++i) {
    ind[0] = members[i];
    solver->addRow(1, ind ,val, newRhs[i], newRhs[i]);
  }
  delete[] ind;
  delete[] val;
  // add columns Ay to rows Ix = b.
  ind = new int[num_cut_rows];
  val = new double[num_cut_rows];
  int num_elem = 0;
  double value;
  for (int i=0; i<num_cut_cols; ++i) {
    // newMatA is col major
    for (int j=0; j<num_cut_rows; ++j) {
      value = newMatA[i*num_cut_rows+j];
      if (value>EPS || value<-EPS) {
        ind[num_elem] = num_orig_rows+j;
        val[num_elem] = value;
        num_elem++;
      }
      else {
        std::cout << "Small coefficient in the linear constraints." << std::endl;
      }
    }
    // add col
    solver->addCol(num_elem, ind, val, -solver->getInfinity(),
                    solver->getInfinity(), 0.0);
    num_elem = 0;
  }
  delete[] ind;
  delete[] val;
  // add cone y in L. Leading variable may not be y_1.
  if (num_cut_cols==1) {
    solver->setColLower(num_orig_cols, 0.0);
  }
  else {
    int * cone_ind = new int[num_cut_cols];
    int var_head = num_orig_cols + cut->getVarHead();
    cone_ind[0] = var_head;
    int k = 1;
    for (int i=0; i<num_cut_cols; ++i) {
      if (i==(var_head-num_orig_cols)) {
        continue;
      }
      else {
        cone_ind[k] = num_orig_cols+i;
        k++;
      }
    }
    OsiLorentzConeType type = OSI_QUAD;
    solver->addConicConstraint(type, num_cut_cols, cone_ind);
    delete[] cone_ind;
  }
}

void CglConicGD1::get_rows(OsiConicSolverInterface const & si,
                           int cut_cone, int & num_eq_rows, int *& rows) const {
  std::vector<int> vrows;
  OsiLorentzConeType type;
  int cone_size;
  int * members;
  si.getConicConstraint(cut_cone, type, cone_size, members);
  int m = si.getNumCols();
  int n = si.getNumRows();
  int * cone_members = new int[m]();
  for (int i=0; i<cone_size; ++i) {
    cone_members[members[i]] = 1;
  }
  CoinPackedMatrix const * mat = si.getMatrixByRow();
  char const * row_sense = si.getRowSense();
  for (int i=0; i<n; ++i) {
    int flag = 0;
    if (row_sense[i]=='E') {
      int first = mat->getVectorFirst(i);
      int last = mat->getVectorLast(i);
      int const * ind = mat->getIndices();
      //int vec_size = mat->getVectorSize(i);
      for (int j=first; j<last; j++) {
        if (cone_members[ind[j]]==0) {
          flag = 1;
          break;
        }
      }
      if (flag==0) {
        vrows.push_back(i);
      }
    }
  }
  num_eq_rows = vrows.size();
  rows = new int[num_eq_rows];
  std::vector<int>::const_iterator it;
  int k=0;
  for (it=vrows.begin(); it!=vrows.end(); ++it) {
    rows[k] = *it;
    k++;
  }
  delete[] members;
  delete[] cone_members;
}

int CglConicGD1::getNumCutsAdded() const {
  return num_cuts_;
}

/// Clone
CglConicCutGenerator * CglConicGD1::clone() const {
  //
  //CglConicGD1 * new_cut = new CglConicGD1(*this);
  std::cerr << "Clone is not implemented yet!" << std::endl;
  throw std::exception();
  CglConicGD1 * new_cut;
  return new_cut;
}

/// Create C++ lines to get to current state
std::string CglConicGD1::generateCpp( FILE * fp) {
  std::cerr << "GenerateCpp is not implemented yet!" << std::endl;
  throw std::exception();
  return std::string();
}

// compute a set of disjunction variables
std::vector<std::pair<int,int> > CglConicGD1::compute_dis_var_cone(
                              OsiConicSolverInterface const & si) const {
  // return the set of variables that are fractional and in a cone.
  std::vector<std::pair<int,int> > candidates;
  double const * sol = si.getColSolution();
  int num_col = si.getNumCols();
  char const * col_type = si.getColType(true);
  int num_cones = si.getNumCones();
  OsiConeType * cone_types = new OsiConeType[num_cones];
  si.getConeType(cone_types);
  for (int i=0; i<num_col; ++i) {
    if (col_type[i]== (char) 0) {
      continue;
    }
    // check wheather it is fractional
    double floor_i = floor(sol[i]);
    double ceil_i = floor_i + 1.0;
    if ((ceil_i-sol[i])>EPS && (sol[i]-floor_i)>EPS) {
      // variable i is fractional, check whether it is in a cone?
      for (int j=0; j<num_cones; ++j) {
        if (cone_types[j]==OSI_LORENTZ) {
          OsiLorentzConeType lctype;
          int cone_size;
          int * members;
          si.getConicConstraint(j, lctype, cone_size, members);
          // skip the cone if it is not binding
          double term1 = sol[members[0]];
          double term2 = 0.0;
          for (int k=0; k<cone_size; ++k) {
            term2 += sol[members[k]]*sol[members[k]];
          }
          double value = term1 - sqrt(term2);
          if (value>1e-5) {
            delete[] members;
            std::cout << "Skipping non-binding cone " << j << " "
                      << value << std::endl;
            continue;
          }
          // todo(aykut) a linear cut would do if it is a leading variable
          // if (members[0]==i) {
          //   std::cout << "leading variable. skipping..." << std::endl;
          // }
          // we are interested in other members
          for (int k=0; k<cone_size; ++k) {
            if (members[k]==i) {
              candidates.push_back(std::make_pair(i,j));
            }
          }
          delete[] members;
        }
        else {
          std::cout << "Implemented for Lorentz cones only!" << std::endl;
        }
      }
    }
  }
  delete[] cone_types;
  if (candidates.empty()) {
    std::cerr << "Could not find a suitable variable to create disjunction!"
              << std::endl;
    //throw std::exception();
    return candidates;
  }
  else {
    return candidates;
  }
}

void CglConicGD1::print_cut(CglConicGD1Cut * cut) const {
}
