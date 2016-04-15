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

// generates cuts and adds them to si.
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
    if (col_type[i]!= (char) 0) {
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
            for (int k=0; k<cone_size; ++k) {
              if (members[k]==i) {
                candidates.push_back(std::make_pair(i,j));
              }
            }
            delete[] members;
          }
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
