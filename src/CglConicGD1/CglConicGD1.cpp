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
#include <numeric>

#define EPS 1e-2

CglConicGD1::CglConicGD1(OsiConicSolverInterface * solver)
  : param_(0) {
  param_ = new CglConicGD1Param();
  num_cuts_ = 0;
  cuts_.clear();
  cuts_cone_ind_.clear();
}

// copy constructor
CglConicGD1::CglConicGD1(const CglConicGD1 & other) {
  // copy param_
  param_ = new CglConicGD1Param(*(other.getParam()));
  num_cuts_ = 0;
  cuts_.clear();
  cuts_cone_ind_.clear();
}

// copy assignment operator
CglConicGD1 & CglConicGD1::operator=(const CglConicGD1 & rhs) {
  // copy param_
  param_ = new CglConicGD1Param(*(rhs.getParam()));
  num_cuts_ = 0;
  cuts_.clear();
  cuts_cone_ind_.clear();
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
  cuts_cone_ind_.clear();
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
  solver->setHintParam(OsiDoReducePrint, true, OsiHintTry);
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
              << " value " << solver->getColSolution()[dis_var]
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
    // get disjunction index relative to the rows Ax = b
    // submatrix A and right hand side
    CoinPackedMatrix * A = NULL;
    double * b = NULL;
    double * sol = NULL;
    int rel_dis_var;
    get_input_set(solver, dis_var, cut_cone, num_eq_rows, rows, A, b,
                  sol, rel_dis_var);
    CglConicGD1Cut * cut = new CglConicGD1Cut(PRIMAL_FORM, A, b, sol);
    double alpha = ceil(sol[rel_dis_var]);
    double beta = floor(sol[rel_dis_var]);
    cut->generateCut(rel_dis_var, alpha, beta);
    if(!cut->success()) {
      std::cout << "tau is " << cut->tau() << std::endl;
      std::cout << "Cut generation did not result any cuts." << std::endl;

      delete A;
      delete[] b;
      delete[] sol;
      delete[] rows;
      delete cut;
      continue;
    }
    else if (cut->infeasible()) {
      // problem is infeasible
      std::cout << "Problem is infeasible!"
                << std::endl;
      delete A;
      delete[] b;
      delete[] sol;
      delete[] rows;
      delete cut;
      continue;
    }
    //cut->print_cut();
    // compute violation
    double viol = 0.0;
    double const * cutA = cut->getCutA();
    double const * cutb = cut->getCutb();
    // viol = (A sol -b)_1 - norm((A sol - b)_{2:m})
    if (cut->getNumRows() > 1) {
      double * axmb = new double[cut->getNumRows()]();
      for (int i=0; i<cut->getNumRows(); ++i) {
        axmb[i] = 0.0;
        for (int j=0; j<cut->getNumCols(); ++j) {
          axmb[i] += cutA[j*cut->getNumRows() + i]*sol[j];
        }
        axmb[i] -= cutb[i];
      }
      double nn = std::inner_product(axmb+1, axmb+cut->getNumRows(), axmb+1, 0.0);
      viol =  sqrt(nn) - axmb[0];
      // print Ax-b
      //std::cout << "==================== A sol - b ====================" << std::endl;
      //for (int i=0; i<cut->getNumRows(); ++i) std::cout << axmb[i] << " ";
      //std::cout << std::endl;
      delete[] axmb;
    }
    else {
      viol = std::inner_product(cutA, cutA+cut->getNumCols(), sol, -cutb[0]);
      viol = -viol;
    }

    std::cout << "tau is       " << cut->tau() << std::endl;
    std::cout << "Violation is " << viol << std::endl;



    delete A;
    delete[] b;
    delete[] sol;
    delete[] rows;


    num_cuts_++;
    cuts_.push_back(cut);
    cuts_cone_ind_.push_back(cut_cone);
  }
  add_cuts(solver);
  // clear cuts after adding them
  clear_cuts();
  //std::copy(cut_row.begin(), cut_row.end(), cut_row_in);
  // once we determined the cut row, we can generate cut using any cone and
  // any cone member we want.
  return solver;
}


// generates cuts for each possible variable and add the best only.
OsiConicSolverInterface * CglConicGD1::generateAndAddBestCut(
                                     OsiConicSolverInterface const & si,
                                     const CglTreeInfo info) {
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
  double initial_obj = si.getObjValue();
  int best_var = -1;
  int best_cone = -1;
  double best_imp = 0.0;
  CglConicGD1Cut * best_cut = NULL;
  // check which (var,cone) pair is the best
  for (it=candidates.begin(); it!=candidates.end(); it++) {
    OsiConicSolverInterface * solver = si.clone();
    solver->setHintParam(OsiDoReducePrint, true, OsiHintTry);
    int dis_var = it->first;
    int cut_cone = it->second;
    std::cout << "Checking cut using var " << dis_var
              << " value " << solver->getColSolution()[dis_var]
              << " cone " << cut_cone << std::endl;
    // get equality constraints the cut cone members are in
    // get Ax=b
    int num_eq_rows;
    int * rows;
    get_rows(si, cut_cone, num_eq_rows, rows);
    if (num_eq_rows==0) {
      delete[] rows;
      delete solver;
      continue;
    }
    CoinPackedMatrix * A = NULL;
    double * b = NULL;
    double * sol = NULL;
    int rel_dis_var;
    get_input_set(solver, dis_var, cut_cone, num_eq_rows, rows, A, b,
                  sol, rel_dis_var);
    CglConicGD1Cut * cut = new CglConicGD1Cut(PRIMAL_FORM, A, b, sol);
    double alpha = ceil(sol[rel_dis_var]);
    double beta = floor(sol[rel_dis_var]);
    cut->generateCut(rel_dis_var, alpha, beta);
    if(!cut->success()) {
      std::cout << "tau is " << cut->tau() << std::endl;
      std::cout << "Cut generation did not result any cuts." << std::endl;

      delete A;
      delete[] b;
      delete[] sol;
      delete[] rows;
      delete cut;
      delete solver;
      continue;
    }
    else if (cut->infeasible()) {
      // problem is infeasible
      std::cout << "Problem is infeasible!"
                << std::endl;
      delete A;
      delete[] b;
      delete[] sol;
      delete[] rows;
      delete cut;
      delete solver;
      continue;
    }

    delete A;
    delete[] b;
    delete[] sol;
    delete[] rows;

    add_cone_from_cut(*solver, cut, cut_cone);
    solver->initialSolve();
    double obj = solver->getObjValue();
    double imp = obj-initial_obj;
    std::cout << "Improvement is " << imp << std::endl;
    if (imp>best_imp) {
      best_var = dis_var;
      best_cone = cut_cone;
      best_imp = imp;
      best_cut = cut;
    }
    else {
      delete cut;
    }
    delete solver;
  }
  // create a copy of the input
  OsiConicSolverInterface * solver = si.clone();
  if (best_cut) {
    num_cuts_++;
    cuts_.push_back(best_cut);
    cuts_cone_ind_.push_back(best_cone);
    // add best cut only
    add_cone_from_cut(*solver, best_cut, best_cone);
  }
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
  cuts_cone_ind_.clear();
}

void CglConicGD1::add_cuts(OsiConicSolverInterface * solver) {
  // std::vector<CglConicGD1Cut*>::const_iterator it;
  // for (it=cuts_.begin(); it!=cuts_.end(); ++it) {
  //   //add_cut(solver, *it);
  //   // add cut to the solver
  //   add_cone_from_cut(*solver, cut, cut_cone);
  // }

  for (int i=0; i<num_cuts_; ++i) {
    CglConicGD1Cut * cut = cuts_[i];
    int cone_ind = cuts_cone_ind_[i];
    add_cone_from_cut(*solver, cut, cone_ind);
  }
}

// when is cut invalid? (1) when no hyperplane intersects quadric.
// we have 3 options, cut_type_ 1, 2, 3.
void CglConicGD1::add_cut(OsiConicSolverInterface * solver,
                          CglConicGD1Cut const * cut) {
  if(!cut->success()) {
    std::cout << "Cut generation did not result any cuts." << std::endl;
    return;
  }
  else if (cut->infeasible()) {
    // cut procedure decided that problem is infeasible.
    std::cout << "Problem is infeasible!" << std::endl;
    throw std::exception();
  }
  // cut is form cutA x - cutb \in L.
  // cut is linear if number of rows of cutA is 1 and then cut is
  // cutA x - cutb >= 0
  //add_cone_from_cut(solver, cut);
}

void CglConicGD1::add_cone_from_cut(OsiConicSolverInterface & si,
                                    CglConicGD1Cut const * cut,
                                    int cut_cone_index) {
  // cut is Ax - b in L
  // or Ax - y = b, y in L in standard form.

  // x is given by the cut_cone_index.
  // x is members of the cone with the index cut_cone_index.
  // y is a new variable that will be introduced to the problem.
  //

  OsiLorentzConeType type;
  int cone_size;
  int * members;
  si.getConicConstraint(cut_cone_index, type, cone_size, members);
  if (type!=OSI_QUAD) {
    std::cout << "Lorentz cones only!." << std::endl;
    throw std::exception();
  }

  int num_orig_rows = si.getNumRows();
  int num_orig_cols = si.getNumCols();

  double const * cutA = cut->getCutA();
  double const * cutb =  cut->getCutb();

  int num_cut_rows = cut->getNumRows();
  // num_cut_cols is same as size of cut cone
  int num_cut_cols = cut->getNumCols();

  if (cone_size != num_cut_cols) {
    std::cerr << "Starting cone size should be same as "
      "number of columns in the cut." << std::endl;
    throw std::exception();
  }

  // check whether cut is linear
  if (num_cut_rows == 1) {
    // add Ax >= b
    si.addRow(num_cut_cols, members, cutA, cutb[0], si.getInfinity());
    delete[] members;
    return;
  }


  // add Ax - y = b to the model
  {
    // add rows Ax = b first.
    double * val = new double[num_cut_cols];
    // for each row of A
    for (int i=0; i<num_cut_rows; ++i) {
      // iterate over columns of A
      for (int j=0; j<num_cut_cols; ++j) {
        // A is col ordered
        val[j] = cutA[j*num_cut_rows+i];
      }
      si.addRow(num_cut_cols, members, val, cutb[i], cutb[i]);
    }
    delete[] val;
    // add columns -y to rows Ax = b.
    int * ind = new int[1];
    val = new double[1];
    val[0] = -1.0;
    for (int i=0; i<num_cut_rows; ++i) {
      ind[0] = num_orig_rows+i;
      // add col
      si.addCol(1, ind, val, -si.getInfinity(),
                     si.getInfinity(), 0.0);
    }
    delete[] ind;
    delete[] val;
  }

  // add y in L
  {
    int * cone_ind = new int[num_cut_rows];
    for (int i=0; i<num_cut_rows; ++i) {
      cone_ind[i] = num_orig_cols+i;
    }
    OsiLorentzConeType type = OSI_QUAD;
    si.addConicConstraint(type, num_cut_rows, cone_ind);
    delete[] cone_ind;
  }
  delete[] members;
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
    std::cout << "Could not find a suitable variable to create disjunction!"
              << std::endl;
    //throw std::exception();
    return candidates;
  }
  else {
    return candidates;
  }
}


void CglConicGD1::print_cut(CglConicGD1Cut const * cut) const {
}





OsiConicSolverInterface * CglConicGD1::generateAndAddCuts(
                              OsiConicSolverInterface const & si,
                              std::vector<CoinPackedMatrix const *> AA,
                              std::vector<double const *> bb) {
  return NULL;
}

OsiConicSolverInterface * CglConicGD1::generateAndAddBestCut(
                              OsiConicSolverInterface const & si,
                              std::vector<CoinPackedMatrix const *> AA,
                              std::vector<double const *> bb) {
  return NULL;
}


// get cone in Ax=b and x in L form
void CglConicGD1::get_input_set(OsiConicSolverInterface const * solver,
                                int dis_var,
                                int cut_cone,
                                int num_eq_rows,
                                int const * rows,
                                CoinPackedMatrix *& A,
                                double *& b,
                                double *& sol,
                                int & rel_dis_var) const {
  CoinPackedMatrix const * matrix = solver->getMatrixByRow();
  OsiLorentzConeType type;
  int csize;
  int * cmembers = NULL;
  solver->getConicConstraint(cut_cone, type, csize, cmembers);
  if (type != OSI_QUAD) {
    delete[] cmembers;
    std::cerr << "Not implemented yet. Only Lorentz cones for now." << std::endl;
    throw std::exception();
  }
  A = new CoinPackedMatrix(*matrix, num_eq_rows, rows, csize, cmembers);
  // get b
  b = new double[num_eq_rows];
  double const * full_b = solver->getRowLower();
  for (int i=0; i<num_eq_rows; ++i) {
    b[i] = full_b[rows[i]];
  }
  sol = new double[csize];
  double const * full_sol = solver->getColSolution();
  for (int i=0; i<csize; ++i) {
    sol[i] = full_sol[cmembers[i]];
  }

  for (int i=0; i<csize; ++i) {
    if (cmembers[i] == dis_var) {
      rel_dis_var = i;
      break;
    }
  }
  delete[] cmembers;
}
