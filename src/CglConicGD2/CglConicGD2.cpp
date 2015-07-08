#include "CglConicGD2.hpp"

CglConicGD2::CglConicGD2(): param_(0),
  cut_generating_cone_(-1),
  beta_(0),
  disjunction_(0),
  cuts_(0) {
  param_ = new CglConicGD2Param();
}

// copy constructor
CglConicGD2::CglConicGD2(const CglConicGD2 & other) {
  // copy param_
  param_ = new CglConicGD2Param(*(other.getParam()));
}

// copy assignment operator
CglConicGD2 & CglConicGD2::operator=(const CglConicGD2 & rhs) {
  // copy param_
  param_ = new CglConicGD2Param(*(rhs.getParam()));
  return *this;
}

CglConicGD2::~CglConicGD2() {
  if (param_)
    delete param_;
}

void CglConicGD2::setParam(const CglConicGD2Param & param) {
  param_ = new CglConicGD2Param(param);
}

// generates MIR cut.
void CglConicGD2::generateCuts(const OsiConicSolverInterface & si,
			       OsiConicCuts & cs,
			       const CglTreeInfo info) {
  // get simplex tableau, ie. Ax=b
  int simplex_interface = si.canDoSimplexInterface();
  if (simplex_interface==1) {
    // has only tableau information
    std::cout << "Conic Cgl: I have tableau access but do not have pivoting"
	      << " implementation, OsiSimplex level 1." << std::endl;
  }
  else if (simplex_interface==0) {
    std::cerr << "Conic Cgl: I can not get simplex tableau information!"
	      << " OsiSimplex level 0." << std::endl;
    std::cerr << "Conic Cgl: Terminating cut generation..." << std::endl;
    return;
  }
  // option 1: pick cut generating cone, a cone in canonical form
  // this is already called in constructor
  // compute_cut_generating_cone();

  // option 2: pick a general disjunction (c1,c10) and (c2,c20) c10 and c20
  // are right-hand-sides, c1 and c2 are coefficient vectors. Disjunction in
  // the following form,
  // C_i := {x \in L: Ax=b, ci^T x >= ci0} i \in {1,2}
  // check 2: General disjunction should satisfy the following two assumptions,
  // Assumption 1:  C_1 does not contain C_2, C_2 does not contain C1.
  // Assumption 2: C_1 and C_2 are strictly feasible sets.
  // This is already called in constructor
  //compute_disjunction();

  double const * c1 = disjunction_->get_c1();
  double c10 = disjunction_->get_c10();
  double const * c2 = disjunction_->get_c2();
  double c20 = disjunction_->get_c20();
  // option 3: pick cut generation parameter beta, a rational number
  // check 3: check if beta satisfies the required conditions.
  // this is already called in constructor
  // compute_beta();

  // generate cut
  // for this we need to add n rows to the model, and
}

/// Clone
CglConicCutGenerator * CglConicGD2::clone() const {
  //
  CglConicGD2 * new_cut = new CglConicGD2(*this);
  return new_cut;
}

/// Create C++ lines to get to current state
std::string CglConicGD2::generateCpp( FILE * fp) {
  return std::string();
}

// compute cut generating cone using param_
void CglConicGD2::compute_cut_generating_cone() {
  cut_generating_cone_=-1;
}

// compute disjunction using param_
void CglConicGD2::compute_disjunction() {
  disjunction_ = new Disjunction(0, 0, 0.0, 0, 0.0);
}

// compute cut parameter beta using param_
void CglConicGD2::compute_beta() {
  beta_=0.0;
}
