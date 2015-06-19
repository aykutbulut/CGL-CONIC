#include <OsiConicSolverInterface.hpp>
#include "CglConicCutGenerator.hpp"

class CglConicGD1Cut {
  // rows to generate cut
  int * rows_;
  // matrix A
  double * matA_;
  // null matrix of A
  double * matH_;
  // matrix V, eigenvectors of Q
  double * matV_;
  // matrix Q of quadratic
  double * matQ_;
  // vector q of quadratic
  double * vecq_;
  // scalar rho of quadratic
  double rho_;
  // point x0
  double * vecx0_;
  // x^0 ^T J
  double * vecx0J_;
  // parameter tau that gives us a cone from quadric
  double tau_;
  double tau1_;
  double tau2_;
  // matrix Q(tau) of quadratic at tau
  double * matQ_tau_;
  // vector q(tau) of quadratic at tau
  double * vecq_tau_;
  // scalar rho(tau) of quadratic at tau
  double rho_tau_;
  // cut generating cone index
  int cindex_;
  // cut generating cone type, for now it works for Lorentz cones only
  OsiConeType ctype_;
  OsiLorentzConeType lctype_;
  // cut generating cone size
  int csize_;
  // cut generating cone members
  int * cmembers_;
  // cone members is 0-1 array
  int * cone_members_;
  // eigenvalues of matrix Q
  double * eigQ_;
  // todo(aykut) I do not see the use of the following, keep it for now.
  double * dirTestV_;
  // todo(aykut) see what valid is for
  bool valid_;
  // disjunction used.
  Disjunction * disjunction_;
  // solver interface
  OsiConicSolverInterface const * solver_;
  // for now we will keep disjunction in w space in u alpha and beta. disjunction_ is the
  // one in x-space.
  double * a_;
  double alpha_;
  double beta_;
  double * new_matA_;
  double * new_rhs_;
  // number of rows of A used for generating cut.
  int num_rows_;
  // disjunction var
  int dis_var_;
  // todo(aykut) get rid of varHead_ variable.
  int varHead_;
  double * dirTestU_;
  // quadric type
  char quad_type_;
  // cut type, 1 for primal, -1 for dual, 2 and 3 for linear cuts, when cut is linear
  // check coef_ and linear_cut_rhs_ to retrieve it.
  // when 2 a^T x >= rhs
  // when 3 a^T x <= rhs
  int cut_type_;
  int * linear_cut_ind_;
  double * linear_cut_coef_;
  double linear_cut_rhs_;
  // whether the cut is active
  bool active_;
  // number of 0 eigenvalues of Q
  int numZEig_;
  void generate_cut();
  void compute_matrixA();
  void compute_matrixH();
  void compute_matrixQ();
  void compute_vectorq();
  void compute_rho();
  void classify_quadric();
  void compute_disjunction();
  // compute tau value that yields a cone
  void compute_tau();
  // compute rho(tau), rho_tau_
  void compute_rho_tau();
  // compute q(tau), q_tau_
  void compute_q_tau();
  // compute Q(tau), Q_tau_
  void compute_Q_tau();
  // compute matrix as a part of the cut to add to model
  void compute_new_A();
  // compute right-hand-side that corresponds to NewA
  void compute_new_rhs();
  void regularize();
  // compute cone at tau
  //void coneAtTau();
  // solve quadratic formula for its largest roots
public:
  CglConicGD1Cut();
  CglConicGD1Cut(OsiConicSolverInterface const * solver,
		 int num_rows, int * rows,
		 int cut_cone, int dis_var);
  bool valid() const;
  // number of rows of the linear part of the cut
  // ie num of rows of the new A matrix.
  int getNumRows() const;
  // number of columns in the cut
  int getNumCols() const;
  // get variables in the cut
  int * getMembers() const;
  int cutType() const;
  // return linear part of the cut, constraint matrix.
  double const * getNewMatA() const;
  // return right hand side of the linear part of the cut.
  double const * getNewRhs() const;
  int getVarHead() const;
  // get the variables in the cut if cut is in linear form
  int const * linear_cut_ind() const;
  // get the coefficients of the cut if cut is in linear form
  double const * linear_cut_coef() const;
  // get the rhs of the cut if the cut is in linear form
  double linear_cut_rhs() const;
  // get the size of the cut if the cut is in the linear form
  int linear_cut_size() const;
  ~CglConicGD1Cut();
};
