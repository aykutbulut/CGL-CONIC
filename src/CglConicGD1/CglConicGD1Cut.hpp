#include <OsiConicSolverInterface.hpp>
#include "CglConicCutGenerator.hpp"

/*!
  Class to represent a conic cut of Belotti et. al. introduced in [1] and [2].
  Typical use of this class:

  Generate an instance with the proper inputs (conic problem, rows and disjunction
  cone and variable index). This constructor will generate the cut and store
  it. Then CglConicGD1 has a method to add this cut to a OsiConicSolverInterface.

  Moreover this class provides functions to querry the generated cut.

  This procedure works for un-rotated Lorentz cones only for now.

  References:
  [1] On families of quadratic surfaces having fixed intersection with two
  hyperplanes
  [2] A conic representation of the convex hull of disjunctive sets and conic
  cuts for integer second order cone optimization
 */

class CglConicGD1Cut {
  ///@name Data obtained from input problem directly.
  //@{
  /// Indices of problem rows that will be used to generate cut.
  int * rows_;
  /// Coefficient matrix obrtained using rows_.
  double * matA_;
  /// vecx0_ is a point that will be obtained from the given solver. It is a
  /// feasible point for the linear constraints and conic constraints for the
  /// given problem.
  double * vecx0_;
  /// Index of cut generating cone.
  int cindex_;
  /// Type of cut generating cone, Lorentz or Scaled.
  OsiConeType ctype_;
  /// Lorentz cone type of cut generating cone, rotated or not.
  OsiLorentzConeType lctype_;
  /// Size of cut generating cone.
  int csize_;
  // Members of cut generating cone.
  int * cmembers_;
  /// cone_members_[i] is 1 if column i is in the cone, 0 otherwise.
  /// length of this array is larger than csize_ when multiple cones
  /// exist in the problem.
  int * cone_members_;
  //@}

  ///@name Quadric data obtained by processing the row and cone input. The
  // quadric gives a quadratic constraint representation of the same set,
  // except nonnegativity of the leading variable.
  //@{
  /// Matrix Q of quadric.
  double * matQ_;
  /// Vector q of quadric.
  double * vecq_;
  /// Scalar rho of quadric.
  double rho_;
  /// matV_ keeps eigenvectors of matQ_.
  double * matV_;
  /// Eigenvalues of matrix Q.
  double * eigQ_;
  /// x^0 is a feasible point (both linear and conic constraints).
  /// vecx0J_ is x^0T J.
  double * vecx0J_;
  /// Columns of matH_ are basis for the null space of matA_. Columns of matH_
  /// are normalized.
  double * matH_;
  //@}

  ///@name Cut data.
  //@{
  /// Parameter tau that gives a cone from convolution of
  /// quadric and disjunction.
  double tau_;
  double tau1_;
  double tau2_;
  /// Matrix Q(tau).
  double * matQ_tau_;
  /// Vector q(tau)/
  double * vecq_tau_;
  /// Scalar rho(tau).
  double rho_tau_;
  // todo(aykut) I do not see the use of the following, keep it for now.
  double * dirTestV_;
  // todo(aykut) see what valid is for
  bool valid_;
  /// Disjunction used to generate the cut, in x-space.
  Disjunction * disjunction_;
  /// Solver interface.
  OsiConicSolverInterface const * solver_;
  /// Disjunction in w space, a_ is coefficient, alpha and beta are right hand
  /// sites.  one in x-space.
  double * a_;
  double alpha_;
  double beta_;
  ///
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
