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

  == Theoretical details in Cut generation ==

  Cut generation happens in 3 different spaces.
  <ul>

    <li> Primal space where the feasible region (i.e., \f$Ax = b\f$, \f$x in
    L\f$) is.

    <li> Dual space, where the feasible region is represented in quadratic
    form, i.e. \f$ w^\top Q w + 2q^\top w + \rho \leq 0\f$, this space is
    referred as w-space. See Belotti et. al. for definitions and evaluation of
    \f$Q\f$, \f$q\f$, and \f$\rho\f$. Sets given by this constraint might be
    referred as \f$(Q, q, \rho)\f$.

    <li> Regularized space, a space where the regularized quadric
    lives. Regularized quadric is in \f$(\tilde{J}, 0, \delta)\f$ form. This
    space will be referred as u-space (u-space since the equations are
    functions of \f$u\f$ in the Belotti et. al. paper.)

    === Generation of the Cut ===

    The cut is generated in the u-space and then represented in (or transformed
    into) the x-space to add into the input problem. To do these
    transformations one must understand the relationship (linear
    transformations) between (1) \f$x\f$ and \f$w\f$ and (2) \f$w\f$ and
    \f$u\f$. Transformation between \f$x\f$ and \f$w\f$ are as follows,

    \f$ x = x^0 + Hw\f$ and \f$w = H^\top x - H^\top x^0\f$

    where \f$x^0\f$ is a point in the affine space given by \f$Ax=b\f$ and
    columns of \f$H\f$ are an orthonormal basis for the null space of \f$A\f$.

    Relationship between \f$w\f$ and \f$u\f$ is as follows,

    \f$u = \frac{\tilde{D}^{\frac{1}{2} V^\top \left( w + Q^{-1}q \right)}
                {\sqrt{\left|q^\top Q^{1} q - \rho\right|}} \f$

    where \f$Q = VDV^\top\f$ and \f$\tilde{D}_{ii} := \left|D_{ii} \right|$.

    Cut is generated in u-space since the quadratic representation is simpler
    and formulas are easier to track. To generate cut in the u-space first the
    disjunction should be represented in u-space. To represent the disjunctions
    in u-space we insert (x^0+Hw) for x in disjunction constraints. This will
    give disjunctions in ters of w (in w-space). We then insert the inverse of
    the complicated formula above (\f$u\f$--\f$w\f$ relationship) for w to
    obtain disjunctions in terms of \f$u\f$.

    Once disjunctions are represented in u-space, we compute $\tau$. This A
    gives the cut in quadratic form. The cut is first represented in Lorentz
    cone form and then transformed into the x-space. The final cut is
    \f$ A' x - y = b' \f$ and \f$y \in \mathbb{L}^m\f$.


  References:
  [1] On families of quadratic surfaces having fixed intersection with two
  hyperplanes

  [2] A conic representation of the convex hull of disjunctive sets and conic
  cuts for integer second order cone optimization

  [3] A Complete Characterization of Disjunctive Conic Cuts for Mixed Integer
  Second Order Cone Optimization
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
  /// disjunctive variable index
  int var_index_;
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

  ///@name regularized cone data
  //@{
  /// regularized matrix Q
  double * Jtilde_;
  // regularized vector q is 0.0, no need to store
  /// regularized scalar rho
  double rho_tilde_;
  //@}

  ///@name Quadric data obtained by processing the row and cone input. The
  // quadric gives a quadratic constraint representation of the same set,
  // except nonnegativity of the leading variable.
  //@{
  /// Columns of matH_ are basis for the null space of matA_. Columns of matH_
  /// are normalized.
  double * matH_;
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
  //@}

  ///@name Cut data.
  //@{
  /// Parameter tau that gives a cone from convolution of
  /// quadric and disjunction.
  double tau_;
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
  // cutting revealed that problem is infeasible
  bool infeasible_;
  /// Disjunction used to generate the cut, in x-space.
  Disjunction * disjunction_;
  /// Solver interface.
  OsiConicSolverInterface const * solver_;
  /// Normalized disjunction in regularized space, a_ is coefficient, alpha and
  /// beta are right hand sites.
  double * a_;
  double alpha_;
  double beta_;
  /// resulting cut, cut is Anew x - y = new_rhs and y in L
  double * new_matA_;
  double * new_rhs_;
  // number of rows of A used for generating cut.
  int num_rows_;
  // disjunction var
  int dis_var_;
  // relative disjunction var
  int rel_dis_var_;
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
  // print matrix, row major 0, col major 1
  void print_matrix(int major, int num_rows, int num_cols,
                    double const * matrix,
                    char const * name = "matrix") const;
  // print vector
  void print_vector(int n, double const * vector,
                    char const * name = "vector") const;
  void print_scalar(double value, char const * name) const;
  // compute cone at tau
  //void coneAtTau();
  // solve quadratic formula for its largest roots
public:
  CglConicGD1Cut();
  CglConicGD1Cut(OsiConicSolverInterface const * solver,
                 int num_rows, int * rows,
                 int cut_cone, int dis_var);
  bool valid() const;
  /// returns true when the problem is infeasible
  bool infeasible() const;
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
  /// print cut
  void print_cut() const;
  ~CglConicGD1Cut();
};
