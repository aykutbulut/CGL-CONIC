#ifndef CglConicGD1Cut_H
#define CglConicGD1Cut_H

#include "CglConicCutGenerator.hpp"


///@file

/*!
 Input set form. Primal form denotes input set is
 \f$ Ax = b, x \in \mathbb{L} \f$ form, whereas dual is
 \f$ x^\top Q x + 2 q^\top x + \rho \leq 0 \f$.
**/

enum CglConicInputType {
  PRIMAL_FORM = 0,
  DUAL_FORM,
  QUAD_FORM
};

/*!

  Class to represent a conic cut of Belotti et. al. introduced in [1] and [2].

  # How to Use #CglConicGD1Cut

  Purpose of this class is to compute disjunctive cuts for convex conic sets
  (sets represented by linear and conic constraints, they are not necessarily
  cones, might be ellipsoids). Inputs to this class can be in three forms,

  <ul>
    <li> Primal form, \f$ Ax = b, x \in \mathbb{L}\f$,
    <li> Dual form, \f$ \|Ax - b\| \leq d^\top x - \gamma \f$,
    <li> Quadratic form, \f$ x^\top Q x + 2q^\top x + \rho \leq 0\f$.
  </ul>

  Primal form can be converted to quadratic form using a basis for the null
  space of matrix \f$A\f$. See the Belotti et. al. papers for details. A linear
  transformation can map the primal form (in x-space) to a dual form in a
  different space (will be referred as w-space). Conversion between dual form
  and quadratic form does not need a mapping. Quadratic form can be obtained by
  algebraic manipulation of the dual form. Finally, note that not every
  quadratic formula leads to conic representable set. This class assumes the
  set corresponding to the input quadratic is conic representable.

  All of the inputs are mapped/converted to a quadratic representation in
  w-space. Cut is computed in w-space in quadratic form then converted to dual
  form for user. User can convert the cut to the primal form or quadratic form
  if she desires (both easy to do, no need for eigenvalue decomposition, null
  space computation etc).

  Typical use of this class:

  // generate cut instance
  // initiate cut generation by providing a disjunction
  // get resulting cut
  // check feasibility of the problem
  // check whether cut is linear
  // generate another cut from a different disjunction

  Intermediate computations (mapping input conic set to its quadratic
  representation in dual space in case of non-quadratic inputs, eigenvalue
  decomposition of matrix \f$Q\f$ of quadratic form) are chached for next cut
  generation calls for a different disjunction.

  # Theoretical details in Cut generation

  Cut generation happens in 3 different spaces.

  <ul>

    <li> Primal space where the feasible region (i.e., \f$Ax = b\f$, \f$x \in
    \mathbb{L}\f$) is in. Cones are in this form in OsiConicSolverInterface instances in
    CglConicGD1 level.

    <li> Dual space, where the feasible region is represented in quadratic
    form, i.e. \f$ w^\top Q w + 2q^\top w + \rho \leq 0\f$, this space is
    referred as w-space. See Belotti et. al. for definitions and evaluation of
    \f$Q\f$, \f$q\f$, and \f$\rho\f$. Sets given by this constraint might be
    referred as \f$(Q, q, \rho)\f$.

    <li> Regularized space, a space where the regularized quadric
    lives. Regularized quadric is in \f$(\tilde{J}, 0, \delta)\f$ form. This
    space will be referred as u-space (u-space since the equations are
    functions of \f$u\f$ in the Belotti et. al. paper.)

  </ul>

  ## Transformations between different spaces

  Conic input in primal form can be represented in w-space as a quadratic as
  follows,

  \f$ w^\top Q w + 2q^\top w + \rho \leq 0\f$,

  where, \f$ Q := H^\top J H\f$, \f$q := H^\top J x^0\f$, \f$\rho := x^{0T} J
  x^{0}\f$, columns of \f$H\f$ form an orthonormal basis for \f$A\f$ and
  \f$J\f$ is identity except \f$J(1,1)\f$ is \f$-1\f$.

  Whereas, conic input set in dual form (\f$Ax-b \in \mathbb{L}\f$) can be
  represented as a quadratic in the same space of input as follows,

  \f$ x^\top Q x + 2q^\top x + \rho \leq 0\f$,

  where, \f$ Q := A^\top J A\f$, \f$q := -A^\top J b \f$, \f$ \rho := b^\top J
  b \f$. Note that this is the same representation in case of \f$H \leftarrow
  A\f$ and \f$x^0 \leftarrow -b\f$. A subtle difference is columns of \f$A\f$ are
  not orthonormal.

  We refer to the primal space as x space and dual space as w space (where
  cones can be represented as quadratics). There is another space that we will
  refer as u-space. This space is where quadratics can be
  represented/regularized as a sphere (when \f$Q\f$ is postive semidefinite and
  there exists an \f$\overline{w}\f$ s.t. \f$Q\overline{w} = q\f$ and \f$ q^\top
  \overline{w} - \rho > 0\f$) or Lorentz cone (when \f$Q\f$ has exactly 1
  negative eigenvalue and there exists an \f$\overline{w}\f$ s.t.
  \f$Q\overline{w} = q\f$ and \f$ q^\top \overline{w} - \rho = 0\f$). We will
  give explicit formula of the linear transformation that maps a quadratic set
  (that satisfies the contions above) to a sphere or Lorentz cone.

  The cut is generated in the u-space and then represented in (or transformed
  into) the x-space for primal case or w-space for dual case (i.e., space of
  the input conic set). To do these transformations one must understand the
  relationship (linear transformations) between (1) \f$x\f$ and \f$w\f$ and (2)
  \f$w\f$ and \f$u\f$. Transformation between \f$x\f$ and \f$w\f$ is as
  follows,

  \f$ x = x^0 + Hw\f$ and \f$w = H^\top x - H^\top x^0\f$

  where \f$x^0\f$ is a point that satisfies \f$Ax=b\f$, and columns of \f$H\f$
  are an orthonormal basis for the null space of \f$A\f$.

  Relationship between \f$w\f$ and \f$u\f$ is not unique and depends on the
  input set. Let $(Q,q,\rho)$ be the quadratic representation of the input
  set. Then we are interested in the following cases only,

  <ul>

    <li> \f$Q\f$ is positive semidefinite and there exists an
         \f$\overline{w}\f$ s.t. \f$Q\overline{w} = q\f$ and \f$ q^\top
         \overline{w} - \rho > 0\f$.

    <li> \f$Q\f$ has exactly 1 negative eigenvalue and there exists an
         \f$\overline{w}\f$ s.t.  \f$Q\overline{w} = q\f$ and \f$ q^\top
         \overline{w} - \rho = 0\f$

  </ul>

  Quadratic set is a feasible ellipsoid in the first case and maps into a
  sphere. It is a cone in the second case and maps into a Lorentz cone in the
  u-space. Cases where the conditions do not hold are edge cases and we are not
  interested in them. Check references for details.

  We use eigenvalue decomposition of \f$Q\f$ for the transformations. Let \f$Q
  := VDV^\top\f$ be the eigenvalue decomposition of \f$Q\f$.

  Then quadratic can be written as follows,

  \f$ \left(w-\overline{w}\right)^\top V \tilde{D}^{\frac{1}{2}}
      \tilde{J}
      \tilde{D}^{\frac{1}{2}} V^\top \left(w - \overline{w}\right)
      \leq
      q^\top \overline{w} - \rho \f$.

  For the first case, linear transformation between w and u spaces is as
  follows,

  \f$u = \frac {\tilde{D}^{\frac{1}{2}} V^\top \left( w - \overline{w} \right)}
               {\sqrt{q^\top \overline{w} - \rho}} \f$.

  For the second case, linear transformation between w and u spaces is as
  follows,

  \f$u = \tilde{D}^{\frac{1}{2}} V^\top \left( w - \overline{w} \right) \f$.

  where \f$ \tilde{D}_{ii} := \left| D_{ii} \right| \f$.

  Inverse transformation between w and u spaces

  Case 1,

  \f$ w = \sqrt{q^\top \overline{w} - \rho} V \tilde{D}^{\frac{-1}{2}} u +
  \overline{w} \f$.

  Case 2,

  \f$ w = V \tilde{D}^{\frac{-1}{2}} u + \overline{w} \f$.

  Cut is generated in u-space since the quadratic representation is simpler and
  formulas are easier to track. To generate cut in the u-space, first the
  disjunction should be represented in the u-space.


  ## Representing disjunctions in the u-space

  ### Primal Form Input

  To represent the disjunctions in u-space we insert \f$ (x^0+Hw) \f$ for x in
  disjunction constraints. This will give disjunctions in terms of w (in
  w-space). We then insert the inverse of the complicated formula above
  (\f$u\f$--\f$w\f$ relationship) for w to obtain disjunctions in terms of
  \f$u\f$.

  Let the disjunction be \f$ a^\top x \geq \alpha \f$ and \f$a^\top x \leq
  \beta \f$. Then this disjunction can be represented in the u-space as
  follows,

  \f$
    \sqrt{q^\top \overline{w} - \rho} a^\top H V \tilde{D}^{\frac{-1}{2}} u
    \geq
    \alpha - a^\top x^0 - a^\top H \overline{w}.
  \f$

  Lower branch can be written similarly.

  ### Dual Form Input

  Let the disjunction be \f$ a^\top w \geq \alpha \f$ and \f$a^\top w \leq
  \beta \f$. Then this disjunction can be represented in the u-space as
  follows,

  \f$
    a^\top V \tilde{D}^{\frac{-1}{2}} u
    \geq
    \alpha - a^\top \overline{w}.
  \f$


  ## Computing cut

  Once disjunctions are represented in u-space, we compute \f$\tau\f$. This
  gives the cut in quadratic form in w space. We need to represent cut in dual
  form in the space of the input set (for users convenience, user is aware of
  the space of her problem only).

  First we represent the cut in the dual
  form in w space. For this we need to compute eigenvalue decomposition of
  $Q(\tau)$. Let \f$ Q(\tau) := V(\tau) D(\tau) V(\tau)^\top \f$. We need to
  compute a \f$\overline{w}(\tau)\f$ such that \f$ Q(\tau) \overline{w}(\tau) =
  q(\tau)\f$ and \f$ q(\tau)^\top \overline{w}(\tau) - \rho(\tau) = 0\f$. Such
  \f$\tau\f$ exists since \f$(Q(\tau), q(\tau), \rho(\tau))\f$ defines a cone and
  unique when \f$Q(\tau)\f$ is not singular.

  Cut can be represented in dual form in w space as follows,

  \f$ \tilde{D}^{\frac{1}{2}} V^\top (w+\overline{w}) \in \mathbb{L} \f$.

  ### Computing cut for inputs in primal form

  We use matrix \f$ x := x^0 + Hw \f$ and its inverse \f$ w = H^\top (x-x^0)
  \f$ to represent cut in the x space. The cut is as follows,

  \f$ \tilde{D}^{\frac{1}{2}} V^\top (H^\top (x-x^0) + \overline{w}) \in
  \mathbb{L} \f$,

  \f$ \tilde{D}^{\frac{1}{2}} V^\top H^\top x - \tilde{D}^{\frac{1}{2}} V^\top
  H^\top x^0 + \tilde{D}^{\frac{1}{2}} V^\top \overline{w}) \in
  \mathbb{L} \f$,

  \f$ cutA \leftarrow \tilde{D}^{\frac{1}{2}} V^\top H^\top \f$,

  \f$ cutb \leftarrow \tilde{D}^{\frac{1}{2}} V^\top
  H^\top x^0 - \tilde{D}^{\frac{1}{2}} V^\top \overline{w} \f$.


  ### Computing cut for inputs in dual form

  Dual form inputs are already in w space. Cut is as follows,

  \f$ cutA \leftarrow \tilde{D}^{\frac{1}{2}} V^\top \f$,

  \f$ cutb \leftarrow - \tilde{D}^{\frac{1}{2}} V^\top \overline{w} \f$.


  # References

  [1] On families of quadratic surfaces having fixed intersection with two
  hyperplanes

  [2] A conic representation of the convex hull of disjunctive sets and conic
  cuts for integer second order cone optimization

  [3] A Complete Characterization of Disjunctive Conic Cuts for Mixed Integer
  Second Order Cone Optimization
**/

class CglConicGD1Cut {
  /// Input cone type, primal or dual, need this to decide whether we have a
  /// null space basis (#matH_) in case of primal.
  CglConicInputType const input_type_;
  /// A solution to Ax = b, used in case of primal input.
  double * vecx0_;

  ///@name Input conic set in primal or dual form (in x-space or w-space)
  /// Primal form is \f$ Ax = b, x \in \mathbb{L} \f$.
  /// Dual form is \f$ Ax - b \in \mathbb{L} \f$.
  /// Assumes rows of A are linearly independent.
  //@{
  /// Number of rows of matrix A.
  int matA_num_rows_;
  /// Number of columns of primal matrix A.
  int matA_num_cols_;
  /// Pointer to primal matrix A, dense and column ordered.
  double * matA_;
  /// Pointer to primal right handside vector.
  double * vecb_;
  //@}

  ///@name Input conic set in dual form (in w-space)
  /// Input cone given as a quadratic
  /// \f$ x^\top Q x + 2q^\top x + \rho \leq 0 \f$.
  //@{
  /// Number of variables in dual space.
  int quad_num_cols_;
  /// Matrix \f$Q\f$ of quadric.
  double * matQ_;
  /// Vector \f$q\f$ of quadric.
  double * vecq_;
  /// Scalar \f$\rho\f$ of quadric.
  double rho_;
  /// Stores a solution to \f$Q w = q \f$.
  double * wbar_;
  //@}

  ///@name Eigenvalue Decomposition of Q.
  //@{
  /// Pointer to dense column ordered matrix \f$V\f$ of eigenvalue
  /// decomposition of #matQ_. Columns of #matV_ are eigenvectors of #matQ_.
  double * matV_;
  /// Pointer to eigenvectors of #matQ_, only diagonal elements are stored.
  double * matD_;
  //@}

  ///@name Disjunction Data
  /// Data of the last disjunction requested by the user.
  //@{
  /// index of disjunctive variable
  int dis_index_;
  /// right hand side of disjunction \f$ a^\top x \leq \alpha \f$
  double alpha_;
  /// right hand side of disjunction \f$ a^\top x \geq \beta \f$
  double beta_;
  /// Disjunction coefficient in w space. For dual input this is same as
  /// e_{dis_index} unit vector. For primal input, \f$e_{dis_index}^\top H\f$.
  double * dis_coef_in_w_;
  /// alpha in w space, same as #alpha_ for dual
  /// and quadratic input. #alpha_-#vecx0_[#dis_index_] for primal input.
  double alpha_in_w_;
  /// beta in w space, same as #beta_ for dual
  /// and quadratic input. #beta_-#vecx0_[#dis_index_] for primal input.
  double beta_in_w_;
  //@}

  ///@name Regularized input conic set
  //@{
  /// Regularized matrix \f$Q\f$.
  double * Jtilde_;
  // regularized vector q is 0.0, no need to store
  /// Regularized scalar \f$\rho\f$.
  double rho_tilde_;
  //@}

  ///@name Resulting cut in quadratic form in w-space
  //@{
  /// Parameter \f$\tau\f$ that gives a cone from convolution of
  /// input conic set and quadratic corresponding to disjunction.
  double tau_;
  /// Matrix \f$Q(\tau)\f$ in col major dense form.
  double * matQ_tau_;
  /// Vector \f$q(\tau)\f$
  double * vecq_tau_;
  /// Scalar \f$\rho(\tau)\f$.
  double rho_tau_;
  /// Stores a solution to \f$Q(\tau) w = q(\tau) \f$.
  double * wbar_tau_;
  //@}

  ///@name Eigenvalue Decomposition of #Qtau_.
  //@{

  /// Pointer to dense column ordered matrix \f$V(\tau)\f$ of eigenvalue
  /// decomposition of #matQtau_. Columns of #matVtau_ are eigenvectors of
  /// #matQtau_.
  double * matV_tau_;
  /// Pointer to eigenvectors of #matQtau_, only diagonal elements are stored.
  double * matD_tau_;
  //@}


  ///@name Resulting cut represented in dual form in the space of input set (x
  /// for primal, w for dual).

  /// User needs the cut represented in the same space as input.  Resulting cut
  /// is \f$ Ax - b in \mathbb{L}\f$ in the input space. Resulting cut is
  /// linear if #cutA_num_rows_ is 1, for this case the cut will be one of
  /// disjunctive constraints. This form is computed from quadratic
  /// representation using eigenvalue decomposition of #Qtau_ (and matrix
  /// #matH_ for input in primal form).
  //@{
  int cutA_num_rows_;
  int cutA_num_cols_;
  double * cutA_;
  double * cutb_;
  //@}

  /// Null space basis of #primalA_ in case of input in primal form. Dense,
  /// column ordered. Columns are orthonormal basis vectors.
  // todo(aykut) H is row ordered in computations, fix it.
  double * matH_;

  /// True if cut process revealed that problem is infeasible.
  bool infeasible_;
  /// True if cut generated successfully, 0 otherwise.
  bool success_;
  ///@name Convert to quadratic form.
  //@{
  /// Converts/represents the input conic set in quadratic form, i.e., #matQ_,
  /// #vecq_, #rho_.
  void compute_quadratic();
  /// Compute matrix #matH_.
  void compute_matrixH();
  /// Compute matrix #matQ_.
  void compute_matrixQ();
  /// Compute vector #vecq_.
  void compute_vectorq();
  /// Compute $rho_.
  void compute_rho();
  //@}

  /// Compute eigenvalue decompostion of matrix #matQ_ and
  /// store results in #matV_ and #matD_.
  void decompose_matrixQ();
  void compute_cut();

  /// Compute disjunctiv hyperplane coefficient and right hand sides in w
  /// space.
  void compute_disjunction_in_w();
  // compute tau value that yields a cone
  void compute_tau();
  //void compute_disjunction();
  // compute rho(tau), rho_tau_
  void compute_rho_tau();
  // compute q(tau), q_tau_
  void compute_q_tau();
  // compute Q(tau), Q_tau_
  void compute_Q_tau();
  /// Compute eigenvalue decomposition of matrix #matQ_tau_.
  void decompose_matrixQtau();

  // compute matrix as a part of the cut to add to model
  //void compute_cutA();
  // compute right-hand-side that corresponds to NewA
  //void compute_cutb();

  ///@name Debug printing functions
  //@{
  /// Print matrix, row major 0, col major 1.
  void print_matrix(int major, int num_rows, int num_cols,
                    double const * matrix,
                    char const * name = "matrix") const;
  /// Print vector.
  void print_vector(int n, double const * vector,
                    char const * name = "vector") const;
  /// Print scalar.
  void print_scalar(double value, char const * name) const;
  //@}


  // compute cone at tau
  // void coneAtTau();
  // solve quadratic formula for its largest roots
public:
  ///@name Constructors and destructor
  //@{
  /// Generate cut from cones given in primal form.
  /// \f$ Ax = b, x in L \f$
  CglConicGD1Cut(CglConicInputType input_type, CoinPackedMatrix const * A,
                 double const * b,
                 double const * x0 = NULL);
  /// Generate cut from cones given in dual form,
  /// \f$ x^\top Q x + 2q^\top x + \rho \leq 0 \f$.
  CglConicGD1Cut(int n, double const * Q, double const * q, double rho);
  /// Destructor.
  ~CglConicGD1Cut();
  //@}

  /// disjunction \f$ x_{dis_index} \leq \alpha,
  /// x_{dis_index} \geq \beta \f$.
  void generateCut(int dis_index, double alpha, double beta);

  /// Returns true when the problem is infeasible.
  bool infeasible() const;
  /// Returns true if cut generated successfuly.
  bool success() const;

  // number of rows of the linear part of the cut
  // ie num of rows of the new A matrix.
  int getNumRows() const;
  // number of columns in the cut
  int getNumCols() const;

  // Get cut.
  // return linear part of the cut, constraint matrix.
  double const * getCutA() const;
  // return right hand side of the linear part of the cut.
  double const * getCutb() const;

  /// Get tau
  double tau() const { return tau_; }

  /// print cut
  void print_cut() const;

private:
  /// Disable default constructor.
  CglConicGD1Cut();
  /// Disable copy constructor.
  CglConicGD1Cut(CglConicGD1Cut const & other);
  /// Disable copy assignment operator.
  CglConicGD1Cut & operator=(CglConicGD1Cut const & other);
};

/// Define struct for eigenvalue--eigenvector pair.
struct EigenPair {
  double value_;
  double * vector_;
};

/// Create function object to help sorting.
struct EigenLess {
  bool operator()(EigenPair const * a, EigenPair const * b) {
    return a->value_ < b->value_;
  }
};


#endif
