#include "CglConicGD1Cut.hpp"
#include <vector>
#include <numeric>

extern "C" {
  // blas routines
  void dcopy_(int*, double*, int*, double*, int*);
  void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*,
              double*, double*, int*);
  void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*,
              double*, int*, double*, double*, int*);
  void dsyrk_(char*, char*, int*, int*, double*, double*, int*, double*,
              double*, int*);
  void dsyr_(char*, int*, double*, double*, int*, double*, int*);
  double ddot_(int*, double*, int*, double*, int*);
  void daxpy_(int*, double*, double*, int*, double*, int*);

  // lapack routines
  void dsysv_(char *uplo, int *n, int *nrhs, double *a,
              int *lda, int * ipiv, double *b, int * ldb,
              double * work, int *lwork, int * info);
  void dsyev_(char *jobz, char *uplo, int *n, double *a,
              int *lda, double *w, double *work, int *lwork,
              int *info);
  void dgesvd_(char *jobu, char *jobvt, int *m, int *n,
               double *a, int *lda, double *S, double *U,
               int *ldu, double *vt, int *ldvt, double *work,
               int *lwork, int *info);
}


#define EPS 1e-4

// helper functions which are not class members
static void svDecompICL(int m, int n, double * A, double * VT);
static void solveSM(int n, double *A, double *b, double *x);
static void eigDecompICL(int n, double * L, double * eig);
// returns largest root of a quadratic formula
static double quad_formula(double quad_term, double lin_term,
                    double ind_term);
// compute roots of a quadratic formula
static void roots(double a, double b, double c, double & root1,
                  double & root2);

CglConicGD1Cut::CglConicGD1Cut(CglConicInputType input_type,
               CoinPackedMatrix const * A, double const * b,
               double const * x0) :
  input_type_(input_type) {
  // if input_type is DUAL_FORM, just ignore x0.
  matA_num_rows_ = A->getNumRows();
  matA_num_cols_ = A->getNumCols();
  int m = matA_num_rows_;
  int n = matA_num_cols_;
  // get A in dense col ordered form.
  matA_ = new double[n*m]();
  bool col_ordered = A->isColOrdered();
  int major_dim = A->getMajorDim();
  int minor_dim = A->getMinorDim();
  for (int i=0; i<major_dim; ++i) {
    int first = A->getVectorFirst(i);
    int size = A->getVectorSize(i);
    int const * ind = A->getIndices() + first;
    double const * val = A->getElements() + first;
    for (int j=0; j<size; ++j) {
      if (col_ordered)
        // column i, row ind[j]
        matA_[ind[j] + i*minor_dim] = val[j];
      else {
        // column ind[j], row i
        matA_[i + ind[j]*major_dim] = val[j];
      }
    }
  }
  // get b
  vecb_ = new double[m];
  std::copy(b, b+m, vecb_);
  // get x0
  vecx0_ = NULL;
  if (input_type_ == PRIMAL_FORM) {
    vecx0_ = new double[n];
    std::copy(x0, x0+n, vecx0_);
  }
  else {
    // in dual form, x0 is -b
    vecx0_ = new double[m];
    for (int j=0; j<m; ++j) vecx0_[j] = -vecb_[j];
  }
  quad_num_cols_ = -1;
  matQ_ = NULL;
  vecq_ = NULL;
  rho_ = -1.0;
  wbar_ = NULL;
  matV_ = NULL;
  matD_ = NULL;
  dis_index_ = -1;
  alpha_ = 0.0;
  beta_ = 0.0;
  dis_coef_in_w_ = NULL;
  alpha_in_w_ = 0.0;
  beta_in_w_ = 0.0;
  Jtilde_ = NULL;
  rho_tilde_ = -1.0;
  tau_ = 0.0;
  matQ_tau_ = NULL;
  vecq_tau_ = NULL;
  rho_tau_ = 0.0;
  wbar_tau_ = NULL;
  matV_tau_ = NULL;
  matD_tau_ = NULL;
  cutA_num_rows_ = -1;
  cutA_num_cols_ = -1;
  cutA_ = NULL;
  cutb_ = NULL;
  infeasible_ = false;
  success_ = false;
  // compute quadratic represtation
  compute_quadratic();
  // compute eigenvalue decomposition
  decompose_matrixQ();
}


CglConicGD1Cut::CglConicGD1Cut(int n, double const * Q,
                               double const * q, double rho)
  : input_type_(QUAD_FORM) {
}

void CglConicGD1Cut::compute_quadratic() {
  if (input_type_ == PRIMAL_FORM) {
    quad_num_cols_ = matA_num_cols_ - matA_num_rows_;
  }
  else {
    quad_num_cols_ = matA_num_cols_;
  }
  compute_matrixH();
  compute_matrixQ();
  compute_vectorq();
  compute_rho();
}

// compute matrix H from matrix A, H is null mat of A.
void CglConicGD1Cut::compute_matrixH() {
  if (input_type_ == DUAL_FORM) {
    matH_ = matA_;
    return;
  }
  int num_cols = matA_num_cols_;
  int num_rows = matA_num_rows_;
  // copy A to a working array
  double * tempA = new double[num_cols*num_rows];
  int blas_nm = num_cols*num_rows;
  int blas_one = 1;
  dcopy_(&blas_nm, matA_, &blas_one, tempA, &blas_one);
  // right hand side singular vectors of A
  double * VT = new double[num_cols*num_cols];
  svDecompICL(num_rows, num_cols, tempA, VT);
  matH_ = new double[(num_cols-num_rows)*num_cols]();
  // Take the last n-m columns of V, lapack returns V^T
  for(int i=0; i<(num_cols-num_rows); ++i) {
    dcopy_(&num_cols, (VT+num_rows+i), &num_cols, (matH_+i*num_cols), &blas_one);
  }
  delete[] tempA;
  delete[] VT;
  //print_matrix(1, num_cols, num_cols-num_rows, matH_, "H");
}

// negate first row of H. Multiply H^T with this matrix.

// todo(aykut) matH_ is stored as row major in the memory. It should be col
// ordered.

// A is col ordered and H is row oredered, in dual form we point H to A, this
// will create problems.

void CglConicGD1Cut::compute_matrixQ() {
  // number of rows of matrix H
  int m;
  // number of columns of matrix H
  int n;
  if (input_type_ == PRIMAL_FORM) {
    m = matA_num_cols_;
    n = matA_num_cols_ - matA_num_rows_;
  }
  else {
    m = matA_num_rows_;
    n = matA_num_cols_;
  }
  matQ_ = new double[n*n]();
  // Temporary array, stores the first row of H
  double * d = new double[n];
  int blas_one = 1;
  int blas_mm1 = m-1;
  dcopy_(&n, matH_, &m, d, &blas_one);
  // Temporary array, exlcudes the first row of H
  double * A = new double[(m-1)*n];
  for(int i=0; i<n; i++) {
    dcopy_(&blas_mm1, matH_+i*m+1, &blas_one, A+i*(m-1), &blas_one);
  }
  // computes Q = A^TA - dd^T =  H^T J H
  char blas_type_c = 'C';
  char blas_type_n = 'N';
  char blas_upper = 'U';
  double blas_double_one = 1.0;
  double blas_double_neg_one = -1.0;
  double blas_zero = 0.0;
  dsyrk_(&blas_upper, &blas_type_c, &n, &blas_mm1,
              &blas_double_one, A, &blas_mm1, &blas_zero, matQ_, &n);
  dsyr_(&blas_upper, &n, &blas_double_neg_one, d,
        &blas_one, matQ_, &n);
  delete[] d;
  delete[] A;
  //print_matrix(1, n, n, matQ_, "Q");
}

// for both primal and dual
// Q is H^\top J H
// q is H^\top J x0
// rho is x0^\top J x0
void CglConicGD1Cut::compute_vectorq() {
  // number of rows of matrix H
  int m;
  // number of columns of matrix H
  int n;
  if (input_type_ == PRIMAL_FORM) {
    m = matA_num_cols_;
    n = matA_num_cols_ - matA_num_rows_;
  }
  else {
    m = matA_num_rows_;
    n = matA_num_cols_;
  }
  vecq_ = new double[n]();
  vecx0_[0] = -1.0*vecx0_[0];

  char blas_type_c = 'C';
  double blas_double_one = 1.0;
  int blas_one = 1;
  double blas_zero = 0.0;
  dgemv_(&blas_type_c, &m, &n, &blas_double_one, matH_, &m,
               vecx0_, &blas_one, &blas_zero, vecq_, &blas_one);
  // reverse negation of x0[0]
  vecx0_[0] = -1.0*vecx0_[0];
  //print_vector(n, vecq_, "q");
}

void CglConicGD1Cut::compute_rho() {
  // number of rows of matrix H
  int m;
  if (input_type_ == PRIMAL_FORM)
    m = matA_num_cols_;
  else {
    m = matA_num_rows_;
  }
  rho_ = - (vecx0_[0]*vecx0_[0]);
  int blas_mm1 = m-1;
  int blas_one = 1;
  rho_ += ddot_(&blas_mm1, vecx0_+1, &blas_one, vecx0_+1, &blas_one);
  //print_scalar(rho_, "rho");
}

void CglConicGD1Cut::decompose_matrixQ() {
  success_ = true;
  int n;
  if (input_type_ == PRIMAL_FORM)
    n = matA_num_cols_ - matA_num_rows_;
  else {
    n = matA_num_cols_;
  }
  // compute q^T Q^-1 q - \rho
  double * Qq = new double[n];
  solveSM(n, matQ_, vecq_, Qq);
  // copy matrix Q
  matV_ = new double[n*n];
  int blas_one = 1;
  int blas_n_sqr = n*n;
  dcopy_(&blas_n_sqr, matQ_, &blas_one, matV_, &blas_one);
  // vector with the eigenvalues of Q
  matD_ = new double[n]();
  // compute the eigenvalue decomposition
  eigDecompICL(n, matV_, matD_);
  // in case nonpositive eigenvalues bail out
  for (int i=0; i<n; ++i) {
    if (matD_[i]<1e-3) {
      std::cout << "Q is not positive definite!" << std::endl;
      success_ = false;
      break;
    }
  }
  delete [] Qq;
  // push eigenvalue--eigenvector pairs into vector
  std::vector<EigenPair*> epair;
  for (int i=0; i<n; ++i) {
    EigenPair * curr = new EigenPair();
    curr->value_ = matD_[i];
    curr->vector_ = matV_ + i*n;
    epair.push_back(curr);
  }
  // sort eigenvalue eigenvector pairs
  std::sort(epair.begin(), epair.end(), EigenLess());

  // restore matD_ and matV_
  double * newMatV = new double[n*n];
  std::vector<EigenPair*>::const_iterator it;
  int k=0;
  for (it=epair.begin(); it!=epair.end(); ++it) {
    matD_[k] = (*it)->value_;
    std::copy((*it)->vector_, (*it)->vector_+n, newMatV+k*n);
    k++;
  }
  delete[] matV_;
  matV_ = newMatV;
  newMatV = NULL;

  // free epair
  std::vector<EigenPair*>::iterator iit;
  for (iit=epair.begin(); iit!=epair.end(); ++iit) {
    delete *iit;
  }
  epair.clear();
  //print_matrix(1, n, n, matV_, "V");
  //print_vector(n, matD_, "D");
}

void CglConicGD1Cut::generateCut(int dis_index, double alpha, double beta) {
  dis_index_ = dis_index;
  alpha_ = alpha;
  beta_ = beta;
  // compute disjunction
  // compute tau
  compute_tau();
  if (infeasible_) {
    // computing tau revealed that problem is infeasible, i.e. intersection of
    // input set with disjunctive half spaces is empty.
    return;
  }
  if (!success_) {
    // computing tau failed. tau might be close to -1.
    return;
  }
  if (cutA_num_rows_==1) {
    // only one of the disjunctive hyperplane intersects input set,
    // cut is linear and already computed.
    return;
  }
  // tau computed successfully
  compute_disjunction_in_w();
  compute_Q_tau();
  compute_q_tau();
  compute_rho_tau();
  decompose_matrixQtau();
  compute_cut();
}

void CglConicGD1Cut::compute_disjunction_in_w() {
  int n;
  int m;
  if (input_type_ == PRIMAL_FORM) {
    m = matA_num_cols_;
    n = matA_num_cols_ - matA_num_rows_;
  }
  else {
    m = matA_num_rows_;
    n = matA_num_cols_;
  }
  dis_coef_in_w_ = new double[n]();
  if (input_type_ == PRIMAL_FORM) {
    for (int i=0; i<n; ++i) {
      dis_coef_in_w_[i] = matH_[i*m + dis_index_];
    }
    alpha_in_w_ = alpha_ - vecx0_[dis_index_];
    beta_in_w_ = beta_ - vecx0_[dis_index_];
  }
  else {
    dis_coef_in_w_[dis_index_] = 1.0;
    alpha_in_w_ = alpha_;
    beta_in_w_ = beta_;
  }

  // normalize
  double norm = 0.0;
  norm = std::inner_product(dis_coef_in_w_, dis_coef_in_w_+n, dis_coef_in_w_, 0.0);
  norm = sqrt(norm);
  for (int i=0; i<n; ++i) {
    dis_coef_in_w_[i] = dis_coef_in_w_[i]/norm;
  }
  alpha_in_w_ = alpha_in_w_ / norm;
  beta_in_w_ = beta_in_w_ / norm;

  //print_vector(n, dis_coef_in_w_, "a_in_w");
  //print_scalar(alpha_in_w_, "alpha_in_w");
  //print_scalar(beta_in_w_, "beta_in_w");
}

void CglConicGD1Cut::compute_tau() {
  success_ = false;
  // number of rows of matrix H
  int m;
  // number of columns of matrix H
  int n;
  if (input_type_ == PRIMAL_FORM) {
    m = matA_num_cols_;
    n = matA_num_cols_ - matA_num_rows_;
  }
  else {
    m = matA_num_rows_;
    n = matA_num_cols_;
  }
  // compute disjunction coef in regularized space.
  // for PRIMAL FORM, disjunction is

  // (1) find \overline{w} by solving Q \overline{w} = q
  // (2) find term1 which is sqrt{q^\top w - rho}
  // (3) compute term1 a^\top H, get corresponding row of H.
  // (4) compute term1 a^\top H V, vector matrix multiplication
  // (5) compute term1 a^\top H V D^{-1/2}, elementwise vector vector multip
  //     this is disjunction coef in u space.
  // (6) find alpha in u space, alpha - a^\top x0 + a^\top H \overline{w}
  // (7) find beta in u space, beta - a^\top x0 + a^\top H \overline{w}
  // (8) check whether disjunctive hyperplanes intersect with the conic
  //     set, if none intersects then problem is infeasible. if only one
  //     intersects then resulting cut is linear. If both intersects the result
  //     is a conic cut.
  // (9) compute quad coef, (alpha-beta)^2 / 4
  // (10) compute linear coef, (1- alpha beta)
  // (11) constant term is 1.
  // (12) solve second degree polynomial to get tau


  // (1) find \overline{w} by solving Q \overline{w} = q
  wbar_ = new double[n]();
  solveSM(n, matQ_, vecq_, wbar_);
  //print_vector(n, wbar_, "wbar");

  // (2) find term1 which is sqrt{q^\top w - rho}
  double term1 = std::inner_product(vecq_, vecq_+n, wbar_, -rho_);
  term1 = sqrt(term1);
  //print_scalar(term1, "term1");

  // (3) compute term1 a^\top H, get corresponding row of H.
  double * aH = new double[n];
  double * dis_coef = new double[n];
  // copy dis_index_ row of H
  for (int i=0; i<n; ++i) {
    aH[i] = matH_[i*m + dis_index_];
  }
  // copy dis_index_ row of H
  for (int i=0; i<n; ++i) {
    dis_coef[i] = term1*aH[i];
  }

  // (4) compute term1 a^\top H V, vector matrix multiplication
  char blas_type_c = 'C';
  double blas_double_one = 1.0;
  int blas_one = 1;
  double blas_zero = 0.0;
  dgemv_(&blas_type_c, &n, &n, &blas_double_one,
              matV_, &n, dis_coef, &blas_one, &blas_zero, dis_coef, &blas_one);

  // (5) compute term1 a^\top H V D^{-1/2}, elementwise vector vector multip
  //     this is disjunction coef in u space.
  for (int i=0; i<n; ++i) {
    dis_coef[i] = dis_coef[i]*(1.0/sqrt(fabs(matD_[i])));
  }

  // (6) find alpha in u space, alpha - a^\top x0 + a^\top H \overline{w}
  double aHwbar = 0.0;
  aHwbar = std::inner_product(aH, aH+n, wbar_, 0.0);
  double alpha = alpha_ - vecx0_[dis_index_] + aHwbar;

  // (7) find beta in u space, beta - a^\top x0 + a^\top H \overline{w}
  double beta = beta_ - vecx0_[dis_index_] + aHwbar;

  // (7.1) normalize disjunction
  double norm_a = 0.0;
  norm_a = std::inner_product(dis_coef, dis_coef+n, dis_coef, 0.0);
  norm_a = sqrt(norm_a);

  //print_vector(n, dis_coef, "unscaled a");

  for (int i=0; i<n; ++i) {
    dis_coef[i] = dis_coef[i]/norm_a;
  }
  alpha = alpha/norm_a;
  beta = beta/norm_a;

  //print_vector(n, dis_coef, "a");
  //print_scalar(alpha, "alpha");
  //print_scalar(beta, "beta");
  // (8) check whether disjunctive hyperplanes intersect with the conic
  //     set, if none intersects then problem is infeasible. if only one
  //     intersects then resulting cut is linear. If both intersects the result
  //     is a conic cut.
  bool alpha_intersects = (alpha*alpha <= 1);
  bool beta_intersects = (beta*beta <= 1);
  if (alpha_intersects && beta_intersects) {
    // both hyperplanes intersects regularized quadric
    // cool, keep working on it.
  }
  else if (alpha_intersects) {
    // only ax=alpha intersects quadric
    // cut is a^top x - alpha in L
    success_ = true;
    infeasible_ = false;
    cutA_num_rows_ = 1;
    cutA_num_cols_ = m;
    cutA_ = new double[m]();
    cutA_[dis_index_] = 1.0;
    cutb_ = new double[1];
    cutb_[0] = alpha_;
    delete[] aH;
    delete[] dis_coef;
    return;
  }
  else if(beta_intersects) {
    // only ax=beta intersects quadric
    // cut is -a^top x + beta in L
    success_ = true;
    infeasible_ = false;
    cutA_num_rows_ = 1;
    cutA_num_cols_ = m;
    cutA_ = new double[m]();
    cutA_[dis_index_] = -1.0;
    cutb_ = new double[1];
    cutb_[0] = -beta_;
    delete[] aH;
    delete[] dis_coef;
    return;
  }
  else {
    // none of the hyperplanes intersect the quadric, problem is infeasible.
    success_ = true;
    infeasible_ = true;
    delete[] aH;
    delete[] dis_coef;
    return;
  }

  // (9) compute quad coef, (alpha-beta)^2 / 4
  double quad_coef = (alpha-beta)*(alpha-beta);

  // (10) compute linear coef, (1- alpha beta)
  double lin_coef = 4.0*(1.0-alpha*beta);

  // (11) constant term is 1.
  double const_term = 4.0;

  // (12) solve second degree polynomial to get tau

  if (lin_coef*lin_coef < 4.0*quad_coef*const_term) {
    //std::cerr << "Imaginary root!" << std::endl;
    success_ = false;
  }
  else {
    // quadformula returns 1 when the roots are imaginary which will
    // not happen if quadric intersects with disjunction hyperplanes
    tau_ = quad_formula(quad_coef, lin_coef, const_term);
    if (tau_ > -1.1) {
      // numerical problems might arise when tau is greater than -1.1.
      // abort.
      success_ = false;
    }
    else {
      success_ = true;
    }
  }
  //print_scalar(tau_, "tau");
  delete[] aH;
  delete[] dis_coef;
}

// compute rho(tau), rho_tau_
void CglConicGD1Cut::compute_rho_tau() {
  rho_tau_ += tau_*alpha_in_w_*beta_in_w_;
  //print_scalar(rho_tau_, "rho(tau)");
}

// compute q(tau), q_tau_
void CglConicGD1Cut::compute_q_tau() {
  // number of columns of matrix H
  int n;
  if (input_type_ == PRIMAL_FORM) {
    n = matA_num_cols_ - matA_num_rows_;
  }
  else {
    n = matA_num_cols_;
  }
  vecq_tau_ = new double[n]();
  std::copy(vecq_, vecq_+n, vecq_tau_);
  double aux = -0.5*tau_*(alpha_in_w_ + beta_in_w_);
  int blas_one = 1;
  daxpy_(&n, &aux, dis_coef_in_w_, &blas_one, vecq_tau_, &blas_one);
  //print_vector(n, vecq_tau_, "q(tau)");
}

 // compute Q(tau), Q_tau_
void CglConicGD1Cut::compute_Q_tau() {
  // number of rows of matrix H
  int m;
  // number of columns of matrix H
  int n;
  if (input_type_ == PRIMAL_FORM) {
    m = matA_num_cols_;
    n = matA_num_cols_ - matA_num_rows_;
  }
  else {
    m = matA_num_rows_;
    n = matA_num_cols_;
  }
  matQ_tau_ = new double[n*n];
  int blas_n_sqr = n*n;
  int blas_one = 1;
  dcopy_(&blas_n_sqr, matQ_, &blas_one, matQ_tau_, &blas_one);
  char blas_upper = 'U';
  dsyr_(&blas_upper, &n, &tau_, dis_coef_in_w_, &blas_one,
             matQ_tau_, &n);
  // copy upper triangular to lower triangular part
  // for each column
  for (int i=0; i<n; ++i) {
    // for each row
    for (int j=0; j<i; ++j) {
      matQ_tau_[j*n+i] = matQ_tau_[i*n+j];
    }
  }
  //print_matrix(1, n, n, matQ_tau_, "Q(tau)");
}


void CglConicGD1Cut::decompose_matrixQtau() {
  int n;
  if (input_type_ == PRIMAL_FORM)
    n = matA_num_cols_ - matA_num_rows_;
  else {
    n = matA_num_cols_;
  }
  wbar_tau_ = new double[n];
  solveSM(n, matQ_tau_, vecq_tau_, wbar_tau_);
  //print_vector(n, wbar_tau_, "wbar_tau");
  // copy matrix Q
  matV_tau_ = new double[n*n];
  int blas_n_sqr = n*n;
  int blas_one = 1;
  dcopy_(&blas_n_sqr, matQ_tau_, &blas_one, matV_tau_, &blas_one);
  // vector with the eigenvalues of Q
  matD_tau_ = new double[n]();
  // compute the eigenvalue decomposition
  eigDecompICL(n, matV_tau_, matD_tau_);

  // check whether eigenvalues are zero.
  for (int i=0; i<n; ++i) {
    if (matD_tau_[i]<0.001 && matD_tau_[i]>-0.001) {
      std::cout << "Zero eigenvalue in $Q(\tau)$." << std::endl;
      success_ = false;
    }
  }

  // at most 1 eigenvalue must be negative
  {
    int num_neg_eigen = 0;
    for (int i=0; i<n; ++i) {
      if (matD_tau_[i]<0.0) {
        num_neg_eigen++;
      }
    }
    if (num_neg_eigen>1) {
      std::cerr << "Number of negative eigenvalues should be at most 1!"
                << std::endl;
      success_ = false;
      return;
    }
  }

  // push eigenvalue--eigenvector pairs into vector
  std::vector<EigenPair*> epair;
  for (int i=0; i<n; ++i) {
    EigenPair * curr = new EigenPair();
    curr->value_ = matD_tau_[i];
    curr->vector_ = matV_tau_ + i*n;
    epair.push_back(curr);
  }
  // sort eigenvalue eigenvector pairs
  std::sort(epair.begin(), epair.end(), EigenLess());

  // restore matD_tau_ and matV_tau_
  double * newMatV = new double[n*n];
  std::vector<EigenPair*>::const_iterator it;
  int k=0;
  for (it=epair.begin(); it!=epair.end(); ++it) {
    matD_tau_[k] = (*it)->value_;
    std::copy((*it)->vector_, (*it)->vector_+n, newMatV+k*n);
    k++;
  }
  delete[] matV_tau_;
  matV_tau_ = newMatV;
  newMatV = NULL;

  //print_matrix(1, n, n, matV_tau_, "V(tau)");
  //print_vector(n, matD_tau_, "D(tau)");

  // free epair
  std::vector<EigenPair*>::iterator iit;
  for (iit=epair.begin(); iit!=epair.end(); ++iit) {
    delete *iit;
  }
  epair.clear();
}

// compute matrix as a part of the cut to add to model
void CglConicGD1Cut::compute_cut() {
  // number of rows of matrix H
  int m;
  // number of columns of matrix H
  int n;
  if (input_type_ == PRIMAL_FORM) {
    m = matA_num_cols_;
    n = matA_num_cols_ - matA_num_rows_;
  }
  else {
    m = matA_num_rows_;
    n = matA_num_cols_;
  }
  cutA_num_rows_ = n;
  cutA_num_cols_ = m;

  // cut is V(tau) Dtilde(tau) ^{1/2} (H^\top (x-x0) + wbar_tau)

  // compute Dtilde^{1/2}
  double * sqrtDtau = new double[n*n]();
  for (int i=0; i<n; ++i) {
    sqrtDtau[i*n+i] = sqrt(fabs(matD_tau_[i]));
  }
  //print_matrix(1, n, n, sqrtDtau, "sqrtDtau");

  // compute sqrtQtau =  V Dtilde^{1/2}
  double * sqrtQtau = new double[n*n]();
  char blas_type_c = 'C';
  char blas_type_n = 'N';
  double blas_double_one = 1.0;
  double blas_double_neg_one = -1.0;
  int blas_one = 1;
  double blas_zero = 0.0;
  dgemm_(&blas_type_c, &blas_type_c, &n, &n, &n, &blas_double_one, matV_tau_,
              &n, sqrtDtau, &n, &blas_zero, sqrtQtau, &n);

  // compute cutA
  if (input_type_ == PRIMAL_FORM) {
    cutA_ = new double[n*m];
    // multiply with H^\top
    dgemm_(&blas_type_c, &blas_type_c, &n, &m, &n,
                &blas_double_one, sqrtQtau, &n, matH_, &m, &blas_zero, cutA_, &n);
    // cutb is cutA x^0 - sqrtQtau wbar_tau

    // print_matrix(1, n, n, sqrtQtau, "sqrtQtau");

    cutb_ = new double[n]();
    dgemv_(&blas_type_n, &n, &m, &blas_double_one,
           cutA_, &n, vecx0_, &blas_one, &blas_zero, cutb_, &blas_one);

    // print_vector(n, cutb_, "cutA x0");

    dgemv_(&blas_type_n, &n, &n, &blas_double_neg_one,
           sqrtQtau, &n, wbar_tau_, &blas_one, &blas_double_one,
           cutb_, &blas_one);

    // print_vector(n, cutb_, "cutb");
  }
  else {
    cutA_ = sqrtQtau;
    sqrtQtau = NULL;
    cutb_ = new double[n];
    // matrix vector multiplication
    // cutb is cutA wbar
  }

  // decide which branch the cut is in
  // branch := ( V(tau) \tilde{D}(tau) ^{1/2} )_1 (-wbar + wbar_tau)
  // if branch > 0.0 then cutA_1 = -cutA_1 and cutb[0] = -cutb[0]
  // first column of sqrtQtau times (-wbar + wbar_tau)
  double branch = 0.0;

  //print_matrix(1, n, n, sqrtQtau, "sqrtQtau");
  //print_vector(n, wbar_, "wbar");
  //print_vector(n, wbar_tau_, "wbar_tau_");

  for (int i=0; i<n; ++i) {
    branch += sqrtQtau[i]*(-wbar_[i]+wbar_tau_[i]);
  }
  //print_scalar(branch, "dec");
  if (branch < 0.0) {
    for (int i=0; i<m; ++i) {
      cutA_[i*n] = -cutA_[i*n];
    }
    cutb_[0] = -cutb_[0];
  }
  delete[] sqrtDtau;
  delete[] sqrtQtau;
}

bool CglConicGD1Cut::success() const {
  return success_;
}

bool CglConicGD1Cut::infeasible() const {
  return infeasible_;
}

// number of rows of the linear part of the cut
// ie num of rows of the new A matrix.
int CglConicGD1Cut::getNumRows() const {
  return cutA_num_rows_;
}

// number of columns in the cut
int CglConicGD1Cut::getNumCols() const {
  return cutA_num_cols_;
}

CglConicGD1Cut::~CglConicGD1Cut() {
  if (vecx0_) {
    delete[] vecx0_;
  }
  if (matA_) {
    delete[] matA_;
  }
  if (vecb_) {
    delete[] vecb_;
  }
  if (matQ_) {
    delete[] matQ_;
  }
  if (vecq_) {
    delete[] vecq_;
  }
  if (wbar_) {
    delete[] wbar_;
  }
  if (matV_) {
    delete[] matV_;
  }
  if (matD_) {
    delete[] matD_;
  }
  if (dis_coef_in_w_) {
    delete[] dis_coef_in_w_;
  }
  if (Jtilde_) {
    delete[] Jtilde_;
  }
  if (matQ_tau_) {
    delete[] matQ_tau_;
  }
  if (vecq_tau_) {
    delete[] vecq_tau_;
  }
  if (wbar_tau_) {
    delete[] wbar_tau_;
  }
  if (matV_tau_) {
    delete[] matV_tau_;
  }
  if (matD_tau_) {
    delete[] matD_tau_;
  }
  if (cutA_) {
    delete[] cutA_;
  }
  if (cutb_) {
    delete[] cutb_;
  }
  if (matH_) {
    delete[] matH_;
  }
}

// return linear part of the cut, constraint matrix.
double const * CglConicGD1Cut::getCutA() const {
  return cutA_;
}

// return right hand side of the linear part of the cut.
double const * CglConicGD1Cut::getCutb() const {
  return cutb_;
}

///////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS WHICH ARE NOT CLASS MEMBERS
///////////////////////////////////////////////////////////////////////////////

static void svDecompICL(int m, int n, double * A, double * VT){
  /* Do not compute U */
  char jobu = 'N';
  /* Compute all the rows of VT */
  char jobvt = 'A';
  int lda = m; //Leading dimension of A
  double * s = (double *) malloc(m*sizeof(double)); //Singular values of A
  //We are not computing U, so not need to assign memory
  double * u = NULL;
  int ldu = 1;
  //Right hand side singular vectors of A in VT
  int ldvt = n;
  // Working space size needed by lapack
  double worksize = 0.0;
  int lwork = -1; // Tells lapack to query only what is the space needed
  int info = 0;
  //  printf("fails Singular query\n");
  // Query what is the optimal size for the computation of the singular values
  dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, u, &ldu,
          VT, &ldvt, &worksize, &lwork, &info);
  // Initilize the working space needed by lapack
  lwork = (int) worksize;
  double * work = (double *) malloc(lwork*sizeof(double));
  //   printf("fails Singular computing\n");
  // Here is where the actually computation is happening
  dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, u, &ldu,
          VT, &ldvt, work, &lwork, &info);
  free(work);
  free(s);
}

// solve symmetric system Ax=b
static void solveSM(int n, double *A, double *b, double *x) {
  //printVector(b, n,"b");
  char uplo = 'U';
  int nrhs = 1;
  double *copyA = new double[n*n];
  memcpy ( copyA, A, (n * n) * sizeof(double) );
  int * ipiv = new int[n];
  std::fill_n(ipiv, n, 0);
  memcpy ( x, b, n * sizeof(double) );
  double worksize = 0.0;
  int lwork = -1;
  int info;
  dsysv_(&uplo, &n, &nrhs, copyA, &n, ipiv, x, &n,
         &worksize, &lwork, &info);
  lwork = (int) worksize;
  double * work = new double[lwork];
  dsysv_(&uplo, &n, &nrhs, copyA, &n, ipiv, x, &n,
         work, &lwork, &info);
  delete[] work;
  delete[] ipiv;
  delete[] copyA;
}

// compute eigenvalues of lower triangular matrix L
static void eigDecompICL(int n, double * L, double * eig) {
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  int lwork = -1;
  double worksize[1];
  dsyev_(&jobz, &uplo, &n, L, &n, eig, worksize, &lwork, &info);
  lwork = (int) worksize[0];
  double * work = new double[lwork];
  dsyev_(&jobz, &uplo, &n, L, &n, eig, work, &lwork, &info);
  delete[] work;
}

// return maximum root of quadratic formula
static double quad_formula(double quad_term, double lin_term,
                           double ind_term) {
  double root1, root2;
  roots(quad_term, lin_term, ind_term, root1, root2);
  /**
   ** This considers the parallel only case where
   ** we need the max root
   **/
  return (std::max(root1, root2));
}

// compute roots of quadratic formula
static void roots(double a, double b, double c,
                  double & root1, double & root2) {
  //Check if it is a linear function
  //Returns the same value in the twoo rots when linear
  if(fabs(a) < EPS) {
    root1 = root2 = - c / b;
  }
  double discr = b*b - 4*a*c;
  //Check if roots are comples
  //Stops when complex roots
  if(discr < -EPS) {
    root1 = root2 = 1.0;
  }
  else {
    double numerator = 0.0;
    if (b > EPS)
      numerator = 0 - b - sqrt(discr);
    else
      numerator = 0 - b + sqrt(discr);
    /**
     ** If numerator is negative one of the roots
     ** will be infinity, this makes no sense in our case
     ** where the too roots are finite.
     ** TODO: Verify that this is not going to generate
     ** false results
     **/
    if (fabs(numerator) < EPS) {
      root1 = root2 = 0;
    }
    else {
      root1 = numerator / (2*a);
      root2 = (2*c) / numerator;
    }
  }
}


// print matrix, row major 0, col major 1
void CglConicGD1Cut::print_matrix(int major, int num_rows, int num_cols,
                                  double const * matrix,
                                  char const * name) const {
  std::cout << "==================== "
            << name
            << " ===================="
            << std::endl;
  for (int i=0; i<num_rows; ++i) {
    for (int j=0; j<num_cols; ++j) {
      if (major) {
        std::cout << matrix[i+j*num_rows] << " ";
      }
      else {
        std::cout << matrix[i*num_cols+j] << " ";
      }
    }
    std::cout << std::endl;
  }
}

// print vector
void CglConicGD1Cut::print_vector(int n, double const * vector,
                                  char const * name) const {
  std::cout << "==================== "
            << name
            << " ===================="
            << std::endl;
  for (int i=0; i<n; ++i) {
    std::cout << vector[i] << " ";
  }
  std::cout << std::endl;
}

void CglConicGD1Cut::print_scalar(double value, char const * name) const {
  std::cout << "==================== "
            << name
            << " ===================="
            << std::endl
            << value
            << std::endl;
}

void CglConicGD1Cut::print_cut() const {
  if (cutA_ && cutb_) {
    print_matrix(1, cutA_num_rows_, cutA_num_cols_, cutA_, "cutA");
    print_vector(cutA_num_rows_, cutb_, "cutb");
  }
}
