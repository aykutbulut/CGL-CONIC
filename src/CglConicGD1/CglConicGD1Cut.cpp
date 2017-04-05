#include "CglConicGD1Cut.hpp"

extern "C" {
  #include <cblas.h>
  // lapack routines
  void dgels_(char *trans, int *m, int *n, int *nrhs,
              double *A, int *lda, double *b, int *ldb,
              double *work, int *lwork, int *info);
  void dsysv_(char *uplo, int *n, int *nrhs, double *a,
              int *lda, int * ipiv, double *b, int * ldb,
              double * work, int *lwork, int * info);
  void dposv_(char *uplo, int *n, int *nrhs, double *a,
              int *lda, double *b, int *ldb, int *info);
  void dpotrs_(char *uplo, int *n, int *nrhs, double *a,
               int *lda, double *b, int *ldb, int *info);
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

CglConicGD1Cut::CglConicGD1Cut(OsiConicSolverInterface const * solver,
                               int num_rows, int * rows,
                               int cut_cone, int dis_var)
: solver_(solver) {
  num_rows_ = num_rows;
  dis_var_ = dis_var;
  rows_ = new int[num_rows];
  std::copy(rows, rows+num_rows, rows_);
  cindex_ = cut_cone;
  // set arrays to 0
  matA_ = 0;
  matH_ = 0;
  matV_ = 0;
  matQ_ = 0;
  vecq_ = 0;
  eigQ_ = 0;
  vecx0_ = 0;
  cmembers_ = 0;
  cone_members_ = 0;
  Jtilde_ = 0;
  rho_tilde_ = 0;
  new_matA_ = 0;
  new_rhs_ = 0;
  cut_type_ = 1;
  dirTestU_ = 0;
  dirTestV_ = 0;
  vecq_tau_ = 0;
  rho_tau_ = 0;
  matQ_tau_ = 0;
  linear_cut_ind_ = 0;
  linear_cut_coef_ = 0;
  // fill csize_, cmembers_, cone_members_, ctype_, lctype_,
  solver_->getConicConstraint(cut_cone, lctype_, csize_, cmembers_);
  ctype_ = OSI_LORENTZ;
  int num_cols = solver->getNumCols();
  cone_members_ = new int[num_cols]();
  for (int i=0; i<csize_; ++i) {
    cone_members_[cmembers_[i]] = 1;
  }
  // todo(aykut) check if csize_ is same as num_rows_
  // I do not know what to do in this case
  if (num_rows_>=csize_) {
    valid_ = false;
    return;
  }
  infeasible_ = false;
  // create disjunction from disj_var_
  compute_disjunction();
  generate_cut();
}

// create disjunction from disj_var_
void CglConicGD1Cut::compute_disjunction() {
  // dis_var should be a cone member
  double * a = new double[csize_]();
  int new_dis_index = 0;
  for (int i=0; i<csize_; ++i) {
    if (cmembers_[i]!=dis_var_)
      new_dis_index++;
    else
      break;
  }
  a[new_dis_index] = 1.0;
  rel_dis_var_ = new_dis_index;
  // define alpha and beta
  double alpha = floor(solver_->getColSolution()[dis_var_]);
  double beta = alpha+1.0;
  // create coefficient matrix
  disjunction_ = new Disjunction(csize_, a, alpha, a, beta);
  delete[] a;
}

void CglConicGD1Cut::generate_cut() {
  valid_ = true;
  compute_matrixA();
  compute_matrixH();
  compute_matrixQ();
  compute_vectorq();
  compute_rho();
  classify_quadric();
  // compute tau value that yields a cone
  compute_tau();
  if (valid_==false) {
    return;
  }
  if (cut_type_==2 || cut_type_==3) {
    // only one disjunction hyperplane intersects with quadric
    // which makes the following computations of this function
    // irrelevant
    return;
  }
  if (valid_ and infeasible_) {
    // problem is infeasible
    return;
  }
  // compute rho(tau), rho_tau_
  compute_rho_tau();
  // compute q(tau), q_tau_
  compute_q_tau();
  // compute Q(tau), Q_tau_
  compute_Q_tau();
  // compute matrix as a part of the cut to add to model
  compute_new_A();
  // compute right-hand-side that corresponds to NewA
  //compute_new_rhs();
}

void CglConicGD1Cut::compute_vectorq() {
  vecx0_ = new double[csize_];
  double const * x0 = solver_->getColSolution();
  for (int i=0; i<csize_; ++i) {
    vecx0_[i] = x0[cmembers_[i]];
  }
  vecx0_[0] = -1.0*vecx0_[0];
  int m = csize_;
  int n = m-num_rows_;
  vecq_ = new double[n]();
  cblas_dgemv (CblasColMajor, CblasTrans, m, n, 1.0, matH_, m, vecx0_, 1,
               0.0, vecq_, 1);
  vecx0_[0] = -1.0*vecx0_[0];
}

void CglConicGD1Cut::compute_rho() {
  int m = csize_;
  rho_ = - (vecx0_[0]*vecx0_[0]);
  rho_ += cblas_ddot(m-1, vecx0_+1, 1, vecx0_+1, 1);
}

void CglConicGD1Cut::classify_quadric() {
  //compute q^T Q^-1 q - \rho
  int m = csize_;
  int n = m-num_rows_;
  double * Qq = new double[n];
  solveSM(n, matQ_, vecq_, Qq);
  // quadricRHS is q^T Q^-1 q - \rho
  double qadricRHS = cblas_ddot(n, vecq_, 1, Qq, 1) - rho_;
  //Copy matrix Q
  matV_ = new double[n*n];
  cblas_dcopy(n*n, matQ_, 1, matV_, 1);
  //Vector with the eigenvalues of Q
  eigQ_ = new double[n]();
  //Compute the eigenvalue decomposition
  eigDecompICL(n, matV_, eigQ_);
  // in case nonpositive eigenvalues bail out
  for (int i=0; i<n; ++i) {
    if (eigQ_[i]<1e-3) {
      std::cout << "Q is not positive definite!" << std::endl;
      valid_ = false;
      break;
    }
  }
  delete [] Qq;
}

void CglConicGD1Cut::compute_matrixA() {
  // we assume cmembers_ and rows_ are ordered.
  // matA_ has num_rows_ many rows and csize_ many columns.
  // matA_ is column ordered
  matA_ = new double[num_rows_*csize_]();
  // rows_ stores the indices of the picked rows of constraint matrix
  // cmembers_ stores the indices of the picked rows of constraint matrix
  int nc = solver_->getNumCols();
  int nr = solver_->getNumRows();
  int * row_members = new int[nr]();
  for (int i=0; i<num_rows_; i++) {
    row_members[rows_[i]] = 1;
  }
  // cone_members_ is filled in constructor
  CoinPackedMatrix const * mat = solver_->getMatrixByCol();
  int const * indices = mat->getIndices();
  double const * elements = mat->getElements();
  // reduced row index
  int rri = 0;
  // reduced column index
  int rci = 0;
  for (int i=0; i<nc; ++i) {
    if (cone_members_[i]==1) {
      for (int j=0; j<nr; ++j) {
        if (row_members[j]==1) {
          matA_[rci*num_rows_+rri] = mat->getCoefficient(j,i);
          rri++;
        }
      }
      rri = 0;
      rci++;
    }
  }
  delete[] row_members;
}

// compute matrix H from matrix A, H is null mat of A.
void CglConicGD1Cut::compute_matrixH() {
  int num_cols = csize_;
  //Copy A to a working array
  double * tempA = new double[num_cols*num_rows_];
  cblas_dcopy(num_cols*num_rows_, matA_, 1, tempA, 1);
  //Right hand side singular vectors of A
  double * VT = new double[num_cols*num_cols];
  svDecompICL(num_rows_, num_cols, tempA, VT);
  matH_ = new double[(num_cols-num_rows_)*num_cols]();
  // Take the last n-m columns of V, lapack returns V^T
  for(int i=0; i<(num_cols-num_rows_); ++i) {
    //cblas_dcopy(num_cols, (VT + num_rows_ + i), num_cols, (matH_+i*num_cols), 1);
    cblas_dcopy(num_cols, (VT+num_rows_+i), num_cols, (matH_+i*num_cols), 1);
  }
  delete[] tempA;
  delete[] VT;

  // Set H hard coded for debugging
  // matH_[1]= 0.3739;
  // matH_[2]= -0.5764;
  // matH_[3]= -0.6203;
  // matH_[4]= 0.3094;
  // matH_[5]= 0.2178;
  // matH_[7]= -0.5067;
  // matH_[8]= -0.6141;
  // matH_[9]= 0.4788;
  // matH_[10]= 0.2179;
  // matH_[11]= 0.2991;

  // change matH_ to col order
  // double * tempH new double[(num_cols-num_rows_)*num_cols];
  // cblas_dcopy((num_cols-num_rows_)*num_cols, matH_, 1, tempH, 1);
  // for(int i=0; i<(num_cols-num_rows_); ++i)
  //   cblas_dcopy(num_cols, (VT + num_rows_ + i), num_cols, (matH_+i*num_cols), 1);


}

void CglConicGD1Cut::compute_matrixQ() {
  // mH_ = m; Number of rows of H, num_cols-num_rows_
  // nH_ = n; Number of columns of H, num_cols
  //int num_cols = solver_->getNumCols();
  // negate first row of H. Multiply H^T with this matrix.
  // matH_ is stored as row major in the memory.
  // it should be col ordered
  int m = csize_;
  int n = m-num_rows_;
  matQ_ = new double[n*n]();
  // Temporary array, stores the first row of H
  double * d = new double[n];
  cblas_dcopy(n, matH_, m, d, 1);
  // Temporary array, exlcudes the first row of H
  double * A = new double[(m-1)*n];
  for(int i=0; i<n; i++)
    cblas_dcopy(m-1, matH_+i*m+1, 1, A+i*(m-1), 1);
  //computes Q = A^TA - dd^T =  H^T J H
  cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, n, m-1,
              1.0, A, m-1, 0.0, matQ_, n);
  cblas_dsyr(CblasColMajor, CblasUpper, n, -1.0, d, 1, matQ_, n);
  delete[] d;
  delete[] A;
}

void CglConicGD1Cut::compute_tau() {
  int m = csize_;
  int n = csize_-num_rows_;

  print_vector(csize_, vecx0_, "x0");

  // === Compute disjunction in u-space ===
  // c is the disjunction coefficient, we compute a from c in u-space

  print_matrix(1, m-n, m, matA_, "A");
  print_matrix(1, m, n, matH_, "H");
  print_matrix(1, n, n, matQ_, "Q");
  print_matrix(1, n, n, matV_, "V");
  print_vector(n, eigQ_, "eigQ");
  print_vector(n, vecq_, "q");
  print_scalar(rho_, "rho");

  a_ = new double[n]();
  double const * c = disjunction_->get_c1();
  double norm_of_q = cblas_dnrm2(n, vecq_, 1);
  cblas_dgemv(CblasColMajor, CblasTrans, m, n, norm_of_q,
              matH_, m, c, 1, 0.0, a_, 1);
  // normalize a
  double norm_of_a = cblas_dnrm2(n, a_, 1);
  cblas_dscal(n, 1.0/norm_of_a, a_, 1);
  print_vector(n, a_, "normalized a");

  // compute alpha_ and beta_ in regularized space
  {
    double cTx0 = vecx0_[rel_dis_var_];
    double * Hq = new double[m]();
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, 1.0,
                matH_, m, vecq_, 1, 0.0, Hq, 1);
    print_vector(m, Hq, "Hq");
    double cTHq = Hq[rel_dis_var_];
    alpha_ = floor(vecx0_[rel_dis_var_]) - cTx0 + cTHq;
    alpha_ = alpha_/norm_of_a;
    beta_ = ceil(vecx0_[rel_dis_var_]) - cTx0 + cTHq;
    beta_ = beta_/norm_of_a;
    delete[] Hq;
    print_scalar(alpha_, "alpha normalized");
    print_scalar(beta_, "beta normalized");
  }

  // == check whether the disjunctive hyperplane (transfered into the
  // regularized space) intersects with the quadric in the regularized space
  {
    bool alpha_intersects = (alpha_*alpha_ <= 1);
    bool beta_intersects = (beta_*beta_ <= 1);
    if (alpha_intersects && beta_intersects) {
      // both hyperplanes intersects regularized quadric
      // cool, keep working on it.
    }
    else if (alpha_intersects) {
      // only ax=alpha intersects quadric
      // hand linear cut ax <= alpha to solver.
      valid_ = true;
      cut_type_ = 3;
      linear_cut_ind_ = new int[1];
      linear_cut_ind_[0] = dis_var_;
      linear_cut_coef_ = new double[1];
      linear_cut_coef_[0] = 1.0;
      linear_cut_rhs_ = ceil(disjunction_->get_c10());
      return;
    }
    else if(beta_intersects) {
      // only ax=beta intersects quadric
      // hand linear cut ax >= beta to solver.

      valid_ = true;
      cut_type_ = 2;
      linear_cut_ind_ = new int[1];
      linear_cut_ind_[0] = dis_var_;
      linear_cut_coef_ = new double[1];
      linear_cut_coef_[0] = 1.0;
      linear_cut_rhs_ = floor(disjunction_->get_c20());
      return;
    }
    else {
      // none of the hyperplanes intersect the quadric, problem is infeasible.
      valid_ = true;
      infeasible_ = true;
    }
  }
  // compute coefficients of polynomial
  double quad_coef = (alpha_-beta_)*(alpha_-beta_);
  double norm_q_square = cblas_ddot(m, vecq_, 1, vecq_, 1);
  double aTq = cblas_ddot(m, a_, 1, vecq_, 1);
  double lin_coef = 4.0*norm_q_square -
    (alpha_+beta_+2*aTq)*(alpha_+beta_+2*aTq) +
    (alpha_-beta_)*(alpha_-beta_);
  double const_term = 4.0*norm_q_square;
  if (lin_coef*lin_coef < 4.0*quad_coef*const_term) {
    std::cerr << "Imaginary root!" << std::endl;
    valid_ = false;
  }
  else {
    // quadformula returns 1 when the roots are imaginary which will
    // not happen if quadric intersects with disjunction hyperplanes
    tau_ = quad_formula(quad_coef, lin_coef, const_term);
  }
  print_scalar(tau_, "tau");
  if (tau_ > -1.3) {
    valid_ = false;
  }
}

// compute rho(tau), rho_tau_
void CglConicGD1Cut::compute_rho_tau() {
  rho_tau_ += tau_*alpha_*beta_;
  print_scalar(rho_tau_, "rho_tau");
}

 // compute q(tau), q_tau_
void CglConicGD1Cut::compute_q_tau() {
  int n = csize_-num_rows_;
  vecq_tau_ = new double[n]();
  std::copy(vecq_, vecq_+n, vecq_tau_);
  double aux = -0.5*tau_*(alpha_ + beta_);
  cblas_daxpy(n, aux, a_, 1, vecq_tau_, 1);
  print_vector(n, vecq_tau_, "q(tau)");
}

 // compute Q(tau), Q_tau_
void CglConicGD1Cut::compute_Q_tau() {
  int m = csize_;
  int n = m-num_rows_;
  matQ_tau_ = new double[n*n];
  cblas_dcopy(n*n, matQ_, 1, matQ_tau_, 1);
  cblas_dsyr(CblasColMajor, CblasUpper, n, tau_, a_, 1,
             matQ_tau_, n);
  // copy upper triangular to lower triangular part
  // for each column
  for (int i=0; i<n; ++i) {
    // for each row
    for (int j=0; j<i; ++j) {
      matQ_tau_[j*n+i] = matQ_tau_[i*n+j];
    }
  }
  print_matrix(1, n, n, matQ_tau_, "Q(tau)");
}

// compute matrix as a part of the cut to add to model
void CglConicGD1Cut::compute_new_A() {
  int m = csize_;
  int n = csize_-num_rows_;
  // compute eigenvalue decomposition of Qtau_
  double * Vtau = new double[n*n];
  cblas_dcopy(n*n, matQ_tau_, 1, Vtau, 1);
  double * eigQtau = new double[n]();
  eigDecompICL(n, Vtau, eigQtau);
  print_matrix(1, n, n, Vtau, "V(tau)");
  print_vector(n, eigQtau, "eig Q(tau)");
  // check whether any of the eigenvalues are zero.
  for (int i=0; i<n; ++i) {
    if (eigQtau[i]<0.001 && eigQtau[i]>-0.001) {
      delete[] Vtau;
      delete[] eigQtau;
      valid_ = false;
      return;
    }
  }

  // at most 1 eigenvalue must be negative
  {
    int num_neg_eigen = 0;
    for (int i=0; i<n; ++i) {
      if (eigQtau[i]<0.0) {
        num_neg_eigen++;
      }
    }
    if (num_neg_eigen>1) {
      std::cerr << "Number of negative eigenvalues should be at most 1!"
                << std::endl;
      delete[] Vtau;
      delete[] eigQtau;
      valid_ = false;
      return;
    }
  }

  // compute Dtilde^{1/2}
  double * sqrtDtau = new double[n*n]();
  for (int i=0; i<n; ++i) {
    sqrtDtau[i*n+i] = sqrt(fabs(eigQtau[i]));
  }
  print_matrix(1, n, n, sqrtDtau, "sqrtDtau");

  // compute W
  double * W = new double[n*n];
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n,
              n, 1.0, Vtau, n, sqrtDtau, n, 0.0, W, n);
  print_matrix(1, n, n, W, "W");

  // compute wbar
  double * wbar = new double[n]();
  solveSM(n, matQ_tau_, vecq_tau_, wbar);
  for (int i=0; i<n; ++i) {
    wbar[i] = -wbar[i];
  }
  print_vector(n, wbar, "wbar");


  // TODO(AYKUT) CAN WE REMOVE THIS SCALING
  double norm_q = cblas_dnrm2(n, vecq_, 1);
  // print_scalar(norm_q, "||q||");
  if (norm_q<0.01) {
    std::cerr << "q is very close to 0, cutting failed."
              << std::endl;
    delete[] Vtau;
    delete[] eigQtau;
    delete[] sqrtDtau;
    delete[] W;
    delete[] wbar;
    valid_ = false;
    return;
  }
  else {
    // Anew is 1/sqrt(q'*q) * W*H'
    new_matA_ = new double[n*m]();
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, m, n,
                1.0/norm_q, W, n, matH_, m, 0.0, new_matA_, n);
    print_matrix(1, n, m, new_matA_, "Anew");

    new_rhs_ = new double[n]();
    // bnew is coef*W*(H'*x0 - q) + W * wbar
    // pick rel_dis_var_ th row of H.
    cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, matH_, m,
                vecx0_, 1,
                0.0, new_rhs_, 1);
    // subtract q
    cblas_daxpy(n, -1.0, vecq_, 1, new_rhs_, 1);
    // W(H'x0-q)
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, W, n,
                new_rhs_, 1,
                0.0, new_rhs_, 1);
    //print_vector(n, new_rhs_, "W*(H'*x0-q)");
    // coef*W(H'x0-q)
    cblas_dscal(n, 1.0/norm_q, new_rhs_, 1);
    // add Wwbar
    double * Wwbar = new double[n];
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, W, n,
                wbar, 1,
                0.0, Wwbar, 1);
    cblas_daxpy(n, 1.0, Wwbar, 1, new_rhs_, 1);
    print_vector(n, new_rhs_, "bnew");

    // negate first row of A and rhs[0] if -W(:,1)'*wbar is less than 0.0
    if (Wwbar[0]>0) {
      // negate first row of A
      for (int i=0; i<m; ++i) {
        new_matA_[i*n] = -new_matA_[i*n];
      }
      // negerate rhs
      new_rhs_[0] = -new_rhs_[0];
    }
    delete[] Wwbar;
  }

  // TODO(AYKUT) EIGENVALUE ORDERING


  delete[] Vtau;
  delete[] eigQtau;
  delete[] sqrtDtau;
  delete[] W;
  delete[] wbar;


  //Initilize the new A matrix section
  // int m = csize_;
  // int n = m-num_rows_;
  // new_matA_ = new double[n*m]();
  // char jobz = 'V';
  // char uplo = 'U';
  // int info = 0;
  // int lwork = -1;
  // double worksize[1];
  // double * L = new double[n*n];
  // cblas_dcopy(n*n, matQ_tau_, 1, L, 1);
  // double * d = new double[n];
  // // eigenvalue decomposition, L holds normalized eigenvectors,
  // // d holds eigenvalues. Q = L^TDL
  // dsyev_(&jobz, &uplo, &n, L, &n, d, worksize, &lwork, &info);
  // lwork = (int) worksize[0];
  // double * work = new double[lwork];
  // dsyev_(&jobz, &uplo, &n, L, &n, d, work, &lwork, &info);
  // if (cut_type_ > 0){
  //   double * LDinv = new double[n*n]();
  //   // todo(aykut) why do we need this? we already know the leading
  //   // variable of the conic constraint.
  //   varHead_ = -1;
  //   // the following loop is for computing LD^-1
  //   for(int i = 0; i < n; i++){
  //     //printf("d[%d] = %f\n", i, d[i]);
  //     std::cout << "d[" << i << "] " << d[i] << std::endl;
  //     if(fabs(d[i]) < 1e-2){
  //       valid_ = false;
  //     }
  //     if (d[i] < 0.0){
  //       //printf("head var = %d\n", i);
  //       varHead_ = i;
  //       cblas_daxpy(n, 1.0/sqrt(fabs(d[i])), L+n*i, 1, LDinv+n*i, 1);
  //       dirTestU_ = new double[n]();
  //       cblas_dcopy(n, L+n*i, 1, dirTestU_, 1);
  //     }
  //     else {
  //       cblas_daxpy(n, 1/sqrt(d[i]), L+n*i, 1, LDinv + n*i, 1);
  //     }
  //   }
  //   cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, csize_, n,
  //               n, -1.0, matH_, csize_, LDinv, n, 0.0, new_matA_, csize_);
  // }
  // else {
  //   varHead_ = -1;
  //   for(int i=0; i<n; i++) {
  //     if(fabs(d[i]) < 1e-2){
  //       valid_ = false;
  //     }
  //     if (d[i] < -EPS) {
  //       varHead_ = i;
  //       cblas_daxpy(n, sqrt(fabs(d[i])), L+n*i, 1, new_matA_+i, n);
  //       dirTestU_ = new double[n]();
  //       cblas_dcopy(n, L+n*i, 1, dirTestU_, 1);
  //     }
  //     else {
  //       cblas_daxpy(n, sqrt(d[i]), L+n*i, 1, new_matA_+i, n);
  //     }
  //   }
  // }
  // delete[] work;
  // delete[] L;
  // delete[] d;
}

// compute right-hand-side that corresponds to NewA
void CglConicGD1Cut::compute_new_rhs() {
  int m = csize_;
  int n = m-num_rows_;
  int * ipiv = new int[n]();
  double * aux = new double[n];
  cblas_dcopy(n, vecq_tau_, 1, aux, 1);
  new_rhs_ = new double[m]();
  char uplo = 'U';
  int nrhs = 1;
  double * L = new double[n*n];
  double worksize;
  int lwork = -1;
  int info;
  cblas_dcopy(n*n, matQ_tau_, 1, L, 1);
  dsysv_(&uplo, &n, &nrhs, L, &n, ipiv, aux, &n,
         &worksize, &lwork, &info);
  lwork = (int) worksize;
  double * work = new double[lwork];
  dsysv_(&uplo, &n, &nrhs, L, &n, ipiv, aux, &n,
         work, &lwork, &info);
  if (cut_type_>0) {
    cblas_daxpy(n, 1.0, aux, 1, dirTestV_, 1);
    double dirtest;
    if (dirTestU_) {
      dirtest = cblas_ddot(n, dirTestV_, 1, dirTestU_, 1);
    }
    else {
      dirtest = cblas_ddot(n, dirTestV_, 1, dirTestV_, 1);
    }
    if(dirtest < -EPS){
      for(int i=0; i<m; i++){
        new_matA_[varHead_*m+i]*= -1.0;
      }
    }
    cblas_dgemv(CblasColMajor, CblasNoTrans, m, n, -1.0, matH_,
                m, aux, 1, 0.0, new_rhs_, 1);
    cblas_daxpy(m, 1.0, vecx0_, 1, new_rhs_, 1);
  }
  else {
    //todo(aykut) understand the direction test. What are we testing?
    // why are we testing, what are dirTestV_ and dirTestU_?
    //Matrix times vector Q^(-1)q
    if (dirTestV_ != 0) {
      delete[] dirTestV_;
    }
    dirTestV_ = new double [n]();
    cblas_dcopy (n, matV_, 1, dirTestV_, 1);
    double dirtest = cblas_ddot(n, dirTestV_, 1, dirTestU_, 1);
    if(dirtest < -EPS){
      for(int i=0; i<n; i++) {
        new_matA_[varHead_ + m*i]*= -1.0;
      }
    }
    //Computes V^T a
    cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, -1.0, new_matA_,
                n, aux, 1, 0.0, new_rhs_, 1);
  }
  delete[] work;
  delete[] ipiv;
  delete[] aux;
  delete[] L;
}


// // this one should go to CglConicGD1 class.
// void add_cut() {
//   //std::cout<<"numVars_ = "<<numVars_ <<std::endl;
//   if(valid_) {
//     IclopsMosek * solver = model_->getSolver();
//     int * vars = ccgenerator_->getVarsSub();
//     solver->addFreeVars(elliDim_);
//     std::cout<<"\n\n\n\n\n\n\n\n\n\nnumVars_ = "<<numVars_ <<std::endl;
//     std::cout<<"Number of rows before = "<<solver->getNumRows()<<std::endl;
//     solver->addRows(numVars_);
//     /* Get the new number of variables */
//     int numVarP = solver->getNumCols();
//     int numConP = solver->getNumRows();
//     std::cout<<"Number of rows after = "<<solver->getNumRows()<<std::endl;
//     //Array to store the variables indices for adding a row
//     int newRowsub[elliDim_+1];
//     //Array to store the values of the coefficients for adding a row
//     double newRowval[elliDim_+1];
//     if(ccgenerator_->getType() > 0){
//       /* Register the indices of the new variables */
//       for(int i = 0; i < elliDim_; i ++) {
//      newRowsub[i+1] = numVarP - elliDim_ + i;
//       }
//       /* In the original variables part we add and identity matrix */
//       newRowval[0] = 1.0;
//       /* Pointer to acces the values
//       of the entries for the new constraints */
//       double * indexA;
//       for(int i = 0; i < numVars_; i++){
//      /* Get the indices for the new variables in the problem */
//      newRowsub[0] = vars[i];
//      /* Get the values of the entries of the new constraints */
//      for(int j = 0; j < elliDim_; j++) {
//        indexA = newAslice_ + i + j * numVars_;
//        newRowval[j+1] = (fabs(*(indexA)) < ICLEPS)?0.0:*(indexA);
//      }
//      /* Add the constraints bounds, i.e the right hand side */
//      solver->setRowBounds(numConP - numVars_ + i,
//                           *(newbslice_+i),
//                           *(newbslice_+i));
//      /* Add the new row to the constraint matrix */
//      solver->setRow(numConP - numVars_ + i, elliDim_+1,
//                     newRowsub, newRowval);
//       }
//     }
//     else {
//       /* Register the indices of the new variables */
//       for(int i = 0; i < elliDim_; i ++) {
//      newRowsub[i] = vars[i];
//       }
//       /* In the new variables part we add and identity matrix */
//       newRowval[elliDim_] = 1.0;
//       double * indexA;
//       for(int i = 0; i < numVars_; i++){
//      /* Get the indices for the new variables in the problem */
//      newRowsub[elliDim_] = numVarP - elliDim_ + i;
//      /* Get the values of the entries of the new constraints */
//      for(int j = 0; j < elliDim_; j++) {
//        indexA = newAslice_ + i + j * numVars_;
//        newRowval[j] = (fabs(*(indexA)) < ICLEPS)?0.0:*(indexA);
//      }
//      /* Add the constraints bounds, i.e the right hand side */
//      solver->setRowBounds(numConP - numVars_ + i,
//                           *(newbslice_+i),
//                           *(newbslice_+i));
//      /* Add the new row to the constraint matrix */
//      solver->setRow(numConP - numVars_ + i, elliDim_+1,
//                     newRowsub, newRowval);
//       }
//     }
//     /* Append the new cone
//     ** VarHead_: This is the index of the leading variable in the
//     ** quadratic cone
//     **/
//     int csub[elliDim_];
//     csub[0] = varHead_ + numVarP - elliDim_;
//     //std::cout<<"csub[0] = "<<csub[0]<<"\n";
//     int l = 1;
//     for(int i = 0; i < elliDim_; i ++) {
//       if ( i != varHead_) {
//      csub[l] = numVarP - elliDim_ + i;
//      //std::cout<<"csub["<<l<<"] = "<<csub[l]<<"\n";
//      l++;
//       }
//     }
//     solver->addCone(elliDim_, csub);
//     //printf("A conic cut added \n");
//   }
// }

bool CglConicGD1Cut::valid() const {
  return valid_;
}

bool CglConicGD1Cut::infeasible() const {
  return infeasible_;
}

// number of rows of the linear part of the cut
// ie num of rows of the new A matrix.
int CglConicGD1Cut::getNumRows() const {
  return csize_ - num_rows_;
}

// number of columns in the cut
int CglConicGD1Cut::getNumCols() const {
  return csize_;
}

// get variables in the cut
int * CglConicGD1Cut::getMembers() const {
  return cmembers_;
}

int CglConicGD1Cut::cutType() const {
  return cut_type_;
}

// todo(aykut) this is irrelevant since the var head is always the first
// cone member
int CglConicGD1Cut::getVarHead() const {
  return varHead_;
}

CglConicGD1Cut::~CglConicGD1Cut() {
  if (rows_) {
    delete[] rows_;
  }
  if (matA_) {
    delete[] matA_;
  }
  if (matH_) {
    delete[] matH_;
  }
  if (matV_) {
    delete[] matV_;
  }
  if (matQ_) {
    delete[] matQ_;
  }
  if (vecq_) {
    delete[] vecq_;
  }
  if (eigQ_) {
    delete[] eigQ_;
  }
  if (cmembers_)
    delete[] cmembers_;
  if (cone_members_)
    delete[] cone_members_;
  if (Jtilde_)
    delete[] Jtilde_;
  if (vecx0_)
    delete[] vecx0_;
  if (matQ_tau_)
    delete[] matQ_tau_;
  if (vecq_tau_)
    delete[] vecq_tau_;
  if (new_matA_)
    delete[] new_matA_;
  if (new_rhs_)
    delete[] new_rhs_;
  if (dirTestU_)
    delete[] dirTestU_;
  if (dirTestV_)
    delete[] dirTestV_;
  if (linear_cut_ind_)
    delete[] linear_cut_ind_;
  if (linear_cut_coef_)
    delete[] linear_cut_coef_;
}

// return linear part of the cut, constraint matrix.
double const * CglConicGD1Cut::getNewMatA() const {
  return new_matA_;
}

// return right hand side of the linear part of the cut.
double const * CglConicGD1Cut::getNewRhs() const {
  return new_rhs_;
}

// get the variables in the cut if cut is in linear form
int const * CglConicGD1Cut::linear_cut_ind() const {
  if (cut_type_!=2 && cut_type_!=3) {
    std::cerr << "Cut is not in linear form!" << std::endl;
    throw std::exception();
  }
  return linear_cut_ind_;
}

// get the coefficients of the cut if cut is in linear form
double const * CglConicGD1Cut::linear_cut_coef() const {
  if (cut_type_!=2 && cut_type_!=3) {
    std::cerr << "Cut is not in linear form!" << std::endl;
    throw std::exception();
  }
  return linear_cut_coef_;
}

// get the rhs of the cut if the cut is in linear form
double CglConicGD1Cut::linear_cut_rhs() const {
  if (cut_type_!=2 && cut_type_!=3) {
    std::cerr << "Cut is not in linear form!" << std::endl;
    throw std::exception();
  }
  return linear_cut_rhs_;
}

// get the size of the cut if the cut is in the linear form
int CglConicGD1Cut::linear_cut_size() const {
  if (cut_type_!=2 && cut_type_!=3) {
    std::cerr << "Cut is not in linear form!" << std::endl;
    throw std::exception();
  }
  return 1;
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
  if (new_matA_==NULL or new_rhs_==NULL) {
    return;
  }
  print_matrix(1, csize_-num_rows_, csize_, new_matA_,
               "cut coefficient matrix");
  print_vector(csize_-num_rows_, new_rhs_, "cut rhs");
}
