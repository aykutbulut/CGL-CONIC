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

#define EPS 1e-10

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
  rows_ = 0;
  matA_ = 0;
  matH_ = 0;
  matV_ = 0;
  matQ_ = 0;
  vecq_ = 0;
  eigQ_ = 0;
  vecx0_ = 0;
  vecx0J_ = 0;
  cmembers_ = 0;
  cone_members_ = 0;
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
  // cone_members_ = new int[csize_]();
  // for (int i=0; i<csize_; ++i) {
  //   cone_members_[cmembers_[i]] = 1;
  // }
  // todo(aykut) check if csize_ is same as num_rows_
  // I do not know what to do in this case
  if (num_rows_>=csize_) {
    valid_ = false;
    return;
  }
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
  // compute rho(tau), rho_tau_
  compute_rho_tau();
  // compute q(tau), q_tau_
  compute_q_tau();
  // compute Q(tau), Q_tau_
  compute_Q_tau();
  // compute matrix as a part of the cut to add to model
  compute_new_A();
  // compute right-hand-side that corresponds to NewA
  compute_new_rhs();
}

void CglConicGD1Cut::compute_vectorq() {
  vecx0_ = new double[csize_];
  double const * x0 = solver_->getColSolution();
  for (int i=0; i<csize_; ++i) {
    vecx0_[i] = x0[cmembers_[i]];
  }
  vecx0J_ = new double[csize_];
  std::copy(vecx0_, vecx0_+csize_, vecx0J_);
  vecx0J_[0] = -1.0*vecx0J_[0];
  int m = csize_;
  int n = m-num_rows_;
  vecq_ = new double[n]();
  cblas_dgemv (CblasColMajor, CblasTrans, m, n, 1.0, matH_, m, vecx0J_, 1,
	       0.0, vecq_, 1);
}

void CglConicGD1Cut::compute_rho() {
  int m = csize_;
  rho_ = - (vecx0J_[0]*vecx0J_[0]);
  rho_ += cblas_ddot(m-1, vecx0J_+1, 1, vecx0J_+1, 1);
}

void CglConicGD1Cut::classify_quadric() {
  //compute q^T Q^-1 q - \rho
  int m = csize_;
  int n = m-num_rows_;
  double * Qq = new double[n];
  solveSM(n, matQ_, vecq_, Qq);
  double qadricRHS = cblas_ddot(n, vecq_, 1, Qq, 1) - rho_;
  //Copy matrix Q
  matV_ = new double[n*n];
  cblas_dcopy(n*n, matQ_, 1, matV_, 1);
  //Vector with the eigenvalues of Q
  eigQ_ = new double[n]();
  //Compute the eigenvalue decomposition
  eigDecompICL(n, matV_, eigQ_);
  quad_type_ = 'e';
  if( eigQ_[0] < -EPS ){
    if( eigQ_[1] < EPS){
      quad_type_ = 'u';
      active_ = false;
    }
    else if (fabs(qadricRHS) < EPS)
      quad_type_ = 'k';
    else
      quad_type_ = 'h';
  }
  else if (fabs(eigQ_[0]) < EPS) {
    numZEig_ = 0;
    for (int i=0; i<n; ++i) {
      if (fabs(eigQ_[i]) < EPS)
	numZEig_++;
    }
    quad_type_ = 'p';
    active_ = false;
    if (eigQ_[1] < EPS) {
      quad_type_ = 'u';
    }
  }
  delete [] Qq;
}

void CglConicGD1Cut::compute_matrixA() {
  // todo(aykut) num_cols is same as cone size. hence cone_members_ becomes irrelevant
  int num_cols = solver_->getNumCols();
  // todo(aykut) we may not need this
  cone_members_ = new int[num_cols]();
  for (int i=0; i<csize_; ++i) {
    cone_members_[cmembers_[i]] = 1;
  }
  matA_ = new double[num_rows_*csize_]();
  // get matrix A in col ordered form
  // cols are ordered as in original matrix
  CoinPackedMatrix const * mat = solver_->getMatrixByCol();
  int const * indices = mat->getIndices();
  double const * elements = mat->getElements();
  for (int i=0; i<num_cols; ++i) {
    // for each column i
    // if column is not a member of the cone, skip
    if (cone_members_[i]==0)
      continue;
    // some columns are all zero, ie first is same as last.
    int first = mat->getVectorFirst(i);
    int last = mat->getVectorLast(i);
    for (int j=first; j<last; ++j) {
      if (indices[j]>num_rows_)
	break;
      matA_[i*num_rows_+indices[j]] = elements[j];
    }
  }
}

// compute matrix H from matrix A, H is null mat of A.
void CglConicGD1Cut::compute_matrixH() {
  int num_cols = csize_;
  //Copy A to a working array
  double * tempA =  new double [num_cols*num_rows_];
  cblas_dcopy(num_cols*num_rows_, matA_, 1, tempA, 1);
  //Right hand side singular vectors of A
  double * VT = new double [num_cols*num_cols];
  svDecompICL(num_rows_, num_cols, tempA, VT);
  matH_ = new double[(num_cols-num_rows_)*num_cols]();
  // Take the last n-m columns of V, lapack returns V^T
  for(int i=0; i<(num_cols-num_rows_); ++i) {
    //cblas_dcopy(num_cols, (VT + num_rows_ + i), num_cols, (matH_+i*num_cols), 1);
    cblas_dcopy(num_cols, (VT+num_rows_+i), num_cols, (matH_+i*num_cols), 1);
  }
  delete [] tempA;
  delete [] VT;
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
  // Temporary array, stores the firs row of H
  //double * d = new double[n];
  double * d = new double[n];
  cblas_dcopy(n, matH_, m, d, 1);
  // Temporary array, exlcudes the first row of H
  //double * A = new double[(m-1)*n];
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
  int n = m-num_rows_;
  // initialize vecq_tau_ to vec_q
  vecq_tau_ = new double[n];
  std::copy(vecq_, vecq_+n, vecq_tau_);
  //Computes a^TQ^(-1)a
  double * chol = new double[n*n];
  cblas_dcopy(n*n, matQ_, 1, chol, 1);
  double * Qinv_a = new double[n];
  // u is the disjunction coefficient, we compute a from u in w-space
  // a = H^T u
  // todo(aykut) what if a is 0?
  double const * u = disjunction_->get_c1();
  a_ = new double[n]();
  cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, matH_, m, u, 1, 0.0, a_, 1);
  // normilize a
  double norm_of_a = cblas_dnrm2(n, a_, 1);
  if (norm_of_a<EPS) {
    // todo(aykut) does this mean that problem is infeasible?
    valid_ = false;
    delete[] chol;
    delete[] Qinv_a;
    return;
  }
  // todo(aykut) do this using blas
  for (int i=0; i<n; ++i) {
    a_[i] = a_[i]/norm_of_a;
  }
  // end of computing a
  cblas_dcopy(n, a_, 1, Qinv_a, 1);
  char U = 'U';
  int nrhs = 1;
  int info = 0;
  //Matrix times vector Q^(-1)a
  dposv_(&U, &n, &nrhs, chol, &n, Qinv_a, &n, &info);
  //Scalar aQinv_a
  double aQinv_a = cblas_ddot(n, a_, 1, Qinv_a, 1);
  //Scalar qQinv_a
  double qQinv_a = cblas_ddot(n, vecq_tau_, 1, Qinv_a, 1);
  //Matrix times vector Q^(-1)q
  dirTestV_ = new double[n]();
  double * Qinv_q = new double[n];
  cblas_dcopy(n, vecq_tau_, 1, Qinv_q, 1);
  dpotrs_(&U, &n, &nrhs, chol, &n, Qinv_q, &n, &info);
  // todo(aykut) dirTestV_ is -Qinv_q
  cblas_daxpy(n, -1.0, Qinv_q, 1, dirTestV_, 1);
  //Scalar qQ^{-1}q
  double qQinv_q = cblas_ddot(n, vecq_tau_, 1, Qinv_q, 1);
  //quadratic equation
  // todo(aykut) see we have right alpha/beta setting
  // todo(aykut) alpha and beta is wrong you should fix them so that they are in w space.
  // alpha <- alpha - u^Tx0
  // beta <- beta - u^Tx0
  double uT_x0 = cblas_ddot(m, u, 1, vecx0_, 1);
  alpha_ = disjunction_->get_c20() - uT_x0;
  beta_ = disjunction_->get_c10() - uT_x0;
  alpha_ = alpha_ / norm_of_a;
  beta_ = beta_ / norm_of_a;
  double quadTerm = (alpha_ - beta_)*(alpha_ - beta_)*aQinv_a;
  double linearTerm = 4*aQinv_a*(qQinv_q - rho_)
    - (alpha_ + beta_ + 2*qQinv_a)*(alpha_ + beta_ + 2*qQinv_a)
    + (alpha_ - beta_)*(alpha_ - beta_);
  double indTerm = 4*(qQinv_q - rho_);
  double intr1 = (beta_ + qQinv_a)*(beta_ + qQinv_a) - aQinv_a*(qQinv_q - rho_);
  double intr2 = (alpha_ + qQinv_a)*(alpha_ + qQinv_a) - aQinv_a*(qQinv_q - rho_);
  if(intr1 > -EPS && intr2 > -EPS)
    // todo(aykut) is it necessary to have it as a class member?
    // this means both hyperplanes (a,alpha) and (a,beta) does not intersect with
    // quadric, I think this may indicate infeasible problem.
    valid_ = false;
  if (intr1 > -EPS) {
    if (intr2 < -EPS) {
      // one piece of disjunction does not intersect with quadric (a,alpha)
      // only (a,beta) intersects with quadric
      cut_type_ = 2;
      linear_cut_ind_ = new int[1];
      linear_cut_ind_[0] = dis_var_;
      linear_cut_coef_ = new double[1];
      linear_cut_coef_[0] = 1.0;
      linear_cut_rhs_ = ceil(disjunction_->get_c20());
    }
  }
  if (intr2 > -EPS){
    if (intr1 < -EPS){
      // one piece of disjunction does not intersect with quadric (a,alpha)
      // only (a,beta) intersects with quadric
      cut_type_ = 3;
      linear_cut_ind_ = new int[1];
      linear_cut_ind_[0] = dis_var_;
      linear_cut_coef_ = new double[1];
      linear_cut_coef_[0] = 1.0;
      linear_cut_rhs_ = floor(disjunction_->get_c10());
    }
  }
  // quadformula returns 1 when the roots are imaginary which will
  // not happen if quadric intersects with disjunction hyperplanes
  tau_ = quad_formula(quadTerm, linearTerm, indTerm);
  // todo(aykut) what is 1e-6? another eps?
  // tau should be negative, check paper, -1 > tau_2 > tau_1
  if (cut_type_!=2 && cut_type_!=3) {
    // todo(aykut) in the essence this check looks same as
    // intr1>-EPS intr2>-EPS to me
    if(fabs(tau_ + (1.0/aQinv_a)) < 1e-6 || tau_ > EPS) {
      // this means tau is very close to zero (in case of unit spheres) or
      // positive
      valid_ = false;
    }
  }
  delete [] chol;
  delete [] Qinv_a;
  delete [] Qinv_q;
}

// compute rho(tau), rho_tau_
void CglConicGD1Cut::compute_rho_tau() {
  rho_tau_ += tau_*alpha_*beta_;
}

 // compute q(tau), q_tau_
void CglConicGD1Cut::compute_q_tau() {
  int m = csize_;
  int n = m-num_rows_;
  double aux = -tau_ * (alpha_ + beta_) / 2;
  cblas_daxpy(n, aux, a_, 1, vecq_tau_, 1);
}

 // compute Q(tau), Q_tau_
void CglConicGD1Cut::compute_Q_tau() {
  int m = csize_;
  int n = m-num_rows_;
  matQ_tau_ = new double[n*n];
  cblas_dcopy(n*n, matQ_, 1, matQ_tau_, 1);
  cblas_dsyr(CblasColMajor, CblasUpper, n, tau_, a_, 1,
	     matQ_tau_, n);
}

// compute matrix as a part of the cut to add to model
void CglConicGD1Cut::compute_new_A() {
  //Initilize the new A matrix section
  int m = csize_;
  int n = m-num_rows_;
  new_matA_ = new double[n*m]();
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  int lwork = -1;
  double worksize[1];
  double * L = new double[n*n];
  cblas_dcopy(n*n, matQ_tau_, 1, L, 1);
  double * d = new double[n];
  // eigenvalue decomposition, L holds normalized eigenvectors,
  // d holds eigenvalues. Q = L^TDL
  dsyev_(&jobz, &uplo, &n, L, &n, d, worksize, &lwork, &info);
  lwork = (int) worksize[0];
  double * work = new double[lwork];
  dsyev_(&jobz, &uplo, &n, L, &n, d, work, &lwork, &info);
  if (cut_type_ > 0){
    double * LDinv = new double[n*n]();
    // todo(aykut) why do we need this? we already know the leading
    // variable of the conic constraint.
    varHead_ = -1;
    // the following loop is for computing LD^-1
    for(int i = 0; i < n; i++){
      //printf("d[%d] = %f\n", i, d[i]);
      if(fabs(d[i]) < 1e-2){
	valid_ = false;
      }
      if (d[i] < -EPS){
	//printf("head var = %d\n", i);
	varHead_ = i;
	cblas_daxpy(n, 1.0/sqrt(fabs(d[i])), L+n*i, 1, LDinv+n*i, 1);
	dirTestU_ = new double[n]();
	cblas_dcopy(n, L+n*i, 1, dirTestU_, 1);
      }
      else {
	cblas_daxpy(n, 1/sqrt(d[i]), L+n*i, 1, LDinv + n*i, 1);
      }
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, csize_, n,
		n, -1.0, matH_, csize_, LDinv, n, 0.0, new_matA_, csize_);
  }
  else {
    varHead_ = -1;
    for(int i=0; i<n; i++) {
      if(fabs(d[i]) < 1e-2){
	valid_ = false;
      }
      if (d[i] < -EPS) {
	varHead_ = i;
	cblas_daxpy(n, sqrt(fabs(d[i])), L+n*i, 1, new_matA_+i, n);
	dirTestU_ = new double[n]();
	cblas_dcopy(n, L+n*i, 1, dirTestU_, 1);
      }
      else {
	cblas_daxpy(n, sqrt(d[i]), L+n*i, 1, new_matA_+i, n);
      }
    }
  }
  delete[] work;
  delete[] L;
  delete[] d;
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
    if (dirTestV_ != 0)
      delete [] dirTestV_;
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


void CglConicGD1Cut::regularize() {
  // double alpha = disjunction_->get_c10();
  // double beta = disjunction_->get_c20();
  double * rega = new double[csize_];
  double tildeD = 0.0;
  double const * a = disjunction_->get_c1();
  //Computes V^T a
  cblas_dgemv(CblasColMajor, CblasTrans, csize_, csize_, 1.0, matV_,
	      csize_, a, 1, 0.0, rega, 1);
  //Computes  tildeD^(-1)v^T a
  for (int i=0; i<csize_; i++){
    tildeD = sqrt(fabs(eigQ_[i]));
    if (tildeD<EPS)
      tildeD = 1.0;
    rega[i] = rega[i]/tildeD;
  };
  double scaling = cblas_ddot(csize_, rega, 1, rega, 1);
  if (quad_type_== 'p') {
    tau_ = -1.0 / scaling;
    valid_ = false;
  }
  else {
    //std::cout<<"Is not paraboloid "<<std::endl;
    double nrma = cblas_dnrm2(csize_, rega, 1);
    cblas_dscal (csize_, 1/nrma, rega, 1);
    double * Qq = new double[csize_];
    solveSM(csize_, matQ_, vecq_, Qq);
    double aQp = cblas_ddot(csize_, a, 1, Qq, 1);
    double regalpha = alpha_ + aQp;
    double regbeta = beta_ + aQp;
    double qadricRHS = cblas_ddot(csize_, vecq_, 1, Qq, 1) - rho_;
    rho_ = 0.0;
    if (qadricRHS < -EPS)
      rho_ = 1.0;
    else if(qadricRHS > EPS)
      rho_ = -1.0;
    double denom = sqrt(fabs(qadricRHS));
    if (denom < EPS)
      denom = nrma;
    else
      denom *= nrma;
    regalpha /= denom;
    regbeta /= denom;
    double firsta2 = rega[0]*rega[0];
    double squarecoeff = (1 - 2*firsta2)*(regalpha - regbeta)*(regalpha - regbeta)/4;
    double linearcoeff = (-1)*( rho_*(1 - 2*firsta2) + regalpha*regbeta);
    double scalar = -rho_;
    tau1_ = tau2_ = 0.0;
    double root1, root2;
    root1 = root2 = 0.0;
    roots(squarecoeff, linearcoeff, scalar, root1, root2);
    tau1_ = std::min(root1, root2);
    tau2_ = std::max(root1, root2);
    if (firsta2 > (0.5 - EPS)) {
      tau_ = tau2_ / scaling;
      quad_type_ = 'c';
      valid_ = false;
      if (tau2_ < (1/(2*firsta2 - 1) - EPS) ){
	quad_type_ = 'k';
	valid_ = true;
      }
    }
    else {
      tau_ = tau1_ / scaling ;
      quad_type_ = 'c';
      valid_ = false;
      if (tau1_ > (1/(2*firsta2 - 1) + EPS) ){
	//std::cout<<"Is a cone\n "<<std::endl;
	quad_type_ = 'k';
	valid_ = true;
      }
    }
    delete [] Qq;
  }
  delete [] rega;
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
// 	newRowsub[i+1] = numVarP - elliDim_ + i;
//       }
//       /* In the original variables part we add and identity matrix */
//       newRowval[0] = 1.0;
//       /* Pointer to acces the values
// 	 of the entries for the new constraints */
//       double * indexA;
//       for(int i = 0; i < numVars_; i++){
// 	/* Get the indices for the new variables in the problem */
// 	newRowsub[0] = vars[i];
// 	/* Get the values of the entries of the new constraints */
// 	for(int j = 0; j < elliDim_; j++) {
// 	  indexA = newAslice_ + i + j * numVars_;
// 	  newRowval[j+1] = (fabs(*(indexA)) < ICLEPS)?0.0:*(indexA);
// 	}
// 	/* Add the constraints bounds, i.e the right hand side */
// 	solver->setRowBounds(numConP - numVars_ + i,
// 			     *(newbslice_+i),
// 			     *(newbslice_+i));
// 	/* Add the new row to the constraint matrix */
// 	solver->setRow(numConP - numVars_ + i, elliDim_+1,
// 		       newRowsub, newRowval);
//       }
//     }
//     else {
//       /* Register the indices of the new variables */
//       for(int i = 0; i < elliDim_; i ++) {
// 	newRowsub[i] = vars[i];
//       }
//       /* In the new variables part we add and identity matrix */
//       newRowval[elliDim_] = 1.0;
//       double * indexA;
//       for(int i = 0; i < numVars_; i++){
// 	/* Get the indices for the new variables in the problem */
// 	newRowsub[elliDim_] = numVarP - elliDim_ + i;
// 	/* Get the values of the entries of the new constraints */
// 	for(int j = 0; j < elliDim_; j++) {
// 	  indexA = newAslice_ + i + j * numVars_;
// 	  newRowval[j] = (fabs(*(indexA)) < ICLEPS)?0.0:*(indexA);
// 	}
// 	/* Add the constraints bounds, i.e the right hand side */
// 	solver->setRowBounds(numConP - numVars_ + i,
// 			     *(newbslice_+i),
// 			     *(newbslice_+i));
// 	/* Add the new row to the constraint matrix */
// 	solver->setRow(numConP - numVars_ + i, elliDim_+1,
// 		       newRowsub, newRowval);
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
// 	csub[l] = numVarP - elliDim_ + i;
// 	//std::cout<<"csub["<<l<<"] = "<<csub[l]<<"\n";
// 	l++;
//       }
//     }
//     solver->addCone(elliDim_, csub);
//     //printf("A conic cut added \n");
//   }
// }

bool CglConicGD1Cut::valid() const {
  return valid_;
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
  if (vecx0_)
    delete[] vecx0_;
  if (vecx0J_)
    delete[] vecx0J_;
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
  double * s = new double[m]; //Singular values of A
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
  double * work = new double[lwork];
  //   printf("fails Singular computing\n");
  // Here is where the actually computation is happening
  dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, u, &ldu,
          VT, &ldvt, work, &lwork, &info);
  delete[] work;
  delete[] s;
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

