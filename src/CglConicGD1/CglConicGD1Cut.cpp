#include "CglConicGD1Cut.hpp"
#include "CglConicGD1Help.hpp"

#define EPS 1e-10

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
  // compute x0 and x0j
  vecx0_ = new double[csize_];
  double const * x0 = solver_->getColSolution();
  for (int i=0; i<csize_; ++i) {
    vecx0_[i] = x0[cmembers_[i]];
  }
  vecx0J_ = new double[csize_];
  std::copy(vecx0_, vecx0_+csize_, vecx0J_);
  vecx0J_[0] = -1.0*vecx0J_[0];
  // create disjunction from disj_var_
  //compute_disjunction();
  generate_cut();
}

CglConicGD1Cut::CglConicGD1Cut(double const * matQ, double const * vecq,
                               double rho,
                               int num_rows, int num_cols, double const * matH,
                               double const * x0,
                               int const * cmembers,
                               int dis_var,
                               Disjunction const & disjunction)
    : solver_(NULL) {
  num_rows_ = num_rows-num_cols;
  csize_ = num_rows;
  // we do not know the following when the input is quadric
  rows_ = NULL;
  cindex_ = -1;
  // set arrays to 0
  matA_ = NULL;
  cmembers_ = NULL;
  cone_members_ = NULL;
  //lctype_= -1;
  //ctype_ = OSI_LORENTZ;
  // compute what we know from input
  matH_ = new double[num_rows*num_cols];
  std::copy(matH, matH+num_rows*num_cols, matH_);
  matV_ = 0;
  matQ_ = new double[num_cols*num_cols];
  std::copy(matQ, matQ+num_cols*num_cols, matQ_);
  vecq_ = new double[num_cols];
  std::copy(vecq, vecq+num_cols, vecq_);
  rho_ = rho;
  eigQ_ = 0;
  vecx0J_ = 0;
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
  // todo(aykut) check if csize_ is same as num_rows_
  // I do not know what to do in this case
  if (num_rows_>=csize_) {
    valid_ = false;
    return;
  }
  dis_var_ = dis_var;
  // create disjunction.
  //disjunction_ = new Disjunction(disjunction);
  a_ = new double[num_cols];
  std::copy(disjunction.get_c1(), disjunction.get_c1()+num_cols, a_);
  alpha_ = disjunction.get_c10();
  beta_ = disjunction.get_c20();
  // set vecx0
  vecx0_ = new double[num_rows];
  std::copy(x0, x0+num_rows, vecx0_);
  // copy cmembers
  cmembers_ = new int[num_rows];
  std::copy(cmembers, cmembers+num_rows, cmembers_);
  generate_cut_dual();
}

// create disjunction from disj_var_
void CglConicGD1Cut::compute_disjunction() {
  int m = csize_;
  int n = m-num_rows_;
  // u is the disjunction coefficient in x-space,
  // a is the disjunction coefficient in w-space,
  // we compute a from u using a <- H^T u
  // todo(aykut) what if a is 0?
  // == compute u first
  double * u = new double[m]();
  int new_dis_index = 0;
  for (int i=0; i<m; ++i) {
    if (cmembers_[i]!=dis_var_)
      new_dis_index++;
    else
      break;
  }
  u[new_dis_index] = 1.0;
  disjunction_ = new Disjunction(csize_, u,
                                 solver_->getColSolution()[dis_var_],
                                 u,
                                 solver_->getColSolution()[dis_var_]+1.0);
  // == compute a from u using a <- H^T u
  a_ = new double[n]();
  cblas_dgemv(CblasColMajor, CblasTrans, m, n, 1.0, matH_, m, u, 1, 0.0, a_, 1);
  // == normalize a
  double norm_of_a = cblas_dnrm2(n, a_, 1);
  if (norm_of_a<EPS) {
    // todo(aykut) does this mean that problem is infeasible?
    std::cerr << "Numerical problems. Disjunction coefficient is very close to 0."
              << std::endl;
    valid_ = false;
  }
  // todo(aykut) do this using blas
  for (int i=0; i<n; ++i) {
    a_[i] = a_[i]/norm_of_a;
  }
  // end of computing a
  // compute alpha and beta
  double uT_x0 = cblas_ddot(m, u, 1, vecx0_, 1);
  alpha_ = floor(solver_->getColSolution()[dis_var_])+1.0 - uT_x0;
  beta_ = floor(solver_->getColSolution()[dis_var_]) - uT_x0;
  alpha_ = alpha_ / norm_of_a;
  beta_ = beta_ / norm_of_a;
  // create coefficient matrix
  //disjunction_ = new Disjunction(csize_, a_, alpha, a_, beta);
  delete[] u;
}

void CglConicGD1Cut::generate_cut() {
  valid_ = true;
  compute_matrixA();
  compute_matrixH();
  compute_disjunction();
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
  int m = csize_;
  int n = m-num_rows_;
  vecq_ = new double[n]();
  cblas_dgemv (CblasColMajor, CblasTrans, m, n, 1.0, matH_, m, vecx0J_, 1,
               0.0, vecq_, 1);
  // print vector q
  std::cout << "==================== vector q ==================== "
            << std::endl;
  for (int i=0; i<n; ++i)
    std::cout << vecq_[i] << " ";
  std::cout << std::endl;
}

void CglConicGD1Cut::compute_rho() {
  int m = csize_;
  rho_ = - (vecx0J_[0]*vecx0J_[0]);
  rho_ += cblas_ddot(m-1, vecx0J_+1, 1, vecx0J_+1, 1);
  // print rho
  std::cout << "==================== rho ==================== "
            << std::endl;
  std::cout << rho_ << std::endl;
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
  // change matH_ to col order
  // double * tempH new double[(num_cols-num_rows_)*num_cols];
  // cblas_dcopy((num_cols-num_rows_)*num_cols, matH_, 1, tempH, 1);
  // for(int i=0; i<(num_cols-num_rows_); ++i)
  //   cblas_dcopy(num_cols, (VT + num_rows_ + i), num_cols, (matH_+i*num_cols), 1);


  // print H
  {
    int m = num_cols;
    int n = num_cols-num_rows_;
    std::cout << "==================== H ===================="
              << std::endl;
    for (int i=0; i<m; ++i) {
      for (int j=0; j<n; ++j) {
        std::cout << matH_[j*m+i] << " ";
      }
      std::cout << std::endl;
    }
  }
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
  // print Q
  std::cout << "==================== Q ===================="
            << std::endl;
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j) {
      std::cout << matQ_[j*m+i] << " ";
    }
    std::cout << std::endl;
  }
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
  double aQinv_q = cblas_ddot(n, a_, 1, Qinv_q, 1);
  //quadratic equation
  double quadTerm = (alpha_ - beta_)*(alpha_ - beta_)*aQinv_a;
  double linearTerm = 4*aQinv_a*(qQinv_q - rho_)
    - (alpha_ + beta_ + 2*aQinv_q)*(alpha_ + beta_ + 2*aQinv_q)
    + (alpha_ - beta_)*(alpha_ - beta_);
  double indTerm = 4*(qQinv_q - rho_);
  double intr1 = (beta_ + qQinv_a)*(beta_ + qQinv_a) - aQinv_a*(qQinv_q - rho_);
  double intr2 = (alpha_ + qQinv_a)*(alpha_ + qQinv_a) - aQinv_a*(qQinv_q - rho_);
  if(intr1 > -EPS || intr2 > -EPS)
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
      linear_cut_rhs_ = floor(vecx0_[dis_var_])+1.0;
    }
  }
  if (intr2 > -EPS){
    if (intr1 < -EPS){
      // one piece of disjunction does not intersect with quadric (a,alpha)
      // only (a,alpha) intersects with quadric
      cut_type_ = 3;
      linear_cut_ind_ = new int[1];
      linear_cut_ind_[0] = dis_var_;
      linear_cut_coef_ = new double[1];
      linear_cut_coef_[0] = 1.0;
      linear_cut_rhs_ = floor(vecx0_[dis_var_]);
    }
  }
  // print coeffcients in the qad formula
  {
    std::cout << "==================== quad formula ==================== "
              << std::endl
              << quadTerm
              << " "
              << linearTerm
              << " "
              << indTerm
              << std::endl;
  }

  // quadformula returns 1 when the roots are imaginary which will
  // not happen if quadric intersects with disjunction hyperplanes
  tau_ = quad_formula(quadTerm, linearTerm, indTerm);
  if(fabs(tau_ + (1.0/aQinv_a)) < 1e-6 || tau_ > EPS) {
    // this means tau is very close to zero (in case of unit spheres) or
    // positive
    valid_ = false;
  }
  delete [] chol;
  delete [] Qinv_a;
  delete [] Qinv_q;
  // print tau
  std::cout << "==================== tau ==================== "
            << std::endl;
  std::cout << tau_ << std::endl;
}

// compute rho(tau), rho_tau_
void CglConicGD1Cut::compute_rho_tau() {
  rho_tau_ += tau_*alpha_*beta_;
  // print rho
  std::cout << "==================== rho(tau) ==================== "
            << std::endl;
  std::cout << rho_tau_ << std::endl;
}

 // compute q(tau), q_tau_
void CglConicGD1Cut::compute_q_tau() {
  int m = csize_;
  int n = m-num_rows_;
  double aux = -tau_ * (alpha_ + beta_) / 2;
  cblas_daxpy(n, aux, a_, 1, vecq_tau_, 1);
  // print vector q
  std::cout << "==================== vector q(tau) ==================== "
            << std::endl;
  for (int i=0; i<n; ++i)
    std::cout << vecq_tau_[i] << " ";
  std::cout << std::endl;
}

 // compute Q(tau), Q_tau_
void CglConicGD1Cut::compute_Q_tau() {
  int m = csize_;
  int n = m-num_rows_;
  matQ_tau_ = new double[n*n];
  cblas_dcopy(n*n, matQ_, 1, matQ_tau_, 1);
  cblas_dsyr(CblasColMajor, CblasUpper, n, tau_, a_, 1,
             matQ_tau_, n);
  // print Q(tau)
  std::cout << "==================== Q(tau) ===================="
            << std::endl;
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j) {
      if (i<=j)
        std::cout << matQ_tau_[j*n+i] << " ";
      else
        std::cout << matQ_tau_[i*n+j] << " ";
    }
    std::cout << std::endl;
  }
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
    delete[] LDinv;
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

  // print new A
  {
    std::cout << "==================== new A ===================="
              << std::endl;
    for (int i=0; i<m; ++i) {
      for (int j=0; j<n; ++j)
        std::cout << new_matA_[j*m+i] << " ";
      std::cout << std::endl;
    }
  }
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
  // solve symmetric system Q(tau)x = q(tau)
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

  // print new_rhs
  {
    std::cout << "==================== new rhs ===================="
              << std::endl;
    for (int i=0; i<m; ++i)
      std::cout << new_rhs_[i] << " ";
    std::cout << std::endl;
  }
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

void CglConicGD1Cut::generate_cut_dual() {
  valid_ = true;
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
