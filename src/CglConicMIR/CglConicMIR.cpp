#include "CglConicMIR.hpp"
#include <string>
#include <sstream>
#include <iomanip>

CglConicMIR::CglConicMIR(): param_(0), solver_(0) {
  param_ = new CglConicMIRParam();
  binding_ = 0;
  num_cuts_ = 0;
}

// copy constructor
CglConicMIR::CglConicMIR(const CglConicMIR & other) {
  // copy param_
  param_ = new CglConicMIRParam(*(other.getParam()));
  binding_ = 0;
  num_cuts_ = 0;
}

// copy assignment operator
CglConicMIR & CglConicMIR::operator=(const CglConicMIR & rhs) {
  // copy param_
  param_ = new CglConicMIRParam(*(rhs.getParam()));
  binding_ = 0;
  num_cuts_ = 0;
  return *this;
}

CglConicMIR::~CglConicMIR() {
  if (param_)
    delete param_;
  if (binding_)
    delete[] binding_;
}

void CglConicMIR::setParam(const CglConicMIRParam & param) {
  param_ = new CglConicMIRParam(param);
}

void CglConicMIR::generateCuts(const OsiConicSolverInterface & si,
			       OsiConicCuts & cs,
			       const CglTreeInfo info) {
}

// generates MIR cut.
void CglConicMIR::generateAndAddCuts(OsiConicSolverInterface & si,
				     const CglTreeInfo info) {
  add_slacks(si);
  solver_ = &si;
  // decide variable to generate cut, it should be a member of cone.
  // decide row to use, it should have a nonzero coef for cut generating
  // variable
  std::set<int> cut_rows = choose_cut_row(si);
  // todo(aykut) pick a variable, then generate cuts for the varaible from cut_rows.
  // and then maybe choose the deepest cut? can we estimate how deep is the cut?
  //si.writeMps("before_cuts");
  double const * sol = si.getColSolution();
  //update_binding_cones();
  int num_cones = si.getNumCones();
  //OsiConeType * types = new OsiConeType[num_cones];
  OsiLorentzConeType * lc_types = new OsiLorentzConeType[num_cones];
  si.getConeType(lc_types);
  int * sizes = new int[num_cones];
  si.getConeSize(sizes);
  std::cout << std::setw(6) << "Cone "
	    << std::setw(10) << "Member"
	    << std::setw(10) << "Var"
	    << std::setw(10) << "Row"
	    << std::endl;
  int ccone = 0;
  int cmember = 0;
  if (lc_types[ccone]==OSI_QUAD)
    cmember = 1;
  else if (lc_types[ccone]==OSI_RQUAD)
    cmember = 2;
  std::set<int>::const_iterator it;
  for (it=cut_rows.begin(); it!=cut_rows.end(); ++it) {
    // add one cut for one row. just decide on the cone and
    // member
    OsiLorentzConeType t;
    int size;
    int * members;
    si.getConicConstraint(ccone, t, size, members);
    std::cout << std::setw(6) << ccone
	      << std::setw(10) << cmember
	      << std::setw(10) << members[cmember]
	      << std::setw(10) << *it
	      << std::endl;
    int success = add_cut(si, ccone, members[cmember], *it);
    //si.writeMps("after_cut");
    delete[] members;
    if (success) {
      num_cuts_++;
      cmember++;
    }
    // update cut_cone and cut_member
    if (cmember==sizes[ccone]) {
      ccone++;
      if (ccone>=num_cones) {
	break;
      }
      if (lc_types[ccone]==OSI_QUAD)
	cmember=1;
      else
	cmember=2;
    }
    if (ccone>=num_cones) {
      break;
    }
  }
  delete[] lc_types;
  delete[] sizes;
}

void CglConicMIR::add_slacks(OsiConicSolverInterface & si) {
  // check number of equality constraints
  int num_rows = si.getNumRows();
  char const * row_sense = si.getRowSense();
  int num_eq_rows = 0;
  for (int i=0; i<num_rows; ++i) {
    if (row_sense[i]=='E')
      num_eq_rows++;
  }
  if (num_eq_rows>100)
    return;
  // modify rows by adding slacks.
  CoinPackedMatrix const * mat = si.getMatrixByRow();
  double const * elements = mat->getElements();
  int const * cols = mat->getIndices();
  int rows[1];
  double vals[1];
  for (int i=0; i<num_rows; ++i) {
    if (row_sense[i]=='L') {
      // add slack column
      rows[0] = i;
      vals[0] = 1.0;
      si.addCol(1, rows, vals, 0.0, si.getInfinity(),
		0.0, std::string("s"));
      // change constraint type
      si.setRowType(i, 'E', si.getRightHandSide()[i], 0.0);
    }
    else if (row_sense[i]=='G') {
      // add slack column
      rows[0] = i;
      vals[0] = -1.0;
      si.addCol(1, rows, vals, 0.0, si.getInfinity(),
		0.0, std::string("s"));
      // change constraint type
      si.setRowType(i, 'E', si.getRightHandSide()[i], 0.0);
    }
  }
  // replace free variables.
  int num_cols = si.getNumCols();
  num_rows = si.getNumRows();
  int old_num_rows = num_rows;
  double const * lb = si.getColLower();
  mat = si.getMatrixByRow();
  elements = mat->getElements();
  cols = mat->getIndices();
  num_cols = si.getNumCols();
  int first;
  int last;
  int * new_cols;
  double * new_vals;
  int xi_cols[1];
  double xi_vals[1];
  int xi_rows[1];
  std::map<int,int> positive;
  std::map<int,int> negative;
  for (int i=0; i<num_cols; ++i) {
    if (lb[i]<0.0) {
      xi_cols[0] = i;
      xi_vals[0] = 1.0;
      // add row for x_i
      si.addRow(1, xi_cols, xi_vals, 0.0, 0.0);
      // add cols for x_i+ and x_i-
      xi_rows[0] = num_rows;
      xi_vals[0] = -1.0;
      positive[i] = num_cols;
      num_cols++;
      si.addCol(1, xi_rows, xi_vals, 0.0, si.getInfinity(),
		0.0, std::string("xi+"));
      xi_vals[0] = 1.0;
      negative[i] = num_cols;
      num_cols++;
      si.addCol(1, xi_rows, xi_vals, 0.0, si.getInfinity(),
		0.0, std::string("xi-"));
      num_rows++;
    }
  }
  mat = si.getMatrixByRow();
  elements = mat->getElements();
  cols = mat->getIndices();
  for (int i=0; i<old_num_rows; ++i) {
    // replace all free variables at row i
    int first = mat->getVectorFirst(i);
    int last = mat->getVectorLast(i);
    int num_elem = last-first;
    // replace the free variable
    new_cols = new int[2*num_elem];
    // add constraints
    new_vals = new double[2*num_elem];
    delete[] new_cols;
    delete[] new_vals;
  }
}

int CglConicMIR::getNumCutsAdded() const {
  return num_cuts_;
}


void CglConicMIR::update_binding_cones() {
  if (solver_==0) {
    std::cerr << "I do not have a solver interface!" << std::endl;
    throw std::exception();
  }
  int num_cones = solver_->getNumCones();
  double const * sol = solver_->getColSolution();
  if (binding_==0)
    binding_ = new int[num_cones]();
  for (int i=0; i<num_cones; ++i) {
    OsiLorentzConeType cut_cone_type;
    int cut_cone_size = -1;
    int * cut_cone_members = 0;
    solver_->getConicConstraint(i, cut_cone_type, cut_cone_size,
			  cut_cone_members);
    // check if it is binding or not?
    double lhs = 0.0;
    // x1 or 2*x1*x2 depending on cone type
    double term1 = 0.0;
    // ||x_{2:n}|| or x_3^2 + ... + x_{cone_size}^2 depending on cone type
    double term2 = 0.0;
    int start = -1;
    if (cut_cone_type==OSI_QUAD) {
      start = 1;
      term1 = sol[cut_cone_members[0]];
    }
    else {
      start = 2;
      term1 = 2*sol[cut_cone_members[0]]*sol[cut_cone_members[1]];
    }
    for (int j=start; j<cut_cone_size; ++j) {
      term2 += sol[cut_cone_members[j]]*sol[cut_cone_members[j]];
    }
    if (cut_cone_type==OSI_QUAD) {
      term2 = sqrt(term2);
    }
    //std::cout << sol[cut_cone_members[0]] - sqrt(lhs) << std::endl;
    if (term1 - term2 < 1e-5) {
      binding_[i] = 1;
    }
    delete[] cut_cone_members;
  }
}

int CglConicMIR::add_cut(OsiConicSolverInterface & si, int cut_cone,
			  int cut_var, int cut_row) {
  // generate cuts
  int num_rows = si.getNumRows();
  int num_cols = si.getNumCols();
  double rhs_cut_row = si.getRightHandSide()[cut_row];
  int first = si.getMatrixByRow()->getVectorFirst(cut_row);
  int last = si.getMatrixByRow()->getVectorLast(cut_row);
  int const * indices = si.getMatrixByRow()->getIndices();
  double const * elements = si.getMatrixByRow()->getElements();
  int num_elem  = last-first;
  // +1 is in case we need when cut_var is not present in the row
  int * cols = new int[num_elem+1];
  double * value = new double[num_elem+1];
  std::copy(indices+first, indices+last, cols);
  std::copy(elements+first, elements+last, value);
  int flag = 0;
  for (int i=first; i<last; ++i) {
    if (indices[i]==cut_var)  {
      value[i-first] = value[i-first] - 1.0;
      flag=1;
      break;
    }
  }
  if (flag==0) {
    cols[num_elem] = cut_var;
    value[num_elem] = -1.0;
    num_elem = num_elem+1;
  }
  double * neg_value = new double[num_elem];
  for (int i=0; i<num_elem; ++i) {
    neg_value[i] = -value[i];
  }
  int t_e[2] = {num_rows, num_rows+1};
  double t_v[2]  = {1.0, 1.0};
  // add cut
  // we have +1 for variable t indexed as num_cols
  int * cut_e = new int[num_elem+1];
  double * cut_v = new double[num_elem+1];
  std::copy(cols, cols+num_elem, cut_e);
  // index of variable t
  cut_e[num_elem] = num_cols;
  // find fractional variable and set alpha
  double alpha = 0.0;
  // from paper:
  // moreover if alpha is chosen such that alpha = a_j and b/a_j>0
  // for some j in {integer var} and a_i<=b for all i in {integer var}\{j}
  // then cut is facet-defining for conv(S).
  // this is the integer var with a fractional value
  // (todo) a better way to choose this.
  int frac_var = -1;
  std::vector<int> fractionals = si.getFractionalIndices();
  std::set<int> frac_candidates;
  for (int i=0; i<last-first; ++i) {
    if (value[i]==0.0)
      continue;
    if (std::find(fractionals.begin(), fractionals.end(), cols[i])!=fractionals.end()) {
      frac_var = i;
      std::cout << "Fractional variable is " << cols[frac_var] << " value is "
		<< si.getColSolution()[cols[frac_var]] << std::endl;
      break;
    }
  }
  if (frac_var==-1) {
    std::cout << "Fractional variable not detected." << std::endl;
    std::cout << "Could not generate cut." << std::endl;
    return 0;
    //throw std::exception();
  }
  alpha = 4.0*value[frac_var];
  //alpha = value[frac_var];
  double f_alpha = rhs_cut_row/alpha - floor(rhs_cut_row/alpha);
  for (int i=0; i<num_elem; ++i) {
    if (si.isInteger(cols[i])) {
      cut_v[i] = phi(value[i]/alpha, f_alpha);
    }
    else {
      cut_v[i] = -1.0/fabs(alpha);
    }
  }
  cut_v[num_elem] = -1.0/fabs(alpha);
  double cut_rhs = phi(rhs_cut_row/alpha, f_alpha);
  // see how much the cut is violated. t<-x_i
  double viol = 0.0;
  double const * sol = solver_->getColSolution();
  double lhs = 0.0;
  for (int i=0; i<num_elem; ++i) {
    lhs += sol[cut_e[i]]*cut_v[i];
  }
  lhs += fabs(sol[cut_var])*cut_v[num_elem];
  if ((lhs-cut_rhs)>1e-4) {
    std::cout << "measure is " << lhs-cut_rhs << ". Adding cut..." << std::endl;
    // add cut
    // add t rows
    si.addRow(num_elem, cols, value, rhs_cut_row, si.getInfinity());
    si.addRow(num_elem, cols, neg_value, -rhs_cut_row, si.getInfinity());
    // add t column
    si.addCol(2, t_e, t_v, 0.0, si.getInfinity(), 0.0);
    // add cut
    si.addRow(num_elem+1, cut_e, cut_v, -si.getInfinity(), cut_rhs);
    // modify cone
    // remove old one and add new one.
    OsiConeType cut_cone_type;
    OsiLorentzConeType lc_cut_cone_type;
    int cut_cone_size = -1;
    int * cut_cone_members = 0;
    si.getConicConstraint(cut_cone, lc_cut_cone_type, cut_cone_size,
			  cut_cone_members);
    // replace t with cut_var in the cone
    for (int i=0; i<cut_cone_size; ++i) {
      if (cut_cone_members[i]==cut_var) {
	cut_cone_members[i] = num_cols;
	break;
      }
    }
    si.modifyConicConstraint(cut_cone, lc_cut_cone_type, cut_cone_size, cut_cone_members);
    delete[] cut_cone_members;
  }
  else {
    std::cout << "measure is " << lhs-cut_rhs << ". Skipping cut..." << std::endl;
  }
  delete[] cols;
  delete[] value;
  delete[] cut_e;
  delete[] cut_v;
  delete[] neg_value;
  if ((lhs-cut_rhs)>0.0) {
    return 1;
  }
  else {
    return 0;
  }
}


std::set<int> CglConicMIR::choose_cut_row(OsiConicSolverInterface const & si) {
  std::set<int> cut_row;
  int num_cones = si.getNumCones();
  int num_cols = si.getNumCols();
  int num_rows = si.getNumRows();
  char const * row_sense = si.getRowSense();
  // get set of integer variables that have fractional value
  std::set<int> cut_var_candidate;
  double const * sol = si.getColSolution();
  for (int i=0; i<num_cols; ++i) {
    if (si.isInteger(i)) {
      double value = sol[i];
      int vfloor = floor(value);
      int vceil = vfloor+1;
      double EPS = 1e-5;
      if ((vfloor+EPS < value) and (value < vceil-EPS)) {
	cut_var_candidate.insert(i);
      }
    }
  }
  // find a row that has a fractional integer variable and a cone variable
  CoinPackedMatrix const * mat;
  mat = si.getMatrixByRow();
  double const * col_lb = si.getColLower();
  for (int i=0; i<num_rows; ++i) {
    if (row_sense[i]!='E')
      continue;
    int first = mat->getVectorFirst(i);
    int last = mat->getVectorLast(i);
    int num_elem  = last-first;
    int const * indices = mat->getIndices();
    double const * elements = mat->getElements();
    // if it has a variable that does not have 0-lower-bound, skip.
    int flag=0;
    for (int j=first; j<last; ++j) {
      if (col_lb[indices[j]]<0.0) {
	flag=1;
	break;
      }
    }
    if (flag)
      continue;
    // check whether row has a fractional integer variable
    for (int j=first; j<last; ++j) {
      if (cut_var_candidate.find(indices[j])!=cut_var_candidate.end()
	  and elements[j]!=0.0) {
	cut_row.insert(i);
	break;
      }
    }
  }
  return cut_row;
  //std::copy(cut_row.begin(), cut_row.end(), cut_row_in);
  // once we determined the cut row, we can generate cut using any cone and
  // any cone member we want.
}

double CglConicMIR::phi(double a, double f) {
  int n= floor(a);
  double value;
  if (a<n+f) {
    value = (1-2*f)*n-(a-n);
  }
  else {
    value = (1-2*f)*n+(a-n)-2*f;
  }
  return value;
}


/// Clone
CglConicCutGenerator * CglConicMIR::clone() const {
  //
  //CglConicMIR * new_cut = new CglConicMIR(*this);
  CglConicMIR * new_cut;
  return new_cut;
}

/// Create C++ lines to get to current state
std::string CglConicMIR::generateCpp( FILE * fp) {
  return std::string();
}
