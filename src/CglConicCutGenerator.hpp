#ifndef CglConicCutGenerator_H
#define CglConicCutGenerator_H

// CGL headers
#include <CglTreeInfo.hpp>
// OsiConic headers
#include <OsiConicSolverInterface.hpp>

// this is an abstract base class for generating conic cuts.
// specific conic cut generators inherit this class.
//
//

// what is common to all conic cut generators?
// for now we assume there is not much in common.
// I will come back here when I implemented all specific cuts
// since I will know then what is common.
// following maybe?
// --cut generating cone
// ----size of cone
// ----members of the cut generating cone
// --disjunction?

class CglConicCutGenerator {
public:
  // default constructor
  CglConicCutGenerator();
  // copy constructor
  CglConicCutGenerator(CglConicCutGenerator const & other);
  // copy assignment operator
  CglConicCutGenerator & operator=(CglConicCutGenerator const & rhs);
  virtual ~CglConicCutGenerator();
  // clone method
  virtual CglConicCutGenerator * clone() const = 0;
  // Create C++ lines so set generator, I do not see why we have this.
  virtual std::string generateCpp( FILE * ) {return "";}

  /// This can be used to refresh any information
  virtual void refreshSolver(OsiConicSolverInterface * si) {}
  // generate cuts
  virtual void generateCuts(const OsiConicSolverInterface & si,
			    OsiConicCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());
  // generate linear/ordinary cuts.
  virtual void generateCuts(const OsiConicSolverInterface & si,
			    OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());

  // generate linear cuts for a solver interface with conic constraints
  // given explicitly
  virtual void generateCuts(OsiSolverInterface const & si, OsiCuts & cuts,
			    int num_cones,
			    OsiLorentzConeType const * cone_type,
			    int const * cone_size,
			    int const * const * members) {}
  // generate linear cuts for a solver interface with conic constraints
  // given explicitly
  virtual void generateCuts(OsiSolverInterface const & si, OsiCuts & cuts,
			    int num_cones,
			    OsiLorentzConeType const * cone_type,
			    int const * cone_size,
			    int const * const * members, int num_points) {}
  int aggressiveness() const;
  void setAggressiveness(int value);
  // set whether can do global cuts
  void setGlobalCuts(bool ability);
  bool canDoGlobalCuts() const;
  /**
     Returns true if may generate cuts in tree (rather than root node).
     Used so know if matrix will change in tree.  Really
     meant so column cut generators can still be active
     without worrying code.
     Default is true
  */
  virtual bool mayGenerateRowCutsInTree() const;
  /// Return true if needs optimal basis to do cuts
  virtual bool needsOptimalBasis() const;
  /// Return maximum length of cut in tree
  virtual int maximumLengthOfCutInTree() const
  { return COIN_INT_MAX;}
private:
  /**
     Aggressiveness - 0 = neutral, 100 is normal root node.
     Really just a hint to cut generator
  */
  int aggressive_;
  /// True if can do global cuts i.e. no general integers
  bool canDoGlobalCuts_;
};

// we define a class to represent a disjunction. This is for conic Cgl
// use only. Not intented for public use.
class Disjunction {
  int size_;
  double * c1_;
  double c10_;
  double * c2_;
  double c20_;
public:
  // defualt constructor
  Disjunction(int size,
	      double const * c1, double c10,
	      double const * c2, double c20): size_(size) {
    c1_ = new double[size];
    std::copy (c1, c1+size, c1_);
    c2_ = new double[size];
    std::copy (c2, c2+size, c2_);
    c10_ = c10;
    c20_ = c20;
  }
  // copy constructor
  Disjunction(Disjunction const & other) {
    // save size_
    size_ = other.size();
    // copy c1_
    if (c1_)
      delete[] c1_;
    c1_ = new double[size_];
    double const * other_c1 = other.get_c1();
    std::copy(other_c1, other_c1+size_, c1_);
    // copy c2_
    if (c2_)
      delete[] c2_;
    c2_ = new double[size_];
    double const * other_c2 = other.get_c2();
    std::copy(other_c2, other_c2+size_, c2_);
    // save c10 and c20
    c10_ = other.get_c10();
    c20_ = other.get_c20();
  }
  // copy assignment operator
  Disjunction & operator=(Disjunction const & rhs) {
    // save size_
    size_ = rhs.size();
    // copy c1_
    if (c1_)
      delete[] c1_;
    c1_ = new double[size_];
    double const * rhs_c1 = rhs.get_c1();
    std::copy(rhs_c1, rhs_c1+size_, c1_);
    // copy c2_
    if (c2_)
      delete[] c2_;
    c2_ = new double[size_];
    double const * rhs_c2 = rhs.get_c2();
    std::copy(rhs_c2, rhs_c2+size_, c2_);
    // save c10 and c20
    c10_ = rhs.get_c10();
    c20_ = rhs.get_c20();
    return *this;
  }
  ~Disjunction() { delete[] c1_; delete[] c2_; }
  double const * get_c1() const { return c1_; }
  double get_c10() const { return c10_; }
  double const * get_c2() const { return c2_; }
  double get_c20() const { return c20_; }
  int size() const {return size_;}
};

#endif
