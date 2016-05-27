#ifndef CglConicOAParam_H
#define CglConicOAParam_H

#include "CglConicParam.hpp"

/*!
  Represents the parameters used when generating cuts with CglConicOA.
  Can specify the following

  <ul>
    <li> How to choose cone to generate cut.
    <li> How to choose the row to generate cut.
    </ul>

  For now we generate cuts for all infeasible cones. The cone is considered
  feasible if x1 - ||x_2:n|| > -cone_tol (similarly for rotated cones).
 */

class CglConicOAParam : public CglConicParam {
  /// cone is feasible if x1 - ||x_2:n|| > -cone_tol
  double coneTol_;
public:
  // constructor
  CglConicOAParam(double coneTol);
  // notes(aykut) copy constructor and copy assignment operator that will be
  // defined by compiler is enough.

  // destructor
  virtual ~CglConicOAParam();
  // Virtual functions inherited.
  /// Clone
  virtual CglParam * clone() const;
  /// Get cone tolerance
  double coneTol() const { return coneTol_; }
  // /** Set INFINIT */
  // virtual void setINFINIT(const double inf);
  // /** Set EPS */
  // virtual void setEPS(const double eps);
  // /** Set EPS_COEFF */
  // virtual void setEPS_COEFF(const double eps_c);
  // /** Set MAX_SUPPORT */
  // virtual void setMAX_SUPPORT(const int max_s);
private:
  /// Disable copy constructor.
  CglConicOAParam();
};

#endif
