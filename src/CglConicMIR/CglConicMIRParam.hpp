#ifndef CglConicMIRParam_H
#define CglConicMIRParam_H

#include "CglConicParam.hpp"

// How to choose cone to generate cut.
// How to choose the row ro generate cut.
// What happens if the cone is in canonical form?
//

class CglConicMIRParam : public CglConicParam {
public:
  // constructor
  CglConicMIRParam();
  // copy constructor
  CglConicMIRParam(const CglConicMIRParam & other);
  // copy assignment operator
  CglConicMIRParam & operator=(const CglConicMIRParam & other);
  // destructor
  virtual ~CglConicMIRParam();
  // Virtual functions inherited.
  /// Clone
  virtual CglParam * clone() const;
  // /** Set INFINIT */
  // virtual void setINFINIT(const double inf);
  // /** Set EPS */
  // virtual void setEPS(const double eps);
  // /** Set EPS_COEFF */
  // virtual void setEPS_COEFF(const double eps_c);
  // /** Set MAX_SUPPORT */
  // virtual void setMAX_SUPPORT(const int max_s);
private:

};

#endif

