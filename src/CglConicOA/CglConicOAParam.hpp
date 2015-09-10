#ifndef CglConicOAParam_H
#define CglConicOAParam_H

#include "CglConicParam.hpp"

// How to choose cone to generate cut.
// How to choose the row to generate cut.
//

class CglConicOAParam : public CglConicParam {
public:
  // constructor
  CglConicOAParam();
  // copy constructor
  CglConicOAParam(const CglConicOAParam & other);
  // copy assignment operator
  CglConicOAParam & operator=(const CglConicOAParam & other);
  // destructor
  virtual ~CglConicOAParam();
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

