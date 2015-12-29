#ifndef CglConicIPMParam_H
#define CglConicIPMParam_H

#include "CglConicParam.hpp"

// How to choose cone to generate cut.
// How to choose the row to generate cut.
//

class CglConicIPMParam : public CglConicParam {
public:
  // constructor
  CglConicIPMParam();
  // copy constructor
  CglConicIPMParam(const CglConicIPMParam & other);
  // copy assignment operator
  CglConicIPMParam & operator=(const CglConicIPMParam & other);
  // destructor
  virtual ~CglConicIPMParam();
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
