#ifndef CglConicGD1Param_H
#define CglConicGD1Param_H

#include "CglConicParam.hpp"

// How to choose cone to generate cut.
// How to choose the row to generate cut.
//

class CglConicGD1Param : public CglConicParam {
public:
  // constructor
  CglConicGD1Param();
  // copy constructor
  CglConicGD1Param(const CglConicGD1Param & other);
  // copy assignment operator
  CglConicGD1Param & operator=(const CglConicGD1Param & other);
  // destructor
  virtual ~CglConicGD1Param();
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

