#ifndef CglConicGD2Param_H
#define CglConicGD2Param_H

#include "CglConicParam.hpp"

class CglConicGD2Param : public CglConicParam {
public:
  // constructor
  CglConicGD2Param();
  // copy constructor
  CglConicGD2Param(const CglConicGD2Param & other);
  // copy assignment operator
  CglConicGD2Param & operator=(const CglConicGD2Param & other);
  // destructor
  virtual ~CglConicGD2Param();
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
};

#endif

