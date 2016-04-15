#ifndef CglConicIPMintParam_H
#define CglConicIPMintParam_H

#include "CglConicParam.hpp"

// How to choose cone to generate cut.
// How to choose the row to generate cut.
//

class CglConicIPMintParam : public CglConicParam {
public:
  // constructor
  CglConicIPMintParam();
  // copy constructor
  CglConicIPMintParam(const CglConicIPMintParam & other);
  // copy assignment operator
  CglConicIPMintParam & operator=(const CglConicIPMintParam & other);
  // destructor
  virtual ~CglConicIPMintParam();
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
