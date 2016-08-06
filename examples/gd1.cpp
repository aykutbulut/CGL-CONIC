//-----------------------------------------------------
// Simple example usage of the cut generation library.
//
// This sample program adds general disjunction cuts from
// Belotti et al. and resolves the relaxation problem to see
// whether the bound is improved.
//
// usage:
//   solve_with_gd1 mpsFileName
// example:
//   solve_with_mir ../../Data/Sample/p0033
//-----------------------------------------------------

// STDLIB headers
#include <cassert>
#include <iostream>
#include <string>
#include <cassert>
// CoinUtils headers
#include <CoinError.hpp>
#include <CoinWarmStartBasis.hpp>
// OSI headers
#include <OsiCuts.hpp>
// OSICONIC header
#include <OsiConicSolverInterface.hpp>
#include <OsiConicCuts.hpp>
// COLA headers
#include <ColaModel.hpp>

// get IPM solver
#if defined(__OSI_MOSEK__)
  // use mosek as IPM solver
  #include <OsiMosekSolverInterface.hpp>
  //#define IPM_SOLVER OsiMosekSolverInterface
  typedef OsiMosekSolverInterface IPM_SOLVER;
#elif defined(__OSI_CPLEX__)
  // use cplex as IPM solver
  #include <OsiCplexSolverInterface.hpp>
  typedef OsiCplexSolverInterface IPM_SOLVER;
#elif defined(__OSI_IPOPT__)
  // use ipopt
  #include <OsiIpoptSolverInterface.hpp>
  typedef OsiIpoptSolverInterface IPM_SOLVER;
#endif

// CGL headers
// #include "CglKnapsackCover.hpp"
// #include "CglSimpleRounding.hpp"
// #include "CglGMI.hpp"
// #include "CglGomory.hpp"
// #include "CglMixedIntegerRounding.hpp"
// #include "CglMixedIntegerRounding2.hpp"
// Conic CGL headers
#include "CglConicGD1.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

int main(int argc, const char *argv[]) {
  // If no parms specified then use these
  string mpsFileName = argv[1];
  try {
    // Instantiate a specific solver interface
    //OsiConicSolverInterface * si = new ColaModel();
    OsiConicSolverInterface * si = new IPM_SOLVER();
    // Read file describing problem
    si->readMps(mpsFileName.c_str(),"mps");
    // Solve continuous problem
    si->initialSolve();
    // Save the orig socp relaxation value for
    // comparisons later
    double origSocpObj = si->getObjValue();
    // Instantiate cut generator
    CglConicGD1 cg(si);
    bool equalObj;
    CoinRelFltEq eq(0.0001);
    int num_cut = 0;
    double obj;
    //---------------------------------------------------
    // Keep applying cuts until no more cuts are generated
    //---------------------------------------------------
    do {
      // Get current solution value
      obj = si->getObjValue();
      // Generate and apply cuts
      OsiConicSolverInterface * nsi = cg.generateAndAddCuts(*si);
      delete si;
      si = nsi;
      //si->writeMps("after_cut");
      si->resolve();
      equalObj = eq(si->getObjValue(), obj);
    } while (!equalObj);
    // double const * sol = si->getColSolution();
    // Print total number of cuts applied,
    // and total improvement in the SOCP objective value
    cout <<endl <<endl;
    cout << "----------------------------------------------------------"
         <<endl;
    cout << "Cut generation phase completed:" <<endl;
    cout << "   " << cg.getNumCutsAdded() << " many cuts added." << endl;
    cout << "   changing the SOCP objective value from " << origSocpObj
         << " to " << si->getObjValue() <<endl;
    cout << "----------------------------------------------------------"
         <<endl;
    cout <<endl <<endl;
  delete si;
  }
  catch (CoinError e) {
    cout << e.className() << "::" << e.methodName() << " - " << e.message() << endl;
  }
  return 0;
}
