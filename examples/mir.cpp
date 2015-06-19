//-----------------------------------------------------
// Simple example usage of the cut generation library.
//
// This sample program iteratively tightens a
// given formulation by adding conic cuts, then calls
// branch-and-bound to solve the tightened
// formulation.
//
// usage:
//   solve_with_mir mpsFileName objectiveSense
// where:
//   mpsFileName: Name of an mps file (without the
//                file extension)
// example:
//   solve_with_mir ../../Data/Sample/p0033
//-----------------------------------------------------

// STDLIB headers
#include <cassert>
#include <iostream>
#include <string>
#include <cassert>
// CoinUtils headers
#include "CoinError.hpp"
#include "CoinWarmStartBasis.hpp"
// OSI headers
#include "OsiCuts.hpp"
#include "OsiClpSolverInterface.hpp"
// OSICONIC header
#include <OsiConicSolverInterface.hpp>
#include <OsiConicCuts.hpp>
// COLA headers
#include <ColaModel.hpp>
// CGL headers
// #include "CglKnapsackCover.hpp"
// #include "CglSimpleRounding.hpp"
// #include "CglGMI.hpp"
// #include "CglGomory.hpp"
// #include "CglMixedIntegerRounding.hpp"
// #include "CglMixedIntegerRounding2.hpp"
// Conic CGL headers
#include "CglConicMIR.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

int main(int argc, const char *argv[])
{
  // If no parms specified then use these
  string mpsFileName = argv[1];
  try {
    // Instantiate a specific solver interface
    OsiConicSolverInterface * si = new ColaModel();
    // Read file describing problem
    si->readMps(mpsFileName.c_str(),"mps");
    // Solve continuous problem
    si->initialSolve();
    // Save the orig lp/soco relaxation value for
    // comparisons later
    double origLpObj = si->getObjValue();
    // Instantiate cut generators
    CglConicMIR cg;
    //---------------------------------------------------
    // Keep applying cuts until
    //   1. no more cuts are generated
    // or
    //   2. the objective function value doesn't change
    //---------------------------------------------------
    bool equalObj;
    CoinRelFltEq eq(0.0001);
    double obj;
    do {
      // Get current solution value
      obj = si->getObjValue();
      // Generate and apply cuts
      cg.generateAndAddCuts(*si);
      // si->writeMps("after_cut");
      si->resolve();
      //si->initialSolve();
      cout <<endl;
      // -----------------------------------------------
      // Set Boolean flag to true if new objective is
      // almost equal to prior value.
      //
      // The test is:
      // abs(oldObj-newObj) <= 0.0001*(CoinMax(abs(oldObj),abs(newObj))+1.);
      // see CoinRelFloatEqual.h
      // -----------------------------------------------
      equalObj = eq(si->getObjValue(), obj);
      break;
    } while (!equalObj);
    double const * sol = si->getColSolution();
    // Print total number of cuts applied,
    // and total improvement in the lp objective value
    cout <<endl <<endl;
    cout << "----------------------------------------------------------"
         <<endl;
    cout << "Cut generation phase completed:" <<endl;
    cout << "   " << cg.getNumCutsAdded() << " many cuts added." << endl;
    cout << "   changing the lp objective value from " << origLpObj
         << " to " << si->getObjValue() <<endl;
    cout << "----------------------------------------------------------"
         <<endl;
    cout <<endl <<endl;
  }
  catch (CoinError e) {
    cout << e.className() << "::" << e.methodName() << " - " << e.message() << endl;
  }
  return 0;
}
