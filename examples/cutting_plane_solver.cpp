// Name:     CglConicOA.cpp
// Author:   Aykut Bulut
//           Lehigh University
//           email: aykut@lehigh.edu, aykutblt@gmail.com
//-----------------------------------------------------------------------------
// Copyright (C) 2015, Lehigh University.  All Rights Reserved.

//-----------------------------------------------------
// Implements simple cutting plane solver for SOCO problems.
// Makes use of CglConicOA cut library.
//
// usage:
//   cutting_plane_solver mpsFileName
// example:
//   cutting_plane_solver ../../Data/Sample/p0033
//-----------------------------------------------------

#include <CglConicOA.hpp>
#include <ColaModel.hpp>

int main(int argc, char ** argv) {
  // create conic solver interface
  ColaModel * conic_solver = new ColaModel();
  // read problem including conic constraints
  conic_solver->readMps(argv[1]);
  // solve initial problem ignoring conic constraints.
  conic_solver->OsiClpSolverInterface::initialSolve();
  CglConicOA cg;
  OsiCuts * cuts;
  // solve problem while we can generate cuts.
  do {
    conic_solver->OsiClpSolverInterface::resolve();
    // generate cuts
    cuts = new OsiCuts();
    cg.generateCuts(*conic_solver, *cuts);
    // add cuts to the problem
    int num_cuts = cuts->sizeRowCuts();
    if (num_cuts==0) {
      break;
    }
    else {
      std::cout << num_cuts << " many cuts produced." << std::endl;
    }
    for (int i=0; i<num_cuts; ++i) {
      OsiRowCut & cut = cuts->rowCut(i);
      int len = cut.row().getNumElements();
      if (len > 0) {
        // Add cut to solver
      }
      else {
        std::cerr << "Problem with cut size." << std::endl;
        throw std::exception();
      }
    }
    delete cuts;
  } while(true);
  delete conic_solver;
  return 0;
}
