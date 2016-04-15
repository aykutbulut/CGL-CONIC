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
#include <OsiMskSolverInterface.hpp>

int main(int argc, char ** argv) {
  // create conic solver interface
  ColaModel * conic_solver = new ColaModel();
  // read problem including conic constraints
  conic_solver->readMps(argv[1]);
  // solve initial problem ignoring conic constraints.
  conic_solver->OsiClpSolverInterface::initialSolve();
  CglConicOA cg;
  OsiCuts * cuts;
  int total_num_cuts = 0;
  clock_t start_time = clock();
  // solve problem while we can generate cuts.
  do {
    // ignore conic constraints and solve LP problem
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
    total_num_cuts += num_cuts;
    conic_solver->OsiSolverInterface::applyCuts(*cuts);
    delete cuts;
  } while(true);
  clock_t duration = clock() - start_time;
  // print solution status
  conic_solver->report_feasibility();
  std::cout << "Total number of cuts: " << total_num_cuts << std::endl;
  std::cout << "Objective value:      " << conic_solver->getObjValue() << std::endl;
  std::cout << "CPU time:             "
            << double(duration)/double(CLOCKS_PER_SEC) << std::endl;
  delete conic_solver;
  return 0;
}
