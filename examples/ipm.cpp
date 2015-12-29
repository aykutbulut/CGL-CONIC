/* Implements a cutting plane solver. Uses CglConicIPM to cut.
 */

#include <CglConicIPM.hpp>
#include <ColaModel.hpp>

void add_nonnegativity(OsiConicSolverInterface * solver);

int main(int argc, char ** argv) {
  // read problem
  ColaModel * cola_solver = new ColaModel();
  // read problem
  cola_solver->readMps(argv[1]);
  // get cone information
  int num_cones = cola_solver->getNumCones();
  OsiLorentzConeType * cone_type = new OsiLorentzConeType[num_cones];
  int * cone_size = new int[num_cones];
  int ** cone_members = new int*[num_cones];
  for (int i=0; i<num_cones; ++i) {
    cola_solver->getConicConstraint(i, cone_type[i], cone_size[i],
				    cone_members[i]);
  }

  // add nonnegativity for leading variables
  add_nonnegativity(cola_solver);
  // solve problem ignoring conic constraints
  cola_solver->OsiClpSolverInterface::initialSolve();
  double const * sol = cola_solver->getColSolution();
  // create cut generator and generate cuts for the problem
  CglConicCutGenerator * cg_ipm = new CglConicIPM();
  OsiCuts * cuts;
  int total_num_cuts = 0;
  // solve problem while we can generate cuts.
  do {
    // ignore conic constraints and solve LP problem
    cola_solver->OsiClpSolverInterface::resolve();
    // generate cuts
    cuts = new OsiCuts();
    cg_ipm->generateCuts(*cola_solver, *cuts, num_cones, cone_type, cone_size,
			 cone_members);
    // add cuts to the problem
    int num_cuts = cuts->sizeRowCuts();
    if (num_cuts==0) {
      break;
    }
    else {
      std::cout << num_cuts << " ipm cuts produced." << std::endl;
    }
    total_num_cuts += num_cuts;
    cola_solver->OsiSolverInterface::applyCuts(*cuts);
    delete cuts;
    // print feasibility of conic constraints
    cola_solver->report_feasibility();
  } while(true);
  // print solution status
  cola_solver->report_feasibility();
  std::cout << total_num_cuts << " ipm cuts produced in total." << std::endl;
  std::cout << "Objective value " << cola_solver->getObjValue() << std::endl;
  delete cola_solver;
  delete[] cone_type;
  delete[] cone_size;
  for (int i=0; i<num_cones; ++i) {
    delete[] cone_members[i];
  }
  delete[] cone_members;
  delete cg_ipm;
  return 0;
}

void add_nonnegativity(OsiConicSolverInterface * solver) {
  int num_cones = solver->getNumCones();
  for (int i=0; i<num_cones; ++i) {
    OsiLorentzConeType type;
    int size;
    int * members;
    solver->getConicConstraint(i, type, size, members);
    if (type==OSI_QUAD) {
      solver->setColLower(members[0], 0.0);
    }
    else if (type==OSI_RQUAD) {
      solver->setColLower(members[0], 0.0);
      solver->setColLower(members[1], 0.0);
    }
    else {
      std::cerr << "Unknown cone type!" << std::endl;
      throw std::exception();
    }
    delete[] members;
  }
}
