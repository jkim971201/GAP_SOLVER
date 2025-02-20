#include <cstdio>
#include <chrono>
#include <vector>
#include <memory>
#include <cassert>

#include "GDSolver.h"
#include "GAPBuilder.h"
#include "ILPSolver.h"

using namespace gapbuilder;
using namespace ilpsolver;
using namespace gdsolver;

int main(int argc, char** argv)
{
  if(argc < 2) 
  {
    printf("No input file\n");
    exit(1);
  }

  const std::string input_file_name = argv[1];

  GAPBuilder gap_builder(input_file_name);
  gap_builder.print();

  auto gap_instance = gap_builder.getGAPInstance();

  auto t1 = std::chrono::high_resolution_clock::now();

  std::unique_ptr<ILPSolver> solver_ilp 
    = std::make_unique<ILPSolver>(gap_instance->profits,
                                  gap_instance->weights,
                                  gap_instance->capacities);

  bool ilp_success = solver_ilp->solve();

  assert(ilp_success == true);

  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> runtime_ilp = t2 - t1;

  printf("ILP Solver Results\n");
  printf("Objective Value : %d\n", solver_ilp->getOptimalValue());
  printf("Runtime : %f\n", runtime_ilp.count());

  auto t3 = std::chrono::high_resolution_clock::now();

  std::unique_ptr<GDSolver> solver_gd
    = std::make_unique<GDSolver>(gap_instance);

  bool gd_success = solver_gd->solve();

  auto t4 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> runtime_gd = t4 - t3;

  printf("GD Solver Results\n");
  printf("Runtime : %f\n", runtime_gd.count());

  return 0;
}
