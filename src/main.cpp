#include <cstdio>
#include <chrono>
#include <vector>
#include <memory>
#include <cassert>

#include "LPSolver.h"
#include "CGSolver.h"

using namespace lpsolver;
using namespace cgsolver;

int main() 
{
  std::vector<std::vector<int>> profits = 
  {{24, 18},
   {16, 21},
   {18, 14},
   {10, 12},
   {17, 26},
   {21, 18}};

  std::vector<std::vector<int>> weights = 
  {{18, 20},
   {21, 16},
   {14,  9},
   {19, 17},
   {17, 12},
   {10, 19}};

  std::vector<int> capacities = {48, 43};

//  std::vector<std::vector<int>> profits = 
//  { {10,  6},
//    { 7,  8},
//    { 5, 11} };
//
//  std::vector<std::vector<int>> weights = 
//  { { 9,  5},
//    { 6,  7},
//    { 3,  9} };
//
//  std::vector<int> capacities = {11, 18};

//  std::vector<std::vector<int>> profits = 
//  { {8, 1},
//    {3, 7},
//    {2, 5},
//    {9, 2} };
//
//  std::vector<std::vector<int>> weights = 
//  { {2, 5},
//    {3, 1},
//    {3, 1},
//    {1, 3} };
//
//  std::vector<int> capacities = {5, 8};

  auto t1 = std::chrono::high_resolution_clock::now();
  std::unique_ptr<LPSolver> solver_ilp = std::make_unique<LPSolver>(profits, weights, capacities);
  bool ilp_success = solver_ilp->solve();
  assert(ilp_success == true);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> runtime_ilp = t2 - t1;

  printf("LP Solver Results\n");
  printf("Objective Value : %f\n", solver_ilp->getOptimalValue());
  printf("Runtime : %f\n", runtime_ilp.count());
  printf("\n");

  auto t3 = std::chrono::high_resolution_clock::now();
  std::unique_ptr<CGSolver> solver_cg = std::make_unique<CGSolver>(profits, weights, capacities);
  bool pd_success = solver_cg->solve();
  auto t4 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> runtime_cg = t4 - t3;

  printf("CG Solver Results\n");
  printf("Objective Value : %f\n", solver_cg->getOptimalValue());
  printf("Runtime : %f\n", runtime_cg.count());

  return 0;
}
