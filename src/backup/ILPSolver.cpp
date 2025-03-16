#include <vector>
#include <unordered_map>
#include "ILPSolver.h"

// Google Ortools
#include <ortools/linear_solver/linear_solver.h>

using operations_research::MPConstraint;
using operations_research::MPObjective;
using operations_research::MPSolver;
using operations_research::MPVariable;

namespace ilpsolver
{

ILPSolver::ILPSolver(std::vector<std::vector<int>>& profits,
                     std::vector<std::vector<int>>& weights,
                     std::vector<int>& capacities)
  : profits_   (profits),
    weights_   (weights),
    capacities_(capacities)
{}

bool
ILPSolver::solve()
{
  int num_items = weights_.size();
  int num_sacks = capacities_.size();

  // Initialize ILP Solver
  std::unique_ptr<MPSolver> ilp_solver(MPSolver::CreateSolver("SCIP"));
  auto set_num_thread_status = ilp_solver->SetNumThreads(1); 

  std::unordered_map<int, std::vector<ILPCandidate*>> item_id2candidates;
  std::unordered_map<int, std::vector<ILPCandidate*>> sack_id2candidates;
  std::unordered_map<ILPCandidate*, MPVariable*> x_table;

  // Make Variables
  for(int sack_id = 0; sack_id < num_sacks; sack_id++)
  {
    for(int item_id = 0; item_id < num_items; item_id++)
    {
      int profit = profits_[item_id][sack_id];
      int weight = weights_[item_id][sack_id];
      ILPCandidate* new_cand = new ILPCandidate(item_id, sack_id, profit, weight);
      MPVariable* const new_x = ilp_solver->MakeIntVar(0, 1, "");
      item_id2candidates[item_id].push_back(new_cand);
      sack_id2candidates[sack_id].push_back(new_cand);
      x_table[new_cand] = new_x;
    }
  }

  // Make Objective
  MPObjective* const objective = ilp_solver->MutableObjective();
  for(auto& [ilp_cand, x_var] : x_table)
    objective->SetCoefficient(x_var, ilp_cand->profit);

  objective->SetMinimization();

  // Make Candidate Selection Constraints (Must select only one row)
  for(auto& [item_id, candidates] : item_id2candidates)
  {
    MPConstraint* const new_select_constraint 
      = ilp_solver->MakeRowConstraint(1, 1);
    for(auto& candidate : candidates)
    {
      auto x_variable = x_table[candidate];
      new_select_constraint->SetCoefficient(x_variable, 1);
    }
  }

  // Make Constraint (lower than bin capacity (== width))
  for(auto& [sack_id, candidates] : sack_id2candidates)
  {
    int sack_capacity = capacities_[sack_id];

    MPConstraint* const new_overlap_constraint 
      = ilp_solver->MakeRowConstraint(0.0, sack_capacity);

    for(auto cand : candidates)
    {
      auto x_var_candi = x_table[cand];
      int weight = cand->weight;
      new_overlap_constraint->SetCoefficient(x_var_candi, weight);
    }
  }

  // Solve
  const MPSolver::ResultStatus result_status = ilp_solver->Solve();

  // When Solver fails, return false to increase bin_cap_margin.
  if(result_status != MPSolver::OPTIMAL)
    return false;

  std::vector<int> solution(num_items);
  // When Solver successes, commit solution to DB.
  for(auto& [candidate, x_var] : x_table)
  {
    if(x_var->solution_value() == 1)
    {
      int item_id = candidate->item_id;
      int sack_id = candidate->sack_id;
      int profit  = candidate->profit;
      solution[item_id] = sack_id;
			printf("move cell %d to bin %d\n", item_id, sack_id);
    }
  }

  //printf("NumConstraint : %d\n", ilp_solver->NumConstraints());
  //printf("NumVariable   : %d\n", ilp_solver->NumVariables());

  obj_value_ = objective->Value();
  solution_ = solution;
  return true;
}

int
ILPSolver::getOptimalValue() const
{
  return obj_value_;
}

std::vector<int>
ILPSolver::getResult() const
{
  return solution_;
}

}
