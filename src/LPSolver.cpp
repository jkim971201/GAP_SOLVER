#include <vector>
#include <unordered_map>
#include "LPSolver.h"

// Google Ortools
#include <ortools/linear_solver/linear_solver.h>

using operations_research::MPConstraint;
using operations_research::MPObjective;
using operations_research::MPSolver;
using operations_research::MPVariable;

namespace lpsolver
{

LPSolver::LPSolver(const std::vector<std::vector<int>>& profits,
                   const std::vector<std::vector<int>>& weights,
                   const std::vector<int>& capacities)
  : profits_   (profits),
    weights_   (weights),
    capacities_(capacities)
{}

bool
LPSolver::solve()
{
  int num_items = weights_.size();
  int num_sacks = capacities_.size();

  // Initialize LP Solver
  std::unique_ptr<MPSolver> lp_solver(MPSolver::CreateSolver("GLOP"));
  std::unordered_map<int, std::vector<LPCandidate*>> item_id2candidates;
  std::unordered_map<int, std::vector<LPCandidate*>> sack_id2candidates;
  std::unordered_map<LPCandidate*, MPVariable*> x_table;

  // Make Variables
  const double infinity = lp_solver->infinity();
  for(int sack_id = 0; sack_id < num_sacks; sack_id++)
  {
    for(int item_id = 0; item_id < num_items; item_id++)
    {
      int profit = profits_[item_id][sack_id];
      int weight = weights_[item_id][sack_id];
      LPCandidate* new_cand = new LPCandidate(item_id, sack_id, profit, weight);
      //MPVariable* const new_x = lp_solver->MakeNumVar(0.0, infinity, "");
      MPVariable* const new_x = lp_solver->MakeNumVar(0.0, 1.0, "");
      item_id2candidates[item_id].push_back(new_cand);
      sack_id2candidates[sack_id].push_back(new_cand);
      x_table[new_cand] = new_x;
    }
  }

  // Make Objective
  MPObjective* const objective = lp_solver->MutableObjective();
  for(auto& [ilp_cand, x_var] : x_table)
    objective->SetCoefficient(x_var, ilp_cand->profit);

  objective->SetMinimization();

  // Make Candidate Selection Constraints (Must select only one row)
  for(auto& [item_id, candidates] : item_id2candidates)
  {
    MPConstraint* const new_select_constraint 
      = lp_solver->MakeRowConstraint(1, 1);
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
      = lp_solver->MakeRowConstraint(0.0, sack_capacity);

    for(auto cand : candidates)
    {
      auto x_var_candi = x_table[cand];
      int weight = cand->weight;
      new_overlap_constraint->SetCoefficient(x_var_candi, weight);
    }
  }

  // Solve
  const MPSolver::ResultStatus result_status = lp_solver->Solve();

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
    }
  }

  obj_value_ = objective->Value();
  solution_ = solution;
  return true;
}

double
LPSolver::getOptimalValue() const
{
  return obj_value_;
}

std::vector<int>
LPSolver::getResult() const
{
  return solution_;
}

}
