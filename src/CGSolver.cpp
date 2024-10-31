#include <vector>
#include <limits>
#include <algorithm>
#include <cassert>
#include <unordered_map>
#include "CGSolver.h"

// Google Ortools
#include <ortools/linear_solver/linear_solver.h>
#include "ortools/algorithms/knapsack_solver.h"

using operations_research::MPConstraint;
using operations_research::MPObjective;
using operations_research::MPSolver;
using operations_research::MPVariable;
using operations_research::KnapsackSolver;

namespace cgsolver
{

int
CGSolver::getIndex(int item_id, int sack_id) const
{
  return item_id + sack_id * num_items_;
}

void
CGSolver::getColRow(int flatten_id, int& col, int& row) const
{
  row = flatten_id / num_items_;
  col = flatten_id % num_items_;
}

CGSolver::CGSolver(const std::vector<std::vector<int>>& profits,
                   const std::vector<std::vector<int>>& weights,
                   const std::vector<int>& capacities)
{
  num_items_ = static_cast<int>(weights.size());
  num_sacks_ = static_cast<int>(capacities.size());

  //printf("NumItem : %d\n", num_items_);
  //printf("NumSack : %d\n", num_sacks_);

  int matrix_size = num_items_ * num_sacks_;
  profits_.resize(matrix_size);
  weights_.resize(matrix_size);

  for(int i = 0; i < num_items_; i++)
  {
    for(int j = 0; j < num_sacks_; j++)
    {
      int profit = profits[i][j];
      int weight = weights[i][j];
      int flatten_id = getIndex(i, j);
      profits_[flatten_id] = profit;
      weights_[flatten_id] = weight;
    }
  }

  capacities_ = capacities;
}

void
CGSolver::makeInitialRMP()
{
  rmp_selection_constraints_.clear();
  rmp_convex_constraints_.clear();

  rmp_selection_constraints_.resize(num_items_);
  rmp_convex_constraints_.resize(num_sacks_);

  // Big-M Method
  // Add Artificial variable to obtain initial basic basis.
  constexpr int BIG_M = 99;
  for(int item_id = 0; item_id < num_items_; item_id++)
  {
		RMPConstraint* new_select_constraint = new RMPConstraint;
    RMPVar* new_artificial_var = new RMPVar(BIG_M);
    new_artificial_var->is_artifical = true;
		new_select_constraint->variables.push_back(new_artificial_var);
    rmp_selection_constraints_[item_id] = new_select_constraint;
    rmp_vars_.push_back(new_artificial_var);
  }

	// Initialize Convex Constraint
  for(int sack_id = 0; sack_id < num_sacks_; sack_id++)
  {
		RMPConstraint* new_convex_constraint = new RMPConstraint;
    rmp_convex_constraints_[sack_id] = new_convex_constraint;
  }
}

int
CGSolver::getProfitSum(int sack_id, 
                       const std::vector<int>& knapsack_solution) const
{
  int sum_profit = 0;
  for(auto item_id : knapsack_solution)
  {
    int flatten_id = getIndex(item_id, sack_id);
    int profit = profits_[flatten_id];
    sum_profit += profit;
  }
  return sum_profit;
}

bool
CGSolver::checkIfFeasible(int sack_id, 
                          const std::vector<int>& knapsack_solution) const
{
	int capacity = capacities_.at(sack_id);
  int sum_weight = 0;
	bool is_feasible = true;
  for(auto item_id : knapsack_solution)
  {
    int flatten_id = getIndex(item_id, sack_id);
    int weight = weights_[flatten_id];
    sum_weight += weight;

		if(sum_weight > capacity)
		{
			is_feasible = false;
			break;
		}
  }

  return is_feasible;
}

bool
CGSolver::solve() 
{
  makeInitialRMP();

	constexpr int max_iter = 1000;
	int column_generation_iter = 0;
	bool rmp_updated = true;
	while(rmp_updated == true)
	{
    rmp_updated = solveRMP();
		column_generation_iter += 1;

		if(column_generation_iter > max_iter)
			break;
	}

  return true;
}

// Solve Restricted Master Problem
bool 
CGSolver::solveRMP()
{
  //printf("Solve Restricted Master Problem\n");
  // Initialize LP Solver
  std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GLOP"));
  std::unordered_map<RMPVar*, MPVariable*> rmp2lp;
  std::vector<MPConstraint*> selection_lp_constraints;

  const double infinity = solver->infinity();
  // Make Variables
  for(auto rmp_var : rmp_vars_)
  {
    MPVariable* new_var = solver->MakeNumVar(0.0, infinity, "");
    rmp2lp[rmp_var] = new_var;
  }

  // Make Objective
  MPObjective* const objective = solver->MutableObjective();
  for(auto& [rmp_var, lp_var] : rmp2lp)
    objective->SetCoefficient(lp_var, rmp_var->cost);

  objective->SetMinimization();

  // Make Constraints (1)
  for(auto constraint : rmp_selection_constraints_) 
  {
    MPConstraint* const new_constraint 
      = solver->MakeRowConstraint(1.0, 1.0);

    for(auto rmp_var : constraint->variables)
    {
      MPVariable* lp_var = rmp2lp.at(rmp_var);
      new_constraint->SetCoefficient(lp_var, 1.0);
    }

    selection_lp_constraints.push_back(new_constraint);
  }

  // Make Constraints (2)
  for(auto constraint : rmp_convex_constraints_) 
  {
    MPConstraint* const new_constraint 
      = solver->MakeRowConstraint(0.0, 1.0);

    for(auto rmp_var : constraint->variables)
    {
      MPVariable* lp_var = rmp2lp.at(rmp_var);
      new_constraint->SetCoefficient(lp_var, 1.0);
    }
  }

  // Solve
  const MPSolver::ResultStatus result_status = solver->Solve();
  if(result_status != MPSolver::OPTIMAL)
    assert(0);

  //printf("RMP Value : %f\n", objective->Value());
	rmp_value_ = objective->Value();

  std::vector<double> dual_values(num_items_);
  for(int item_id = 0; item_id < num_items_; item_id++)
  {
    MPConstraint* lp_constraint = selection_lp_constraints[item_id];
    double dual_val = lp_constraint->dual_value();
    dual_values[item_id] = dual_val;
    //printf("Dual Value : %f\n", dual_val);
  }

  // kp_add_new_var : if knapsack finds new basic variable
  bool kp_add_new_var = solveKP(dual_values);
  return kp_add_new_var;
}

bool 
CGSolver::solveKP(std::vector<double>& dual_value)
{
  assert(num_items_ == dual_value.size());

  bool rmp_updated = false;
  for(int sack_id = 0; sack_id < num_sacks_; sack_id++)
  {
    // printf("Solve Knapsack for sack %d\n", sack_id);
    int capacity = capacities_[sack_id];

    std::vector<std::vector<int64_t>> multi_weights(1);
    std::vector<int64_t>& weights = multi_weights[0];
    std::vector<int64_t> capacities(1, capacity);
    std::vector<int64_t> values;

    std::vector<double> dual_costs;
    for(int item_id = 0; item_id < num_items_; item_id++)
    {
      int flatten_id = getIndex(item_id, sack_id);
      int profit = profits_[flatten_id];
      int weight = weights_[flatten_id];
      double dual_cost = -profit + dual_value[item_id];
      dual_costs.push_back(dual_cost);
      weights.push_back(weight);
    }

    // This solver maximizes the sum of profit
    KnapsackSolver solver(KnapsackSolver::KNAPSACK_DYNAMIC_PROGRAMMING_SOLVER, "");

    for(auto cost : dual_costs)
      values.push_back(static_cast<int64_t>(cost));

    solver.Init(values, multi_weights, capacities);
    int64_t computed_value = solver.Solve();

    // printf("Knapsack Value : %ld\n", computed_value);

    std::vector<int> packed_items;
    for(std::size_t i = 0; i < values.size(); ++i) 
    {
      if(solver.BestSolutionContains(i)) 
        packed_items.push_back(i);
    }

    bool is_feasible = checkIfFeasible(sack_id, packed_items);
    if(is_feasible == false)
      continue;

    std::string string_for_hash = std::to_string(sack_id);
    std::vector<int> packed_items_flag(num_items_, 0);
    for(int item_id : packed_items)
      packed_items_flag[item_id] = 1;

    for(int i = 0; i < num_items_; i++)
    {
      if(packed_items_flag[i] == 1)
        string_for_hash += "1";
      else
        string_for_hash += "0";
    }

    std::size_t rmp_var_hash = std::hash<std::string>()(string_for_hash);

    auto rmp_var_itr = rmp_var_table_.find(rmp_var_hash);
    if(rmp_var_itr == rmp_var_table_.end())
    {
      rmp_updated = true;
			int sum_profit = getProfitSum(sack_id, packed_items);
      RMPVar* new_rmp_var = new RMPVar(sum_profit);
      new_rmp_var->hash_val = rmp_var_hash;
      rmp_var_table_[rmp_var_hash] = new_rmp_var;
      rmp_vars_.push_back(new_rmp_var);

      // 1. Add this variable to selection constraints
      for(auto item_id : packed_items)
      {
        // Selection Constraint for this item
        RMPConstraint* sel_constr = rmp_selection_constraints_.at(item_id);
        sel_constr->variables.push_back(new_rmp_var);
      }

      // 2. Add this variable to convex constraints
      // Convex Constraint for this item
      RMPConstraint* cvx_constr = rmp_convex_constraints_.at(sack_id);
      cvx_constr->variables.push_back(new_rmp_var);
    }
  }

  // return new variable is added to restricted master problem
  return rmp_updated;
}

double
CGSolver::getOptimalValue() const
{
  return rmp_value_;
}

std::vector<int>
CGSolver::getResult() const
{
  return solution_;
}

}
