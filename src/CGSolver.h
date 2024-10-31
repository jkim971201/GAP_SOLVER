#ifndef CG_SOLVER_H
#define CG_SOLVER_H

#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <string>

namespace cgsolver
{

// Variable of Restriced Master Problem
struct RMPVar
{
  int cost;
  bool is_artifical = false;
  int sack_id = -1;
  std::vector<int> items;
  std::size_t hash_val;
  RMPVar(int _cost) : cost(_cost) {}
};

struct RMPConstraint
{
  std::vector<RMPVar*> variables;
};

// CG := Column Generation
class CGSolver
{
  public:

    CGSolver(const std::vector<std::vector<int>>& profits,
             const std::vector<std::vector<int>>& weights,
             const std::vector<int>& capacities);

    bool solve();
    std::vector<int> getResult() const;
    double getOptimalValue() const;

  private:

    int getIndex(int item_id, int sack_id) const;
    void getColRow(int flatten_id, int& col, int& row) const;

    void makeInitialRMP();

    // Solve Restricted Master Problem
    bool solveRMP();
    bool solveKP(std::vector<double>& dual_val); 
    // vector<double> dual_val : item_id -> dual value

    int getProfitSum(int sack_id, const std::vector<int>& kp_sol) const;
    bool checkIfFeasible(int sack_id, const std::vector<int>& knapsack_solution) const;

    int num_items_;
    int num_sacks_;
    double rmp_value_;
    std::vector<int> solution_;
    std::vector<int> profits_;
    std::vector<int> weights_;
    std::vector<int> capacities_;

    std::map<int, std::vector<RMPVar*>> item2mp_vars_;
    std::map<int, std::vector<RMPVar*>> sack2mp_vars_;
    std::unordered_map<std::size_t, RMPVar*> rmp_var_table_;

    // this will be updated while cg iteration
    std::vector<RMPVar*>        rmp_vars_;
    std::vector<RMPConstraint*> rmp_selection_constraints_;
    std::vector<RMPConstraint*> rmp_convex_constraints_;
};

}

#endif
