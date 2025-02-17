#ifndef ILP_SOLVER_H
#define ILP_SOLVER_H

#include <vector>

namespace ilpsolver
{

struct ILPCandidate
{
  int item_id;
  int sack_id;
  int profit;
  int weight;

  ILPCandidate(int _item, int _sack, int _profit, int _weight)
    : item_id(_item), sack_id(_sack), profit(_profit), weight(_weight)
  {}
};

class ILPSolver
{
  public:

    ILPSolver(std::vector<std::vector<int>>& profits,
              std::vector<std::vector<int>>& weights,
              std::vector<int>& capacities);

    bool solve();
    std::vector<int> getResult() const;
    int getOptimalValue() const;

  private:

    int obj_value_;
    std::vector<int> solution_;

    std::vector<std::vector<int>> profits_;
    std::vector<std::vector<int>> weights_;
    std::vector<int> capacities_;
};

}

#endif
