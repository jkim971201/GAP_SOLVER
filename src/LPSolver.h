
#ifndef LP_SOLVER_H
#define LP_SOLVER_H

#include <vector>

namespace lpsolver
{

struct LPCandidate
{
  int item_id;
  int sack_id;
  int profit;
  int weight;

  LPCandidate(int _item, int _sack, int _profit, int _weight)
    : item_id(_item), sack_id(_sack), profit(_profit), weight(_weight)
  {}
};

class LPSolver
{
  public:

    LPSolver(const std::vector<std::vector<int>>& profits,
              const std::vector<std::vector<int>>& weights,
              const std::vector<int>& capacities);

    bool solve();
    std::vector<int> getResult() const;
    double getOptimalValue() const;

  private:

    int obj_value_;
    std::vector<int> solution_;

    std::vector<std::vector<int>> profits_;
    std::vector<std::vector<int>> weights_;
    std::vector<int> capacities_;
};

}

#endif
