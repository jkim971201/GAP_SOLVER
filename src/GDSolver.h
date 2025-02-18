#ifndef GD_SOLVER_H
#define GD_SOLVER_H

#include "GAPBuilder.h"

#include <vector>

namespace gdsolver
{

struct GAPInstance;

struct ILPCandidate
{
  int cand_id = -1;
  int cell_id = -1;
  int bin_id = -1;
  float val_y = 0.0; // probability
  float val_x = 0.0; // values to be optimized
};

class GDSolver
{
  public:

    GDSolver(std::shared_ptr<gapbuilder::GAPInstance> instance);

    bool solve();

  private:

    float softmax_tmpr_;

    void updateBinPenalty();
    void computeBinSlack();
    void computeSubGradient();
    void computeFlattenInfo();
		void computeSumExpX();

    std::shared_ptr<gapbuilder::GAPInstance> instance_;

    int num_cells_;
    int num_bins_;
    int num_candidates_;

    std::vector<float> disps_;
    std::vector<float> widths_;
    std::vector<float> capacities_;

    std::vector<ILPCandidate> candidates_;

    std::vector<float> bin_slack_;
    std::vector<float> penalty_;

    std::vector<float> vector_x_;
    std::vector<float> vector_y_;

    std::vector<float> vector_sum_exp_;
		std::vector<float> vector_min_x_;

    std::vector<int> num_cands_each_cell_;
    std::vector<int> cell_id_to_cand_start_;
    std::vector<int> cands_in_cell_;
};

}

#endif
