#ifndef ADMM_SOLVER_H
#define ADMM_SOLVER_H

#include "GAPBuilder.h"

#include <vector>

namespace admmsolver
{

struct GAPInstance;

struct MoveCandidate
{
  int cand_id = -1;
  int cell_id = -1;
  int bin_id = -1;
};

class ADMMSolver
{
  public:

    ADMMSolver(std::shared_ptr<gapbuilder::GAPInstance> instance);

    bool solve();

  private:

    float rho_;

    void initializeX(const int    num_cells, 
                     const float* disps,
                           int*   cell_to_best_cand);

    void updatePrimalX(const int    num_cells,
                       const float  rho, 
                       const float* widths,
                       const float* disps,
                       const float* capacities,
                       const float* y_cur,
                       const float* l_cur,
                             float* bin_slack,
                             int*   cell_to_best_cand); /* return vector */

    void updatePrimalY(const int    num_bins,
                       const float  rho, 
                       const float* bin_slack,
                       const float* l_cur,
                             float* y_next); /* return vector */

    void updateDual(const int    num_bins,
                    const float  rho, 
                    const float* bin_slack,
                    const float* y_next,
                    const float* l_cur,
                          float* l_next); /* return vector */

    float computeDisplacement(const int    num_cells,
                              const int*   cell_to_best_cand,
                              const float* disp);

    void computeBinSlack(const int    num_bins,
                         const int*   cell_to_best_cand,
                         const float* widths,
                         const float* capacities,
                               float* bin_slack); /* return vector */

    bool checkTermination();

    void updateNextIter(int iter);
    void setHyperParmeter();
    void computeFlattenInfo();

    std::shared_ptr<gapbuilder::GAPInstance> instance_;

    int num_cells_;
    int num_bins_;
    int num_candidates_;

    std::vector<float> disps_;
    std::vector<float> widths_;
    std::vector<float> capacities_;

    std::vector<MoveCandidate> candidates_;

    std::vector<float> bin_slack_;

    std::vector<float> vector_x_cur_;
    std::vector<float> vector_y_cur_;
    std::vector<float> vector_l_cur_; // l : lambda

    std::vector<float> vector_x_next_;
    std::vector<float> vector_y_next_;
    std::vector<float> vector_l_next_;

    std::vector<float> cell_to_demand_;

    std::vector<int> cand_id_to_cell_id_;
    std::vector<int> cand_id_to_bin_id_;

    std::vector<int> cell_to_best_cand_;

    std::vector<int> bin_ids_;
    std::vector<int> cell_ids_;

    std::vector<std::vector<int>> bin_id_to_cands_;
    std::vector<std::vector<int>> cell_id_to_cands_;
};

}

#endif
