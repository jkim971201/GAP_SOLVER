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
    float lambda_;

    void initializeX(const int    num_cells, 
                     const int*   cell_id_to_cand_id_start,
                     const float* disps,
                           float* vector_x);

    void updatePrimalX(const int    num_candidates,
                       const int    max_pgd_iter,
                       const float  rho, 
                       const float  lambda,
                       const float* widths,
                       const float* disps,
                       const float* capacities,
                       const float* x_cur,
                       const float* y_cur,
                       const float* u_cur,
                             float* x_next); /* return vector */

    void simplexProjection(const int    num_cells,
                           const int*   cell_id_to_num_cand,
                           const float* vector_input, 
                                 float* vector_sorted,
                                 float* vector_output); /* return vector */ 

    void updatePrimalY(const int    num_bins,
                       const float  rho, 
                       const float* bin_usages,
                       const float* u_cur,
                             float* y_next); /* return vector */

    void updateDual(const int    num_bins,
                    const float  rho, 
                    const float* bin_usages,
                    const float* y_next,
                    const float* u_cur,
                          float* u_next); /* return vector */

    float computeDisplacement(const int    num_candidates,
                              const float* disp,
                              const float* x_vector);

    float computeFractionalCost(const int    num_candidates,
                                const float  lambda, 
                                const float* x_vector);

    float computeOverflowCost(const int    num_bins,
                              const float  rho,
                              const float* bin_usage,
                              const float* vector_y,
                              const float* vector_u);

    void computeBinUsage(const int    num_bins,
                         const int*   bin_id_to_cand_id_start,
                         const int*   cand_id_to_cell_id,
                         const float* widths,
                         const float* capacities,
                         const float* vector_x,
                               float* bin_usages); /* return vector */

    void computeGradient(const int    num_candidates,
                         const float  rho,
                         const int*   cand_id_to_cell_id,
                         const int*   cand_id_to_bin_id,
                         const float* disp,
                         const float* widths,
                         const float* bin_usage,
                         const float* vector_y,
                         const float* vector_u,
                               float* grad); /* return vector */

		void makeIntegerSolution(const float* vector_in, float* vector_out); /* return vector */

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

    std::vector<float> bin_usage_;

    std::vector<float> vector_x_cur_;
    std::vector<float> vector_y_cur_;
    std::vector<float> vector_u_cur_;

    std::vector<float> vector_x_next_;
    std::vector<float> vector_y_next_;
    std::vector<float> vector_u_next_;

    std::vector<float> vector_grad_;

    std::vector<int> cand_id_to_cell_id_;
    std::vector<int> cand_id_to_bin_id_;

    std::vector<int> cell_id_to_num_cand_;
    std::vector<int> cell_id_to_cand_id_start_;
    std::vector<int> cand_id_foreach_cells_;

    std::vector<int> bin_id_to_num_cand_;
    std::vector<int> bin_id_to_cand_id_start_;
    std::vector<int> cand_id_foreach_bins_;
};

}

#endif
