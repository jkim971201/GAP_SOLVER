#include <algorithm>
#include <numeric>
#include <cmath>

#include "ADMMSolver.h"

namespace admmsolver
{

ADMMSolver::ADMMSolver(std::shared_ptr<gapbuilder::GAPInstance> instance)
{
  setHyperParmeter();

  instance_ = instance;

  num_cells_ = instance_->num_items;
  num_bins_  = instance_->num_sacks;

  const auto& profits_2d = instance_->profits;
  const auto& weights_2d = instance_->weights;
  const auto& capacities_int = instance_->capacities;

  disps_.resize(num_cells_ * num_bins_);
  widths_.resize(num_cells_);
  capacities_.resize(num_bins_);

  std::transform(capacities_int.begin(), 
                 capacities_int.end(),
                 capacities_.begin(),
                 [] (int cap) { return static_cast<float>(cap); });

  int num_cand_made = 0;
  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
  {
    for(int cell_id = 0; cell_id < num_cells_; cell_id++)
    {
      int cand_id = num_cand_made;
      MoveCandidate new_cand = {cand_id, cell_id, bin_id};
      candidates_.push_back(new_cand);

      disps_[cand_id] = static_cast<float>(profits_2d.at(cell_id).at(bin_id));
      widths_[cell_id] = static_cast<float>(weights_2d.at(cell_id).at(bin_id));

      num_cand_made++;
    }
  }

  num_candidates_ = static_cast<int>(candidates_.size());

  bin_usage_.resize(num_bins_, 0.0);

  vector_x_cur_.resize(num_candidates_, 0.0);
  vector_y_cur_.resize(num_candidates_, 0.0);
  vector_u_cur_.resize(num_candidates_, 0.0);

  vector_x_next_.resize(num_candidates_, 0.0);
  vector_y_next_.resize(num_candidates_, 0.0);
  vector_u_next_.resize(num_candidates_, 0.0);

  vector_grad_.resize(num_candidates_, 0.0);

  cand_id_to_cell_id_.resize(num_candidates_, 0);
  cand_id_to_bin_id_.resize(num_candidates_, 0);

  cell_id_to_num_cand_.resize(num_cells_, 0);
  cell_id_to_cand_id_start_.resize(num_cells_ + 1, 0);
  cand_id_foreach_cells_.resize(num_candidates_, 0);

  bin_id_to_num_cand_.resize(num_bins_, 0);
  bin_id_to_cand_id_start_.resize(num_bins_ + 1, 0);
  cand_id_foreach_bins_.resize(num_candidates_, 0);

  computeFlattenInfo();
}

void
ADMMSolver::setHyperParmeter()
{
}

void
ADMMSolver::computeFlattenInfo()
{
  // 1. Cell to candidates
  std::vector<std::vector<int>> cell2candidates(num_cells_);
  for(int cand_id = 0; cand_id < num_candidates_; cand_id++)
  {
    auto& candidate = candidates_.at(cand_id);
    cand_id_to_cell_id_[cand_id] = candidate.cell_id;
    cand_id_to_bin_id_[cand_id] = candidate.bin_id;
    cell2candidates[candidate.cell_id].push_back(cand_id);
  }

  int flatten_idx_cell = 0;
  for(int cell_id = 0; cell_id < num_cells_; cell_id++)
  {
    cell_id_to_cand_id_start_[cell_id] = flatten_idx_cell;

    const auto& cands = cell2candidates.at(cell_id);
    int num_cand_this_cell = static_cast<int>(cands.size());
    cell_id_to_num_cand_[cell_id] = num_cand_this_cell;

    for(int cand_id : cands)
    {
      cand_id_foreach_cells_[flatten_idx_cell] = cand_id;
      flatten_idx_cell++;
      // printf("cell_id : %d cand_id : %d\n", cell_id, cand_id);
    }
  }

  cell_id_to_cand_id_start_.back() = num_candidates_;

  // 2. Bin to candidates
  std::vector<std::vector<int>> bin2candidates(num_bins_);
  for(int cand_id = 0; cand_id < num_candidates_; cand_id++)
  {
    auto& candidate = candidates_.at(cand_id);
    bin2candidates[candidate.bin_id].push_back(cand_id);
  }

  int flatten_idx_bin = 0;
  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
  {
    bin_id_to_cand_id_start_[bin_id] = flatten_idx_bin;

    const auto& cands = bin2candidates.at(bin_id);
    int num_cand_this_bin = static_cast<int>(cands.size());
    bin_id_to_num_cand_[bin_id] = num_cand_this_bin;

    for(int cand_id : cands)
    {
      cand_id_foreach_bins_[flatten_idx_bin] = cand_id;
      flatten_idx_bin++;
      // printf("cell_id : %d cand_id : %d\n", cell_id, cand_id);
    }
  }

  bin_id_to_cand_id_start_.back() = num_candidates_;
}

void
ADMMSolver::updateParameter()
{
}

float
ADMMSolver::computeDisplacement(const int    num_candidates,
                                const float* disp,
                                const float* x_val)
{
  float sum_disp = 0.0;
  for(int cand_id = 0; cand_id < num_candidates; cand_id++)
    sum_disp += disp[cand_id] * x_val[cand_id];

  return sum_disp;
}

void
ADMMSolver::computeBinUsage(const int    num_bins,
                            const int*   bin_id_to_cand_id_start,
                            const int*   cand_id_to_cell_id,
                            const float* widths,
                            const float* capacities,
                            const float* vector_x,
                                  float* bin_usages) /* return vector */
{
  for(int bin_id = 0; bin_id < num_bins; bin_id++)
  {
    float bin_usage = 0.0;
    float bin_capacity = capacities[bin_id]; 
    
    int cand_id_start = bin_id_to_cand_id_start[bin_id];
    int cand_id_end = bin_id_to_cand_id_start[bin_id + 1];
    for(int cand_id = cand_id_start; cand_id < cand_id_end; cand_id++)
    {
      int cell_id = cand_id_to_cell_id_.at(cand_id);
      float cell_width = widths_.at(cell_id);
      float x_val = vector_x[cand_id];

      bin_usage += cell_width * x_val;
    }

    bin_usages[bin_id] = bin_usage - bin_capacity;
  }
}

void
ADMMSolver::updatePrimalX(const int    num_candidates,
                          const float  rho, 
                          const float  lambda,
                          const float* widths,
                          const float* disps,
                          const float* capacities,
                          const float* x_cur,
                          const float* y_cur,
                          const float* u_cur,
                                float* x_next)
{
  const int max_pgd_iter = 200;
  const float step_size = 7.8e-9;

  std::vecotr<float> x_vector_temp(num_candidates, 0.0);
  std::vector<float> v_vector(num_candidates, 0.0);
  std::vector<float> x_vector_k_minus_1(num_candidates, 0.0);
  std::vector<float> x_vector_k_minus_2(num_candidates, 0.0);

  for(int cand_id = 0; cand_id < num_candidates_; cand_id++)
  {
    x_vector_k_minus_1.at(cand_id) = x_cur[cand_id];
    x_vector_k_minus_2.at(cand_id) = x_cur[cand_id];
  }

  for(int pgd_iter = 0; pgd_iter < max_pgd_iter; pgd_iter++)
  {
    // 1. update v_vector
    float pgd_iter_float = static_cast<float>(pgd_iter);
    for(int v_idx = 0; v_idx < num_candidates; v_idx++)
    {
      float x_k_minus_1 = x_vector_k_minus_1.at(v_idx);
      float x_k_minus_2 = x_vector_k_minus_2.at(v_idx);

      v_vector[v_idx] 
        = (2 * pgd_iter_float - 1) / (pgd_iter_float + 1) * x_k_minus_1
        + (pgd_iter_float - 2) / (pgd_iter_float + 1) * x_k_minus_2;
    }

    // 2. compute gradient
    for(int grad_idx = 0; grad_idx < num_candidates; grad_idx++)
    {

    }

    // 3. simplex projection
  }
}

void
ADMMSolver::computeGradient(const int    num_candidates,
                            const float  rho,
                            const float* disp,
                            const float* widths,
                            const float* bin_usage,
                            const float* vector_y,
                            const float* vector_u,
                                  float* grad)
{

}

void bubbleSort(float* arr, int n)
{
  int i,j;
  float temp;
  for(j = n; j > 0; j--)
  {
    for(i = 0; i < j; i++)
    {
      if(arr[i] < arr[i + 1])
      {
        // swap
        temp = arr[i + 1];
        arr[i + 1] = arr[i];
        arr[i] = temp; 
      }
    }
  }
}

void
ADMMSolver::simplexProjection(float* vector_to_project,
                              float* vector_projection_workspace)
{
  for(int cell_id = 0; cell_id < num_cells_; cell_id++)
  {
    int cand_id_start = cell_id_to_cand_id_start_.at(cell_id);
    int num_candidates_this_cell = cell_id_to_num_cand_.at(cell_id);
    float* vector_this_cell = vector_to_project + 
  }
}

bool
ADMMSolver::solve()
{
  bool admm_success = false;

  int max_admm_iter = 200;
  int admm_iter = 0;

  while(admm_iter++ < max_admm_iter)
  {
    int max_pgd_iter = 200;

    /* ========================= Main ADMM Iteration ========================== */
    /* 1. Solve subproblem w.r.t. x_k (Projected Gradient) */
    updatePrimalX(num_candidates_,
                  rho_,
                  lambda_,
                  max_pgd_iter,
                  widths_.data(), 
                  disps_.data(), 
                  capacities_.data(), 
                  vector_x_cur_.data(), 
                  vector_y_cur_.data(),
                  vector_u_cur_.data(),
                  vector_x_next_.data()); /* return new x_vector */

    /* 2. Solve subproblem w.r.t. y_k */
    updatePrimalY(rho_, 
                  widths_.data(), 
                  capacities_.data(), 
                  vector_x_next_.data(), 
                  vector_y_cur_.data(), 
                  vector_u_cur_.data(),
                  vector_y_next_.data()); /* return new y_vector */

    /* 3. Update Dual Variable u */
    updateDual(rho_, 
               widths_.data(),
               capacities_.data(), 
               vector_x_next_.data(), 
               vector_y_next_().data(), 
               vector_u_cur_.data(),
               vector_u_next_.data()); /* return new u_vector */
  
    updateNextIter(admm_iter);
    /* ======================================================================== */
  }
  
  return admm_success;
}

}
