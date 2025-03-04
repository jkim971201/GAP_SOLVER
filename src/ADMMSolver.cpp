#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <numeric>

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
  for(int cell_id = 0; cell_id < num_cells_; cell_id++)
  {
    for(int bin_id = 0; bin_id < num_bins_; bin_id++)
    {
      int cand_id = num_cand_made;
      MoveCandidate new_cand = {cand_id, cell_id, bin_id};
      // printf("cand_id : %d cell_id : %d bin_id : %d\n", cand_id, cell_id, bin_id);
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

  initializeX(num_cells_, 
              cell_id_to_cand_id_start_.data(), 
              disps_.data(), 
              vector_x_cur_.data());
}

void
ADMMSolver::initializeX(const int    num_cells,
                        const int*   cell_id_to_cand_id_start,
                        const float* disps,
                              float* vector_x)
{
  for(int cell_id = 0; cell_id < num_cells; cell_id++)
  {
    int min_cand_id = -1;
    float min_disp = std::numeric_limits<float>::max();
    const int cand_id_start = cell_id_to_cand_id_start[cell_id];
    const int cand_id_end = cell_id_to_cand_id_start[cell_id + 1];
  
    for(int index = cand_id_start; index < cand_id_end; index++)
    {
      int cand_id = cand_id_foreach_cells_.at(index);
      int disp_this_cand = disps[cand_id];
      if(disp_this_cand < min_disp)
      {
        min_disp = disp_this_cand;
        min_cand_id = cand_id;
      }
    }
    
    // printf("Cell : %d Min Disp : %f Min Index : %d\n", cell_id, min_disp, min_cand_id);
    vector_x[min_cand_id] = 1.0;
  }

  float initial_disp = computeDisplacement(num_candidates_, disps_.data(), vector_x);
  printf("Initial Displacement : %f\n", initial_disp);

  computeBinUsage(num_bins_,
                  bin_id_to_cand_id_start_.data(),
                  cand_id_to_cell_id_.data(),
                  widths_.data(),
                  capacities_.data(),
                  vector_x,
                  bin_usage_.data());

  printf("Bin Usage\n");
  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
    printf("  %f\n", bin_usage_.at(bin_id) + capacities_.at(bin_id));

  float initial_ovf 
    = computeOverflowCost(num_bins_, 
                          rho_, 
                          bin_usage_.data(), 
                          vector_y_cur_.data(), 
                          vector_u_cur_.data());
  printf("Initial OverflowCost : %f\n", initial_ovf);
}

void
ADMMSolver::updateNextIter(int iter)
{
  std::swap(vector_x_cur_, vector_x_next_);
  std::swap(vector_y_cur_, vector_y_next_);
  std::swap(vector_u_cur_, vector_u_next_);

  lambda_ = lambda_ * 1.05;

  float disp_cost 
    = computeDisplacement(num_candidates_, disps_.data(), vector_x_cur_.data());

  float ovf_cost
    = computeOverflowCost(num_bins_, rho_, bin_usage_.data(), vector_y_cur_.data(), vector_u_cur_.data());

  float frac_cost
    = computeFractionalCost(num_candidates_, lambda_, vector_x_cur_.data());

  float total_cost = disp_cost + ovf_cost + frac_cost;

  printf("ADMM Iter : %3d Lambda: %f DisplacementCost : %f OverflowCost : %f FractionalCost : %f\n", 
          iter, lambda_, disp_cost, ovf_cost, frac_cost);
}

void
ADMMSolver::setHyperParmeter()
{
  rho_ = 4.0;
  lambda_ = 1.0;
}

void
ADMMSolver::computeFlattenInfo()
{
  // 1. Cell to candidates
  std::vector<std::vector<int>> cell2candidates(num_cells_);
  for(int cand_id = 0; cand_id < num_candidates_; cand_id++)
  {
    auto& candidate = candidates_.at(cand_id);
    assert(cand_id == candidate.cand_id);
    cand_id_to_cell_id_[cand_id] = candidate.cell_id;
    cand_id_to_bin_id_[cand_id] = candidate.bin_id;
    cell2candidates[candidate.cell_id].push_back(cand_id);
  }

//  float sum_min_disp = 0.0;
//  for(int cell_id = 0; cell_id < num_cells_; cell_id++)
//  {
//    auto& cand_vec = cell2candidates.at(cell_id);
//    int min_disp_id = -1;
//    float min_disp = 99999;
//    for(auto cand_id : cand_vec)
//    {
//      float disp_test = disps_.at(cand_id);
//      if(disp_test < min_disp)
//      {
//        min_disp = disp_test;
//        min_disp_id = cand_id;
//      }
//    }
//    sum_min_disp += min_disp;
//    printf("Cell ID : %d Min Disp : %f Min Disp ID : %d\n", cell_id, min_disp, min_disp_id);
//  }
//  printf("Sum Min Disp : %f\n", sum_min_disp);

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

  assert(flatten_idx_cell == num_candidates_);
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

  assert(flatten_idx_bin == num_candidates_);
  bin_id_to_cand_id_start_.back() = num_candidates_;

//  for(int cand_id = 0; cand_id < num_candidates_; cand_id++)
//  {
//    int cell_id = cand_id_to_cell_id_.at(cand_id);
//    int bin_id = cand_id_to_bin_id_.at(cand_id);
//    printf("cand_id %d -> cell : %d bin : %d\n", cand_id, cell_id, bin_id);
//  }
//
//  for(int cell_id = 0; cell_id < num_cells_; cell_id++)
//    printf("cell_id %d cell_id_to_cand_id_start : %d\n", cell_id, cell_id_to_cand_id_start_.at(cell_id));
//
//  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
//    printf("bin_id %d bin_id_to_cand_id_start : %d\n", bin_id, bin_id_to_cand_id_start_.at(bin_id));
//
//  for(int i = 0; i < num_candidates_; i++)
//  {
//    int cell_id = cand_id_to_cell_id_.at(cand_id_foreach_cells_.at(i));
//    printf("i : %d cand_id_foreach_cells : %d cell_id : %d\n", i, cand_id_foreach_cells_.at(i), cell_id);
//  }
//
//  for(int i = 0; i < num_candidates_; i++)
//  {
//    int bin_id = cand_id_to_bin_id_.at(cand_id_foreach_bins_.at(i));
//    printf("i : %d cand_id_foreach_bins : %d bin_id : %d\n", i, cand_id_foreach_bins_.at(i), bin_id);
//  }
}

float
ADMMSolver::computeDisplacement(const int    num_candidates,
                                const float* disp,
                                const float* vector_x)
{
  float sum_disp = 0.0;
  for(int cand_id = 0; cand_id < num_candidates; cand_id++)
    sum_disp += disp[cand_id] * vector_x[cand_id];
  return sum_disp;
}

float
ADMMSolver::computeFractionalCost(const int    num_candidates,
                                  const float  lambda,
                                  const float* vector_x)
{
  float sum_frac_cost = 0.0;
  for(int cand_id = 0; cand_id < num_candidates; cand_id++)
  {
    float x_val = vector_x[cand_id];
    sum_frac_cost += lambda * x_val * (1.0 - x_val);
  }
  return sum_frac_cost;
}

float
ADMMSolver::computeOverflowCost(const int    num_bins,
                                const float  rho,
                                const float* bin_usage,
                                const float* vector_y,
                                const float* vector_u)
{
  float sum_ovf_cost = 0.0;
  for(int bin_id = 0; bin_id < num_bins; bin_id++)
  {
    sum_ovf_cost += rho / 2.0 * std::pow(bin_usage[bin_id] 
                                       + vector_y[bin_id] 
                                       + vector_u[bin_id] / rho, 2);
  }
  return sum_ovf_cost;
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
    for(int index = cand_id_start; index < cand_id_end; index++)
    {
      int cand_id = cand_id_foreach_bins_.at(index);
      int cell_id = cand_id_to_cell_id[cand_id];
      float cell_width = widths_[cell_id];
      float x_val = vector_x[cand_id];

      bin_usage += cell_width * x_val;
    }

    bin_usages[bin_id] = bin_usage - bin_capacity;
  }
}

void
ADMMSolver::updatePrimalX(const int    num_candidates,
                          const int    max_pgd_iter,
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
  const float step_size = 7.8e-9 / rho;

  std::vector<float> vector_x_temp(num_candidates, 0.0);
  std::vector<float> vector_workspace(num_candidates, 0.0);
  std::vector<float> vector_v(num_candidates, 0.0);
  std::vector<float> vector_x_k_minus_1(num_candidates, 0.0);
  std::vector<float> vector_x_k_minus_2(num_candidates, 0.0);

  for(int cand_id = 0; cand_id < num_candidates_; cand_id++)
  {
    vector_x_k_minus_1.at(cand_id) = x_cur[cand_id];
    vector_x_k_minus_2.at(cand_id) = x_cur[cand_id];
  }

  for(int pgd_iter = 0; pgd_iter < max_pgd_iter; pgd_iter++)
  {
    // 1. update v_vector
    float pgd_iter_float = static_cast<float>(pgd_iter);
    for(int v_idx = 0; v_idx < num_candidates; v_idx++)
    {
      float x_k_minus_1 = vector_x_k_minus_1.at(v_idx);
      float x_k_minus_2 = vector_x_k_minus_2.at(v_idx);

      vector_v[v_idx] = 
          (2 * pgd_iter_float - 1) / (pgd_iter_float + 1) * x_k_minus_1
        - (    pgd_iter_float - 2) / (pgd_iter_float + 1) * x_k_minus_2;
    }

    // 2. compute gradient
    computeBinUsage(num_bins_,
                    bin_id_to_cand_id_start_.data(),
                    cand_id_to_cell_id_.data(),
                    widths,
                    capacities,
                    vector_v.data(),
                    bin_usage_.data());

    computeGradient(num_candidates_,
                    rho_,
                    cand_id_to_cell_id_.data(),
                    cand_id_to_bin_id_.data(),
                    disps,
                    widths,
                    bin_usage_.data(),
                    vector_v.data(),
                    y_cur,
                    u_cur,
                    vector_grad_.data());

    std::transform(vector_v.begin(), 
                   vector_v.end(),
                   vector_grad_.begin(),
                   vector_x_temp.begin(),
                   [&] (float v, float grad) { return v - step_size * grad; });

    // 3. simplex projection
    simplexProjection(num_cells_,
                      cell_id_to_num_cand_.data(),
                      vector_x_temp.data(),
                      vector_workspace.data(),
                      vector_x_temp.data());

    if(pgd_iter != max_pgd_iter - 1)
    {
      std::swap(vector_x_k_minus_2, vector_x_k_minus_1);
      std::swap(vector_x_k_minus_1, vector_x_temp);
    }
  }

  for(int cand_id = 0; cand_id < num_candidates_; cand_id++)
    x_next[cand_id] = vector_x_temp.at(cand_id);
}

void 
ADMMSolver::updatePrimalY(const int    num_bins,
                          const float  rho, 
                          const float* bin_usages,
                          const float* u_cur,
                                float* y_next)
{
  for(int bin_id = 0; bin_id < num_bins; bin_id++)
  {
    y_next[bin_id] = 
      std::max(static_cast<float>(0.0), -bin_usages[bin_id] - u_cur[bin_id] / rho);
  }
}

void 
ADMMSolver::updateDual(const int    num_bins,
                       const float  rho, 
                       const float* bin_usages,
                       const float* y_next,
                       const float* u_cur,
                             float* u_next)
{
  for(int bin_id = 0; bin_id < num_bins; bin_id++)
  {
    u_next[bin_id] 
      = u_cur[bin_id] + rho * (bin_usages[bin_id] + y_next[bin_id]);
  }
}

void
ADMMSolver::computeGradient(const int    num_candidates,
                            const float  rho,
                            const int*   cand_id_to_cell_id,
                            const int*   cand_id_to_bin_id,
                            const float* disp,
                            const float* widths,
                            const float* bin_usage,
                            const float* vector_x,
                            const float* vector_y,
                            const float* vector_u,
                                  float* grad)
{
  for(int cand_id = 0; cand_id < num_candidates_; cand_id++)
  {
    int cell_id = cand_id_to_cell_id[cand_id];
    int bin_id = cand_id_to_bin_id[cand_id];
    grad[cand_id] = disp[cand_id] 
                  + rho * widths[cell_id] * (bin_usage[bin_id] 
                                            + vector_y[bin_id] 
                                            + vector_u[bin_id] / rho)
                  + lambda_ * (1.0 - 2.0 * vector_x[cand_id]);
  }
}

void bubbleSort(float* arr, int n)
{
  //printf("Input Array : ");
  //for(int k = 0; k < n; k++)
  //  printf(" %f ", arr[k]);
  //printf("\n");

  int i,j;
  float temp;
  for(j = n - 1; j > 0; j--)
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

  //printf("Sorted Array : ");
  //for(int k = 0; k < n; k++)
  //  printf(" %f ", arr[k]);
  //printf("\n");
}

void 
ADMMSolver::simplexProjection(const int    num_cells,
                              const int*   cell_id_to_num_cand,
                              const float* const vector_input, 
                                    float* vector_workspace,
                                    float* vector_output)
{
//  printf("Input\n");
//  for(int i = 0; i < num_cells_; i++)
//  {
//    for(int j = 0; j < num_bins_; j++)
//      printf(" %f ", vector_input[i * num_bins_ + j]);
//    printf("\n");
//  }

  // Copy vector_input to vector_workspace
  for(int cand_id = 0; cand_id < num_candidates_; cand_id++)
    vector_workspace[cand_id] = vector_input[cand_id];

  for(int cell_id = 0; cell_id < num_cells_; cell_id++)
  {
    int cand_id_start = cell_id_to_cand_id_start_.at(cell_id);
    int num_cand_this_cell = cell_id_to_num_cand_.at(cell_id);

    // printf("cell_id : %d cand_id_start : %d\n", cell_id, cand_id_start);

    float* vector_sorted = vector_workspace + cand_id_start;
    bubbleSort(vector_sorted, num_cand_this_cell);

    int alpha = 0;
    float sum_sorted = 0.0;

    for(int i = 0; i < num_cand_this_cell; i++)
    {
      sum_sorted += vector_sorted[i];
      float i_float = static_cast<float>(i);
      if(vector_sorted[i] + 1.0 / (i_float + 1.0) * (1.0 - sum_sorted) > 0)
        alpha = i;
      else
      {
        assert(i != 0);
        sum_sorted -= vector_sorted[i];
        break;
      }
    }

    float alpha_float = static_cast<float>(alpha);
    float beta = 1.0 / (alpha_float + 1.0) * (1.0 - sum_sorted);

    const float* const vector_this_cell_input = vector_input + cand_id_start;
    float* vector_this_cell_output = vector_output + cand_id_start;
    for(int i = 0; i < num_cand_this_cell; i++)
      vector_this_cell_output[i] = std::max(vector_this_cell_input[i] + beta, static_cast<float>(0.0));
  }

//  printf("Output\n");
//  for(int i = 0; i < num_cells_; i++)
//  {
//    float sum = 0.0;
//    for(int j = 0; j < num_bins_; j++)
//    {
//      sum += vector_output[i * num_bins_ + j];
//      printf(" %f ", vector_output[i * num_bins_ + j]);
//    }
//    printf(" Sum : %f\n", sum);
//  }
}

void
ADMMSolver::makeIntegerSolution(const float* vector_in,
                                      float* vector_out)
{
  for(int cell_id = 0; cell_id < num_cells_; cell_id++)
  {
    int cand_id_start = cell_id_to_cand_id_start_[cell_id];
    int cand_id_end = cell_id_to_cand_id_start_[cell_id + 1];

    int max_cand_id = -1;
    float max_x_val = 0.0;
    for(int index = cand_id_start; index < cand_id_end; index++)
    {
      int cand_id = cand_id_foreach_cells_.at(index);
      float x_val = vector_in[cand_id];

      // printf("x_val : %f\n", x_val);
      if(x_val >= max_x_val)
      {
        max_x_val = x_val;
        max_cand_id = cand_id;
      }
    }

    assert(max_cand_id != -1);
    vector_out[max_cand_id] = 1.0;
  }
}

bool
ADMMSolver::solve()
{
  bool admm_success = false;

  int max_admm_iter = 200;
  for(int admm_iter = 0; admm_iter < max_admm_iter; admm_iter++)
  {
    int max_pgd_iter = 200;

    /* ========================= Main ADMM Iteration ========================== */
    /* 1. Solve subproblem w.r.t. x_k (Projected Gradient) */
    updatePrimalX(num_candidates_,
                  max_pgd_iter,
                  rho_,
                  lambda_,
                  widths_.data(), 
                  disps_.data(), 
                  capacities_.data(), 
                  vector_x_cur_.data(), 
                  vector_y_cur_.data(),
                  vector_u_cur_.data(),
                  vector_x_next_.data()); /* return new x_vector */

    /* 2. Solve subproblem w.r.t. y_k */
    updatePrimalY(num_bins_,
                  rho_, 
                  bin_usage_.data(), 
                  vector_u_cur_.data(),
                  vector_y_next_.data()); /* return new y_vector */

    /* 3. Update Dual Variable u */
    updateDual(num_bins_,
               rho_, 
               bin_usage_.data(), 
               vector_y_next_.data(), 
               vector_u_cur_.data(),
               vector_u_next_.data()); /* return new u_vector */
  
    updateNextIter(admm_iter);
    /* ======================================================================== */
  }
  
  float final_disp = computeDisplacement(num_candidates_, disps_.data(), vector_x_cur_.data());
  printf("Final Displacement : %f\n", final_disp);

  computeBinUsage(num_bins_,
                  bin_id_to_cand_id_start_.data(),
                  cand_id_to_cell_id_.data(),
                  widths_.data(),
                  capacities_.data(),
                  vector_x_cur_.data(),
                  bin_usage_.data());

  printf("Bin Usage\n");
  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
    printf("  %f\n", bin_usage_.at(bin_id) + capacities_.at(bin_id));

  std::vector<float> vector_x_int(num_candidates_, 0.0);
  makeIntegerSolution(vector_x_cur_.data(), vector_x_int.data());

  float final_disp_int = computeDisplacement(num_candidates_, disps_.data(), vector_x_int.data());
  printf("Final Displacement (int) : %f\n", final_disp_int);

  computeBinUsage(num_bins_,
                  bin_id_to_cand_id_start_.data(),
                  cand_id_to_cell_id_.data(),
                  widths_.data(),
                  capacities_.data(),
                  vector_x_int.data(),
                  bin_usage_.data());

  printf("Bin Usage (int)\n");
  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
    printf("  %f\n", bin_usage_.at(bin_id) + capacities_.at(bin_id));

  return admm_success;
}

}
