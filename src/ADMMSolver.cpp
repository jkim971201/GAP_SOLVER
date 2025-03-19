#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <numeric>

#include "ADMMSolver.h"

#define NUM_THREAD 1

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
  widths_.resize(num_cells_ * num_bins_);
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
      widths_[cand_id] = static_cast<float>(weights_2d.at(cell_id).at(bin_id));

      num_cand_made++;
    }
  }

  num_candidates_ = static_cast<int>(candidates_.size());

  bin_slack_.resize(num_bins_, 0.0);

  vector_y_cur_.resize(num_candidates_, 0.0);
  vector_l_cur_.resize(num_candidates_, 0.0);

  vector_y_next_.resize(num_candidates_, 0.0);
  vector_l_next_.resize(num_candidates_, 0.0);

  cand_id_to_cell_id_.resize(num_candidates_, 0);
  cand_id_to_bin_id_.resize(num_candidates_, 0);

  cell_to_best_cand_.resize(num_cells_, 0);

  computeFlattenInfo();

  initializeX(num_cells_, 
              disps_.data(), 
              cell_to_best_cand_.data());
}

void
ADMMSolver::initializeX(const int    num_cells,
                        const float* disps,
                              int*   cell_to_best_cand)
{
  for(int cell_id = 0; cell_id < num_cells; cell_id++)
  {
    int min_cand_id = -1;
    float min_disp = std::numeric_limits<float>::max();
  
    const auto& candidate_this_cell = cell_id_to_cands_.at(cell_id);
    for(int cand_id : candidate_this_cell)
    {
      int disp_this_cand = disps[cand_id];
      if(disp_this_cand < min_disp)
      {
        min_disp = disp_this_cand;
        min_cand_id = cand_id;
      }
    }

    assert(min_cand_id != -1);
    cell_to_best_cand[cell_id] = min_cand_id;
  }

  float initial_disp 
    = computeDisplacement(num_cells_, 
                          cell_to_best_cand_.data(),
                          disps_.data());

  printf("Initial Displacement : %f\n", initial_disp);

  computeBinSlack(num_cells_,
                  cell_to_best_cand_.data(),
                  widths_.data(),
                  capacities_.data(),
                  bin_slack_.data());

  printf("Bin Usage\n");
  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
    printf("  %f\n", bin_slack_.at(bin_id) + capacities_.at(bin_id));
}

void
ADMMSolver::updateNextIter(int iter)
{
  std::swap(vector_y_cur_, vector_y_next_);
  std::swap(vector_l_cur_, vector_l_next_);

  rho_ = rho_ * 1.01;

  int disp 
    = computeDisplacement(num_cells_, 
                          cell_to_best_cand_.data(),
                          disps_.data());

  printf("Iter: %03d ", iter);
  double sum_ovf = 0.0;
  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
  {
    double slack = bin_slack_.at(bin_id);
    sum_ovf += std::max(0.0, slack);
  }
  printf("Disp : %d SumOvf : %f\n", disp, sum_ovf);

  //printf("Bin Usage\n");
  //for(int bin_id = 0; bin_id < num_bins_; bin_id++)
  //  printf("  %f\n", bin_slack_.at(bin_id) + capacities_.at(bin_id));
}

void
ADMMSolver::setHyperParmeter()
{
  rho_ = 0.0001;
}

void
ADMMSolver::computeFlattenInfo()
{
  // 1. Cell to candidates
  cell_id_to_cands_.resize(num_cells_);
  for(int cand_id = 0; cand_id < num_candidates_; cand_id++)
  {
    auto& candidate = candidates_.at(cand_id);
    int cell_id = candidate.cell_id;
    cell_ids_.push_back(cell_id);
    cell_id_to_cands_[cell_id].push_back(cand_id);
    cand_id_to_cell_id_[cand_id] = cell_id;
  }

  // 2. Bin to candidates
  bin_id_to_cands_.resize(num_bins_);
  for(int cand_id = 0; cand_id < num_candidates_; cand_id++)
  {
    auto& candidate = candidates_.at(cand_id);
    int bin_id = candidate.bin_id;
    bin_ids_.push_back(bin_id);
    bin_id_to_cands_[bin_id].push_back(cand_id);
    cand_id_to_bin_id_[cand_id] = bin_id;
  }
}

float
ADMMSolver::computeDisplacement(const int    num_cells,
                                const int*   cell_to_best_cand,
                                const float* disp)
{
  float sum_disp = 0.0;
  for(int cell_id = 0; cell_id < num_cells; cell_id++)
  {
    int best_cand_this_cell = cell_to_best_cand[cell_id];
    sum_disp += disp[best_cand_this_cell];
  }
  return sum_disp;
}

void
ADMMSolver::computeBinSlack(const int    num_cells,
                            const int*   cell_to_best_cand,
                            const float* widths,
                            const float* capacities,
                                  float* bin_slack) /* return vector */
{
  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
    bin_slack[bin_id] = -capacities[bin_id];

  for(int cell_id = 0; cell_id < num_cells; cell_id++)
  {
    int best_cand_id = cell_to_best_cand[cell_id];
    int bin_id = cand_id_to_bin_id_.at(best_cand_id);
    bin_slack[bin_id] += widths[best_cand_id];
  }
}

void
ADMMSolver::updatePrimalX(const int    num_cells,
                          const float  rho, 
                          const float* widths,
                          const float* disps,
                          const float* capacities,
                          const float* y_cur,
                          const float* l_cur,
                                float* bin_slack,
                                int*   cell_to_best_cand)
{
  std::sort(cell_ids_.begin(), cell_ids_.end(),
            [&] (int id1, int id2) 
            {
              int cand1 = cell_to_best_cand[id1];
              int cand2 = cell_to_best_cand[id2];
              int bin1 = cand_id_to_bin_id_.at(cand1);
              int bin2 = cand_id_to_bin_id_.at(cand2);
              return bin_slack[bin1] > bin_slack[bin2];
            });

  for(int cell_id : cell_ids_)
  {
    int cand_id_cur = cell_to_best_cand[cell_id];
    int bin_id_cur_cand = cand_id_to_bin_id_.at(cand_id_cur);
    float width_cur_cand = widths[cand_id_cur];
    bin_slack[bin_id_cur_cand] -= width_cur_cand;

    int min_cand_id = -1;
    float min_cost = std::numeric_limits<float>::max();

    auto& cands_this_cell = cell_id_to_cands_.at(cell_id);
    for(int cand_id : cands_this_cell)
    {
      int bin_id = cand_id_to_bin_id_.at(cand_id);
      float cost = disps[cand_id] 
        + widths[cand_id] * std::max(float(0), rho * bin_slack[bin_id] + rho * y_cur[bin_id] + l_cur[bin_id]);

      assert(cost >= 0.0);
      if(cost < min_cost)
      {
        min_cost = cost;
        min_cand_id = cand_id;
      }
    }

    assert(min_cand_id != -1);
    
    int cand_id_new = min_cand_id;
    int bin_id_new_cand = cand_id_to_bin_id_.at(cand_id_new);
    float width_new_cand = widths[cand_id_new];
    bin_slack[bin_id_new_cand] += width_new_cand;

    cell_to_best_cand[cell_id] = cand_id_new;
  }
}

void 
ADMMSolver::updatePrimalY(const int    num_bins,
                          const float  rho, 
                          const float* bin_slack,
                          const float* l_cur,
                                float* y_next)
{
  for(int bin_id = 0; bin_id < num_bins; bin_id++)
  {
    y_next[bin_id] 
      = std::max(float(0.0), -bin_slack[bin_id] - l_cur[bin_id] / rho);
  }
}

void 
ADMMSolver::updateDual(const int    num_bins,
                       const float  rho, 
                       const float* bin_slack,
                       const float* y_next,
                       const float* l_cur,
                             float* l_next)
{
  for(int bin_id = 0; bin_id < num_bins; bin_id++)
  {
    float s = bin_slack[bin_id];
    float y = y_next[bin_id];
    float l = l_cur[bin_id];
    l_next[bin_id] = std::max(float(0.0), l + rho * (s + y));
  }
}

void bubbleSort(float* arr, int n)
{
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
}

bool
ADMMSolver::checkTermination()
{
  double sum_ovf = 0.0;
  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
  {
    double slack = bin_slack_.at(bin_id);
    sum_ovf += std::max(0.0, slack);
  }

  if(sum_ovf < 1000.0)
    return true;
  else
    return false;
}

bool
ADMMSolver::solve()
{
  bool admm_success = false;

  int max_admm_iter = 200;
  for(int admm_iter = 0; admm_iter < max_admm_iter; admm_iter++)
  {
    /* ===================== Main ADMM Iteration ====================== */
    /* 1. Solve subproblem w.r.t. x_k */
    updatePrimalX(num_cells_,
                  rho_,
                  widths_.data(), 
                  disps_.data(), 
                  capacities_.data(), 
                  vector_y_cur_.data(),
                  vector_l_cur_.data(),
                  bin_slack_.data(),
                  cell_to_best_cand_.data()); /* return new x_vector */

    /* 2. Solve subproblem w.r.t. y_k */
    updatePrimalY(num_bins_,
                  rho_, 
                  bin_slack_.data(), 
                  vector_l_cur_.data(),
                  vector_y_next_.data()); /* return new y_vector */

    /* 3. Update Dual Variable lambda */
    updateDual(num_bins_,
               rho_, 
               bin_slack_.data(), 
               vector_y_next_.data(), 
               vector_l_cur_.data(),
               vector_l_next_.data()); /* return new l_vector */
  
    updateNextIter(admm_iter);

    bool is_termination = checkTermination();
    if(is_termination == true)
      break;
    /* ================================================================ */
  }
  
  float final_disp = computeDisplacement(num_cells_, cell_to_best_cand_.data(), disps_.data());
  printf("Final Displacement : %f\n", final_disp);

  computeBinSlack(num_cells_,
                  cell_to_best_cand_.data(),
                  widths_.data(),
                  capacities_.data(),
                  bin_slack_.data());

  printf("Bin Usage\n");
  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
    printf("  %f\n", bin_slack_.at(bin_id) + capacities_.at(bin_id));

  return admm_success;
}

}
