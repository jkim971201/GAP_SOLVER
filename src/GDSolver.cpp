#include <algorithm>
#include <numeric>
#include <cmath>

#include "GDSolver.h"

namespace gdsolver
{

GDSolver::GDSolver(std::shared_ptr<gapbuilder::GAPInstance> instance)
{
	softmax_tmpr_ = 1.0;

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
      ILPCandidate new_cand = {cand_id, cell_id, bin_id, 0.0, 0.0};
      candidates_.push_back(new_cand);

      disps_[cand_id] = static_cast<float>(profits_2d.at(cell_id).at(bin_id));
      widths_[cell_id] = static_cast<float>(weights_2d.at(cell_id).at(bin_id));

      num_cand_made++;
    }
  }

  num_candidates_ = static_cast<int>(candidates_.size());

  bin_slack_.resize(num_bins_, 0.0);
  penalty_.resize(num_bins_, 1.0);

  vector_x_.resize(num_candidates_, 0.0);
  vector_y_.resize(num_candidates_, 0.0);

  vector_sum_exp_.resize(num_cells_, 0.0);
	vector_min_x_.resize(num_cells_, 0.0);

  computeFlattenInfo();
}

void
GDSolver::computeFlattenInfo()
{
  num_cands_each_cell_.resize(num_cells_, 0);
  cell_id_to_cand_start_.resize(num_cells_ + 1, 0);
  cands_in_cell_.resize(num_cells_, -1);

  std::vector<std::vector<int>> cell2candidates(num_cells_);
  for(auto& cand : candidates_)
    cell2candidates[cand.cell_id].push_back(cand.cand_id);

  int flatten_idx = 0;
  for(int cell_id = 0; cell_id < num_cells_; cell_id++)
  {
    cell_id_to_cand_start_[cell_id] = flatten_idx;

    const auto& cands = cell2candidates.at(cell_id);
    int num_cand_this_cell = static_cast<int>(cands.size());
    num_cands_each_cell_[cell_id] = num_cand_this_cell;

    for(int cand_id : cands)
    {
      cands_in_cell_[flatten_idx] = cand_id;
      flatten_idx++;
    }
  }
}

void
GDSolver::updateBinPenalty()
{
}

void
GDSolver::computeBinSlack()
{
  // Initialize as Bin Capacity
  for(int i = 0; i < num_bins_; i++)
    bin_slack_[i] = capacities_.at(i);

  for(auto& cand : candidates_)
  {
    int cell_id = cand.cell_id;
    int bin_id  = cand.bin_id;
   
    // subtract cell width
    bin_slack_[bin_id] -= widths_.at(cell_id) * cand.val_y;
  }

  for(auto slack : bin_slack_)
    printf("Bin Slack : %f\n", slack);
}

void
GDSolver::computeSumExpX()
{
	// Step #1. Find Min X
	// To avoid overflow of exponential, softmax takes x_min as offset.
	for(int cell_id = 0; cell_id < num_cells_; cell_id++)
	{
    int cand_id_start = cell_id_to_cand_start_.at(cell_id);
		int cand_id_start_next = cell_id_to_cand_start_.at(cell_id + 1);
    int num_cands_this_cell = num_cands_each_cell_.at(cell_id);

    float min_x_this_cell = std::numeric_limits<float>::max();
		for(int cand_id = cand_id_start; cand_id < cand_id_start_next; cand_id++)
      min_x_this_cell = std::min(min_x_this_cell, vector_x_.at(cand_id));

		vector_min_x_[cell_id] = min_x_this_cell;
	}

  // Step #2. Compute Sum of Exponential
	for(int cell_id = 0; cell_id < num_cells_; cell_id++)
	{
    float sum_exp = 0.0;
		float min_x_this_cell = static_cast<float>(vector_min_x_.at(cell_id));

    int cand_id_start = cell_id_to_cand_start_.at(cell_id);
		int cand_id_start_next = cell_id_to_cand_start_.at(cell_id + 1);
		for(int cand_id = cand_id_start; cand_id < cand_id_start_next; cand_id++)
			sum_exp += std::exp((vector_x_.at(cand_id) - min_x_this_cell) / softmax_tmpr_);

		vector_sum_exp_[cell_id] = sum_exp;
	}
}

void
GDSolver::computeSubGradient()
{
  std::vector<float> df_dy(num_cells_ * num_bins_, 0.0);
  std::vector<float> dy_dx(num_cells_ * num_bins_, 0.0);
  std::vector<float> df_dx(num_cells_ * num_bins_, 0.0);

  computeSumExpX();

  for(auto& cand : candidates_)
  {
    int i = cand.cand_id; // index of current ilp candidate

    // Step #1. Compute partial f / partial y regard to i
    int bin_id  = cand.bin_id;
    int cell_id = cand.cell_id;

    float slack_b = bin_slack_.at(bin_id);
    float penalty_b = penalty_.at(bin_id);

    float disp_i = disps_.at(i);
    float width_i = widths_.at(cell_id);

    df_dy[i] = slack_b < 0.0 ? disp_i - 2.0 * penalty_b * width_i * slack_b
                             : disp_i;

    // Step #2. Compute partial y / partial x regard to i
		float exp_x = std::exp((vector_x_.at(i) - vector_min_x_.at(cell_id)) / softmax_tmpr_);
		float sum_exp_x = vector_sum_exp_.at(cell_id);
		dy_dx[i] = (exp_x * sum_exp_x - exp_x * exp_x) / (softmax_tmpr_ * sum_exp * sum_exp_x);

		// Step #3. Apply chain rule
    df_dx[i] = df_dy[i] * dy_dx[i];
  }
}

bool
GDSolver::solve()
{
  bool gd_success = false;

  computeBinSlack();
  computeSubGradient();

  return gd_success;
}

}
