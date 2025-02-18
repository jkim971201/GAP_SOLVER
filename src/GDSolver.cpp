#include <algorithm>

#include "GDSolver.h"

namespace gdsolver
{

GDSolver::GDSolver(std::shared_ptr<gapbuilder::GAPInstance> instance)
{
  instance_ = instance;

  num_cells_ = instance_->num_items;
  num_bins_  = instance_->num_sacks;

  const auto& profits_2d = instance_->profits;
  const auto& weights_2d = instance_->weights;
  const auto& capacities_int = instance_->capacities;

  disps_.resize(num_cells_ * num_bins_);
  widths_.resize(num_cells_);
  capacities_.resize(num_bins_);

  printf("1\n");
  std::transform(capacities_int.begin(), 
                 capacities_int.end(),
                 capacities_.begin(),
                 [] (int cap) { return static_cast<float>(cap); });

  printf("2\n");
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

  printf("3\n");
  num_candidates_ = static_cast<int>(candidates_.size());

  bin_slack_.resize(num_bins_, 0.0);
  penalty_.resize(num_bins_, 0.0);

  vector_x_.resize(num_candidates_, 0.0);
  vector_y_.resize(num_candidates_, 0.0);

  vector_s_c.resize(num_cells_, 0.0);

  printf("4\n");
  computeFlattenInfo();

  printf("5\n");
}

void
GDSolver::computeFlattenInfo()
{
  num_bin_cell_.resize(num_bins_, 0);
  bin_id_start_.resize(num_bins_, 0);
  cells_in_bin_.resize(num_candidates_, -1);

  std::vector<std::vector<int>> bin2cells(num_bins_);
  for(auto& cand : candidates_)
  {
    int cell_id = cand.cell_id;
    int bin_id  = cand.bin_id;
    bin2cells[bin_id].push_back(cell_id);
  }

  int flatten_idx = 0;
  for(int bin_id = 0; bin_id < num_bins_; bin_id++)
  {
    bin_id_start_[bin_id] = flatten_idx;

    const auto& cells = bin2cells.at(bin_id);
    int num_cell_in_this_bin = static_cast<int>(cells.size());
    num_bin_cell_[bin_id] = num_cell_in_this_bin;

    for(int cell_id : cells)
    {
      cells_in_bin_[flatten_idx] = cell_id;
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
GDSolver::computeSubGradient()
{
  std::vector<float> df_dy(num_cells_ * num_bins_, 0.0);
  std::vector<float> dy_dx(num_cells_ * num_bins_, 0.0);
  std::vector<float> df_dx(num_cells_ * num_bins_, 0.0);

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
