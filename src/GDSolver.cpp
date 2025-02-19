#include <algorithm>
#include <numeric>
#include <cmath>

#include "GDSolver.h"

namespace gdsolver
{

GDSolver::GDSolver(std::shared_ptr<gapbuilder::GAPInstance> instance)
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
      ILPCandidate new_cand = {cand_id, cell_id, bin_id, 0.0, 0.0};
      candidates_.push_back(new_cand);

      disps_[cand_id] = static_cast<float>(profits_2d.at(cell_id).at(bin_id));
      widths_[cell_id] = static_cast<float>(weights_2d.at(cell_id).at(bin_id));

      num_cand_made++;
    }
  }

  num_candidates_ = static_cast<int>(candidates_.size());

  bin_overflow_.resize(num_bins_, 0.0);
  penalty_.resize(num_bins_, overflow_penalty_);

  vector_x_.resize(num_candidates_, 0.0);
  vector_y_.resize(num_candidates_, 0.0);

  vector_sum_exp_.resize(num_cells_, 0.0);
  vector_exp_.resize(num_candidates_, 0.0);
  vector_min_x_.resize(num_cells_, 0.0);

  df_dy_.resize(num_cells_ * num_bins_, 0.0);
  dy_dx_.resize(num_cells_ * num_bins_, 0.0);
  df_dx_.resize(num_cells_ * num_bins_, 0.0);

  computeFlattenInfo();
}

void
GDSolver::setHyperParmeter()
{
  overflow_penalty_ = 1.0;
  softmax_tmpr_ = 2.0;
  step_size_ = 0.001;
  max_gd_iter_ = 50;
}

void
GDSolver::computeFlattenInfo()
{
  num_cands_each_cell_.resize(num_cells_, 0);
  cell_id_to_cand_start_.resize(num_cells_ + 1, 0);
  cands_in_cell_.resize(num_candidates_, -1);

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

      printf("cell_id : %d cand_id : %d\n", cell_id, cand_id);
    }
  }

  cell_id_to_cand_start_.back() = num_candidates_;
}

void
GDSolver::updateParameter()
{
	softmax_tmpr_ *= 0.95;

	for(int bin_id = 0; bin_id < num_bins_; bin_id++)
	{
		float ovf_b = bin_overflow_.at(bin_id);
    float lambda_b = penalty_.at(bin_id);
    
		if(ovf_b > 0.0)
			penalty_[bin_id] = lambda_b * 1.05;
	}
}

void
GDSolver::computeBinOverflow()
{
  // Initialize as minus of bin capacity
  for(int i = 0; i < num_bins_; i++)
    bin_overflow_[i] = -1.0 * capacities_.at(i);

  for(auto& cand : candidates_)
  {
    int i = cand.cand_id;
    int cell_id = cand.cell_id;
    int bin_id  = cand.bin_id;
   
    // Add cell width * probability of the candidate
    bin_overflow_[bin_id] += widths_.at(cell_id) * vector_y_.at(i);
  }

  //for(auto ovf : bin_overflow_)
  //  printf("Bin Overflow : %f\n", ovf);
}

void
GDSolver::computeExpX()
{
  // Step #1. Find Min X
  // To avoid overflow of exponential, softmax takes x_min as offset.
  for(int cell_id = 0; cell_id < num_cells_; cell_id++)
  {
    int iter_start = cell_id_to_cand_start_.at(cell_id);
    int iter_end = cell_id_to_cand_start_.at(cell_id + 1);

    float min_x_this_cell = std::numeric_limits<float>::max();
    for(int iter = iter_start; iter < iter_end; iter++)
    {
      int cand_id = cands_in_cell_.at(iter);
      min_x_this_cell = std::min(min_x_this_cell, vector_x_.at(cand_id));
    }

    vector_min_x_[cell_id] = min_x_this_cell;
  }

  // Step #2. Compute Sum of Exponential
  for(int cell_id = 0; cell_id < num_cells_; cell_id++)
  {
    float sum_exp = 0.0;
    float min_x_this_cell = static_cast<float>(vector_min_x_.at(cell_id));

    int iter_start = cell_id_to_cand_start_.at(cell_id);
    int iter_end = cell_id_to_cand_start_.at(cell_id + 1);

    for(int iter = iter_start; iter < iter_end; iter++)
    {
      int cand_id = cands_in_cell_.at(iter);
      sum_exp += std::exp((vector_x_.at(cand_id) - min_x_this_cell) / softmax_tmpr_);
    }

    vector_sum_exp_[cell_id] = sum_exp;
  }

  // Step #3. Compute Exponential of x_i
  for(auto& cand : candidates_)
  {
    int i = cand.cand_id;
    int cell_id = cand.cell_id;
    float exp_x = std::exp((vector_x_.at(i) - vector_min_x_.at(cell_id)) / softmax_tmpr_);
    vector_exp_[i] = exp_x;
  }
}

void
GDSolver::computeSubGradient()
{
  for(auto& cand : candidates_)
  {
    int i = cand.cand_id; // index of current ilp candidate

    // Step #1. Compute partial f / partial y regard to i
    int bin_id  = cand.bin_id;
    int cell_id = cand.cell_id;

    float ovf_b = bin_overflow_.at(bin_id);
    float penalty_b = penalty_.at(bin_id);

    float disp_i = disps_.at(i);
    float width_i = widths_.at(cell_id);

    df_dy_[i] = ovf_b > 0.0 ? disp_i + 2.0 * penalty_b * width_i * ovf_b : disp_i;

    // Step #2. Compute partial y / partial x regard to i
    float exp_x = vector_exp_.at(i);
    float sum_exp_x = vector_sum_exp_.at(cell_id);
    dy_dx_[i] = (exp_x * sum_exp_x - exp_x * exp_x) / (softmax_tmpr_ * sum_exp_x * sum_exp_x);

    // Step #3. Apply chain rule
    df_dx_[i] = df_dy_[i] * dy_dx_[i];

    //printf("cand[%2d] : df_dy = :%f\n", i, df_dy_.at(i));
    //printf("cand[%2d] : dy_dx = :%f\n", i, dy_dx_.at(i));
    //printf("cand[%2d] : df_dx = :%f\n", i, df_dx_.at(i));
  }
}

void
GDSolver::updateXYVector()
{
  std::transform(vector_x_.begin(), 
                 vector_x_.end(),
                 df_dx_.begin(),
                 vector_x_.begin(),
                 [&] (float x, float sub_g) { return x - step_size_ * sub_g; });

  for(auto& cand : candidates_)
  {
    int i = cand.cand_id;
    int cell_id = cand.cell_id;
    float sum_exp_x = vector_sum_exp_.at(cell_id);
    float exp_x = vector_exp_.at(i);
    vector_y_[i] = exp_x / sum_exp_x;
  }

  // To debug
  for(int cell_id = 0; cell_id < num_cells_; cell_id++)
  {
    float sum_y = 0.0;

    int iter_start = cell_id_to_cand_start_.at(cell_id);
    int iter_end = cell_id_to_cand_start_.at(cell_id + 1);

    for(int iter = iter_start; iter < iter_end; iter++)
    {
      int cand_id = cands_in_cell_.at(iter);
      // printf("cand_id : %d cell_id : %d\n", cand_id, cell_id);
      sum_y += vector_y_.at(cand_id);
    }

    // printf("Cell : %d Sum_y = %f\n", cell_id, sum_y);
  }
}

bool
GDSolver::solve()
{
  bool gd_success = false;

  computeExpX();
  computeBinOverflow();
  updateXYVector();

  for(int gd_iter = 0; gd_iter < max_gd_iter_; gd_iter++)
  {
	  int disp_int = 0;
    int ovf_int = 0;
	  float disp_float = 0.0;
	  float ovf_float = 0.0;
    for(int i = 0; i < num_candidates_; i++)
		{
			float y_float = vector_y_.at(i);
			int y_int = y_float > 0.5 ? 1 : 0;

			disp_float += disps_.at(i) * y_float;
			ovf_float += bin_overflow_.at(candidates_.at(i).bin_id) * y_float;

			disp_int += disps_.at(i) * y_int;
			ovf_int += bin_overflow_.at(candidates_.at(i).bin_id) * y_int;
		}

		printf("Iter : %d ovf_float : %f ovf_int : %d disp_float : %f disp_int : %d\n",
				gd_iter, ovf_float, ovf_int, disp_float, disp_int);

//    for(int i = 0; i < num_candidates_; i++)
//    {
//      printf("x[%02d] : %f ", i, vector_x_.at(i));
//      printf("y[%02d] : %f\n", i, vector_y_.at(i));
//    }

    // 1. compute subgrad df/dx
    computeSubGradient();

    // 2. x <- x - step * df/dx, y = softmax(x)
    updateXYVector();

    // 3. precompute softmax
    computeExpX();

    // 4. update bin overflow
    computeBinOverflow();

    // 5. update overflow penalty and stepsize
    updateParameter();
  }

	std::vector<int> width_sum(num_bins_, 0);
  for(auto& cand : candidates_)
	{
    int i = cand.cand_id;
		int cell_id = cand.cell_id;
		int bin_id = cand.bin_id;

		if(vector_y_.at(i) >= 0.5)
		{
		  printf("move cell %d to bin %d\n", cell_id, bin_id);
			width_sum[bin_id] += widths_.at(cell_id);
		}
	}

	for(auto w_s : width_sum)
		printf("WidthSum : %d\n", w_s);

  return gd_success;
}

}
