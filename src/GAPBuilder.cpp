#include "GAPBuilder.h"

namespace gapbuilder
{

GAPBuilder::GAPBuilder(std::string_view file_name)
{
  instance_ = std::make_shared<GAPInstance>();
  buildInstance(file_name);
}

std::shared_ptr<GAPInstance> 
GAPBuilder::getGAPInstance()
{
  return instance_;
}

void
GAPBuilder::buildInstance(std::string_view file_name)
{
  std::vector<std::vector<int>>& profits = instance_->profits;
  std::vector<std::vector<int>>& weights = instance_->weights;
  std::vector<int>& capacities = instance_->capacities;

  int& num_items = instance_->num_items;
  int& num_sacks = instance_->num_sacks;
  std::ifstream input_file(file_name.data());

  input_file >> num_sacks >> num_items;

  profits.resize(num_items);
  weights.resize(num_items);

  for(int sack_id = 0; sack_id < num_sacks; sack_id++)
  {
    for(int item_id = 0; item_id < num_items; item_id++)
    {
      int new_profit = 0;
      input_file >> new_profit;
      profits[item_id].push_back(new_profit);
    }
  }

  for(int sack_id = 0; sack_id < num_sacks; sack_id++)
  {
    for(int item_id = 0; item_id < num_items; item_id++)
    {
      int new_weight = 0;
      input_file >> new_weight;
      weights[item_id].push_back(new_weight);
    }
  }

  for(int sack_id = 0; sack_id < num_sacks; sack_id++)
  {
    int new_cap = 0;
    input_file >> new_cap;
    capacities.push_back(new_cap);
  }

}

void 
GAPBuilder::print() const
{
  const int num_items = instance_->num_items;
  const int num_sacks = instance_->num_sacks;

  auto& profits = instance_->profits;
  auto& weights = instance_->weights;
  auto& capacities = instance_->capacities;

  printf("(NumItems, NumSacks) : (%d, %d)\n", 
      num_items, num_sacks);

  // Print Objective
  printf("\n");
  printf("Minimize\n");
  printf("\n");
  for(int item_id = 0; item_id < num_items; item_id++)
  {
    printf(" ");
    for(int sack_id = 0; sack_id < num_sacks; sack_id++)
      printf(" + %2d * x_{%d,%d}", profits.at(item_id).at(sack_id), item_id, sack_id);
    printf("\n");
  }

  // Print Selection Constraint
  printf("\n");
  printf("Satisfying\n");
  printf("\n");
  for(int item_id = 0; item_id < num_items; item_id++)
  {
    printf(" ");
    for(int sack_id = 0; sack_id < num_sacks; sack_id++)
      printf(" + x_{%d,%d}", item_id, sack_id);
    printf(" = 1\n");
  }

  // Print Capacity Constraint
  printf("\n");
  for(int sack_id = 0; sack_id < num_sacks; sack_id++)
  {
    printf(" ");
    for(int item_id = 0; item_id < num_items; item_id++)
      printf(" + %2d * x_{%d,%d}", weights.at(item_id).at(sack_id), item_id, sack_id);
    printf(" <= %d\n", capacities.at(sack_id));
  }

  printf("\n");
}

};
