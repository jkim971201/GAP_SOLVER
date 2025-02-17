#ifndef GAP_BUILDER_H
#define GAP_BUILDER_H

#include <vector>
#include <fstream>
#include <memory>
#include <string_view>

namespace gapbuilder
{

struct GAPInstance
{
  int num_items = 0;
  int num_sacks = 0;
  std::vector<std::vector<int>> profits;
  std::vector<std::vector<int>> weights;
  std::vector<int> capacities;
};

class GAPBuilder
{
  public:

    GAPBuilder(std::string_view file_name);

    std::shared_ptr<GAPInstance> getGAPInstance();

    void print() const;

  private:

    void buildInstance(std::string_view file_name);
    std::shared_ptr<GAPInstance> instance_;
};

}

#endif
