#ifndef GD_SOLVER_H
#define GD_SOLVER_H

#include "GAPBuilder.h"
#include <vector>

namespace gdsolver
{

struct GAPInstance;

class GDSolver
{
  public:

    GDSolver(std::shared_ptr<GAPInstance> instance);

  private:

    std::shared_ptr<GAPInstance> instance_;
};

}

#endif
