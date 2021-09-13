#pragma once

#include "Fluid/FluidSolver.h"
#include "Geometry/Geometry.h"
#include "Util/ColorMap.h"

#include <vector>

class FluidVoxelCutPlane : public Geometry {
public:
  FluidVoxelCutPlane(const std::string &name, const FluidSolver *solver);

  virtual const char *GetTypeName() {
    return typeid(FluidVoxelCutPlane).name();
  }

  virtual void Render();

  virtual void Initialize() {}

  virtual void Update();

  virtual void SetTransform(const glm::mat4&transform);

protected:
  const FluidSolver *mSolver;

  float mDx;

  std::vector<bool> mFluidVoxels;
  std::vector<glm::vec3> mFluidVoxelsPositions;
};
