#pragma once

#include <Geometry/Bbox.h>
#include <Geometry/CSG.h>
#include <Geometry/Geometry.h>
#include <Geometry/Implicit.h>
#include <Geometry/Sphere.h>
#include <vector>

/*! \brief Fractal of spheres class*/
class SphereFractal : public Implicit {
public:
  //! Constructor. Level denotes the level of the fractal. Level can be between
  //! 2 and 5
  SphereFractal(int level = 2);

  ~SphereFractal();

  virtual float GetValue(float x, float y, float z) const;

  //! Builds the fractal. Returns a pointer to an implicit geometry object.
  Implicit *buildFractal();

private:
  std::vector<Implicit *> mSpheres;
  Implicit *mFractal;
};
