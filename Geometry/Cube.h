/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 * Acknowledgements for original code base:
 * - Gunnar Johansson
 * - Ken Museth
 * - Michael Bang Nielsen
 * - Ola Nilsson
 * - Andreas Soderstrom
 *
 * Code updated in the period 2017-2018 by Jochen Jankowai
 *
 *************************************************************************************************/
#ifndef __CUBE_H__
#define __CUBE_H__

#include <Geometry/Quadric.h>
#include <vector>

class Cube : public Implicit {
public:
  Cube();
  ~Cube();

  virtual float GetValue(float x, float y, float z) const;

private:
  std::vector<Quadric *> mPlanes;
  bool mEuclideanDistance;
};

#endif
