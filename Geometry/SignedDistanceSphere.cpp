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
#include <Geometry/SignedDistanceSphere.h>

SignedDistanceSphere::SignedDistanceSphere(float r) {
  this->radius = r;
  this->mBox = Bbox(glm::vec3(-r, -r, -r), glm::vec3(r, r, r));
}

SignedDistanceSphere::~SignedDistanceSphere() {}

float SignedDistanceSphere::getValue(float x, float y, float z) const {
  TransformW2O(x, y, z);
  return std::sqrt(x * x + y * y + z * z) - radius;
}
