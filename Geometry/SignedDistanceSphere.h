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
#ifndef __signed_distance_sphere_h__
#define __signed_distance_sphere_h__

#include <Geometry/Implicit.h>

/*!  \brief Sphere base class */
class SignedDistanceSphere : public Implicit {
public:
  SignedDistanceSphere(float r);
  virtual ~SignedDistanceSphere();
  virtual float getValue(float x, float y, float z) const;

protected:
  float radius;
};

#endif
