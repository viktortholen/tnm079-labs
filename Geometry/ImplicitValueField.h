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
#ifndef __implicit_value_field_h__
#define __implicit_value_field_h__

#include <Geometry/Implicit.h>
#include <Math/Function3D.h>

class ImplicitValueField : public Function3D<float> {
protected:
  const Implicit *mImplicit;

public:
  ImplicitValueField(const Implicit *implicit) : mImplicit(implicit) {}
  virtual ~ImplicitValueField() {}

  //! Evaluate the function at x,y,z
  virtual float GetValue(float x, float y, float z) const {
    return mImplicit->GetValue(x, y, z);
  }

  //! Return a bound on the maximum value of the function
  virtual float GetMaxValue() const { return 1; }
};

#endif
