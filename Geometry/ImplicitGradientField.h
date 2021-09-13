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
#ifndef __implicit_gradient_field_h__
#define __implicit_gradient_field_h__

#include <Geometry/Implicit.h>
#include <Math/Function3D.h>

class ImplicitGradientField : public Function3D<glm::vec3> {
protected:
    const Implicit *mImplicit;

public:
    ImplicitGradientField(const Implicit *implicit) : mImplicit(implicit){};
    virtual ~ImplicitGradientField() {}

    //! Evaluate the function at x,y,z
    virtual glm::vec3 GetValue(float x, float y, float z) const {
        return mImplicit->GetGradient(x, y, z);
    }

    //! Return a bound on the maximum value of the function
    virtual glm::vec3 GetMaxValue() const { return glm::vec3(1.0f, 1.0f, 1.0f); }
};

#endif
