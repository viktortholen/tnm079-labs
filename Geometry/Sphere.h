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
#pragma once

#include <Geometry/Implicit.h>

/*!  \brief Sphere base class */
class Sphere : public virtual Implicit {
public:
    Sphere(float r, bool euclideanDistance = false);
    virtual ~Sphere();
    virtual float GetValue(float x, float y, float z) const;

protected:
    float radius2;
    bool mEuclideanDistance;
};
