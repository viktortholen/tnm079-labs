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

/*!  \brief Quadric base class */
class Quadric : public Implicit {
public:
    //! Initialize the quadric from matrix q
    Quadric(const glm::mat4 &q);
    virtual ~Quadric();
    //! evaluate the quadric at world coordinates x y z
    virtual float GetValue(float x, float y, float z) const;
    //! calculate the gradient at world coordinates x y z
    virtual glm::vec3 GetGradient(float x, float y, float z) const;

protected:
    //! The quadrics coefficent matrix
    glm::mat4 mQuadric;
};
