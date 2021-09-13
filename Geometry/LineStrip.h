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

#include <Geometry/Geometry.h>
#include <iostream>
#include <vector>

class LineStrip : public Geometry {
protected:
    std::vector<glm::vec3> mJoints;

    // display information
    glm::vec3 mJointColor;
    glm::vec3 mLineColor;
    float mJointSize;
    float mLineWidth;

public:
    LineStrip(const std::vector<glm::vec3> &joints);

    virtual void Update() {}

    virtual void Initialize() {}

    virtual void Render();

    virtual const char *GetTypeName() { return typeid(LineStrip).name(); }
};
