/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#ifndef __VortexVectorField_h__
#define __VortexVectorField_h__

#include "Math/Function3D.h"
#include <cmath>

class VortexVectorField : public Function3D<glm::vec3> {
public:
  VortexVectorField() {
    rScale = 0.1f;
    thetaScale = 0.3f;
    phiScale = 0.3f;
  }
  ~VortexVectorField(void) {}

  glm::vec3 GetMaxValue() const {
    return glm::vec3(1.0f, 1.0f, 1.0f);
  }

  glm::vec3 GetValue(float xf, float yf, float zf) const {
    glm::vec3 v;
    // 3D VORTEX FIELD
    float r = sqrt(xf * xf + yf * yf + zf * zf);

    if (r != 0) {
      float theta = atan2(yf, xf);
      float phi = acos(zf / r);

      v[0] = rScale * cos(theta) * sin(phi) + thetaScale * (-sin(theta)) +
             phiScale * cos(theta) * cos(phi);

      v[1] = rScale * sin(theta) * sin(phi) + thetaScale * (cos(theta)) +
             phiScale * sin(theta) * cos(phi);

      v[2] = rScale * cos(phi) + thetaScale * (0) + phiScale * (-sin(phi));
    } else {
      v[0] = 0;
      v[1] = 0;
      v[2] = 0;
    }
    return v;
  }

protected:
  float t;
  float rScale, thetaScale, phiScale;
};

#endif
