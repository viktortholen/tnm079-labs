/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sˆderstrˆm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#ifndef __uniform_cubic_spline_h__
#define __uniform_cubic_spline_h__

#include <Geometry/Geometry.h>
#include <Geometry/LineStrip.h>
#include <cmath>
#include <iostream>
#include <vector>

class UniformCubicSpline : public Geometry {
protected:
  //! The coefficients dictating the shape
  std::vector<glm::vec3> mCoefficients;
  //! The control polygon is simply a LineStrip
  LineStrip mControlPolygon;

  //! Decides the length of the linear approximating segments used when drawing
  float mDt;

  int mBSplineEvaluations;

  // display information
  glm::vec3 mLineColor;
  float mLineWidth;

public:
  UniformCubicSpline(const std::vector<glm::vec3> &joints,
                     glm::vec3 lineColor = glm::vec3(0.0f, 1.0f,
                                                               0.0f),
                     float lineWidth = 2.0f, float segmentLength = 0.01f);

  /*! The BSpline value is calculated from one of the four cardinal BSpline
   * segments
   */
  float GetBSplineValue(size_t i, float t);

  /*! Evaluate the spline as the sum of the coefficients times the bsplines */
  glm::vec3 GetValue(float t);

  virtual void Update() {}
  virtual void Initialize() {}

  virtual void Render();

  virtual const char *GetTypeName() {
    return typeid(UniformCubicSpline).name();
  }
};
#endif
