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
#include <Util/ColorMap.h>
#include <glm.hpp>

//-----------------------------------------------------------------------------
glm::vec3 ColorMap::Map(float val, float low, float high) const {
  val = (val != val) ? 0 : val;
  float h = Switch1(val, low, high);
  h = glm::clamp(h, 0.f, .99999f);
  float pos = h * (mColors.size() - 1);
  float t = pos - floorf(pos);
  size_t index = (size_t)floorf(pos);
  // linear interpolation
  return mColors.at(index) * (1 - t) + mColors.at(index + 1) * t;
}
//-----------------------------------------------------------------------------
glm::vec3 ColorMap::Map(const glm::vec3 &vec, float low,
                             float high) const {
  // Map data in vec componentwise to [0,1]
  return glm::vec3(Switch1(vec[0], low, high), Switch1(vec[1], low, high),
                        Switch1(vec[2], low, high));
}
//-----------------------------------------------------------------------------
