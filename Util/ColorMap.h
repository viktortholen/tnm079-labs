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
#ifndef __COLOR_MAP_H__
#define __COLOR_MAP_H__

#include "Util/Util.h"
#include <vector>

class ColorMap {
public:
  ColorMap() {}
  virtual ~ColorMap() {}

  virtual glm::vec3 Map(float val, float low, float high) const;
  virtual glm::vec3 Map(const glm::vec3 &val, float low,
                             float high) const;

protected:
  // A vector containing the colors to be interpolated
  std::vector<glm::vec3 > mColors;
};

#endif
