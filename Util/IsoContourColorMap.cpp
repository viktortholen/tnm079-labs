/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Söderström (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/

#include "Util/IsoContourColorMap.h"

ColorMapFactory::FactoryRegistration
    IsoContourColorMap::mFactoryRegistration("Iso contour",
                                             new IsoContourColorMap());

IsoContourColorMap::IsoContourColorMap() {
  mColors.push_back(glm::vec3(0.0f, 1.0f, 0.0f));
  mColors.push_back(glm::vec3(1.0f, 0.0f, 0.0f));
}

glm::vec3 IsoContourColorMap::Map(float val, float low, float high) const {
  // Take absolute value
  if (val < 0.0f)
    val = -val;

  // Compute fraction to do wrap-around
  float fraction = val * 10.0f - (size_t)(val * 10.0f);

  // Do color mapping
  return ColorMap::Map(fraction, 0.0f, 1.0f);
}
