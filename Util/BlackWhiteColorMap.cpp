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

#include "Util/BlackWhiteColorMap.h"
#include <glm.hpp>

ColorMapFactory::FactoryRegistration
    BlackWhiteColorMap::mFactoryRegistration("Black-White",
                                             new BlackWhiteColorMap());

BlackWhiteColorMap::BlackWhiteColorMap() {
  mColors.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
  mColors.push_back(glm::vec3(1.0f, 1.0f, 1.0f));
}
