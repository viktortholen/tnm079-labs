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

#include <Util/HotColorMap.h>
#include <glm.hpp>

ColorMapFactory::FactoryRegistration
    HotColorMap::mFactoryRegistration("Hot", new HotColorMap());

HotColorMap::HotColorMap() {
  mColors.push_back(glm::vec3(0.041667f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.083333f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.125000f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.166667f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.208333f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.250000f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.291667f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.333333f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.375000f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.416667f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.458333f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.500000f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.541667f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.583333f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.625000f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.666667f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.708333f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.750000f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.791667f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.833333f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.875000f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.916667f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.958333f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.041667f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.083333f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.125000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.166667f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.208333f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.250000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.291667f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.333333f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.375000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.416667f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.458333f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.500000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.541667f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.583333f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.625000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.666667f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.708333f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.750000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.791667f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.833333f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.875000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.916667f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.958333f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.062500f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.125000f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.187500f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.250000f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.312500f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.375000f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.437500f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.500000f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.562500f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.625000f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.687500f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.750000f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.812500f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.875000f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 0.937500f));
  mColors.push_back(glm::vec3(1.000000f, 1.000000f, 1.000000f));
}
