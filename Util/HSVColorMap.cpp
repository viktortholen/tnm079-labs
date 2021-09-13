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

#include "Util/HSVColorMap.h"
#include <glm.hpp>

ColorMapFactory::FactoryRegistration
    HSVColorMap::mFactoryRegistration("HSV", new HSVColorMap());

HSVColorMap::HSVColorMap() {
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.093750f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.187500f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.281250f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.375000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.468750f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.562500f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.656250f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.750000f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.843750f, 0.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.937500f, 0.000000f));
  mColors.push_back(glm::vec3(0.968750f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.875000f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.781250f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.687500f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.593750f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.500000f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.406250f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.312500f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.218750f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.125000f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.031250f, 1.000000f, 0.000000f));
  mColors.push_back(glm::vec3(0.000000f, 1.000000f, 0.062500f));
  mColors.push_back(glm::vec3(0.000000f, 1.000000f, 0.156250f));
  mColors.push_back(glm::vec3(0.000000f, 1.000000f, 0.250000f));
  mColors.push_back(glm::vec3(0.000000f, 1.000000f, 0.343750f));
  mColors.push_back(glm::vec3(0.000000f, 1.000000f, 0.437500f));
  mColors.push_back(glm::vec3(0.000000f, 1.000000f, 0.531250f));
  mColors.push_back(glm::vec3(0.000000f, 1.000000f, 0.625000f));
  mColors.push_back(glm::vec3(0.000000f, 1.000000f, 0.718750f));
  mColors.push_back(glm::vec3(0.000000f, 1.000000f, 0.812500f));
  mColors.push_back(glm::vec3(0.000000f, 1.000000f, 0.906250f));
  mColors.push_back(glm::vec3(0.000000f, 1.000000f, 1.000000f));
  mColors.push_back(glm::vec3(0.000000f, 0.906250f, 1.000000f));
  mColors.push_back(glm::vec3(0.000000f, 0.812500f, 1.000000f));
  mColors.push_back(glm::vec3(0.000000f, 0.718750f, 1.000000f));
  mColors.push_back(glm::vec3(0.000000f, 0.625000f, 1.000000f));
  mColors.push_back(glm::vec3(0.000000f, 0.531250f, 1.000000f));
  mColors.push_back(glm::vec3(0.000000f, 0.437500f, 1.000000f));
  mColors.push_back(glm::vec3(0.000000f, 0.343750f, 1.000000f));
  mColors.push_back(glm::vec3(0.000000f, 0.250000f, 1.000000f));
  mColors.push_back(glm::vec3(0.000000f, 0.156250f, 1.000000f));
  mColors.push_back(glm::vec3(0.000000f, 0.062500f, 1.000000f));
  mColors.push_back(glm::vec3(0.031250f, 0.000000f, 1.000000f));
  mColors.push_back(glm::vec3(0.125000f, 0.000000f, 1.000000f));
  mColors.push_back(glm::vec3(0.218750f, 0.000000f, 1.000000f));
  mColors.push_back(glm::vec3(0.312500f, 0.000000f, 1.000000f));
  mColors.push_back(glm::vec3(0.406250f, 0.000000f, 1.000000f));
  mColors.push_back(glm::vec3(0.500000f, 0.000000f, 1.000000f));
  mColors.push_back(glm::vec3(0.593750f, 0.000000f, 1.000000f));
  mColors.push_back(glm::vec3(0.687500f, 0.000000f, 1.000000f));
  mColors.push_back(glm::vec3(0.781250f, 0.000000f, 1.000000f));
  mColors.push_back(glm::vec3(0.875000f, 0.000000f, 1.000000f));
  mColors.push_back(glm::vec3(0.968750f, 0.000000f, 1.000000f));
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.937500f));
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.843750f));
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.750000f));
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.656250f));
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.562500f));
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.468750f));
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.375000f));
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.281250f));
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.187500f));
  mColors.push_back(glm::vec3(1.000000f, 0.000000f, 0.093750f));
}
