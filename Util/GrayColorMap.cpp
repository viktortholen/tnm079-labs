#include "Util/GrayColorMap.h"
#include <glm.hpp>

ColorMapFactory::FactoryRegistration
GrayColorMap::mFactoryRegistration("Gray",
    new GrayColorMap());

GrayColorMap::GrayColorMap() {
    mColors.push_back(glm::vec3(0.5f, 0.5f, 0.5f));
    mColors.push_back(glm::vec3(0.5f, 0.5f, 0.5f));
}
