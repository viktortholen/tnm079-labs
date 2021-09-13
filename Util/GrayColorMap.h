#pragma once

#include "Util/ColorMapFactory.h"

class GrayColorMap : public ColorMap {
public:
    GrayColorMap();

protected:
    static ColorMapFactory::FactoryRegistration mFactoryRegistration;
};
