#pragma once

#include "Geometry/Geometry.h"
#include "Math/Function3D.h"
#include "Util/ColorMap.h"

class ScalarCutPlane : public Geometry {
protected:
  float mDx;
  const Function3D<float> *mFunction;

  GLuint mTextureID;
  GLuint mWidth, mHeight;

public:
  ScalarCutPlane(const std::string &name, float dx,
                 const Function3D<float> *function);

  virtual const char *GetTypeName() { return typeid(ScalarCutPlane).name(); }

  virtual void Render();

  virtual void Initialize() {}

  virtual void Update();

  virtual void SetTransform(const glm::mat4&transform);

  virtual void SetColorMap(ColorMap *colormap) {
    GLObject::SetColorMap(colormap);
    Update();
  }
};
