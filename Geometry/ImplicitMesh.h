/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 * Acknowledgements for original code base:
 * - Gunnar Johansson
 * - Ken Museth
 * - Michael Bang Nielsen
 * - Ola Nilsson
 * - Andreas Soderstrom
 *
 * Code updated in the period 2017-2018 by Jochen Jankowai
 *
 *************************************************************************************************/
#ifndef __implicitmesh_h__
#define __implicitmesh_h__

#include <Geometry/Implicit.h>
#include <Geometry/SimpleMesh.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <Math/Volume.h>
#include <Util/JetColorMap.h>

class ImplicitMesh : public Implicit {

public:
  ImplicitMesh(SimpleMesh *mesh);
  virtual ~ImplicitMesh();

  virtual void Initialize();

  virtual float GetValue(float x, float y, float z) const {
    // Transform (x,y,z) to grid coordinates
    TransformW2O(x, y, z);
    x = (x - mBox.pMin[0]) / mMeshSampling;
    y = (y - mBox.pMin[1]) / mMeshSampling;
    z = (z - mBox.pMin[2]) / mMeshSampling;
    return mData->GetValue(x, y, z);
  }

  virtual void SetMeshSampling(float sampling) {
    Implicit::SetMeshSampling(sampling);
    Initialize();
  }

protected:
  SimpleMesh *mSourceMesh;
  Volume<float> *mData;

  //! Computes the closest distance from a point in space to the mesh
  float DistanceToPoint(float x, float y, float z,
                        const SimpleMesh &mesh) const;

  // Computes the closest distance between a triangle and a point in space.
  static std::pair<float, bool> DistanceSquared(const glm::vec3 &p,
                                                const glm::vec3 &v1,
                                                const glm::vec3 &v2,
                                                const glm::vec3 &v3);
};

#endif
