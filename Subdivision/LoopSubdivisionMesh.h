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
#ifndef _LOOP_SUBDIVISION_MESH_
#define _LOOP_SUBDIVISION_MESH_

#include "Geometry/HalfEdgeMesh.h"
#include "Subdivision.h"

/*! \brief Subdivision mesh that implements the Loop scheme
 */
class LoopSubdivisionMesh : public Subdivision, public HalfEdgeMesh {
public:
  LoopSubdivisionMesh(const HalfEdgeMesh &m, size_t s)
      : HalfEdgeMesh(m), mNumSubDivs(s) {}
  LoopSubdivisionMesh() : mNumSubDivs(0) {}

  virtual ~LoopSubdivisionMesh() {}

  //! Subdivides the mesh uniformly one step
  virtual void Subdivide();

  virtual const char *GetTypeName() {
    return typeid(LoopSubdivisionMesh).name();
  }

  //! Return weights for interior verts
  static float Beta(size_t valence);

protected:
  //! The number of accumulated subdivisions
  size_t mNumSubDivs;

  //! Subdivides the face at faceIndex and returns a vector of faces
  virtual std::vector<std::vector<glm::vec3>> Subdivide(size_t faceIndex);

  //! Computes a new vertex, replacing a vertex in the old mesh
  virtual glm::vec3 VertexRule(size_t vertexIndex);

  //! Computes a new vertex, placed along an edge in the old mesh
  virtual glm::vec3 EdgeRule(size_t edgeIndex);
};

#endif
