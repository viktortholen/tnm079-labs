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
#ifndef _A_LOOP_SUBDIVISION_MESH_
#define _A_LOOP_SUBDIVISION_MESH_

#include "Subdivision/LoopSubdivisionMesh.h"

/*! \brief Abstract base class for Adaptive Subdivision that implements the Loop
 * scheme
 */
class AdaptiveLoopSubdivisionMesh : public LoopSubdivisionMesh {
public:
  AdaptiveLoopSubdivisionMesh(const HalfEdgeMesh &m, size_t s)
      : LoopSubdivisionMesh(m, s) {}
  AdaptiveLoopSubdivisionMesh() {}

  virtual ~AdaptiveLoopSubdivisionMesh() {}

  //! Subdivides the mesh uniformly one step
  virtual void Subdivide();

protected:
  virtual bool Subdividable(size_t faceIndex) { return true; }

  //! Subdivides the face at faceIndex given 1 not subdividable neighbor
  virtual std::vector<std::vector<glm::vec3>> Subdivide1(size_t faceIndex);

  //! Subdivides the face at faceIndex given 2 not subdividable neighbors
  virtual std::vector<std::vector<glm::vec3>> Subdivide2(size_t faceIndex);

  //! Subdivides the face at faceIndex given 3 not subdividable neighbors
  //! (trivial)
  virtual std::vector<std::vector<glm::vec3>> Subdivide3(size_t faceIndex);

  //! Computes a new vertex, replacing a vertex in the old mesh
  virtual glm::vec3 VertexRule(size_t vertexIndex);
};

#endif
