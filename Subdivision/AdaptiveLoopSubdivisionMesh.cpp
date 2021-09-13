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

#include "Subdivision/AdaptiveLoopSubdivisionMesh.h"

/*! Subdivides the mesh one step, depending on subdividability
 */
void AdaptiveLoopSubdivisionMesh::Subdivide() {
  // Create new mesh and copy all the attributes
  HalfEdgeMesh subDivMesh;
  subDivMesh.SetTransform(GetTransform());
  subDivMesh.SetName(GetName());
  subDivMesh.SetColorMap(GetColorMap());
  subDivMesh.SetWireframe(GetWireframe());
  subDivMesh.SetShowNormals(GetShowNormals());
  subDivMesh.SetOpacity(GetOpacity());
  if (IsHovering())
    subDivMesh.Hover();
  if (IsSelected())
    subDivMesh.Select();
  subDivMesh.mMinCMap = mMinCMap;
  subDivMesh.mMaxCMap = mMaxCMap;
  subDivMesh.mAutoMinMax = mAutoMinMax;

  // loop over each face and create new ones
  for (size_t i = 0; i < GetNumFaces(); i++) {

    // find neighbor faces
    size_t f1, f2, f3;
    EdgeIterator eit = GetEdgeIterator(f(i).edge);
    f1 = eit.Pair().GetEdgeFaceIndex();
    eit.Pair();
    f2 = eit.Next().Pair().GetEdgeFaceIndex();
    eit.Pair();
    f3 = eit.Next().Pair().GetEdgeFaceIndex();

    size_t numNotSubdividable =
        !Subdividable(f1) + !Subdividable(f2) + !Subdividable(f3);

    // Do not subdivide if "self" is not subdividable
    if (!Subdividable(i)) {
      numNotSubdividable = 3;
    }

    std::vector<std::vector<glm::vec3>> faces;
    switch (numNotSubdividable) {
    case 0:
      // normal subdivision (from LoopSubdivisionMesh)
      faces = LoopSubdivisionMesh::Subdivide(i);
      break;
    case 1:
      // special case 1
      faces = Subdivide1(i);
      break;
    case 2:
      // special case 2
      faces = Subdivide2(i);
      break;
    case 3:
      // trivial case, no subdivision, same as if subdividable(fi) == false
      faces = Subdivide3(i);
      break;
    }

    // add the faces (if any) to subDivMesh
    for (size_t j = 0; j < faces.size(); j++) {
      subDivMesh.AddFace(faces.at(j));
    }
  }

  // Assign the new mesh
  *this = AdaptiveLoopSubdivisionMesh(subDivMesh, ++mNumSubDivs);
  Update();
}

/*! Computes a new vertex, replacing a vertex in the old mesh
  If any of the neighboring faces have subdividability == false
  then we should not move this vertex, else use the rules from loop subdivision
*/
glm::vec3 AdaptiveLoopSubdivisionMesh::VertexRule(size_t vertexIndex) {
  // Get the current vertex
  glm::vec3 vtx = v(vertexIndex).pos;

  // Get face neighborhood
  std::vector<size_t> nb = FindNeighborFaces(vertexIndex);

  for (size_t i = 0; i < nb.size(); i++) {
    if (Subdividable(nb.at(i)) == false) {
      // don't move position
      return vtx;
    }
  }

  // else return VertexRule for LSMesh
  return LoopSubdivisionMesh::VertexRule(vertexIndex);
}

/*! Subdivides the face at faceIndex given 1 not subdividable neighbor
 */
std::vector<std::vector<glm::vec3>>
AdaptiveLoopSubdivisionMesh::Subdivide1(size_t faceIndex) {

  // We know we have one false face

  // 1. Start by orienting the triangle so that we always have the same case
  // We find the start edge as the edge who shares face with the
  // _not_subdividable_ neighbor
  HalfEdgeMesh::EdgeIterator eit = GetEdgeIterator(f(faceIndex).edge);
  size_t num = 0;
  while (Subdividable(eit.Pair().GetEdgeFaceIndex())) {
    eit.Pair().Next();
    assert(num++ < 3);
  }
  // go back to inner edge
  eit.Pair();

  // 2. Now find the vertices
  std::vector<std::vector<glm::vec3>> faces;

  // corner vertices, labeled from current edge's origin vertex
  glm::vec3 v1 = VertexRule(eit.GetEdgeVertexIndex());
  glm::vec3 v2 = VertexRule(eit.Next().GetEdgeVertexIndex());
  glm::vec3 v3 = VertexRule(eit.Next().GetEdgeVertexIndex());

  // edge vertices, labeled from start edge
  glm::vec3 e2v = EdgeRule(eit.Next().Next().GetEdgeIndex());
  glm::vec3 e3v = EdgeRule(eit.Next().GetEdgeIndex());

  // 3. Create the 3 faces and push them on the vector
  std::vector<glm::vec3> face;
  face.push_back(v1);
  face.push_back(v2);
  face.push_back(e2v);
  faces.push_back(face);
  face.clear();

  face.push_back(e2v);
  face.push_back(v3);
  face.push_back(e3v);
  faces.push_back(face);
  face.clear();

  face.push_back(e3v);
  face.push_back(v1);
  face.push_back(e2v);
  faces.push_back(face);
  face.clear();

  // 4. Return
  return faces;
}

/*! Subdivides the face at faceIndex given 2 not subdividable neighbors
 */
std::vector<std::vector<glm::vec3>>
AdaptiveLoopSubdivisionMesh::Subdivide2(size_t faceIndex) {

  // We know we have otwo false faces

  // 1. Start by orienting the triangle so that we always have the same case
  // We find the start edge as the edge who shares face with the _subdividable_
  // neighbor
  EdgeIterator eit = GetEdgeIterator(f(faceIndex).edge);
  size_t num = 0;
  while (!Subdividable(eit.Pair().GetEdgeFaceIndex())) {
    eit.Pair().Next();
    assert(num++ < 3);
  }
  // go back to inner edge
  eit.Pair();

  // 2. Now find the vertices
  std::vector<std::vector<glm::vec3>> faces;

  // corner vertices, labeled from current edge's origin vertex
  glm::vec3 v1 = VertexRule(eit.GetEdgeVertexIndex());
  glm::vec3 v2 = VertexRule(eit.Next().GetEdgeVertexIndex());
  glm::vec3 v3 = VertexRule(eit.Next().GetEdgeVertexIndex());

  // edge vertices, labeled from start edge
  glm::vec3 e1v = EdgeRule(eit.Next().GetEdgeIndex());

  // 3. Create the 2 faces and push them on the vector
  std::vector<glm::vec3> face;
  face.push_back(v1);
  face.push_back(e1v);
  face.push_back(v3);
  faces.push_back(face);
  face.clear();

  face.push_back(e1v);
  face.push_back(v2);
  face.push_back(v3);
  faces.push_back(face);
  face.clear();

  // 4. Return
  return faces;
}

/*! Subdivides the face at faceIndex given 3 not subdividable neighbors
  or if subdividable for this face is false
*/
std::vector<std::vector<glm::vec3>>
AdaptiveLoopSubdivisionMesh::Subdivide3(size_t faceIndex) {

  EdgeIterator eit = GetEdgeIterator(f(faceIndex).edge);
  glm::vec3 v1 = VertexRule(eit.GetEdgeVertexIndex());
  glm::vec3 v2 = VertexRule(eit.Next().GetEdgeVertexIndex());
  glm::vec3 v3 = VertexRule(eit.Next().GetEdgeVertexIndex());
  std::vector<glm::vec3> face;
  face.push_back(v1);
  face.push_back(v2);
  face.push_back(v3);
  std::vector<std::vector<glm::vec3>> faces;
  faces.push_back(face);
  return faces;
}
