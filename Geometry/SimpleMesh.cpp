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
#include <Geometry/SimpleMesh.h>
#include <Util/ColorMap.h>
#include <glm.hpp>
#include <gtc/type_ptr.hpp>

//-----------------------------------------------------------------------------
SimpleMesh::SimpleMesh() {}

//-----------------------------------------------------------------------------
SimpleMesh::~SimpleMesh() {}

//-----------------------------------------------------------------------------
bool SimpleMesh::AddFace(const std::vector<glm::vec3 > &verts) {
  const auto ind1 = AddVertex(verts.at(0));
  const auto ind2 = AddVertex(verts.at(1));
  const auto ind3 = AddVertex(verts.at(2));

  Face tri(ind1, ind2, ind3);
  mFaces.push_back(tri);
  // Compute and assign a normal
  mFaces.back().normal = FaceNormal(mFaces.size() - 1);

  return true;
}

//-----------------------------------------------------------------------------
size_t SimpleMesh::AddVertex(const glm::vec3 &v) {
  std::map<glm::vec3, size_t>::iterator it = mUniqueVerts.find(v);
  if (it != mUniqueVerts.end()) {
    const auto indx = (*it).second; // (*it).second == it->second, which to me was not as
                         // clear...
    return indx;
  }

  const auto indx = mVerts.size();
  mUniqueVerts[v] = indx; // op. [ ] constructs a new entry in map
  Vertex vert;
  vert.pos = v;
  mVerts.push_back(vert);

  return indx;
}

//-----------------------------------------------------------------------------
glm::vec3 SimpleMesh::FaceNormal(size_t faceIndex) const {
  const Face &tri = mFaces.at(faceIndex);
  glm::vec3 e1 = mVerts.at(tri.v2).pos - mVerts.at(tri.v1).pos;
  glm::vec3 e2 = mVerts.at(tri.v3).pos - mVerts.at(tri.v1).pos;
  return glm::normalize(glm::cross(e1, e2));
}
//-----------------------------------------------------------------------------
glm::vec3 SimpleMesh::VertexNormal(size_t vertexIndex) const {

  std::vector<size_t> neighborFaces = FindNeighborFaces(vertexIndex);

  glm::vec3 n(0, 0, 0);

  const size_t numCandidates = neighborFaces.size();
  for (size_t i = 0; i < numCandidates; i++) {
    const Face &triangle = mFaces.at(neighborFaces.at(i));

    // NB Assumes face normals already calculated
    n += triangle.normal;
  }
  n=glm::normalize(n);
  return n;
}
//-----------------------------------------------------------------------------
float SimpleMesh::VertexCurvature(size_t vertexIndex) const {

  std::vector<size_t> oneRing = FindNeighborVertices(vertexIndex);
  assert(oneRing.size() != 0);

  size_t curr, next;
  const glm::vec3 &vi = mVerts.at(vertexIndex).pos;
  float angleSum = 0;
  float area = 0;
  for (size_t i = 0; i < oneRing.size(); i++) {
    // connections
    curr = oneRing.at(i);
    if (i < oneRing.size() - 1)
      next = oneRing.at(i + 1);
    else
      next = oneRing.front();

    // find vertices in 1-ring according to figure 5 in lab text
    // next - beta
    const glm::vec3 &nextPos = mVerts.at(next).pos;
    const glm::vec3 &vj = mVerts.at(curr).pos;

    // compute angle and area
    angleSum += acos(glm::dot(vj - vi , nextPos - vi) /
                     (glm::length(vj - vi) * glm::length(nextPos - vi)));
    area += glm::length(glm::cross(vi - vj, nextPos - vj)) * 0.5f;
  }
  return (2.0f * static_cast<float>(M_PI) - angleSum) / area;
}

float SimpleMesh::FaceCurvature(size_t faceIndex) const {
  // NB Assumes vertex curvature already computed
  const Face &tri = mFaces.at(faceIndex);
  return (mVerts.at(tri.v1).curvature + mVerts.at(tri.v2).curvature +
          mVerts.at(tri.v3).curvature) /
         3.0f;
}

/*! Loops over the neighborhood of a vertex and collects all the vertices sorted
 * counter clockwise. \param [in] vertexIndex  the index to vertex, size_t
 * \return a vector containing the indices to all the found vertices.
 */
std::vector<size_t> SimpleMesh::FindNeighborVertices(size_t vertexIndex) const {
  std::vector<size_t> neighborFaces = FindNeighborFaces(vertexIndex);
  std::vector<size_t> oneRing;
  size_t currVert;

  // pick next counter clock wise vert
  if (mFaces.at(neighborFaces.at(0)).v1 == vertexIndex)
    currVert = mFaces.at(neighborFaces.at(0)).v2;
  if (mFaces.at(neighborFaces.at(0)).v2 == vertexIndex)
    currVert = mFaces.at(neighborFaces.at(0)).v3;
  if (mFaces.at(neighborFaces.at(0)).v3 == vertexIndex)
    currVert = mFaces.at(neighborFaces.at(0)).v1;
  oneRing.push_back(currVert);

  // collect one ring vertices
  for (size_t i = 0; i < neighborFaces.size() - 1; i++) {
    if (mFaces.at(neighborFaces.at(i)).v1 == currVert)
      currVert = mFaces.at(neighborFaces.at(i)).v2;
    else if (mFaces.at(neighborFaces.at(i)).v2 == currVert)
      currVert = mFaces.at(neighborFaces.at(i)).v3;
    else if (mFaces.at(neighborFaces.at(i)).v3 == currVert)
      currVert = mFaces.at(neighborFaces.at(i)).v1;
    oneRing.push_back(currVert);
  }
  return oneRing;
}

/*! Loops over the neighborhood of a vertex and collects all the faces sorted
 * counter clockwise. \param [in] vertexIndex  the index to vertex, size_t
 * \return a vector containing the indices to all the found faces.
 */
std::vector<size_t> SimpleMesh::FindNeighborFaces(size_t vertexIndex) const {
  std::vector<size_t> foundFaces;

  // Find other triangles that include this vertex
  std::vector<Face>::const_iterator iter = mFaces.begin();
  std::vector<Face>::const_iterator iend = mFaces.end();
  while (iter != iend) {
    if ((*iter).v1 == vertexIndex || (*iter).v2 == vertexIndex ||
        (*iter).v3 == vertexIndex) {
      foundFaces.push_back(iter - mFaces.begin());
    }
    iter++;
  }

  // Pick prev vertex
  size_t currVertex = 0;
  const size_t v1 = mFaces.at(foundFaces.at(0)).v1;
  const size_t v2 = mFaces.at(foundFaces.at(0)).v2;
  const size_t v3 = mFaces.at(foundFaces.at(0)).v3;
  if (vertexIndex == v1)
    currVertex = v3;
  else if (vertexIndex == v2)
    currVertex = v1;
  else if (vertexIndex == v3)
    currVertex = v2;

  for (size_t i = 1; i < foundFaces.size() - 1; i++) {
    for (size_t j = i; j < foundFaces.size(); j++) {
      if (mFaces.at(foundFaces.at(j)).v1 == currVertex) {
        // pick the next vert
        currVertex = mFaces.at(foundFaces.at(j)).v2;
        // and swap
        std::swap(foundFaces.at(i), foundFaces.at(j));
        break;
      }
      if (mFaces.at(foundFaces.at(j)).v2 == currVertex) {
        // pick the next vert
        currVertex = mFaces.at(foundFaces.at(j)).v3;
        // and swap
        std::swap(foundFaces.at(i), foundFaces.at(j));
        break;
      }
      if (mFaces.at(foundFaces.at(j)).v3 == currVertex) {
        // pick the next vert
        currVertex = mFaces.at(foundFaces.at(j)).v1;
        // and swap
        std::swap(foundFaces.at(i), foundFaces.at(j));
        break;
      }
    }
  }

  return foundFaces;
}

//-----------------------------------------------------------------------------
void SimpleMesh::Initialize() {
  // Calculate and store all differentials and area

  // First update all face normals and triangle areas
  for (size_t i = 0; i < mFaces.size(); i++) {
    mFaces.at(i).normal = FaceNormal(i);
  }
  // Then update all vertex normals and curvature
  for (size_t i = 0; i < mVerts.size(); i++) {
    // Vertex normals are just weighted averages
    mVerts.at(i).normal = VertexNormal(i);
  }

  // Then update vertex curvature
  for (size_t i = 0; i < mVerts.size(); i++) {
    mVerts.at(i).curvature = VertexCurvature(i);
    // std::cerr <<   mVerts.at(i).curvature << "\n";
  }

  // Finally update face curvature
  for (size_t i = 0; i < mFaces.size(); i++) {
    mFaces.at(i).curvature = FaceCurvature(i);
  }
}

//-----------------------------------------------------------------------------
void SimpleMesh::Update() {

  // Update vertex and face colors
  if (mVisualizationMode == CurvatureVertex) {
    std::vector<Vertex>::iterator iter = mVerts.begin();
    std::vector<Vertex>::iterator iend = mVerts.end();
    float minCurvature = (std::numeric_limits<float>::max)();
    float maxCurvature = -(std::numeric_limits<float>::max)();
    while (iter != iend) {
      if (minCurvature > (*iter).curvature)
        minCurvature = (*iter).curvature;
      if (maxCurvature < (*iter).curvature)
        maxCurvature = (*iter).curvature;
      iter++;
    }
    if (mAutoMinMax == false) {
      std::cerr << "Mapping color based on vertex curvature with range ["
                << mMinCMap << "," << mMaxCMap << "]" << std::endl;
      minCurvature = mMinCMap;
      maxCurvature = mMaxCMap;
    } else {
      std::cerr
          << "Automatic mapping of color based on vertex curvature with range ["
          << minCurvature << "," << maxCurvature << "]" << std::endl;
      mMinCMap = minCurvature;
      mMaxCMap = maxCurvature;
    }
    iter = mVerts.begin();
    while (iter != iend) {
      (*iter).color =
          mColorMap->Map((*iter).curvature, minCurvature, maxCurvature);
      iter++;
    }
  } else if (mVisualizationMode == CurvatureFace) {
    std::vector<Face>::iterator iter = mFaces.begin();
    std::vector<Face>::iterator iend = mFaces.end();
    float minCurvature = (std::numeric_limits<float>::max)();
    float maxCurvature = -(std::numeric_limits<float>::max)();
    while (iter != iend) {
      if (minCurvature > (*iter).curvature)
        minCurvature = (*iter).curvature;
      if (maxCurvature < (*iter).curvature)
        maxCurvature = (*iter).curvature;
      iter++;
    }
    std::cerr << "Mapping color based on face curvature with range ["
              << minCurvature << "," << maxCurvature << "]" << std::endl;
    iter = mFaces.begin();
    while (iter != iend) {
      (*iter).color =
          mColorMap->Map((*iter).curvature, minCurvature, maxCurvature);
      iter++;
    }
  }
}

size_t SimpleMesh::Genus() const {
  std::set<MyEdge> uniqueEdges;
  for (size_t i = 0; i < mFaces.size(); i++) {
    uniqueEdges.insert(MyEdge(mFaces[i].v1, mFaces[i].v2));
    uniqueEdges.insert(MyEdge(mFaces[i].v1, mFaces[i].v3));
    uniqueEdges.insert(MyEdge(mFaces[i].v2, mFaces[i].v3));
  }
  auto E = uniqueEdges.size();
  auto V = mVerts.size();
  auto F = mFaces.size();

  std::cerr << "Number of edges: " << E << ", F: " << F << ", V: " << V << "\n";
  return -(static_cast<int64_t>(V) - static_cast<int64_t>(E) +
           static_cast<int64_t>(F) - 2) /
         2;
}

void SimpleMesh::Dilate(float amount) {
  std::vector<Vertex>::iterator iter = mVerts.begin();
  std::vector<Vertex>::iterator iend = mVerts.end();
  while (iter != iend) {
    (*iter).pos += amount * (*iter).normal;
    iter++;
  }

  Initialize();
  Update();
}

void SimpleMesh::Erode(float amount) {
  std::vector<Vertex>::iterator iter = mVerts.begin();
  std::vector<Vertex>::iterator iend = mVerts.end();
  while (iter != iend) {
    (*iter).pos -= amount * (*iter).normal;
    iter++;
  }

  Initialize();
  Update();
}

void SimpleMesh::Smooth(float amount) {
  std::vector<Vertex>::iterator iter = mVerts.begin();
  std::vector<Vertex>::iterator iend = mVerts.end();
  while (iter != iend) {
    (*iter).pos -= amount * (*iter).normal * (*iter).curvature;
    iter++;
  }

  Initialize();
  Update();
}

void SimpleMesh::Render() {
  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // Apply transform
  glPushMatrix(); // Push modelview matrix onto stack

  // Convert transform-matrix to format matching GL matrix format
  // Load transform into modelview matrix
  glMultMatrixf(glm::value_ptr(mTransform));

  if (mWireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // Draw geometry
  glBegin(GL_TRIANGLES);
  const auto numFaces = mFaces.size();
  for (size_t i = 0; i < numFaces; i++) {
    Face &triangle = mFaces[i];

    glm::vec3 &p0 = mVerts[triangle.v1].pos;
    glm::vec3 &p1 = mVerts[triangle.v2].pos;
    glm::vec3 &p2 = mVerts[triangle.v3].pos;

    if (mVisualizationMode == CurvatureVertex) {
      glm::vec3 &c1 = mVerts.at(triangle.v1).color;
      glColor4f(c1[0], c1[1], c1[2], mOpacity);
      glNormal3fv(glm::value_ptr(mVerts.at(triangle.v1).normal));
      glVertex3fv(glm::value_ptr(p0));

      glm::vec3 &c2 = mVerts.at(triangle.v2).color;
      glColor4f(c2[0], c2[1], c2[2], mOpacity);
      glNormal3fv(glm::value_ptr(mVerts.at(triangle.v2).normal));
      glVertex3fv(glm::value_ptr(p1));

      glm::vec3 &c3 = mVerts.at(triangle.v3).color;
      glColor4f(c3[0], c3[1], c3[2], mOpacity);
      glNormal3fv(glm::value_ptr(mVerts.at(triangle.v3).normal));
      glVertex3fv(glm::value_ptr(p2));
    } else {
      glm::vec3 &color = triangle.color;
      glColor4f(color[0], color[1], color[2], mOpacity);
      glNormal3fv(glm::value_ptr(triangle.normal));

      glVertex3fv(glm::value_ptr(p0));
      glVertex3fv(glm::value_ptr(p1));
      glVertex3fv(glm::value_ptr(p2));
    }
  }
  glEnd();

  // Mesh normals by courtesy of Richard Khoury
  if (mShowNormals) {
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    for (size_t i = 0; i < mFaces.size(); i++) {

      Face &face = mFaces.at(i);

      Vertex &v1 = mVerts.at(face.v1);
      Vertex &v2 = mVerts.at(face.v2);
      Vertex &v3 = mVerts.at(face.v3);

      glm::vec3 faceStart = (v1.pos + v2.pos + v3.pos) / 3.0f;
      glm::vec3 faceEnd = faceStart + face.normal * 0.1f;

      glColor3f(1, 0, 0); // Red for face normal
      glVertex3fv(glm::value_ptr(faceStart));
      glVertex3fv(glm::value_ptr(faceEnd));

      glColor3f(0, 1, 0); // Vertex normals in Green
      glVertex3fv(glm::value_ptr(v1.pos));
      glVertex3fv(glm::value_ptr(v1.pos + v1.normal * 0.1f));
      glVertex3fv(glm::value_ptr(v2.pos));
      glVertex3fv(glm::value_ptr(v2.pos + v2.normal * 0.1f));
      glVertex3fv(glm::value_ptr(v3.pos));
      glVertex3fv(glm::value_ptr(v3.pos + v3.normal * 0.1f));
    }
    glEnd();
    glEnable(GL_LIGHTING);
  }

  if (mWireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // Restore modelview matrix
  glPopMatrix();

  GLObject::Render();
}
