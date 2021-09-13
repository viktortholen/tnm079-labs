#include <Geometry/HalfEdgeMesh.h>
#include <gtc/type_ptr.hpp>
#include <iterator>

HalfEdgeMesh::HalfEdgeMesh() {}

HalfEdgeMesh::~HalfEdgeMesh() {}

/*! \lab1 Implement the addFace */
/*!
 * \param[in] v1 vertex 1, glm::vec3
 * \param[in] v2 vertex 2, glm::vec3
 * \param[in] v3 vertex 3, glm::vec3
 */
bool HalfEdgeMesh::AddFace(const std::vector<glm::vec3>& verts) {
    // Add your code here

    // Add the vertices of the face/triangle
    const auto ind1 = AddVertex(verts.at(0));
    const auto ind2 = AddVertex(verts.at(1));
    const auto ind3 = AddVertex(verts.at(2));

    // Add all half-edge pairs
    const auto edge1 = AddHalfEdgePair(ind1, ind2);
    const auto edge2 = AddHalfEdgePair(ind2, ind3);
    const auto edge3 = AddHalfEdgePair(ind3, ind1);

    // Connect inner ring
    e(edge1.first).next = edge2.first;
    e(edge2.first).next = edge3.first;
    e(edge3.first).next = edge1.first;

    e(edge1.first).prev = edge3.first;
    e(edge2.first).prev = edge1.first;
    e(edge3.first).prev = edge2.first;

    // Finally, create the face, don't forget to set the normal (which should be
    // normalized)
    Face face;
    face.edge = edge1.first;
    mFaces.push_back(face);

    size_t faceIndex = mFaces.size() - 1;
    f(faceIndex).normal = FaceNormal(faceIndex);

    // All half-edges share the same left face (previously added)
    e(edge1.first).face = faceIndex;
    e(edge2.first).face = faceIndex;
    e(edge3.first).face = faceIndex;

    // Optionally, track the (outer) boundary half-edges
    // to represent non-closed surfaces
    return true;
}

/*!
 * \param [in] v the vertex to add, glm::vec3
 * \return the index to the vertex
 */
size_t HalfEdgeMesh::AddVertex(const glm::vec3& v) {
    std::map<glm::vec3, size_t>::iterator it = mUniqueVerts.find(v);
    if (it != mUniqueVerts.end()) {
        return (*it).second;  // get the index of the already existing vertex
    }

    const auto indx = GetNumVerts();
    mUniqueVerts[v] = indx;  // op. [ ] constructs a new entry in map
    Vertex vert;
    vert.pos = v;
    mVerts.push_back(vert);  // add it to the vertex list

    return indx;
}

/*!
 * Inserts a half edge pair between HalfEdgeMesh::Vertex pointed to by v1 and
 * v2. The first HalfEdgeMesh::HalfEdge (v1->v2) is the inner one, and the
 * second (v2->v1) is the outer.
 * \param [in] v1 size_t index of vertex 1
 * \param [in] v2 size_t index of vertex 2
 * \return a pair the indices to the half-edges
 */
std::pair<size_t, size_t> HalfEdgeMesh::AddHalfEdgePair(size_t v1, size_t v2) {
    std::map<OrderedPair, size_t>::iterator it = mUniqueEdgePairs.find(OrderedPair(v1, v2));
    if (it != mUniqueEdgePairs.end()) {
        auto indx1 = it->second;
        auto indx2 = e(it->second).pair;
        if (v1 != e(indx1).vert) {
            std::swap(indx1, indx2);  // sort correctly
        }
        return {indx1, indx2};
    }

    // If not found, calculate both half-edges indices
    const auto indx1 = mEdges.size();
    const auto indx2 = indx1 + 1;

    // Create edges and set pair index
    HalfEdge edge1, edge2;
    edge1.pair = indx2;
    edge2.pair = indx1;

    // Connect the edges to the verts
    edge1.vert = v1;
    edge2.vert = v2;

    // Connect the verts to the edges
    v(v1).edge = indx1;
    v(v2).edge = indx2;

    // Store the edges in mEdges
    mEdges.push_back(edge1);
    mEdges.push_back(edge2);

    // Store the first edge in the map as an OrderedPair
    OrderedPair op(v1, v2);
    mUniqueEdgePairs[op] = indx1;  // op. [ ] constructs a new entry in map, ordering not important
    // sorting done when retrieving

    return {indx1, indx2};
}

/*! \lab1 HalfEdgeMesh Implement the MergeAdjacentBoundaryEdge */
/*!
 * Merges the outer UNINITIALIZED/BORDER to an already set inner half-edge.
 * \param [in] indx the index of the INNER half-edge, size_t
 */
void HalfEdgeMesh::MergeOuterBoundaryEdge(size_t innerEdge) {
    // Add your code here
    // 1. Merge first loop (around innerEdge->vert)
    // 2. Find leftmost edge, last edge counter clock-wise
    // 3. Test if there's anything to merge
    // 3a. If so merge the gap
    // 3b. And set border flags
    // 4. Merge second loop (around innerEdge->pair->vert)
}

/*! Proceeds to check if the mesh is valid. All indices are inspected and
 * checked to see that they are initialized. The method checks: mEdges, mFaces
 * and mVerts. Also checks to see if all verts have a neighborhood using the
 * findNeighbourFaces method.
 */
void HalfEdgeMesh::Validate() {
    std::vector<HalfEdge>::iterator iterEdge = mEdges.begin();
    std::vector<HalfEdge>::iterator iterEdgeEnd = mEdges.end();
    while (iterEdge != iterEdgeEnd) {
        if ((*iterEdge).face == EdgeState::Uninitialized ||
            (*iterEdge).next == EdgeState::Uninitialized ||
            (*iterEdge).pair == EdgeState::Uninitialized ||
            (*iterEdge).prev == EdgeState::Uninitialized ||
            (*iterEdge).vert == EdgeState::Uninitialized)
            std::cerr << "HalfEdge " << iterEdge - mEdges.begin() << " not properly initialized"
                      << std::endl;

        iterEdge++;
    }
    std::cerr << "Done with edge check (checked " << GetNumEdges() << " edges)" << std::endl;

    std::vector<Face>::iterator iterTri = mFaces.begin();
    std::vector<Face>::iterator iterTriEnd = mFaces.end();
    while (iterTri != iterTriEnd) {
        if ((*iterTri).edge == EdgeState::Uninitialized)
            std::cerr << "Tri " << iterTri - mFaces.begin() << " not properly initialized"
                      << std::endl;

        iterTri++;
    }
    std::cerr << "Done with face check (checked " << GetNumFaces() << " faces)" << std::endl;

    std::vector<Vertex>::iterator iterVertex = mVerts.begin();
    std::vector<Vertex>::iterator iterVertexEnd = mVerts.end();
    while (iterVertex != iterVertexEnd) {
        if ((*iterVertex).edge == EdgeState::Uninitialized)
            std::cerr << "Vertex " << iterVertex - mVerts.begin() << " not properly initialized"
                      << std::endl;

        iterVertex++;
    }
    std::cerr << "Done with vertex check (checked " << GetNumVerts() << " vertices)" << std::endl;

    std::cerr << "Looping through triangle neighborhood of each vertex... ";
    iterVertex = mVerts.begin();
    iterVertexEnd = mVerts.end();
    int emptyCount = 0;
    std::vector<size_t> problemVerts;
    while (iterVertex != iterVertexEnd) {
        std::vector<size_t> foundFaces = FindNeighborFaces(iterVertex - mVerts.begin());
        std::vector<size_t> foundVerts = FindNeighborVertices(iterVertex - mVerts.begin());
        if (foundFaces.empty() || foundVerts.empty()) emptyCount++;
        std::set<size_t> uniqueFaces(foundFaces.begin(), foundFaces.end());
        std::set<size_t> uniqueVerts(foundVerts.begin(), foundVerts.end());
        if (foundFaces.size() != uniqueFaces.size() || foundVerts.size() != uniqueVerts.size())
            problemVerts.push_back(iterVertex - mVerts.begin());
        iterVertex++;
    }
    std::cerr << std::endl << "Done: " << emptyCount << " isolated vertices found" << std::endl;
    if (problemVerts.size()) {
        std::cerr << std::endl
                  << "Found " << problemVerts.size() << " duplicate faces in vertices: ";
        std::copy(problemVerts.begin(), problemVerts.end(),
                  std::ostream_iterator<size_t>(std::cerr, ", "));
        std::cerr << "\n";
    }
    std::cerr << std::endl
              << "The mesh has genus " << Genus() << ", and consists of " << Shells()
              << " shells.\n";

    std::cerr << "# Faces: " << std::to_string(mFaces.size()) << std::endl;
    std::cerr << "# Edges: " << std::to_string(mEdges.size() / 2) << std::endl;
    std::cerr << "# Vertices: " << std::to_string(mVerts.size()) << std::endl;
}

/*! \lab1 Implement the FindNeighborVertices */
/*! Loops over the neighborhood of a vertex and collects all the vertices sorted
 * counter clockwise. \param [in] vertexIndex  the index to vertex, size_t
 * \return a vector containing the indices to all the found vertices.
 */
std::vector<size_t> HalfEdgeMesh::FindNeighborVertices(size_t vertexIndex) const {
    // Collected vertices, sorted counter clockwise!
    std::vector<size_t> oneRing;

    // Add your code here
    size_t startIndex = v(vertexIndex).edge;
    size_t index = startIndex;
    size_t prevEdge, nextEdge;

    do {
        nextEdge = e(index).next;
        prevEdge = e(index).prev;
        oneRing.push_back(e(nextEdge).vert);
        index = e(prevEdge).pair;
    } while(index != startIndex);

    return oneRing;
}

/*! \lab1 Implement the FindNeighborFaces */
/*! Loops over the neighborhood of a vertex and collects all the faces sorted
 * counter clockwise. \param [in] vertexIndex  the index to vertex, size_t
 * \return a vector containing the indices to all the found faces.
 */
std::vector<size_t> HalfEdgeMesh::FindNeighborFaces(size_t vertexIndex) const {
    // Collected faces, sorted counter clockwise!
    std::vector<size_t> foundFaces;

    // Add your code here
    size_t startIndex = v(vertexIndex).edge;
    size_t index = startIndex;

    do {
        foundFaces.push_back(e(index).face);
        index = e(e(index).prev).pair;
    } while(index != startIndex);

    /*HalfEdge halfEdge = e(v(vertexIndex).edge);

    size_t startIndex = halfEdge.face;
    size_t index = startIndex;

    do {
        foundFaces.push_back(index);
        halfEdge = e(halfEdge.pair);
        halfEdge = e(halfEdge.prev);
        index = halfEdge.face;

    } while(index != startIndex);*/

    return foundFaces;
}

/*! \lab1 Implement the curvature */
float HalfEdgeMesh::VertexCurvature(size_t vertexIndex) const {
    // Copy code from SimpleMesh or compute more accurate estimate
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

float HalfEdgeMesh::FaceCurvature(size_t faceIndex) const {
    // NB Assumes vertex curvature already computed
    size_t indx = f(faceIndex).edge;
    const EdgeIterator it = GetEdgeIterator(indx);

    const auto& v1 = v(it.GetEdgeVertexIndex());
    const auto& v2 = v(it.Next().GetEdgeVertexIndex());
    const auto& v3 = v(it.Next().GetEdgeVertexIndex());

    return (v1.curvature + v2.curvature + v3.curvature) / 3.0f;
}

glm::vec3 HalfEdgeMesh::FaceNormal(size_t faceIndex) const {
    size_t indx = f(faceIndex).edge;
    const EdgeIterator it = GetEdgeIterator(indx);

    const auto& p1 = v(it.GetEdgeVertexIndex()).pos;
    const auto& p2 = v(it.Next().GetEdgeVertexIndex()).pos;
    const auto& p3 = v(it.Next().GetEdgeVertexIndex()).pos;

    const auto e1 = p2 - p1;
    const auto e2 = p3 - p1;
    return glm::normalize(glm::cross(e1, e2));
}

glm::vec3 HalfEdgeMesh::VertexNormal(size_t vertexIndex) const {

    glm::vec3 n(0.0f, 0.0f, 0.0f);

    // Add your code here

    auto faceVec = FindNeighborFaces(vertexIndex);

    // Calculate normal for each face
    // Sum all normals
    for (auto face : faceVec) {
        n += FaceNormal(face);
    }

    return glm::normalize(n);
}

void HalfEdgeMesh::Initialize() {
    Validate();
    Update();
}

void HalfEdgeMesh::Update() {
    // Calculate and store all differentials and area

    // First update all face normals and triangle areas
    for (size_t i = 0; i < GetNumFaces(); i++) {
        f(i).normal = FaceNormal(i);
    }
    // Then update all vertex normals and curvature
    for (size_t i = 0; i < GetNumVerts(); i++) {
        // Vertex normals are just weighted averages
        mVerts.at(i).normal = VertexNormal(i);
    }

    // Then update vertex curvature
    for (size_t i = 0; i < GetNumVerts(); i++) {
        mVerts.at(i).curvature = VertexCurvature(i);
        //    std::cerr <<   mVerts.at(i).curvature << "\n";
    }

    // Finally update face curvature
    for (size_t i = 0; i < GetNumFaces(); i++) {
        f(i).curvature = FaceCurvature(i);
    }

    std::cerr << "Area: " << Area() << ".\n";
    std::cerr << "Volume: " << Volume() << ".\n";

    // Update vertex and face colors
    if (mVisualizationMode == CurvatureVertex) {
        std::vector<Vertex>::iterator iter = mVerts.begin();
        std::vector<Vertex>::iterator iend = mVerts.end();
        float minCurvature = (std::numeric_limits<float>::max)();
        float maxCurvature = -(std::numeric_limits<float>::max)();
        while (iter != iend) {
            if (minCurvature > (*iter).curvature) minCurvature = (*iter).curvature;
            if (maxCurvature < (*iter).curvature) maxCurvature = (*iter).curvature;
            iter++;
        }
        std::cerr << "Mapping color based on vertex curvature with range [" << minCurvature << ","
                  << maxCurvature << "]" << std::endl;
        iter = mVerts.begin();
        while (iter != iend) {
            if (mColorMap != nullptr)
                (*iter).color = mColorMap->Map((*iter).curvature, minCurvature, maxCurvature);
            iter++;
        }
    } else if (mVisualizationMode == CurvatureFace) {
        std::vector<Face>::iterator iter = mFaces.begin();
        std::vector<Face>::iterator iend = mFaces.end();
        float minCurvature = (std::numeric_limits<float>::max)();
        float maxCurvature = -(std::numeric_limits<float>::max)();
        while (iter != iend) {
            if (minCurvature > (*iter).curvature) minCurvature = (*iter).curvature;
            if (maxCurvature < (*iter).curvature) maxCurvature = (*iter).curvature;
            iter++;
        }
        std::cerr << "Mapping color based on face curvature with range [" << minCurvature << ","
                  << maxCurvature << "]" << std::endl;
        iter = mFaces.begin();
        while (iter != iend) {
            if (mColorMap != nullptr)
                (*iter).color = mColorMap->Map((*iter).curvature, minCurvature, maxCurvature);
            iter++;
        }
    }
}

/*! \lab1 Implement the area */
float HalfEdgeMesh::Area() const {
    float area = 0.f;
    // Add code here
    for (auto face : mFaces) {
        auto E0 = e(face.edge);
        auto E1 = e(E0.next);
        auto E2 = e(E0.prev);

        area += 0.5f * glm::length(glm::cross((v(E1.vert).pos - v(E0.vert).pos), (v(E2.vert).pos - v(E1.vert).pos)));
    }
    return area;
}

/*! \lab1 Implement the volume */
float HalfEdgeMesh::Volume() const {
    float volume = 0;
    // Add code here
    for (auto face : mFaces) {
        auto v1 = v(e(face.edge).vert).pos;
        auto v2 = v(e(e(face.edge).next).vert).pos;
        auto v3 = v(e(e(face.edge).prev).vert).pos;

        auto vec1 = (v2 - v1);
        auto vec2 = (v3 - v1);

        float area = glm::length(glm::cross(vec1, vec2)) * 0.5f;
        volume += glm::dot(((v1 + v2 + v3) / glm::vec3(3.0f, 3.0f, 3.0f)), face.normal * area) / 3.0f;
    }
    std::cerr << "Volume calculation not implemented for half-edge mesh!\n";
    return volume;
}

/*! \lab1 Calculate the number of shells  */
size_t HalfEdgeMesh::Shells() const { return 1; }

/*! \lab1 Implement the genus */
size_t HalfEdgeMesh::Genus() const {
    // Add code here
    std::cerr << "Genus calculation not implemented for half-edge mesh!\n";
    return 0;
}

void HalfEdgeMesh::Dilate(float amount) {
    std::vector<Vertex>::iterator iter = mVerts.begin();
    std::vector<Vertex>::iterator iend = mVerts.end();
    while (iter != iend) {
        (*iter).pos += amount * (*iter).normal;
        iter++;
    }

    Initialize();
    Update();
}

void HalfEdgeMesh::Erode(float amount) {
    std::vector<Vertex>::iterator iter = mVerts.begin();
    std::vector<Vertex>::iterator iend = mVerts.end();
    while (iter != iend) {
        (*iter).pos -= amount * (*iter).normal;
        iter++;
    }

    Initialize();
    Update();
}

void HalfEdgeMesh::Smooth(float amount) {
    std::vector<Vertex>::iterator iter = mVerts.begin();
    std::vector<Vertex>::iterator iend = mVerts.end();
    while (iter != iend) {
        (*iter).pos -= amount * (*iter).normal * (*iter).curvature;
        iter++;
    }

    Initialize();
    Update();
}

void HalfEdgeMesh::Render() {
    glEnable(GL_LIGHTING);
    glMatrixMode(GL_MODELVIEW);

    // Apply transform
    glPushMatrix();  // Push modelview matrix onto stack

    // Convert transform-matrix to format matching GL matrix format
    // Load transform into modelview matrix
    glMultMatrixf(glm::value_ptr(mTransform));

    if (mWireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // Draw geometry
    glBegin(GL_TRIANGLES);
    const auto numTriangles = GetNumFaces();
    for (size_t i = 0; i < numTriangles; i++) {

        auto& face = f(i);

        auto* edge = &e(face.edge);

        auto& v1 = v(edge->vert);
        edge = &e(edge->next);

        auto& v2 = v(edge->vert);
        edge = &e(edge->next);

        auto& v3 = v(edge->vert);

        if (mVisualizationMode == CurvatureVertex) {
            glColor3fv(glm::value_ptr(v1.color));
            glNormal3fv(glm::value_ptr(v1.normal));
            glVertex3fv(glm::value_ptr(v1.pos));

            glColor3fv(glm::value_ptr(v2.color));
            glNormal3fv(glm::value_ptr(v2.normal));
            glVertex3fv(glm::value_ptr(v2.pos));

            glColor3fv(glm::value_ptr(v3.color));
            glNormal3fv(glm::value_ptr(v3.normal));
            glVertex3fv(glm::value_ptr(v3.pos));
        } else {
            glColor3fv(glm::value_ptr(face.color));
            glNormal3fv(glm::value_ptr(face.normal));

            glVertex3fv(glm::value_ptr(v1.pos));
            glVertex3fv(glm::value_ptr(v2.pos));
            glVertex3fv(glm::value_ptr(v3.pos));
        }
    }
    glEnd();

    // Mesh normals by courtesy of Richard Khoury
    if (mShowNormals) {
        glDisable(GL_LIGHTING);
        glBegin(GL_LINES);
        const auto numTriangles = GetNumFaces();
        for (size_t i = 0; i < numTriangles; i++) {

            auto& face = f(i);

            auto* edge = &e(face.edge);

            auto& v1 = v(edge->vert);
            edge = &e(edge->next);

            auto& v2 = v(edge->vert);
            edge = &e(edge->next);

            auto& v3 = v(edge->vert);

            auto faceStart = (v1.pos + v2.pos + v3.pos) / 3.0f;
            auto faceEnd = faceStart + face.normal * 0.1f;

            glColor3f(1.0f, 0.0f, 0.0f);  // Red for face normal
            glVertex3fv(glm::value_ptr(faceStart));
            glVertex3fv(glm::value_ptr(faceEnd));

            glColor3f(0.0f, 1.0f, 0.0f);  // Vertex normals in Green
            glVertex3fv(glm::value_ptr(v1.pos));
            glVertex3fv(glm::value_ptr((v1.pos + v1.normal * 0.1f)));
            glVertex3fv(glm::value_ptr(v2.pos));
            glVertex3fv(glm::value_ptr((v2.pos + v2.normal * 0.1f)));
            glVertex3fv(glm::value_ptr(v3.pos));
            glVertex3fv(glm::value_ptr((v3.pos + v3.normal * 0.1f)));
        }
        glEnd();
        glEnable(GL_LIGHTING);
    }

    if (mWireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // Restore modelview matrix
    glPopMatrix();

    GLObject::Render();
}
