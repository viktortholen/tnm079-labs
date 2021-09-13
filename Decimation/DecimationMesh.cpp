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
#include <cassert>
#include <gtc/type_ptr.hpp>
#include <Decimation/DecimationMesh.h>
#include <GUI/GUI.h>

const DecimationMesh::VisualizationMode DecimationMesh::CollapseCost =
    NewVisualizationMode("Collapse cost");

void DecimationMesh::Initialize() {
    // Allocate memory for the 'collapsed flags'
    mCollapsedVerts.reserve(mVerts.size());
    mCollapsedEdges.reserve(mEdges.size());
    mCollapsedFaces.reserve(mFaces.size());

    // Set all flags to false
    mCollapsedVerts.assign(mVerts.size(), false);
    mCollapsedEdges.assign(mEdges.size(), false);
    mCollapsedFaces.assign(mFaces.size(), false);

    // Allocate memory for the references from half-edge
    // to edge collapses
    mHalfEdge2EdgeCollapse.reserve(mEdges.size());
    mHalfEdge2EdgeCollapse.assign(mEdges.size(), NULL);

    // Loop through the half-edges (we know they are stored
    // sequentially) and create an edge collapse operation
    // for each pair
    auto numCollapses = mEdges.size() / 2;
    for (size_t i = 0; i < numCollapses; i++) {
        EdgeCollapse *collapse = new EdgeCollapse();

        // Connect the edge collapse with the half-edge pair
        collapse->halfEdge = i * 2;

        // Check if the collapse is valid
        if (!isValidCollapse(collapse))
            delete collapse;

        else {
            mHalfEdge2EdgeCollapse[i * 2] = collapse;
            mHalfEdge2EdgeCollapse[i * 2 + 1] = collapse;

            // Compute the cost and push it to the heap
            computeCollapse(collapse);
            mHeap.push(collapse);
        }
    }
    // mHeap.print(std::cout);

    HalfEdgeMesh::Initialize();
}

void DecimationMesh::Update() {
    // Calculate and store all differentials and area

    // First update all face normals and triangle areas
    for (size_t i = 0; i < GetNumFaces(); i++) {
        if (!isFaceCollapsed(i)) f(i).normal = FaceNormal(i);
    }
    // Then update all vertex normals and curvature
    for (size_t i = 0; i < GetNumVerts(); i++) {
        // Vertex normals are just weighted averages
        if (!isVertexCollapsed(i)) mVerts.at(i).normal = VertexNormal(i);
    }

    // Then update vertex curvature
    for (size_t i = 0; i < GetNumVerts(); i++) {
        if (!isVertexCollapsed(i)) mVerts.at(i).curvature = VertexCurvature(i);
        //    std::cerr <<   mVerts.at(i).curvature << "\n";
    }

    // Finally update face curvature
    for (size_t i = 0; i < GetNumFaces(); i++) {
        if (!isFaceCollapsed(i)) f(i).curvature = FaceCurvature(i);
    }

    //  std::cerr << "Area: " << Area() << ".\n";
    //  std::cerr << "Volume: " << Volume() << ".\n";

    // Update vertex and face colors
    if (mVisualizationMode == CurvatureVertex) {
        std::vector<Vertex>::iterator iter = mVerts.begin();
        std::vector<Vertex>::iterator iend = mVerts.end();
        float minCurvature = (std::numeric_limits<float>::max)();
        float maxCurvature = -(std::numeric_limits<float>::max)();
        while (iter != iend) {
            if (!isVertexCollapsed(iter - mVerts.begin())) {
                if (minCurvature > (*iter).curvature) minCurvature = (*iter).curvature;
                if (maxCurvature < (*iter).curvature) maxCurvature = (*iter).curvature;
            }
            iter++;
        }
        std::cerr << "Mapping color based on vertex curvature with range [" << minCurvature << ","
                  << maxCurvature << "]" << std::endl;
        iter = mVerts.begin();
        while (iter != iend) {
            if (!isVertexCollapsed(iter - mVerts.begin()))
                (*iter).color = mColorMap->Map((*iter).curvature, minCurvature, maxCurvature);
            iter++;
        }
    } else if (mVisualizationMode == CurvatureFace) {
        std::vector<Face>::iterator iter = mFaces.begin();
        std::vector<Face>::iterator iend = mFaces.end();
        float minCurvature = (std::numeric_limits<float>::max)();
        float maxCurvature = -(std::numeric_limits<float>::max)();
        while (iter != iend) {
            if (!isFaceCollapsed(iter - mFaces.begin())) {
                if (minCurvature > (*iter).curvature) minCurvature = (*iter).curvature;
                if (maxCurvature < (*iter).curvature) maxCurvature = (*iter).curvature;
            }
            iter++;
        }
        std::cerr << "Mapping color based on face curvature with range [" << minCurvature << ","
                  << maxCurvature << "]" << std::endl;
        iter = mFaces.begin();
        while (iter != iend) {
            if (!isFaceCollapsed(iter - mFaces.begin()))
                (*iter).color = mColorMap->Map((*iter).curvature, minCurvature, maxCurvature);
            iter++;
        }
    }
}

bool DecimationMesh::decimate(size_t targetFaces) {
    // We can't collapse down to less than two faces
    if (targetFaces < 2) targetFaces = 2;

    // Keep collapsing one edge at a time until the target is reached
    // or the heap is empty (when we have no possible collapses left)
    while (mFaces.size() - mNumCollapsedFaces > targetFaces && !mHeap.isEmpty()) decimate();

    // Return true if target is reached
    std::cout << "Collapsed mesh to " << mFaces.size() - mNumCollapsedFaces << " faces"
              << std::endl;
    return mFaces.size() - mNumCollapsedFaces == targetFaces;
}

bool DecimationMesh::decimate() {
    EdgeCollapse *collapse = static_cast<EdgeCollapse *>(mHeap.pop());
    if (collapse == NULL) return false;

    // Stop the collapse when we only have two triangles left
    // (the smallest entity representable)
    if (mFaces.size() - mNumCollapsedFaces == 2) return false;

    size_t e1 = collapse->halfEdge;
    size_t e2 = mEdges[e1].pair;

    size_t v1 = mEdges[e1].vert;
    size_t v2 = mEdges[e2].vert;
    size_t v3 = mEdges[mEdges[e1].prev].vert;
    size_t v4 = mEdges[mEdges[e2].prev].vert;

    size_t f1 = mEdges[e1].face;
    size_t f2 = mEdges[e2].face;

#ifndef NDEBUG
    std::cout << "Collapsing faces " << f1 << " and " << f2 << std::endl;
    std::cout << "Collapsing edges " << e1 << ", " << mEdges[e1].next << ", " << mEdges[e1].prev;
    std::cout << ", " << e2 << ", " << mEdges[e2].next << " and " << mEdges[e2].prev << std::endl;
    std::cout << "Collapsing vertex " << v1 << std::endl;
#endif

    // Verify that the collapse is valid, exit if not so
    if (!isValidCollapse(collapse)) {
        delete collapse;
        mHalfEdge2EdgeCollapse[e1] = NULL;
        mHalfEdge2EdgeCollapse[e2] = NULL;
        std::cout << "failed..." << std::endl;
        return false;
    }

    // We want to remove v1, so we need to connect all of v1's half-edges to v2
    size_t edge = mVerts[v1].edge;
    do {
        mEdges[edge].vert = v2;
        edge = mEdges[mEdges[edge].pair].next;
    } while (edge != mVerts[v1].edge);

    // Make sure v2 points to a valid edge
    while (mEdges[mVerts[v2].edge].face == f1 || mEdges[mVerts[v2].edge].face == f2)
        mVerts[v2].edge = mEdges[mEdges[mVerts[v2].edge].pair].next;

    // Make sure v3 points to a valid edge
    while (mEdges[mVerts[v3].edge].face == f1)
        mVerts[v3].edge = mEdges[mEdges[mVerts[v3].edge].pair].next;

    // Make sure v4 points to a valid edge
    while (mEdges[mVerts[v4].edge].face == f2)
        mVerts[v4].edge = mEdges[mEdges[mVerts[v4].edge].pair].next;

    // Redirect pair pointers
    mEdges[mEdges[mEdges[e1].next].pair].pair = mEdges[mEdges[e1].prev].pair;
    mEdges[mEdges[mEdges[e1].prev].pair].pair = mEdges[mEdges[e1].next].pair;

    mEdges[mEdges[mEdges[e2].next].pair].pair = mEdges[mEdges[e2].prev].pair;
    mEdges[mEdges[mEdges[e2].prev].pair].pair = mEdges[mEdges[e2].next].pair;

    // Move v2 to its new position
    mVerts[v2].pos = collapse->position;

    // One edge collapse further removes 2 additional collapse
    // candidates from the heap
    if (mHalfEdge2EdgeCollapse[mEdges[e1].prev] != NULL)
        delete mHeap.remove(mHalfEdge2EdgeCollapse[mEdges[e1].prev]);
    mHalfEdge2EdgeCollapse[mEdges[mEdges[e1].prev].pair] = mHalfEdge2EdgeCollapse[mEdges[e1].next];

    if (mHalfEdge2EdgeCollapse[mEdges[e2].next] != NULL)
        delete mHeap.remove(mHalfEdge2EdgeCollapse[mEdges[e2].next]);
    mHalfEdge2EdgeCollapse[mEdges[mEdges[e2].next].pair] = mHalfEdge2EdgeCollapse[mEdges[e2].prev];

    // Make sure the edge collapses point to valid edges
    if (mHalfEdge2EdgeCollapse[mEdges[e1].next] != NULL)
        mHalfEdge2EdgeCollapse[mEdges[e1].next]->halfEdge = mEdges[mEdges[e1].prev].pair;
    if (mHalfEdge2EdgeCollapse[mEdges[e2].prev] != NULL)
        mHalfEdge2EdgeCollapse[mEdges[e2].prev]->halfEdge = mEdges[mEdges[e2].next].pair;

    delete collapse;

    // Collapse the neighborhood
    collapseFace(f1);
    collapseFace(f2);

    collapseEdge(e1);
    collapseEdge(mEdges[e1].next);
    collapseEdge(mEdges[e1].prev);

    collapseEdge(e2);
    collapseEdge(mEdges[e2].next);
    collapseEdge(mEdges[e2].prev);

    collapseVertex(v1);

    // Finally, loop through neighborhood of v2 and update all edge collapses
    // (and remove possible invalid cases)
    updateVertexProperties(v2);
    edge = mVerts[v2].edge;
    do {
        size_t face = mEdges[edge].face;
        size_t vert = mEdges[mEdges[edge].pair].vert;
        if (!isFaceCollapsed(face)) updateFaceProperties(face);
        if (!isVertexCollapsed(vert)) updateVertexProperties(vert);

        collapse = mHalfEdge2EdgeCollapse[edge];
        if (collapse != NULL) {
            if (!isValidCollapse(collapse)) {
                delete mHeap.remove(collapse);
                mHalfEdge2EdgeCollapse[edge] = NULL;
                mHalfEdge2EdgeCollapse[mEdges[edge].pair] = NULL;
#ifndef NDEBUG
                std::cout << "Removed one invalid edge collapse" << std::endl;
#endif
            } else {
                computeCollapse(collapse);
                mHeap.update(collapse);
            }
        }

        edge = mEdges[mEdges[edge].pair].next;
    } while (edge != mVerts[v2].edge);

    // mHeap.print(std::cout);

    return true;
}

void DecimationMesh::updateVertexProperties(size_t ind) {
    std::vector<size_t> neighborFaces = FindNeighborFaces(ind);

    // Approximate vertex normal
    glm::vec3 n(0, 0, 0);

    const auto numCandidates = neighborFaces.size();
    for (size_t i = 0; i < numCandidates; i++) {
        Face &triangle = mFaces[neighborFaces[i]];

        // Calculate face normal
        HalfEdge *edge = &mEdges[triangle.edge];

        glm::vec3 &p0 = mVerts[edge->vert].pos;
        edge = &mEdges[edge->next];

        glm::vec3 &p1 = mVerts[edge->vert].pos;
        edge = &mEdges[edge->next];

        glm::vec3 &p2 = mVerts[edge->vert].pos;

        glm::vec3 v1 = p1 - p0;
        glm::vec3 v2 = p2 - p0;
        glm::vec3 faceNormal = glm::cross(v1, v2);
        faceNormal = glm::normalize(faceNormal);

        n += faceNormal;
    }

    n = glm::normalize(n);
    mVerts[ind].normal = n;
}

void DecimationMesh::updateFaceProperties(size_t ind) {
    HalfEdge *edge = &mEdges[mFaces[ind].edge];

    glm::vec3 &p0 = mVerts[edge->vert].pos;
    edge = &mEdges[edge->next];

    glm::vec3 &p1 = mVerts[edge->vert].pos;
    edge = &mEdges[edge->next];

    glm::vec3 &p2 = mVerts[edge->vert].pos;

    // Calculate face normal
    glm::vec3 v1 = p1 - p0;
    glm::vec3 v2 = p2 - p0;
    glm::vec3 n = glm::cross(v1, v2);
    n = glm::normalize(n);

    mFaces[ind].normal = n;
}

bool DecimationMesh::isValidCollapse(EdgeCollapse *collapse) {
    size_t e1 = collapse->halfEdge;
    size_t e2 = mEdges[e1].pair;

    size_t v1 = mEdges[e1].vert;
    size_t v2 = mEdges[e2].vert;
    size_t v3 = mEdges[mEdges[e1].prev].vert;
    size_t v4 = mEdges[mEdges[e2].prev].vert;

    // Do a dummy check
    if (isEdgeCollapsed(e1) || isEdgeCollapsed(e1) || isVertexCollapsed(v1) ||
        isVertexCollapsed(v2))
        return false;

    size_t edge = mVerts[v2].edge;
    std::vector<size_t> neighbors;
    do {
        size_t ind = mEdges[mEdges[edge].pair].vert;
        if (ind != v3 && ind != v4) neighbors.push_back(ind);
        edge = mEdges[mEdges[edge].pair].next;
    } while (edge != mVerts[v2].edge);

    edge = mVerts[v1].edge;
    do {
        size_t ind = mEdges[mEdges[edge].pair].vert;
        if (find(neighbors.begin(), neighbors.end(), ind) != neighbors.end()) return false;

        edge = mEdges[mEdges[edge].pair].next;
    } while (edge != mVerts[v1].edge);

    return true;
}

void DecimationMesh::Cleanup() {
    //  HalfEdgeMesh mesh;
    //  *this = mesh;
}

void DecimationMesh::Render() {
    glEnable(GL_LIGHTING);
    glMatrixMode(GL_MODELVIEW);

    // Apply transform
    glPushMatrix();  // Push modelview matrix onto stack

    // Convert transform-matrix to format matching GL matrix format
    // Load transform into modelview matrix
    glMultMatrixf(glm::value_ptr(mTransform));

    if (mWireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // Draw geometry
    const auto numTriangles = mFaces.size();
    for (size_t i = 0; i < numTriangles; i++) {

        if (isFaceCollapsed(i)) continue;

        // Render without notations
        Face &f = mFaces[i];

        HalfEdge *edge = &mEdges[f.edge];

        Vertex &v1 = mVerts[edge->vert];
        edge = &mEdges[edge->next];

        Vertex &v2 = mVerts[edge->vert];
        edge = &mEdges[edge->next];

        Vertex &v3 = mVerts[edge->vert];

        // Render with notations
        //  Uncomment this block, and comment the block above
        //    to render with notations
        /*
            Face & f = mFaces[i];

            char buffer[10];
            glColor3f(1.0, 0.0, 0.0);
            HalfEdge* edge = &mEdges[mFaces[i].edge];

            // draw face
            sprintf(buffer, "f%i\n", i);
            glm::vec3 vec = (mVerts[edge->vert].pos +
           mVerts[mEdges[edge->pair].vert].pos)*0.5; vec += 0.5 *
           (mVerts[mEdges[edge->prev].vert].pos - vec); drawText(vec, buffer);

            // draw e1
            sprintf(buffer, "e%i\n", mFaces[i].edge);
            vec = (mVerts[edge->vert].pos +
           mVerts[mEdges[edge->pair].vert].pos)*0.5; vec += 0.1 *
           (mVerts[mEdges[edge->prev].vert].pos - vec); drawText(vec, buffer);

            // draw v1
            Vertex& v1 = mVerts[edge->vert];
            sprintf(buffer, "v%i\n", edge->vert);
            drawText(vec, buffer);

            sprintf(buffer, "e%i\n", edge->next);
            edge = &mEdges[edge->next];

            // draw e2
            vec = (mVerts[edge->vert].pos +
           mVerts[mEdges[edge->pair].vert].pos)*0.5; vec += 0.1 *
           (mVerts[mEdges[edge->prev].vert].pos - vec); drawText(vec, buffer);

            // draw v2
            Vertex& v2 = mVerts[edge->vert];
            sprintf(buffer, "v%i\n", edge->vert);
            drawText(vec, buffer);

            sprintf(buffer, "e%i\n", edge->next);
            edge = &mEdges[edge->next];

            // draw e3
            vec = (mVerts[edge->vert].pos +
           mVerts[mEdges[edge->pair].vert].pos)*0.5; vec += 0.1 *
           (mVerts[mEdges[edge->prev].vert].pos - vec); drawText(vec, buffer);

            // draw v3
            Vertex& v3 = mVerts[edge->vert];
            sprintf(buffer, "v%i\n", edge->vert);
            drawText(vec, buffer);
         */

        glBegin(GL_TRIANGLES);
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
            glColor3fv(glm::value_ptr(f.color));
            glNormal3fv(glm::value_ptr(f.normal));

            glVertex3fv(glm::value_ptr(v1.pos));
            glVertex3fv(glm::value_ptr(v2.pos));
            glVertex3fv(glm::value_ptr(v3.pos));
        }
        glEnd();
    }

    if (mWireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    if (mVisualizationMode == CollapseCost) {

        float minCost = std::numeric_limits<float>::max();
        float maxCost = -std::numeric_limits<float>::max();
        std::vector<HalfEdge>::const_iterator iter = mEdges.begin();
        std::vector<HalfEdge>::const_iterator iend = mEdges.end();
        while (iter != iend) {
            EdgeCollapse *collapse = mHalfEdge2EdgeCollapse.at(iter - mEdges.begin());
            if (collapse != NULL) {
                if (minCost > collapse->cost) minCost = collapse->cost;
                if (maxCost < collapse->cost) maxCost = collapse->cost;
            }

            iter++;
        }
        std::cout << "Colormapping collapse cost with range: [" << minCost << "," << maxCost << "]"
                  << std::endl;

        GLfloat lineWidth;
        glGetFloatv(GL_LINE_WIDTH, &lineWidth);
        glLineWidth(5.0f);
        glDisable(GL_LIGHTING);
        glEnable(GL_LINE_SMOOTH);
        glBegin(GL_LINES);
        iter = mEdges.begin();
        while (iter != iend) {
            size_t ind = iter - mEdges.begin();
            if (!isEdgeCollapsed(ind)) {

                const Vertex &v1 = v(e(ind).vert);
                const Vertex &v2 = v(e(e(ind).pair).vert);

                EdgeCollapse *collapse = mHalfEdge2EdgeCollapse.at(ind);
                if (collapse == NULL)
                    glColor3f(1, 1, 1);
                else
                    glColor3fv(glm::value_ptr(mColorMap->Map(collapse->cost, minCost, maxCost)));

                glVertex3fv(glm::value_ptr(v1.pos));
                glVertex3fv(glm::value_ptr(v2.pos));
            }

            iter++;
        }
        glEnd();
        glLineWidth(lineWidth);
        glEnable(GL_LIGHTING);
        glDisable(GL_LINE_SMOOTH);
    }

    // Restore modelview matrix
    glPopMatrix();

    GLObject::Render();
}

void DecimationMesh::drawText(const glm::vec3 &pos, const char *str) {
    glRasterPos3f(pos[0], pos[1], pos[2]);
    for (size_t i = 0; str[i] != '\n'; i++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
}
