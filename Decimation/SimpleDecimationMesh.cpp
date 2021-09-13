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
#include "SimpleDecimationMesh.h"

void SimpleDecimationMesh::computeCollapse(EdgeCollapse *collapse) {
    // The new vertex position is implicitly stored as the
    // position halfway along the edge. The cost is computed as
    // the vertex-to-vertex distance between the new vertex
    // and the old vertices at the edge's endpoints
    const glm::vec3 &v0 = mVerts[mEdges[collapse->halfEdge].vert].pos;
    const glm::vec3 &v1 = mVerts[mEdges[mEdges[collapse->halfEdge].pair].vert].pos;

    collapse->position = 0.5f * (v0 + v1);
    collapse->cost = glm::distance(collapse->position, v0);
}
