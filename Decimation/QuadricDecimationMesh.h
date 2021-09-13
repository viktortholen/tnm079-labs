#pragma once

#include "Decimation/DecimationMesh.h"
#include <iomanip>

#ifdef __APPLE__
#include "GLUT/glut.h"
#else
#include "GL/glut.h"
#endif

class QuadricDecimationMesh : public virtual DecimationMesh {
public:
    static const VisualizationMode QuadricIsoSurfaces;

    virtual std::list<VisualizationMode> GetVisualizationModes() {
        std::list<VisualizationMode> L = DecimationMesh::GetVisualizationModes();
        L.push_back(QuadricIsoSurfaces);
        return L;
    }

    QuadricDecimationMesh() {}
    virtual ~QuadricDecimationMesh() {}

    //! Initialize member data (error quadrics)
    virtual void Initialize();

protected:
    //! Compute the cost and new position for an edge collapse
    virtual void computeCollapse(EdgeCollapse* collapse);
    //! Update vertex properties. Used after an edge collapse
    virtual void updateVertexProperties(size_t ind);
    //! Compute the quadric for a vertex
    glm::mat4 createQuadricForVert(size_t indx) const;
    //! Copmute the quadric for a face
    glm::mat4 createQuadricForFace(size_t indx) const;
    //! Render (redefined)
    virtual void Render();

    //! The quadrics used in the decimation
    std::vector<glm::mat4> mQuadrics;
};
