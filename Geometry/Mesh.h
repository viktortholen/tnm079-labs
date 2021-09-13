#pragma once

#include <vector>
#include <Geometry/Geometry.h>

class Mesh : public Geometry {
public:
    static const VisualizationMode CurvatureVertex;
    static const VisualizationMode CurvatureFace;

    virtual std::list<VisualizationMode> GetVisualizationModes() {
        std::list<VisualizationMode> L = Geometry::GetVisualizationModes();
        L.push_back(CurvatureVertex);
        L.push_back(CurvatureFace);
        return L;
    }

    virtual bool save(std::ostream &os) = 0;

protected:
    //! Adds a vertex to the mesh
    virtual size_t AddVertex(const glm::vec3 &v) = 0;

    //! Given a vertex, find all triangles that includes this vertex (sorted
    //! counter clockwise)
    virtual std::vector<size_t> FindNeighborFaces(size_t vertexIndex) const = 0;

    //! Given a vertex, find the one-ring neighborhood. Ie all vertices that are
    //! connected to this vertex (sorted counter clockwise)
    virtual std::vector<size_t> FindNeighborVertices(size_t vertexIndex) const = 0;

    //! Compute and return the curvature at vertex at vertexIndex
    virtual float VertexCurvature(size_t vertexIndex) const = 0;

    //! Compute and return the curvature at face
    virtual float FaceCurvature(size_t faceIndex) const = 0;

    //! Compute and return the normal at face at faceIndex
    virtual glm::vec3 FaceNormal(size_t faceIndex) const = 0;

    //! Compute and return the normal at vertex at vertexIndex
    virtual glm::vec3 VertexNormal(size_t vertexIndex) const = 0;

    bool mVisualizeNormals;

public:
    //! Minimal requirements for all meshes, inherited
    struct Face {
        Face(const glm::vec3 &n = glm::vec3(0.0f, 0.0f, 0.0f),
             const glm::vec3 &c = glm::vec3(0.5f, 0.1f, 0.7f), float u = 0)
            : normal(n), color(c), curvature(u) {}
        glm::vec3 normal;
        glm::vec3 color;
        float curvature;
    };
    //! Minimal requirements for all meshes, inherited
    struct Vertex {
        Vertex(const glm::vec3 &p = glm::vec3(0.0f, 0.0f, 0.0f),
               const glm::vec3 &n = glm::vec3(0.0f, 0.0f, 0.0f),
               const glm::vec3 &c = glm::vec3(0.5f, 0.1f, 0.7f), float u = 0)
            : pos(p), normal(n), color(c), curvature(u) {}
        glm::vec3 pos;
        glm::vec3 normal;
        glm::vec3 color;
        float curvature;
    };

    Mesh() : mVisualizeNormals(false) { mVisualizationMode = CurvatureFace; }
    virtual ~Mesh() {}

    //! Adds a face to the mesh.
    virtual bool AddFace(const std::vector<glm::vec3> &verts) = 0;

    //! Compute area of mesh
    virtual float Area() const;
    //! Compute volume of mesh
    virtual float Volume() const;

    //! Compute genus of mesh
    virtual size_t Genus() const;

    virtual void Dilate(float epsilon) {}

    virtual void SetColorMap(ColorMap *colormap) {
        GLObject::SetColorMap(colormap);
        Update();
    }

    virtual void VisualizeNormals(bool flag = true) { mVisualizeNormals = flag; }

    virtual void SetVisualizationMode(const VisualizationMode &mode) {
        GLObject::SetVisualizationMode(mode);
        Update();
    }

    virtual const char *GetTypeName() { return typeid(Mesh).name(); }
};
