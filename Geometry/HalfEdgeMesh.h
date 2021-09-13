#pragma once

#include <Geometry/Mesh.h>
#include <Util/ObjIO.h>
#include <Util/Util.h>
#include <cassert>
#include <limits>
#include <map>
#include <set>

/*! \brief A half edge triangle mesh class.
 * A mesh data structure with fast access to neighborhood entities.
 * The neighborhood links are stored through indices.
 */
class HalfEdgeMesh : public Mesh {
public:
    HalfEdgeMesh();
    ~HalfEdgeMesh();
    virtual void Update();
    virtual void Initialize();

    class EdgeIterator {
        friend class HalfEdgeMesh;  // why friend?
    protected:
        // should not change what mesh the pointer points to, can change mesh though
        HalfEdgeMesh const *mHem;
        mutable size_t mIndex;
        EdgeIterator(HalfEdgeMesh const *he, size_t index) {
            mHem = he;
            // Guard against outside access
            assert(index < mHem->GetNumEdges());
            mIndex = index;
        }

    public:
        EdgeIterator &Next() {
            mIndex = mHem->e(mIndex).next;
            return *this;
        }
        EdgeIterator &Prev() {
            mIndex = mHem->e(mIndex).prev;
            return *this;
        }
        EdgeIterator &Pair() {
            mIndex = mHem->e(mIndex).pair;
            return *this;
        }
        const EdgeIterator &Next() const {
            mIndex = mHem->e(mIndex).next;
            return *this;
        }
        const EdgeIterator &Prev() const {
            mIndex = mHem->e(mIndex).prev;
            return *this;
        }
        const EdgeIterator &Pair() const {
            mIndex = mHem->e(mIndex).pair;
            return *this;
        }
        size_t GetEdgeIndex() const { return mIndex; }
        size_t GetEdgeVertexIndex() const { return mHem->e(mIndex).vert; }
        size_t GetEdgeFaceIndex() const { return mHem->e(mIndex).face; }
        const bool operator==(const EdgeIterator &eit) { return this->mIndex == eit.mIndex; }
        const bool operator!=(const EdgeIterator &eit) { return this->mIndex != eit.mIndex; }
    };

    EdgeIterator GetEdgeIterator(size_t i) { return EdgeIterator(this, i); }
    const EdgeIterator GetEdgeIterator(size_t i) const { return EdgeIterator(this, i); }

    //! Adds a triangle to the mesh.
    virtual bool AddFace(const std::vector<glm::vec3> &verts);

    //! Calculates the area of the mesh
    virtual float Area() const;

    //! Calculates the volume of the mesh
    virtual float Volume() const;

    //! Calculates the genus of a mesh
    virtual size_t Genus() const;

    //! Calculates the number of shells
    virtual size_t Shells() const;

    //! Calculates the curvature at a vertex
    virtual float VertexCurvature(size_t vertexIndex) const;

    //! Calculates the curvature at a vertex
    virtual float FaceCurvature(size_t faceIndex) const;

    //! Calculates the normal at a face
    virtual glm::vec3 FaceNormal(size_t faceIndex) const;

    //! Calculates the normal at a vertex
    virtual glm::vec3 VertexNormal(size_t vertexIndex) const;

    //! Checks to see if the mesh is valid
    void Validate();

    virtual void Render() override;

protected:

    enum EdgeState : size_t {
        Uninitialized = std::numeric_limits<size_t>::max() - 1,
        Border = std::numeric_limits<size_t>::max()
    };

    /*! \brief The core half-edge struct
     *  Implements the linked data structure edge type
     */
    struct HalfEdge {
        HalfEdge()
            : vert(EdgeState::Uninitialized)
            , face(EdgeState::Uninitialized)
            , next(EdgeState::Uninitialized)
            , prev(EdgeState::Uninitialized)
            , pair(EdgeState::Uninitialized) {}
        size_t vert;  //!< index into mVerts (the origin vertex)
        size_t face;  //!< index into mFaces
        size_t next;  //!< index into mEdges
        size_t prev;  //!< index into mEdges
        size_t pair;  //!< index into mEdges
    };

    /*! \brief A vertex is a point and an edge index
     */
    struct Vertex : public Mesh::Vertex {
        Vertex() : edge(EdgeState::Uninitialized) {}
        size_t edge;  //!< index into mEdges
    };

    /*! \brief A face has a normal and an index to an edge
     */
    struct Face : public Mesh::Face {
        Face() : edge(EdgeState::Uninitialized) {}
        size_t edge;  //!< index into mEdges
    };

    //! The edges of the mesh
    std::vector<HalfEdge> mEdges;
    //! The vertices in the mesh
    std::vector<Vertex> mVerts;
    //! The faces in the mesh
    std::vector<Face> mFaces;

    //! A utility data structure to speed up removal of redundant vertices
    std::map<glm::vec3, size_t, glm_vec3_comparator<float>> mUniqueVerts;

    struct OrderedPair {
        size_t p1, p2;
        OrderedPair(size_t i1, size_t i2) {
            p1 = std::min(i1, i2);
            p2 = std::max(i1, i2);
        }
        bool operator<(const OrderedPair &rhs) const {
            if (this->p1 < rhs.p1) {
                return true;
            } else if (this->p1 == rhs.p1) {
                if (this->p2 < rhs.p2) {
                    return true;
                }
            }
            return false;
        }
    };
    //! A utility data structure to speed up removal of redundant edges
    std::map<OrderedPair, size_t> mUniqueEdgePairs;

    //! Adds a vertex to the mesh
    virtual size_t AddVertex(const glm::vec3 &v) override;

    //! Adds a half edge pair, from vertex 1 to vertex 2, to the mesh
    std::pair<size_t, size_t> AddHalfEdgePair(size_t v1, size_t v2);

    //! Merges the outer (uninitialized) edge for a newly created triangle
    void MergeOuterBoundaryEdge(size_t innerEdge);

    //! Finds all triangles that includes a given vertex.
    virtual std::vector<size_t> FindNeighborFaces(size_t vertexIndex) const;

    //! Finds all vertices that includes a given vertex.
    virtual std::vector<size_t> FindNeighborVertices(size_t vertexIndex) const;

    //! Return the edge at index i
    HalfEdge &e(size_t i) { return mEdges.at(i); }
    const HalfEdge &e(size_t i) const { return mEdges.at(i); }
    //! Return the face at index i
    Face &f(size_t i) { return mFaces.at(i); }
    const Face &f(size_t i) const { return mFaces.at(i); }
    //! Return the Vertex at index i
    Vertex &v(size_t i) { return mVerts.at(i); }
    const Vertex v(size_t i) const { return mVerts.at(i); }
    //! Return number of vertices
    size_t GetNumVerts() const { return mVerts.size(); }
    //! Return number of faces
    size_t GetNumFaces() const { return mFaces.size(); }
    //! Return number of edges
    size_t GetNumEdges() const { return mEdges.size(); }

    virtual void Dilate(float amount);
    virtual void Erode(float amount);
    virtual void Smooth(float amount);

    virtual bool save(std::ostream &os) {
        os << "# HalfEdgeMesh obj streamer\n# M&A 2008\n\n";
        os << "# Vertices\n";
        for (size_t i = 0; i < GetNumVerts(); i++) {
            os << "v " << v(i).pos[0] << " " << v(i).pos[1] << " " << v(i).pos[2] << "\n";
        }
        os << "\n# Faces\n";
        for (size_t i = 0; i < GetNumFaces(); i++) {
            size_t ind = f(i).edge;
            os << "f ";
            os << e(ind).vert + 1 << " ";
            ind = e(ind).next;
            os << e(ind).vert + 1 << " ";
            ind = e(ind).next;
            os << e(ind).vert + 1 << "\n";
        }
        return os.good();
    }
};
