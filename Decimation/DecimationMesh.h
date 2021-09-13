#pragma once

#include <Decimation/DecimationInterface.h>
#include <Geometry/HalfEdgeMesh.h>
#include <Util/ColorMap.h>
#include <Util/Heap.h>

class DecimationMesh : public DecimationInterface, public HalfEdgeMesh {
public:
  static const VisualizationMode CollapseCost;
  virtual std::list<VisualizationMode> GetVisualizationModes() {
    std::list<VisualizationMode> L = Mesh::GetVisualizationModes();
    L.push_back(CollapseCost);
    return L;
  }

  DecimationMesh()
      : mNumCollapsedVerts(0), mNumCollapsedEdges(0), mNumCollapsedFaces(0) {}
  virtual ~DecimationMesh() {}

  /*! The EdgeCollapse is a Heapable type */
  struct EdgeCollapse : public Heap::Heapable {
    size_t halfEdge;
    glm::vec3 position;
  };

  virtual void Initialize();

  virtual void Update();

  virtual bool decimate();

  virtual bool decimate(size_t targetFaces);

  virtual void Render() override;

  virtual const char *GetTypeName() { return typeid(DecimationMesh).name(); }

protected:
  virtual void updateVertexProperties(size_t ind);

  virtual void updateFaceProperties(size_t ind);

  virtual void computeCollapse(EdgeCollapse *collapse) = 0;

  virtual void Cleanup();

  bool isValidCollapse(EdgeCollapse *collapse);

  inline bool isVertexCollapsed(size_t ind) { return mCollapsedVerts[ind]; }
  inline bool isEdgeCollapsed(size_t ind) { return mCollapsedEdges[ind]; }
  inline bool isFaceCollapsed(size_t ind) { return mCollapsedFaces[ind]; }

  inline void collapseVertex(size_t ind) {
    mCollapsedVerts[ind] = true;
    mNumCollapsedVerts++;
  }
  inline void collapseEdge(size_t ind) {
    mHalfEdge2EdgeCollapse[ind] = NULL;
    mCollapsedEdges[ind] = true;
    mNumCollapsedEdges++;
  }
  inline void collapseFace(size_t ind) {
    mCollapsedFaces[ind] = true;
    mNumCollapsedFaces++;
  }

  //! State array of 'active' verts
  std::vector<bool> mCollapsedVerts;
  //! State array of 'active' edges
  std::vector<bool> mCollapsedEdges;
  //! State array of 'active' faces
  std::vector<bool> mCollapsedFaces;

  //! Number of collapsed verts
  size_t mNumCollapsedVerts;
  //! Number of collapsed edges
  size_t mNumCollapsedEdges;
  //! Number of collapsed faces
  size_t mNumCollapsedFaces;

  //! Utility mapping between half edges and collapses
  std::vector<EdgeCollapse *> mHalfEdge2EdgeCollapse;

  //! The heap that stores the edge collapses
  Heap mHeap;

  void drawText(const glm::vec3 &pos, const char *str);

  virtual bool save(std::ostream &os) {
    os << "# DecimationMesh obj streamer\n# M&A 2008\n\n";
    os << "# Vertices\n";
    for (size_t i = 0; i < GetNumVerts(); i++) {
      os << "v " << v(i).pos[0] << " " << v(i).pos[1] << " " << v(i).pos[2]
         << "\n";
    }
    os << "\n# Faces\n";
    for (size_t i = 0; i < GetNumFaces(); i++) {
      if (!isFaceCollapsed(i)) {
        size_t ind = f(i).edge;
        os << "f ";
        os << e(ind).vert + 1 << " ";
        ind = e(ind).next;
        os << e(ind).vert + 1 << " ";
        ind = e(ind).next;
        os << e(ind).vert + 1 << "\n";
      }
    }
    return os.good();
  }
};
