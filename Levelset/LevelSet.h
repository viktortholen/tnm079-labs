#pragma once

#include <Geometry/Implicit.h>
#include <Levelset/LevelSetGrid.h>
#include <iostream>
#include <gtx/string_cast.hpp>

class LevelSet : public Implicit {
    friend class LevelSetOperator;

protected:
    LevelSetGrid mGrid;
    //! Grid spacing
    float mDx;

public:
    static const VisualizationMode NarrowBand;

    virtual std::list<VisualizationMode> GetVisualizationModes() {
        std::list<VisualizationMode> L = Implicit::GetVisualizationModes();
        L.push_back(NarrowBand);
        return L;
    }

    LevelSet(float dx);
    LevelSet(float dx, const Implicit &impl);
    LevelSet(float dx, const Implicit &impl, const Bbox &box);
    LevelSet(float dx, const Volume<float> &vol);

    virtual ~LevelSet() {}

    virtual const char *GetTypeName() { return typeid(LevelSet).name(); }

    //! Evaluate the implicit at world coordinates x y z
    virtual float GetValue(float x, float y, float z) const;

    //! Calculate the gradient of the implicit at world coordinates x y z
    virtual glm::vec3 GetGradient(float x, float y, float z) const;

    //! calculate the curvature of the implicit at world coordinates x y z
    virtual float GetCurvature(float x, float y, float z) const;

    //! Sets the bounding box in current frame coordinates
    virtual void SetBoundingBox(const Bbox &b);

    //! Set narrow band width (in number of grid points)
    virtual void SetNarrowBandWidth(int width);

    //! Get narrow band width (in number of grid points)
    int GetNarrowBandWidth();

    inline const float GetDx() const { return mDx; }

    //! First order negative differential in x
    float DiffXm(size_t i, size_t j, size_t k) const;
    //! First order positive differential in x
    float DiffXp(size_t i, size_t j, size_t k) const;
    //! Second order central differential in x
    float DiffXpm(size_t i, size_t j, size_t k) const;
    //! Second order second central differential in x
    float Diff2Xpm(size_t i, size_t j, size_t k) const;

    //! First order negative differential in y
    float DiffYm(size_t i, size_t j, size_t k) const;
    //! First order positive differential in y
    float DiffYp(size_t i, size_t j, size_t k) const;
    //! Second order central differential in y
    float DiffYpm(size_t i, size_t j, size_t k) const;
    //! Second order second central differential in y
    float Diff2Ypm(size_t i, size_t j, size_t k) const;

    //! First order negative differential in z
    float DiffZm(size_t i, size_t j, size_t k) const;
    //! First order positive differential in z
    float DiffZp(size_t i, size_t j, size_t k) const;
    //! Second order central differential in z
    float DiffZpm(size_t i, size_t j, size_t k) const;
    //! Second order second central differential in z
    float Diff2Zpm(size_t i, size_t j, size_t k) const;

    //! First order mixed differential in x,y
    float Diff2XYpm(size_t i, size_t j, size_t k) const;
    //! First order mixed differential in y,z
    float Diff2YZpm(size_t i, size_t j, size_t k) const;
    //! First order mixed differential in z,x
    float Diff2ZXpm(size_t i, size_t j, size_t k) const;

    float WENO(float v1, float v2, float v3, float v4, float v5) const;

    friend std::ostream &operator<<(std::ostream &os, const LevelSet &LS) {
        os << "Level set bounding box: " << glm::to_string(LS.mBox.pMin) << " -> "
           << glm::to_string(LS.mBox.pMax) << std::endl;
        os << "Level set dx: " << std::to_string(LS.mDx) << std::endl;
        os << LS.mGrid << std::endl;
        return os;
    }

    virtual void Render();

    void TransformWorldToGrid(float &i, float &j, float &k) const;

    void TransformGridToWorld(float &x, float &y, float &z) const;
};
