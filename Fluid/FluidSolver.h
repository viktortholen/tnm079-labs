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
#ifndef _FLUID_SOLVER_H
#define _FLUID_SOLVER_H

#include <Geometry/Implicit.h>
#include <Levelset/LevelSet.h>
#include <Math/CoordMatrix.h>
#include <Math/Function3D.h>
#include <Math/Volume.h>
#include <Util/Util.h>
#include <cassert>
#include <set>

class FluidSolver : public Function3D<glm::vec3> {
public:
    friend class FluidVoxelCutPlane;

    FluidSolver(float dx) : mDx(dx), mInitialVolume(0), mCurrentVolume(0), mExternalForces(NULL) {}

    void AddSolid(Implicit *impl);
    void AddFluid(LevelSet *LS);

    std::set<Implicit *> &GetSolids() { return mSolids; }
    std::set<LevelSet *> &GetFluids() { return mFluids; }

    void SetExternalForces(Function3D<glm::vec3> *forces) { mExternalForces = forces; }

    int Solve(float time);

    const Bbox &GetBoundingBox() const { return mBox; }

    float ComputeTimestep();

    float ComputePotentialEnergy();

    float ComputeKineticEnergy();

    virtual glm::vec3 GetValue(float x, float y, float z) const;

    virtual glm::vec3 GetMaxValue() const;

protected:
    std::set<Implicit *> mSolids;
    std::set<LevelSet *> mFluids;

    Bbox mBox;
    float mDx;
    float mInitialVolume, mCurrentVolume;

    Volume<glm::vec3> mVelocityField;
    Volume<float> mVoxels;
    Volume<bool> mSolidMask;
    Function3D<glm::vec3> *mExternalForces;

    void ExternalForces(float dt);
    void SelfAdvection(float dt, int steps);
    void EnforceDirichletBoundaryCondition();
    void Projection();
    void VelocityExtension();

    bool IsSolid(size_t i, size_t j, size_t k) const;
    bool IsFluid(size_t i, size_t j, size_t k) const;

    void ClassifyVoxels();
    void ClassifyVoxel(int i, int j, int k);

    inline void TransformGridToWorld(int i, int j, int k, float &x, float &y, float &z) {
        x = mBox.pMin[0] + i * mDx;
        y = mBox.pMin[1] + j * mDx;
        z = mBox.pMin[2] + k * mDx;
    }
};

/*
 * Below follows a set of operators used in the conjugate gradient solver
 */

template <typename Real>
Real dot(const std::vector<Real> &v1, const std::vector<Real> &v2) {
    assert(v1.size() == v2.size());

    Real sum = 0;
    const size_t size = v1.size();
    for (size_t i = 0; i < size; i++) sum += v1[i] * v2[i];

    return sum;
}

template <typename Real>
Real norm(const std::vector<Real> &v) {
    return std::sqrt(dot(v, v));
}

template <typename Real>
std::vector<Real> operator-(const std::vector<Real> &v1, const std::vector<Real> &v2) {
    assert(v1.size() == v2.size());

    const size_t size = v1.size();
    std::vector<Real> v(size);
    for (size_t i = 0; i < size; i++) v[i] = v1[i] - v2[i];

    return v;
}

template <typename Real>
std::vector<Real> operator+(const std::vector<Real> &v1, const std::vector<Real> &v2) {
    assert(v1.size() == v2.size());

    const size_t size = v1.size();
    std::vector<Real> v(size);
    for (size_t i = 0; i < size; i++) v[i] = v1[i] + v2[i];

    return v;
}

template <typename Real>
std::vector<Real> operator*(Real r, const std::vector<Real> &v1) {
    const size_t size = v1.size();
    std::vector<Real> v(size);
    for (size_t i = 0; i < size; i++) v[i] = v1[i] * r;

    return v;
}

#include "Math/ConjugateGradient.h"

#endif
