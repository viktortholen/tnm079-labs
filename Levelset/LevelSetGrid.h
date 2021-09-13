/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sˆderstrˆm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#ifndef __levelset_grid_h__
#define __levelset_grid_h__

#include <Math/Volume.h>
#include <Util/BitMask3D.h>
#include <iostream>
#include <limits>

class LevelSetGrid {
protected:
    Volume<float> mPhi;
    BitMask3D mMask;
    float mInsideConstant, mOutsideConstant;

public:
    LevelSetGrid(int dimX = 0, int dimY = 0, int dimZ = 0,
                 float insideConstant = -(std::numeric_limits<float>::max)(),
                 float outsideConstant = (std::numeric_limits<float>::max)())
        : mPhi(dimX, dimY, dimZ, outsideConstant)
        , mMask(dimX, dimY, dimZ)
        , mInsideConstant(insideConstant)
        , mOutsideConstant(outsideConstant) {}

    ~LevelSetGrid() {}

    class Iterator {
        friend class LevelSetGrid;

    protected:
        const BitMask3D *mask;
        size_t i, j, k;
        size_t iMax, jMax, kMax;
        bool endState;

        Iterator(const BitMask3D *mask, size_t i = 0, size_t j = 0, size_t k = 0);

    public:
        inline Iterator &operator++(int) {
            endState = false;
            if (endState != true) {
                do {
                    ++k;
                    if (k >= kMax) {
                        k = 0;
                        ++j;
                        if (j >= jMax) {
                            j = 0;
                            ++i;
                            if (i >= iMax) {
                                k = kMax;
                                j = jMax;
                                i = iMax;
                                endState = true;
                            }
                        }
                    }
                } while (!endState && (mask->GetValue(i, j, k) != true));
            }
            return *this;
        }

        bool operator!=(const Iterator &b) const {
            return (this->i != b.i || this->j != b.j || this->k != b.k);
        }

        auto GetI() const { return i; }
        auto GetJ() const { return j; }
        auto GetK() const { return k; }
    };

    Iterator BeginNarrowBand() { return Iterator(&mMask); }
    const Iterator BeginNarrowBand() const { return Iterator(&mMask); }

    Iterator EndNarrowBand() {
        return Iterator(&mMask, mMask.GetDimX(), mMask.GetDimY(), mMask.GetDimZ());
    }
    const Iterator EndNarrowBand() const {
        return Iterator(&mMask, mMask.GetDimX(), mMask.GetDimY(), mMask.GetDimZ());
    }

    inline auto GetDimX() const { return mPhi.GetDimX(); }
    inline auto GetDimY() const { return mPhi.GetDimY(); }
    inline auto GetDimZ() const { return mPhi.GetDimZ(); }

    //! Return grid dimensions as measured in number of grid cells
    glm::ivec3 GetDimensions();

    inline float GetValue(size_t i, size_t j, size_t k) const { return mPhi.GetValue(i, j, k); }
    inline void SetValue(size_t  i, size_t j, size_t k, float f) {
        SetMask(i, j, k, true);
        mPhi.SetValue(i, j, k, f);
    }

    inline bool GetMask(size_t  i, size_t  j, size_t  k) const { return mMask.GetValue(i, j, k); }
    inline void SetMask(size_t  i, size_t j, size_t k, bool b) { mMask.SetValue(i, j, k, b); }

    void SetInsideConstant(float insideConstant) { mInsideConstant = insideConstant; }
    inline const float GetInsideConstant() const { return mInsideConstant; }

    void SetOutsideConstant(float outsideConstant) { mOutsideConstant = outsideConstant; }
    inline const float GetOutsideConstant() const { return mOutsideConstant; }

    //! Dilates the narrow band with 6 connectivity
    void Dilate();

    //! Rebuild the narrow band by culling too large values from mask
    void Rebuild();

    friend std::ostream &operator<<(std::ostream &os, const LevelSetGrid &grid) {
        os << "Grid dimensions: "
           << "(" << grid.GetDimX() << "x" << grid.GetDimY() << "x" << grid.GetDimZ() << ")"
           << std::endl;
        Iterator iter = grid.BeginNarrowBand();
        Iterator iend = grid.EndNarrowBand();
        while (iter != iend) {
            auto i = iter.GetI();
            auto j = iter.GetJ();
            auto k = iter.GetK();
            os << "(" << i << "," << j << "," << k << "): " << grid.GetValue(i, j, k) << std::endl;
            iter++;
        }
        return os;
    }
};

#endif
