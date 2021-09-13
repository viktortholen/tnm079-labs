#pragma once

#include "Util/Util.h"
#include <cassert>
#include <cstdio>
#include <vector>
#include <glm.hpp>

/*!
 * A 3D volume of templated type T.
 * Stores values in an stl vector<T>. Uses row major storage, i.e. dimZ's is
 * changing most rapidly.
 */
template <class T>
class Volume {
protected:
    //! An stl vector to hold actual data
    std::vector<T> mData;
    //! The dimensions of the volume
    size_t mDimX;
    size_t mDimY;
    size_t mDimZ;
    //! premult = dimY*dimZ, avoids this multiplication for each getValue
    size_t mPremult;

public:
    //! Default constructor initializes to zero volume
    Volume() : mData(0), mDimX(0), mDimY(0), mDimZ(0), mPremult(0) {}
    //! Sized constructor initializes volume of size dimXxdimYxdimZ to default
    //! value for type T
    Volume(size_t dimX, size_t dimY, size_t dimZ)
        : mData(dimX * dimY * dimZ, T())
        , mDimX(dimX)
        , mDimY(dimY)
        , mDimZ(dimZ)
        , mPremult(dimY * dimZ) {}

    Volume(size_t dimX, size_t dimY, size_t dimZ, T defaultVal)
        : mData(dimX * dimY * dimZ, defaultVal)
        , mDimX(dimX)
        , mDimY(dimY)
        , mDimZ(dimZ)
        , mPremult(dimY * dimZ) {}

    inline auto GetDimX() const { return mDimX; }
    inline auto GetDimY() const { return mDimY; }
    inline auto GetDimZ() const { return mDimZ; }
    //! Returns the value at i,j,k
    inline T GetValue(size_t i, size_t j, size_t k) const {
        i = glm::clamp(i, size_t{0}, mDimX - 1);
        j = glm::clamp(j, size_t{0}, mDimY - 1);
        k = glm::clamp(k, size_t{0}, mDimZ - 1);
        // return mData.at(i*mPremult + j*mDimZ + k); // .at() does bound checking,
        // throws exception
        return mData[i * mPremult + j * mDimZ + k];  // op [] does no bound checking
    }

    //! Returns the value at x,y,z (uses trilinear interpolation
    T GetValue(float x, float y, float z) const {
        auto i = static_cast<size_t>(x);
        auto j = static_cast<size_t>(y);
        auto k = static_cast<size_t>(z);

        float bx = x - static_cast<float>(i);
        float by = y - static_cast<float>(j);
        float bz = z - static_cast<float>(k);

        T val = GetValue(i, j, k) * (1 - bx) * (1 - by) * (1 - bz) +
                GetValue(i + 1, j, k) * bx * (1 - by) * (1 - bz) +
                GetValue(i + 1, j + 1, k) * bx * by * (1 - bz) +
                GetValue(i, j + 1, k) * (1 - bx) * by * (1 - bz) +
                GetValue(i, j, k + 1) * (1 - bx) * (1 - by) * bz +
                GetValue(i + 1, j, k + 1) * bx * (1 - by) * bz +
                GetValue(i + 1, j + 1, k + 1) * bx * by * bz +
                GetValue(i, j + 1, k + 1) * (1 - bx) * by * bz;

        return val;
    }

    size_t ComputeLinearIndex(size_t i, size_t j, size_t k) const {
        i = glm::clamp(i, size_t{0}, mDimX - 1);
        j = glm::clamp(j, size_t{0}, mDimY - 1);
        k = glm::clamp(k, size_t{0}, mDimZ - 1);
        return i * mPremult + j * mDimZ + k;
    }

    void ComputeIndices(size_t ind, int &i, int &j, int &k) const {
        i = i / mPremult;
        j = (i % mPremult) / mDimZ;
        k = (i % mPremult) % mDimZ;
    }

    //! Sets the value at i,j,k to val
    void SetValue(size_t i, size_t j, size_t k, const T &val) {
        assert(i < mDimX && i >= 0 && j < mDimY && j >= 0 && k < mDimZ && k >= 0);
        mData.at(i * mPremult + j * mDimZ + k) = val;  // .at() does bound checking, throws exception
                                                       // mData[i*premult + j*dimZ + k] = val;
    }

    //! Load a volume from binary stream is
    void Load(std::istream &is) {
        is.read((char *)&mDimX, sizeof(mDimX));
        is.read((char *)&mDimY, sizeof(mDimY));
        is.read((char *)&mDimZ, sizeof(mDimZ));

        if (IsBigEndian()) {
            mDimX = EndianSwap(mDimX);
            mDimY = EndianSwap(mDimY);
            mDimZ = EndianSwap(mDimZ);
        }

        std::cerr << "Loading dims: " << mDimX << ", " << mDimY << ", " << mDimZ << std::endl;
        mData.resize(mDimX * mDimY * mDimZ);  // resize (only upwards)
        std::vector<T>(mData).swap(mData);    // shrink to fit (if needed)
        mPremult = mDimY * mDimZ;

        is.read((char *)&mData.at(0), sizeof(T) * mDimX * mDimY * mDimZ);
        if (IsBigEndian()) {
            EndianSwap(&mData.at(0), mData.size());
        }
    }

    void Load(const std::string &path) {
        FILE *file;
        if ((file = fopen(path.c_str(), "rb")) != NULL) {
            int w;
            fscanf(file, "%d \n", &w);
            if (w == 4) {
                fscanf(file, "%u %u %u \n", &mDimZ, &mDimY, &mDimX);
            } else {
                fclose(file);
                std::cerr << "Error: Only float volumes supported, file bit depth = " << w
                          << std::endl;
                return;
            }

            std::cerr << "Loading dims: " << mDimX << ", " << mDimY << ", " << mDimZ << std::endl;
            mData.resize(mDimX * mDimY * mDimZ);  // resize (only upwards)
            std::vector<T>(mData).swap(mData);    // shrink to fit (if needed)
            mPremult = mDimY * mDimZ;

            fread((void*)mData.data(), sizeof(T), mDimX * mDimY * mDimZ, file);
            fclose(file);
            if (IsBigEndian()) {
                EndianSwap(&mData.at(0), mData.size());
            }
        }
    }

    //! Save volume to binary stream os
    void Save(std::ostream &os) const {
        int x = static_cast<int>(mDimX);
        int y = static_cast<int>(mDimX);
        int z = static_cast<int>(mDimX);

        os.write((const char *)&x, sizeof(x));
        os.write((const char *)&y, sizeof(y));
        os.write((const char *)&z, sizeof(z));

        os.write((const char *)mData.data(), sizeof(T) * x * y * z);
    }
};
