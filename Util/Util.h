#pragma once

#include <algorithm>  //for transform
#include <cmath>
#include <functional>
#include <glm.hpp>
#include <iostream>
#include <string>  // for string
#include <sstream>

//! Load a vector from an istream
template <typename T>
void Load(std::istream &is, glm::tvec3<T> &v) {
    is >> v[0] >> v[1] >> v[2];
}

//! Constructor from std::string object
template <typename T>
glm::tvec3<T> vec3(std::string &str) {
    glm::tvec3<T> v;

    std::stringstream s;
    s << str;
    Load(s, v);

    return v;
}

template <typename T>
struct glm_vec3_comparator {
    bool operator()(const glm::tvec3<T> &a, const glm::tvec3<T> &b) const {
        if (a[0] < b[0])
            return true;
        else if (a[0] > b[0])
            return false;
        else if (a[1] < b[1])
            return true;
        else if (a[1] > b[1])
            return false;
        else if (a[2] < b[2])
            return true;
        else
            return false;
    }
};

/*! Stable computation of cotangents */
inline float Cotangent(const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c) {
    const glm::vec3 ba = a - b;
    const glm::vec3 bc = c - b;

    return glm::dot(bc, ba) / glm::length(glm::cross(bc, ba));
}

/*! \brief Return the sign of a value
 *
 *  \f{eqnarray*}
 *  \begin{cases} -1& x<0\\ 1& x\geq 0 \end{cases}
 *  \f}
 *
 * \param[in] value the value to get the sign of
 * \return the sign as a T
 */
template <typename Scalar>
inline Scalar Sign(Scalar value) {
    Scalar mysign = -(Scalar)1.0;
    if (value >= (Scalar)0.0) {
        mysign = (Scalar)1.0;
    }
    return mysign;
}

/*! \brief First order switch function
 *
 * Maps the input data [x1, x2] to [0,1]. Outside values are clamped to fit in
 * [0,1]. \param[in] val the value to switch \param[in] x1 the lower map border
 * \param[in] x2 the higher map border
 * \return a value in [0,1] as a Scalar
 */
template <class Scalar>
inline Scalar Switch1(Scalar val, Scalar x1, Scalar x2) {
    return (val <= x1 ? 0 : (val >= x2 ? 1 : (val - x1) / (x2 - x1)));
}

/*! \brief First order root finding
 *
 *  \f{eqnarray*}
 *  f(0) = y_0 \\
 *  f(1) = y_1 \\
 *  f(x) = kx + m\\
 *  f(x) = 0, \quad \rightarrow \quad x = -m/k
 *  \f}

 * \param[in] y0 the value of \f$f(0)\f$
 * \param[in] y1 the value of \f$f(1)\f$
 * \return the value for which \f$f(x)=0\f$
 */
inline float Root(float y0, float y1) {
    const float m = y0;
    const float k = y1 - y0;
    return -m / k;
}

inline float Round(float d) { return std::floor(d + 0.5f); }

//! Test to see if host has big endian byte order
inline bool IsBigEndian() {
    short int word = 0x001;
    char *byte = (char *)&word;
    return (byte[0] ? false : true);
}

//! Byte swap 32bit uint
inline uint32_t EndianSwap(uint32_t x) {
    return ((x >> 24) | ((x << 8) & 0x00FF0000) | ((x >> 8) & 0x0000FF00) | (x << 24));
}

//! Byte swap array, only works for 32bit types
template <typename T>
inline void EndianSwap(T *data, size_t num) {
    if (sizeof(T) == sizeof(uint32_t)) {
        auto tmpPtr = reinterpret_cast<uint32_t*>(data);
        std::transform(tmpPtr, tmpPtr + num, tmpPtr, [](const auto x) {
            return ((x >> 24) | ((x << 8) & 0x00FF0000) | ((x >> 8) & 0x0000FF00) | (x << 24));
        });
    } else {
        // badness
    }
}