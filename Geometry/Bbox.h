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

#pragma once

#include <algorithm>
#include <iostream>
#include <glm.hpp>
#include <gtx/string_cast.hpp>

struct Bbox {
public:
    glm::vec3 pMin;
    glm::vec3 pMax;

    //! Default constructor initializes to 'zero' cube
    Bbox(const glm::vec3 &pmin = glm::vec3(0, 0, 0), const glm::vec3 &pmax = glm::vec3(0, 0, 0))
        : pMin(pmin), pMax(pmax) {}

    Bbox(float pmin, float pmax)
        : pMin(glm::vec3(pmin, pmin, pmin)), pMax(glm::vec3(pmax, pmax, pmax)) {}

    friend std::ostream &operator<<(std::ostream &os, const Bbox &b) {
        os << glm::to_string(b.pMin) << " -> " << glm::to_string(b.pMax);
        return os;
    }

    //! Union of two Bbox :es
    static Bbox BoxUnion(const Bbox &b1, const Bbox &b2) {
        Bbox b;
        b.pMin[0] = std::min(b1.pMin[0], b2.pMin[0]);
        b.pMin[1] = std::min(b1.pMin[1], b2.pMin[1]);
        b.pMin[2] = std::min(b1.pMin[2], b2.pMin[2]);

        b.pMax[0] = std::max(b1.pMax[0], b2.pMax[0]);
        b.pMax[1] = std::max(b1.pMax[1], b2.pMax[1]);
        b.pMax[2] = std::max(b1.pMax[2], b2.pMax[2]);
        return b;
    }

    static Bbox PointUnion(const Bbox &b, const glm::vec4 &v) {
        Bbox b2;
        b2.pMin[0] = std::min(b.pMin[0], v[0]);
        b2.pMin[1] = std::min(b.pMin[1], v[1]);
        b2.pMin[2] = std::min(b.pMin[2], v[2]);

        b2.pMax[0] = std::max(b.pMax[0], v[0]);
        b2.pMax[1] = std::max(b.pMax[1], v[1]);
        b2.pMax[2] = std::max(b.pMax[2], v[2]);
        return b2;
    }

    //! Intersection of two Bbox :es
    static Bbox BoxIntersection(const Bbox &b1, const Bbox &b2) {
        Bbox b;
        b.pMin[0] = std::max(b1.pMin[0], b2.pMin[0]);
        b.pMin[1] = std::max(b1.pMin[1], b2.pMin[1]);
        b.pMin[2] = std::max(b1.pMin[2], b2.pMin[2]);

        b.pMax[0] = std::min(b1.pMax[0], b2.pMax[0]);
        b.pMax[1] = std::min(b1.pMax[1], b2.pMax[1]);
        b.pMax[2] = std::min(b1.pMax[2], b2.pMax[2]);
        return b;
    }

    Bbox Transform(const glm::mat4 &t) const {
        Bbox b{};
        glm::vec4 v{};

        v = t * glm::vec4(pMin[0], pMin[1], pMin[2], 1.f);
        b.pMax = b.pMin = glm::vec3(v[0], v[1], v[2]);

        v = t * glm::vec4(pMax[0], pMin[1], pMin[2], 1.f);
        b = PointUnion(b, v);
        v = t * glm::vec4(pMax[0], pMax[1], pMin[2], 1.f);
        b = PointUnion(b, v);
        v = t * glm::vec4(pMin[0], pMax[1], pMin[2], 1.f);
        b = PointUnion(b, v);

        v = t * glm::vec4(pMin[0], pMin[1], pMax[2], 1.f);
        b = PointUnion(b, v);
        v = t * glm::vec4(pMax[0], pMin[1], pMax[2], 1.f);
        b = PointUnion(b, v);
        v = t * glm::vec4(pMax[0], pMax[1], pMax[2], 1.f);
        b = PointUnion(b, v);
        v = t * glm::vec4(pMin[0], pMax[1], pMax[2], 1.f);
        b = PointUnion(b, v);

        return b;
    }
};
