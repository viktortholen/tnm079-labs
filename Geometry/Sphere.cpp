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
#include <Geometry/Sphere.h>

Sphere::Sphere(float r, bool euclideanDistance) : mEuclideanDistance(euclideanDistance) {
    this->radius2 = r * r;
    this->mBox = Bbox(glm::vec3(-r, -r, -r), glm::vec3(r, r, r));
}

Sphere::~Sphere() {}

float Sphere::GetValue(float x, float y, float z) const {
    TransformW2O(x, y, z);

    if (mEuclideanDistance)
        return sqrt(x * x + y * y + z * z) - sqrt(radius2);
    else
        return (x * x + y * y + z * z - radius2);
}
