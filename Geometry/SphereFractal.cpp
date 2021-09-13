#include <Geometry/SphereFractal.h>
#include <fstream>

SphereFractal::SphereFractal(int level) {

    // Reading fractal from file..
    std::ifstream inFile;

    if (level < 2) {
        level = 2;
    }
    if (level > 5) {
        level = 5;
    }

    if (level == 2) {
        inFile.open("SupportCode/balls2levels.txt");
    } else if (level == 3) {
        inFile.open("SupportCode/balls3levels.txt");
    } else if (level == 4) {
        inFile.open("SupportCode/balls4levels.txt");
    } else if (level == 5) {
        inFile.open("SupportCode/balls5levels.txt");
    }

    while (inFile.good()) {
        float tx;
        inFile >> tx;

        if (!inFile.good()) break;

        float ty;
        inFile >> ty;

        if (!inFile.good()) break;

        float tz;
        inFile >> tz;

        if (!inFile.good()) break;

        float r;
        inFile >> r;

        if (!inFile.good()) break;

        // Building transform (position and scale the unit sphere)
        glm::mat4 transform{1};
        transform = transform * glm::translate(transform, glm::vec3{tx, ty, tz});
        transform = transform * glm::scale(transform, glm::vec3{r, r, r});

        // Add sphere to array of spheres
        Implicit *sphere = new Sphere(1.0);
        sphere->SetTransform(transform);
        mSpheres.push_back(sphere);
    }

    // Building fractal
    mFractal = buildFractal();

    // Setting bounding box for sphere fractal object
    Bbox box = mFractal->GetBoundingBox();

    float maxDist = std::abs(box.pMax[0] - box.pMin[0]);
    maxDist = std::max(maxDist, std::abs(box.pMax[1] - box.pMin[1]));
    maxDist = std::max(maxDist, std::abs(box.pMax[2] - box.pMin[2]));

    glm::vec3 pad(0.05f, 0.05f, 0.05f);
    box.pMin -= maxDist * pad;
    box.pMax += maxDist * pad;
    SetBoundingBox(box);
}

SphereFractal::~SphereFractal() {
    for (size_t i = 0; i < mSpheres.size(); i++) {
        delete mSpheres[i];
    }
}

Implicit *SphereFractal::buildFractal() {
    Union *u = new Union(mSpheres[0], mSpheres[1]);

    for (size_t i = 2; i < mSpheres.size(); i++) {
        u = new Union(u, mSpheres[i]);
    }

    return u;
}

float SphereFractal::GetValue(float x, float y, float z) const {
    return mFractal->GetValue(x, y, z);
}
