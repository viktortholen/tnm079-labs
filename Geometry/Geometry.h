#pragma once

#include <GUI/GLObject.h>
#include <gtx/transform.hpp>

/*!  \brief Geometry base class */
class Geometry : public GLObject {
public:
    Geometry() {}
    virtual ~Geometry(){};

    //! Update geometry, e.g. differential quantities after modifications
    virtual void Update() = 0;

    //! Initialize geometry
    virtual void Initialize() = 0;

    //! Translate geometry
    void Translate(float x, float y, float z) {
        SetTransform(mTransform * glm::translate(glm::vec3{x, y, z}));
    }

    //! Scale geometry
    void Scale(float s) { SetTransform(mTransform * glm::scale(glm::vec3{s, s, s})); }

    //! Scale geometry
    void Scale(float x, float y, float z) {
        SetTransform(mTransform * glm::scale(glm::vec3{x, y, z}));
    }

    //! Rotate geometry
    void Rotate(float rx, float ry, float rz) {
        SetTransform(mTransform * glm::rotate(rx, glm::vec3{1, 0, 0}));
        SetTransform(mTransform * glm::rotate(ry, glm::vec3{0, 1, 0}));
        SetTransform(mTransform * glm::rotate(rz, glm::vec3{0, 0, 1}));
    }

    //! Dilate geometry
    virtual void Dilate(float amount) {
        std::cout << "Dilate() not implemented for this type of geometry" << std::endl;
    }

    //! Erode geometry
    virtual void Erode(float amount) {
        std::cout << "Erode() not implemented for this type of geometry" << std::endl;
    }

    //! Smooth geometry
    virtual void Smooth(float amount) {
        std::cout << "Smooth() not implemented for this type of geometry" << std::endl;
    }

    //! Set transformation
    virtual void SetTransform(const glm::mat4& transform) { mTransform = transform; }

    //! Get transform
    auto& GetTransform() { return mTransform; }
    const auto& GetTransform() const { return mTransform; }

protected:
    glm::mat4 mTransform;
};
