#include <Subdivision/UniformCubicSpline.h>
#include <glm.hpp>
#include <gtc/type_ptr.hpp>

UniformCubicSpline::UniformCubicSpline(const std::vector<glm::vec3> &joints, glm::vec3 lineColor,
                                       float lineWidth, float segmentLength)
    : mCoefficients(joints), mControlPolygon(joints) {
    this->mLineColor = lineColor;
    this->mLineWidth = lineWidth;
    this->mDt = segmentLength;
}

/*! The BSpline value is calculated from one of the four cardinal BSpline
 * segments
 */
float UniformCubicSpline::GetBSplineValue(size_t i, float t) {
    mBSplineEvaluations++;

    // Find the offset from spline i
    t = t - (float)i;
    if (std::abs(t) >= 2.0f) {
        // if outside of support return zero
        return 0.0f;
    } else if (0.0f <= t && t < 1.0f) {
        // get the fractional part
        float ft = t;
        return 1.0f / 6.0f * (3.0f * ft * ft * ft - 6.0f * ft * ft + 4.0f);
    } else if (1.0f <= t && t < 2.0f) {
        // get the fractional part
        float ft = t - 1.0f;
        return 1.0f / 6.0f * (-ft * ft * ft + 3.0f * ft * ft - 3.0f * ft + 1.0f);
    } else if (-1.0f <= t && t < 0.0f) {
        // get the fractional part
        float ft = t + 1.0f;
        return 1.0f / 6.0f * (-3.0f * ft * ft * ft + 3.0f * ft * ft + 3.0f * ft + 1.0f);
    } else if (-2 < t && t < -1) {
        // get the fractional part
        float ft = t + 2.0f;
        return 1.0f / 6.0f * (ft * ft * ft);
    }
    return 0.0f;
}

/*! Evaluate the spline as the sum of the coefficients times the bsplines */
glm::vec3 UniformCubicSpline::GetValue(float t) {
    glm::vec3 val;
    float sum = 0;
    for (size_t i = 0; i < mCoefficients.size(); i++) {
        float bval = GetBSplineValue(i, t);
        val += mCoefficients.at(i) * bval;
        sum += bval;
    }
    return val;
}

void UniformCubicSpline::Render() {
    // Apply transform
    glPushMatrix();  // Push modelview matrix onto stack

    // Convert transform-matrix to format matching GL matrix format
    // Load transform into modelview matrix
    glMultMatrixf(glm::value_ptr(mTransform));

    mControlPolygon.Render();

    // save line point and color states
    glPushAttrib(GL_POINT_BIT | GL_LINE_BIT | GL_CURRENT_BIT);

    // draw segments
    glLineWidth(mLineWidth);
    glColor3fv(glm::value_ptr(mLineColor));
    glBegin(GL_LINE_STRIP);

    // We only have full BSpline support from spline at index 1, thus we begin
    // evaluating at 1.0
    mBSplineEvaluations = 0;
    for (float i = 1; i < mCoefficients.size() - 2; i += mDt) {
        glVertex3fv(glm::value_ptr(this->GetValue(i)));
    }
    glEnd();

    std::cout << "Number of B-spline evaluations " << mBSplineEvaluations << std::endl;

    // restore attribs
    glPopAttrib();

    glPopMatrix();

    GLObject::Render();
}
