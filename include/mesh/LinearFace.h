#ifndef LINEAR_FACE_H
#define LINEAR_FACE_H

#include "Face.h"
class LinearFace: public Face{
    public:
    friend class TriangularMesh;
    LinearFace(const Eigen::Vector2i&, double, std::string="");

    bool isCurvedFace() const noexcept override{ return false; }
    Eigen::Vector2d normal(std::size_t) const noexcept override{ return _n; }
    // double detJ(std::size_t) const noexcept override{ return _length; }

    protected:
    Eigen::Vector2d _n;
};

#endif