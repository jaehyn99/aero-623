#ifndef CURVED_FACE_H
#define CURVED_FACE_H

#include "Face.h"
class CurvedFace: public Face{
    public:
    friend class TriangularMesh;
    CurvedFace(const Eigen::Vector2i&, double, int, std::string="");

    bool isCurvedFace() const noexcept override{ return true; }
    Eigen::Vector2d normal(std::size_t q) const noexcept override{ return _n.col(q); }
    double detJ(std::size_t q) const noexcept override{ return _detJ[q]; }

    protected:
    Eigen::Matrix2Xd _xL; // internal Lagrange nodes for geometry approximation
    Eigen::Matrix2Xd _xq; // position of each quadrature points
    Eigen::Matrix2Xd _n;  // normal vectors at each quadrature points
    Eigen::VectorXd _detJ; // Jacobian determinant at each quadrature points
};

#endif