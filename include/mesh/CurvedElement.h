#ifndef CURVED_ELEMENT_H
#define CURVED_ELEMENT_H

#include "Element.h"
class CurvedElement: public Element{
    public:
    friend class TriangularMesh;
    CurvedElement(const Eigen::Vector3i&, const Eigen::Vector3i&, double);

    bool isCurvedElement() const noexcept override{ return true; }
    Eigen::Matrix2d jacobian(std::size_t q) const noexcept override{ return _J[q]; }
    double detJacobian(std::size_t q) const noexcept override {return _detJ[q]; }
    auto MLLT() const noexcept{ return _MLLT; }

    protected:
    Eigen::Vector2d _internal;                    // internal Lagrange node (q = 3)
    std::vector<Eigen::Matrix2d> _J;              // (nQi) det(J) at internal quad points
    Eigen::VectorXd _detJ;                        // (nQi) det(J) at internal quad points
    Eigen::LLT<Eigen::MatrixXd> _MLLT;            // LLT factorization of the mass matrix
};

#endif