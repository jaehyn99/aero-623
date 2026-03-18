#ifndef LINEAR_ELEMENT_H
#define LINEAR_ELEMENT_H

#include "Element.h"
class LinearElement: public Element{
    public:
    friend class TriangularMesh;
    LinearElement(const Eigen::Vector3i&, const Eigen::Vector3i&, double, const Eigen::Matrix2d&);
    virtual ~LinearElement() = default;

    bool isCurvedElement() const noexcept override{ return false; };
    Eigen::Matrix2d jacobian(std::size_t) const noexcept override{ return _J; };
    double detJacobian(std::size_t) const noexcept override{ return _detJ; };

    protected:
    Eigen::Matrix2d _J;
    double _detJ;
};

#endif