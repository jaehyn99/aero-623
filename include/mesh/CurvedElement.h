#ifndef CURVED_ELEMENT_H
#define CURVED_ELEMENT_H

#include "Element.h"
class CurvedElement: public Element{
    public:
    friend class TriangularMesh;
    CurvedElement(const Eigen::Vector3i&, const Eigen::Vector3i&, double);

    bool isCurvedElement() const noexcept override{ return true; }
    Eigen::Matrix2d edgeJacobianMatrix(std::size_t q) const noexcept override{ return _Je[q]; }
    double edgeJacobianDeterminant(std::size_t q) const noexcept override{ return _detJe[q]; }
    Eigen::Matrix2d internalJacobianMatrix(std::size_t q) const noexcept override{ return _Ji[q]; }
    double internalJacobianDeterminant(std::size_t q) const noexcept override {return _detJi[q]; }

    protected:
    // Eigen::MatrixXd normals;                    // (nQe x 2) unit outward normals at edge quad points
    std::vector<Eigen::Matrix2d> _Je;              // (nQe) arc length Jacobian at edge quad points
    Eigen::VectorXd _detJe;                        // (nQe) det(J) at edge quad points
    std::vector<Eigen::Matrix2d> _Ji;              // (nQi) det(J) at internal quad points
    Eigen::VectorXd _detJi;                        // (nQi) det(J) at internal quad points
};

#endif