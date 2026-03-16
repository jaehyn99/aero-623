#ifndef ELEMENT_H
#define ELEMENT_H

#include "Eigen/Dense"
class Element{
    public:
    friend class TriangularMesh;
    Element(const Eigen::Vector3i&, const Eigen::Vector3i&, double); //, const Eigen::Vector2d&);
    virtual ~Element() = default;

    Eigen::Vector3i pointID() const noexcept { return _pointID; }
    Eigen::Vector3i faceID() const noexcept { return _faceID; }
    double area() const noexcept{ return _area; }
    // Eigen::Vector2d centroid() const noexcept{ return _centroid; }

    virtual bool isCurvedElement() const noexcept = 0;
    virtual Eigen::Matrix2d edgeJacobianMatrix(std::size_t) const noexcept = 0;
    virtual double edgeJacobianDeterminant(std::size_t) const noexcept = 0;
    virtual Eigen::Matrix2d internalJacobianMatrix(std::size_t) const noexcept = 0;
    virtual double internalJacobianDeterminant(std::size_t) const noexcept = 0;

    protected:
    Eigen::Vector3i _pointID; // point indices, arranged in counter-clockwise order
    Eigen::Vector3i _faceID; // face indices, such that the first face is opposite the first point, etc
    double _area;
    Eigen::Vector2d _centroid;
};

#endif