#ifndef FACE_H
#define FACE_H

#include "Eigen/Dense"
class Face{
    public:
    friend class TriangularMesh;
    Face(const Eigen::Vector2i&, double, std::string="");

    bool isBoundaryFace() const noexcept{ return _title == "Curve1" || _title == "Curve3" || _title == "Curve5" || _title == "Curve7"; }
    bool isPeriodicFace() const noexcept{ return _title == "Curve2" || _title == "Curve4" || _title == "Curve6" || _title == "Curve8"; }
    virtual bool isCurvedFace() const noexcept = 0;
    bool operator==(const Face& other) const noexcept;

    Eigen::Vector2i pointID() const noexcept{ return _pointID; }
    Eigen::Vector2i elemID() const noexcept{ return _elemID; }
    int pointID(Eigen::Index i) const noexcept{ return _pointID[i]; }
    int elemID(Eigen::Index i) const noexcept{ return _elemID[i]; }
    double length() const noexcept{ return _length; }
    std::string title() const noexcept{ return _title; }
    virtual Eigen::Vector2d normal(std::size_t=0) const noexcept = 0; // get the normal vector at the specified quadrature point
    virtual double detJ(std::size_t=0) const noexcept = 0; // get the Jacobian determinant (ds/dsigma) at the specified quadrature point
    
    protected:
    Eigen::Vector2i _pointID; // point indices
    Eigen::Vector2i _elemID = Eigen::Vector2i::Constant(-1); // neighboring element indices, fill in during construction of TriangularMesh
    double _length;
    std::string _title;
};

#endif