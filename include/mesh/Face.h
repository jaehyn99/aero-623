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
    double length() const noexcept{ return _length; }
    // std::size_t Q() const noexcept{ return _Q; }
    std::string title() const noexcept{ return _title; }
    virtual Eigen::Vector2d normal(std::size_t) const noexcept = 0; // get the normal vector at the specified quadrature point
    
    protected:
    Eigen::Vector2i _pointID; // point indices
    Eigen::Vector2i _elemID; // neighboring element indices
    double _length;
    // std::size_t _Q; // number of Lagrange nodes to approximate geometry
    std::string _title;
};

#endif