#ifndef FACE_H
#define FACE_H

#include "Eigen/Dense"
class Face{
    public:
    friend class TriangularMesh;
    Face(const Eigen::Vector2i&, double, std::size_t=2, std::string="");

    bool isBoundaryFace() const noexcept{ return _title == "Curve1" || _title == "Curve3" || _title == "Curve5" || _title == "Curve7"; }
    bool isPeriodicFace() const noexcept{ return _title == "Curve2" || _title == "Curve4" || _title == "Curve6" || _title == "Curve8"; }
    bool isCurvedFace() const noexcept{ return _Q > 2; }
    bool operator==(const Face& other) const noexcept;

    Eigen::Vector2i pointID() const noexcept{ return _pointID; }
    double length() const noexcept{ return _length; }
    std::size_t Q() const noexcept{ return _Q; }
    std::string title() const noexcept{ return _title; }
    
    protected:
    Eigen::Vector2i _pointID; // point indices
    Eigen::Vector2i _elemID; // neighboring element indices
    double _length;
    std::size_t _Q; // number of Lagrange nodes to approximate geometry
    std::string _title;
};

#endif