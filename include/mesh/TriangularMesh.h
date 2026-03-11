#ifndef TRIANGULAR_MESH_H
#define TRIANGULAR_MESH_H

#include "Eigen/Dense"
#include <string>
#include <vector>

class TriangularMesh{
    public:
    struct Face{
        Eigen::Vector2i _pointID; // point indices
        double _length;
        std::size_t _Q; // number of linear nodes
        std::string _title;

        Face(const Eigen::Vector2i&, double, std::size_t=2, std::string="");

        bool isBoundaryFace() const noexcept{ return _title == "Curve1" || _title == "Curve5"; }
        bool isPeriodicFace() const noexcept{ return _title == "Curve2" || _title == "Curve4" || _title == "Curve6" || _title == "Curve8"; }
        bool isCurvedFace() const noexcept{ return _Q > 2; }
        bool operator==(const Face& other) const noexcept;
    };

    struct Element{
        Eigen::Vector3i _pointID; // point indices, arranged in counter-clockwise order
        Eigen::Vector3i _faceID; // face indices, such that the first face is opposite the first point, etc
        std::size_t _order;
        std::string _basis;
        double _area;
        Eigen::Vector2d _centroid;

        Element(const Eigen::Vector3i&, const Eigen::Vector3i&, double, const Eigen::Vector2d&, std::size_t ord=1, const std::string& basis="TriLagrange");
    };

    /*
        p = number of Lagrange nodes for solution estimation
        q = number of Lagrange nodes for curved edge construction
        r = quadrature order of accuracy
    */
    TriangularMesh(const std::string& file_name, std::size_t p, std::size_t q, std::size_t Nq);

    std::size_t numNodes() const noexcept{ return _nodes.size(); }
    std::size_t numFaces() const noexcept{ return _faces.size(); }
    std::size_t numElems() const noexcept{ return _elems.size(); }

    Eigen::Vector2d& node(std::size_t nodeID) noexcept{ return _nodes[nodeID]; }
    const Eigen::Vector2d& node(std::size_t nodeID) const noexcept{ return _nodes[nodeID]; }
    Face& face(std::size_t faceID) noexcept{ return _faces[faceID]; }
    const Face& face(std::size_t faceID) const noexcept{ return _faces[faceID]; }
    Element& elem(std::size_t elemID) noexcept{ return _elems[elemID]; }
    const Element& elem(std::size_t elemID) const noexcept{ return _elems[elemID]; }

    const std::vector<Eigen::Vector2d>& getNodes() const noexcept{ return _nodes; }
    const std::vector<Face>& getFaces() const noexcept{ return _faces; }
    const std::vector<Element>& getElements() const noexcept{ return _elems; }

    double length(std::size_t faceID) const noexcept{ return _faces[faceID]._length; }
    double area(std::size_t elemID) const noexcept{ return _elems[elemID]._area; }
    Eigen::Vector2d centroid(std::size_t elemID){ return _elems[elemID]._centroid; }

    // I2E and B2E are both zero-based, add one if they need to be one-based
    const auto& I2E() const noexcept{ return _I2E; }
    const auto& B2E() const noexcept{ return _B2E; }
    const auto& In() const noexcept{ return _In; }
    const auto& Bn() const noexcept{ return _Bn; }

    protected:
    std::vector<Eigen::Vector2d> _nodes;
    std::vector<Face> _faces;
    std::vector<Element> _elems;

    Eigen::Array<int, Eigen::Dynamic, 4, Eigen::RowMajor> _I2E;
    Eigen::Array<int, Eigen::Dynamic, 3, Eigen::RowMajor> _B2E;
    Eigen::Array<double, Eigen::Dynamic, 2, Eigen::RowMajor> _In;
    Eigen::Array<double, Eigen::Dynamic, 2, Eigen::RowMajor> _Bn;

    std::vector<std::string> split(std::string& str) const noexcept;
};

#endif