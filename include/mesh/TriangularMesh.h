#ifndef TRIANGULAR_MESH_H
#define TRIANGULAR_MESH_H

#include "Eigen/Dense"
#include <string>

class CubicSpline;
class Element;
class Face;
class TriangularMesh{
    public:
    /*
        p = number of Lagrange nodes for solution estimation
        q = number of Lagrange nodes for curved edge construction
        r = quadrature order of accuracy
    */
    TriangularMesh(const std::string& file_name, std::size_t p, std::size_t q, std::size_t r);
    TriangularMesh(TriangularMesh&&);
    ~TriangularMesh();
    TriangularMesh& operator=(TriangularMesh&&);

    std::size_t numNodes() const noexcept{ return _nodes.size(); }
    std::size_t numFaces() const noexcept{ return _faces.size(); }
    std::size_t numElems() const noexcept{ return _elems.size(); }

    Eigen::Vector2d& node(std::size_t nodeID) noexcept{ return _nodes[nodeID]; }
    const Eigen::Vector2d& node(std::size_t nodeID) const noexcept{ return _nodes[nodeID]; }
    Face& face(std::size_t faceID) noexcept{ return *_faces[faceID]; }
    const Face& face(std::size_t faceID) const noexcept{ return *_faces[faceID]; }
    Element& elem(std::size_t elemID) noexcept{ return *_elems[elemID]; }
    const Element& elem(std::size_t elemID) const noexcept{ return *_elems[elemID]; }

    const auto& getNodes() const noexcept{ return _nodes; }
    const auto& getFaces() const noexcept{ return _faces; }
    const auto& getElements() const noexcept{ return _elems; }

    double length(std::size_t faceID) const noexcept;
    double length(std::size_t elemID, std::size_t localFaceID) const noexcept;
    double area(std::size_t elemID) const noexcept;
    // Eigen::Vector2d centroid(std::size_t elemID){ return _elems[elemID]._centroid; }

    // I2E and B2E are both zero-based, add one if they need to be one-based
    const auto& I2E() const noexcept{ return _I2E; }
    const auto& B2E() const noexcept{ return _B2E; }
    const auto& In() const noexcept{ return _In; }
    const auto& Bn() const noexcept{ return _Bn; }

    protected:
    std::vector<Eigen::Vector2d> _nodes;
    std::vector<std::unique_ptr<Face>> _faces;
    std::vector<std::unique_ptr<Element>> _elems;
    std::unique_ptr<CubicSpline> _lower, _upper;

    Eigen::Array<int, Eigen::Dynamic, 4, Eigen::RowMajor> _I2E;
    Eigen::Array<int, Eigen::Dynamic, 3, Eigen::RowMajor> _B2E;
    Eigen::Array<double, Eigen::Dynamic, 2, Eigen::RowMajor> _In;
    Eigen::Array<double, Eigen::Dynamic, 2, Eigen::RowMajor> _Bn;

    std::vector<std::string> split(std::string& str) const noexcept;
};

#endif