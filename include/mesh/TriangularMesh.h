#ifndef TRIANGULAR_MESH_H
#define TRIANGULAR_MESH_H

#include "Eigen/Dense"
#include <string>
#include <vector>

class TriangularMesh{
    public:
    struct Face{
        std::array<int, 2> _pointID;
        std::array<int, 2> _elemID{-1, -1};
        double _length;
        std::size_t _nf; // number of linear nodes
        std::string _title;
        int _periodicFaceID = -1;
        int _periodicElemID = -1;
        Eigen::Vector2d _normal;

        Face(const TriangularMesh& mesh, std::size_t pointID1, std::size_t pointID2, std::size_t nf=2, const std::string& title="");

        bool isBoundaryFace() const noexcept{ return _title != ""; }
        bool isPeriodicFace() const noexcept{ return _periodicFaceID != -1; }
        bool operator==(const Face& other) const noexcept;
    };

    struct Element{
        Eigen::Vector3i _pointID;
        Eigen::Vector3i _faceID;
        std::size_t _order;
        std::string _basis;
        double _area;
        Eigen::Vector2d _centroid;

        Element(TriangularMesh& mesh, std::size_t pointID1, std::size_t pointID2, std::size_t pointID3, std::size_t ord=1, const std::string& basis="TriLagrange");
    };

    TriangularMesh(const std::string& file_name);

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
    Eigen::Vector2d normal(std::size_t elemID, std::size_t localFaceID) const noexcept;
    Eigen::Vector2d centroid(std::size_t elemID){ return _elems[elemID]._centroid; }

    const auto& I2E() const noexcept{ return _I2E; }
    const auto& B2E() const noexcept{ return _B2E; }
    const auto& In() const noexcept{ return _In; }
    const auto& Bn() const noexcept{ return _Bn; }

    protected:
    std::vector<Eigen::Vector2d> _nodes;
    std::vector<Face> _faces;
    std::vector<Element> _elems;

    Eigen::Matrix<int, Eigen::Dynamic, 4, Eigen::RowMajor> _I2E;
    Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> _B2E;
    Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> _In;
    Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> _Bn;

    std::vector<std::string> split(std::string& str) const noexcept;
    void fillMatrices() noexcept;
};

#endif