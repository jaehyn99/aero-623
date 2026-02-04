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

        bool isBoundaryFace() const noexcept{ return _title != ""; };
        bool operator==(const Face& other) const noexcept;
    };

    struct Element{
        std::array<int, 3> _pointID;
        std::array<int, 3> _faceID;
        std::size_t _order;
        std::string _basis;
        double _area;

        Element(TriangularMesh& mesh, std::size_t pointID1, std::size_t pointID2, std::size_t pointID3, std::size_t ord=1, const std::string& basis="TriLagrange");
    };

    TriangularMesh(const std::string& file_name);

    double numNodes() const noexcept{ return _nodes.size(); }
    double numFaces() const noexcept{ return _faces.size(); }
    double numElems() const noexcept{ return _elems.size(); }

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

    void writeGri(const std::string& fileName) const noexcept;

    protected:
    std::vector<Eigen::Vector2d> _nodes;
    std::vector<Face> _faces;
    std::vector<Element> _elems;    

    std::vector<std::string> split(std::string& str) const noexcept;
};

#endif