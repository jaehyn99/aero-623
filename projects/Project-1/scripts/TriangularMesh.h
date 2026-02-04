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
        std::size_t _nf; // number of linear nodes
        std::string _title;
        int _periodicFaceID = -1;
        int _periodicElemID = -1;

        Face(std::size_t pointID1, std::size_t pointID2, std::size_t nf=2, const std::string& title=""):
            _pointID{int(pointID1), int(pointID2)}, _nf(nf), _title(title)
        {}

        bool isBoundaryFace() const noexcept{ return _title != ""; };
        bool operator==(const Face& other) const noexcept;
    };

    struct Element{
        std::array<int, 3> _pointID;
        std::array<int, 3> _faceID;
        std::size_t _order;
        std::string _basis;

        Element(std::size_t pointID1, std::size_t pointID2, std::size_t pointID3, std::size_t ord=1, const std::string& basis="TriLagrange"):
            _pointID{int(pointID1), int(pointID2), int(pointID3)}, _order(ord), _basis(basis) {}
    };

    TriangularMesh(const std::string& file_name);

    const std::vector<Eigen::Vector2d>& getNodes() const noexcept{ return _nodes; }
    const std::vector<Face>& getFaces() const noexcept{ return _faces; }
    const std::vector<Element>& getElements() const noexcept{ return _elems; }

    Eigen::Vector2d vect(std::size_t faceID) const noexcept;
    double length(std::size_t faceID) const noexcept;
    double area(std::size_t elemID) const noexcept;
    Eigen::Vector2d normal(std::size_t elemID, std::size_t localFaceID) const noexcept;

    void writeGri(const std::string& fileName) const noexcept;

    protected:
    std::vector<Eigen::Vector2d> _nodes;
    std::vector<Face> _faces;
    std::vector<Element> _elems;    

    std::vector<std::string> split(std::string& str) const noexcept;
};

#endif