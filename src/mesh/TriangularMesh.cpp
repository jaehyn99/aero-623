#include "mesh/TriangularMesh.h"
#include <algorithm>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>

TriangularMesh::Face::Face(const Eigen::Vector2i& pointID, double length, std::size_t Q, std::string title):
    _pointID(pointID),
    _length(length),
    _Q(Q),
    _title(std::move(title))
{}

bool TriangularMesh::Face::operator==(const Face& other) const noexcept{
    return (_pointID[0] == other._pointID[0] && _pointID[1] == other._pointID[1]) || (_pointID[0] == other._pointID[1] && _pointID[1] == other._pointID[0]);
}

TriangularMesh::Element::Element(const Eigen::Vector3i& pointID, const Eigen::Vector3i& faceID, double area, const Eigen::Vector2d& centroid, std::size_t ord, const std::string& basis):
    _pointID(pointID),
    _faceID(faceID),
    _order(ord),
    _basis(std::move(basis)),
    _area(area),
    _centroid(centroid)
{
    // 


}

TriangularMesh::TriangularMesh(const std::string& fileName, std::size_t p, std::size_t q, std::size_t r){
    std::ifstream f(fileName);
    std::string line;
    std::vector<std::string> v;

    auto splitNextLine = [this, &f, &line, &v]()
        {
            std::getline(f, line);
            v = split(line);
        };

    splitNextLine();
    std::size_t nNode = std::stoi(v[0]);
    std::size_t nElemTot = std::stoi(v[1]);

    // Fill in the nodes
    _nodes.reserve(nNode);
    for (std::size_t i = 0; i < nNode; i++){
        splitNextLine();
        _nodes.emplace_back(std::stod(v[0]), std::stod(v[1]));
    }

    // Fill in the boundary faces
    _faces.reserve(nElemTot); // not exact size, only an approximate
    splitNextLine();
    std::size_t nBGroup = std::stoi(v[0]); // number of boundary groups
    for (std::size_t i = 0; i < nBGroup; i++){
        splitNextLine();
        std::size_t nBFace = std::stoi(v[0]);
        std::size_t Q = 2;
        std::string title = v[2];
        if (title == "Curve1" || title == "Curve5") Q = q+1; // curved boundaries may have a q larger than 2
        for (std::size_t j = 0; j < nBFace; j++){
            splitNextLine();
            std::size_t ind1 = std::stoi(v[0]) - 1;
            std::size_t ind2 = std::stoi(v[1]) - 1;
            Eigen::Vector2i pointID{ind1, ind2};
            double length = (_nodes[ind1] - _nodes[ind2]).lpNorm<2>();
            _faces.emplace_back(pointID, length, Q, title);
        }
    }

    // Fill in the elements, interior nodes, and connectivity matrices
    Eigen::Index icount, bcount;
    _I2E = decltype(_I2E)::Constant(nElemTot, 4, -1);
    _B2E.resize(nElemTot, Eigen::NoChange);
    _elems.reserve(nElemTot);

    while (nElemTot > 0){
        splitNextLine();
        std::size_t nElem = std::stoi(v[0]);
        std::size_t ord = std::stoi(v[1]);
        std::string basis = v[2];
        nElemTot -= nElem;
        for (std::size_t i = 0; i < nElem; i++){
            splitNextLine();
            std::size_t ind1 = std::stoi(v[0]) - 1;
            std::size_t ind2 = std::stoi(v[1]) - 1;
            std::size_t ind3 = std::stoi(v[2]) - 1;
            Eigen::Vector3i pointID{ind1, ind2, ind3};
            Eigen::Vector3i faceID;
            Eigen::Vector3d length;
            Eigen::Vector2d centroid{ (_nodes[ind1].x() + _nodes[ind2].x() + _nodes[ind3].x())/3, (_nodes[ind1].y() + _nodes[ind2].y() + _nodes[ind3].y())/3 };

            for (std::size_t j = 0; j < 3; j++){
                // Look if its faces have already been added to the list of faces
                Eigen::Vector2i facePointID{pointID[(j+1)%3], pointID[(j+2)%3]};
                Face iface(facePointID, 0.0);
                auto it = std::find(_faces.cbegin(), _faces.cend(), iface);
                std::size_t ind = it - _faces.cend();
                faceID[j] = ind;

                // If not found, create a new face and add it to the list
                if (it == _faces.cend()){
                    iface._length = (_nodes[pointID[(j+1)%3]] - _nodes[pointID[(j+2)%3]]).lpNorm<2>();
                    length[j] = iface._length;
                    _faces.push_back(std::move(iface));
                } else length[j] = _faces[ind]._length;

                // Determine if the face is on the boundary
                if (_faces[ind].isBoundaryFace()){
                    _B2E(bcount, 0) = i;
                    _B2E(bcount, 1) = j;
                    _B2E(bcount, 2) = (_faces[ind]._title == "Curve1") ? 0 : 4;
                    bcount++;
                } else{
                    if (_I2E(icount, 0) == -1){
                        // This row on _I2E has not been viisted
                        _I2E(icount, 0) = i;
                        _I2E(icount, 1) = j;
                    } else{
                        _I2E(icount, 2) = i;
                        _I2E(icount, 3) = j;
                    }
                    icount++;
                }
            }

            double s = (length[0]+length[1]+length[2])/2;
            double area = std::sqrt(s*(s-length[0])*(s-length[1])*(s-length[2]));
            _elems.emplace_back(pointID, faceID, area, centroid, ord, basis);
        }
    }
    _I2E.conservativeResize(icount, Eigen::NoChange);
    _B2E.conservativeResize(bcount, Eigen::NoChange);

    // Periodic boundaries, manually, only works on this one mesh
    std::deque<std::reference_wrapper<Face>> curve2Faces, curve4Faces, curve6Faces, curve8Faces;
    for (Face& face: _faces){
        if (face._title == "Curve2") curve2Faces.push_back(std::ref(face));
        else if (face._title == "Curve4") curve4Faces.push_front(std::ref(face));
        else if (face._title == "Curve6") curve6Faces.push_back(std::ref(face));
        else if (face._title == "Curve8") curve8Faces.push_front(std::ref(face));
    }

    // Matching curve2 and curve4
    for (std::size_t i = 0; i < curve2Faces.size(); i++){
        Face& curve2Face = curve2Faces[i].get();
        Face& curve4Face = curve4Faces[i].get();
        std::size_t curve2ID = std::find(_faces.cbegin(), _faces.cend(), curve2Face) - _faces.cbegin();
        std::size_t curve4ID = std::find(_faces.cbegin(), _faces.cend(), curve4Face) - _faces.cbegin();

        if (_I2E(curve2ID, 0) < _I2E(curve4ID, 0)){
            // Left element is on row curve2ID
            _I2E(curve2ID, 2) = _I2E(curve4ID, 0);
            _I2E(curve2ID, 3) = _I2E(curve4ID, 1);
            _I2E.row(curve4ID) = _I2E.row(curve2ID);
        } else{
            // Left element is on row curve4ID
            _I2E(curve4ID, 2) = _I2E(curve2ID, 0);
            _I2E(curve4ID, 3) = _I2E(curve2ID, 1);
            _I2E.row(curve2ID) = _I2E.row(curve4ID);            
        }
    }

    // Matching curve6 and curve8
    for (std::size_t i = 0; i < curve6Faces.size(); i++){
        Face& curve6Face = curve6Faces[i].get();
        Face& curve8Face = curve8Faces[i].get();
        std::size_t curve6ID = std::find(_faces.cbegin(), _faces.cend(), curve6Face) - _faces.cbegin();
        std::size_t curve8ID = std::find(_faces.cbegin(), _faces.cend(), curve8Face) - _faces.cbegin();

        if (_I2E(curve6ID, 0) < _I2E(curve8ID, 0)){
            // Left element is on row curve2ID
            _I2E(curve6ID, 2) = _I2E(curve8ID, 0);
            _I2E(curve6ID, 3) = _I2E(curve8ID, 1);
            _I2E.row(curve8ID) = _I2E.row(curve6ID);
        } else{
            // Left element is on row curve4ID
            _I2E(curve8ID, 2) = _I2E(curve6ID, 0);
            _I2E(curve8ID, 3) = _I2E(curve6ID, 1);
            _I2E.row(curve6ID) = _I2E.row(curve8ID);            
        }
    }

    // Update normal vectors on each edge, always pointing from L to R
    _In.resize(icount, Eigen::NoChange);
    for (Eigen::Index i = 0; i < icount; i++){
        // Consturct a unit normal vector
        const Element& elem = _elems[_I2E(i,0)];
        std::size_t localFaceID = _I2E(i,1);
        const Face& face = _faces[elem._faceID[localFaceID]];
        Eigen::Vector2d edge = node(face._pointID[1]) - node(face._pointID[0]);
        Eigen::Vector2d normal = Eigen::Vector2d{-edge[1], edge[0]};
        normal.normalize();

        // Find the remaining point on the "left" element and direct normal away from it
        const Eigen::Vector2d& p1 = node(elem._pointID[localFaceID]); // Point not on this edge
        const Eigen::Vector2d& p2 = node(elem._pointID[(localFaceID+1)%3]); // One of the points on this edge
        if ((p1-p2).dot(normal) > 0) normal *= -1;
        _In.row(i) = normal;
    }

    _Bn.resize(icount, Eigen::NoChange);
    for (Eigen::Index i = 0; i < bcount; i++){
        // Consturct a unit normal vector
        const Element& elem = _elems[_B2E(i,0)];
        std::size_t localFaceID = _B2E(i,1);
        const Face& face = _faces[elem._faceID[localFaceID]];
        Eigen::Vector2d edge = node(face._pointID[1]) - node(face._pointID[0]);
        Eigen::Vector2d normal = Eigen::Vector2d{-edge[1], edge[0]};
        normal.normalize();

        // Find the remaining point on the "left" element and direct normal away from it
        const Eigen::Vector2d& p1 = node(elem._pointID[localFaceID]); // Point not on this edge
        const Eigen::Vector2d& p2 = node(elem._pointID[(localFaceID+1)%3]); // One of the points on this edge
        if ((p1-p2).dot(normal) > 0) normal *= -1;
        _Bn.row(i) = normal;
    }
}

// Eigen::Vector2d TriangularMesh::normal(std::size_t elemID, std::size_t localFaceID) const noexcept{
//     const Element& elem = _elems[elemID];
//     int faceID = elem._faceID[localFaceID];
//     const Face& face = _faces[faceID];

//     Eigen::Vector2d n = face._normal;
//     if (face._elemID[1] == int(elemID)) return -n; // on element R, revert the normal vector
//     return n; // on element L or on a boundary, periodic or not
// }

std::vector<std::string> TriangularMesh::split(std::string& str) const noexcept{
    std::vector<std::string> v;
    int count = 0;
    while (!str.empty()){
        count++;
        std::size_t ind = str.find(' ');
        if (ind == std::string::npos){
            v.push_back(str);
            str.erase(str.begin(), str.end());
        } else{
            v.push_back(str.substr(0, ind));
            str.erase(str.begin(), str.begin() + ind+1);
        }
    }
    return v;
}