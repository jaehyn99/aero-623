#include "TriangularMesh.h"
#include "LinearElement.h"

#include <algorithm>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>

// TriangularMesh::Face::Face(const Eigen::Vector2i& pointID, double length, std::size_t Q, std::string title):
//     _pointID(pointID),
//     _length(length),
//     _Q(Q),
//     _title(std::move(title))
// {}

// TriangularMesh::Element::Element(const Eigen::Vector3i& pointID, const Eigen::Vector3i& faceID, double area, const Eigen::Vector2d& centroid, std::size_t ord, const std::string& basis):
//     _pointID(pointID),
//     _faceID(faceID),
//     _order(ord),
//     _basis(std::move(basis)),
//     _area(area),
//     _centroid(centroid)
// {}

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
        std::string title = v[2];
        std::size_t Q = 2;
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
    Eigen::Index bcount = 0;
    _elems.reserve(nElemTot);
    _B2E.resize(_faces.size(), Eigen::NoChange); // _faces only contains boundary and periodic bases, for now.

    std::vector<int> elemL(_faces.size(), -1);
    std::vector<int> faceL(_faces.size(), -1);
    std::vector<int> elemR(_faces.size(), -1);
    std::vector<int> faceR(_faces.size(), -1);

    elemL.reserve(nElemTot);
    faceL.reserve(nElemTot);
    elemR.reserve(nElemTot);
    faceR.reserve(nElemTot);    

    while (nElemTot > 0){
        splitNextLine();
        std::size_t nElem = std::stoi(v[0]);
        // std::size_t ord = std::stoi(v[1]);
        // std::string basis = v[2];
        nElemTot -= nElem;
        for (std::size_t i = 0; i < nElem; i++){
            splitNextLine();
            std::size_t ind1 = std::stoi(v[0]) - 1;
            std::size_t ind2 = std::stoi(v[1]) - 1;
            std::size_t ind3 = std::stoi(v[2]) - 1;
            Eigen::Vector3i pointID{int(ind1), int(ind2), int(ind3)};
            Eigen::Vector3i faceID;
            Eigen::Vector3d length;
            //Eigen::Vector2d centroid{ (_nodes[ind1].x() + _nodes[ind2].x() + _nodes[ind3].x())/3, (_nodes[ind1].y() + _nodes[ind2].y() + _nodes[ind3].y())/3 };

            bool isCurvedElement = false;
            for (std::size_t j = 0; j < 3; j++){
                // Look if its faces have already been added to the list of faces
                Eigen::Vector2i facePointID{pointID[(j+1)%3], pointID[(j+2)%3]};
                Face iface(facePointID, 0.0);
                auto it = std::find(_faces.cbegin(), _faces.cend(), iface);
                std::size_t ind = it - _faces.cbegin();
                faceID[j] = ind;

                // If not found, create a new face and add it to the list
                if (it == _faces.cend()){
                    iface._length = (_nodes[pointID[(j+1)%3]] - _nodes[pointID[(j+2)%3]]).lpNorm<2>();
                    length[j] = iface._length;
                    _faces.push_back(std::move(iface));

                    elemL.push_back(-1);
                    faceL.push_back(-1);
                    elemR.push_back(-1);
                    faceR.push_back(-1);
                } else length[j] = _faces[ind]._length;
                isCurvedElement = isCurvedElement || _faces[ind].isCurvedFace();

                // Determine if the face is on the boundary, then update B2E because it's easier
                if (_faces[ind].isBoundaryFace()){
                    _B2E(bcount, 0) = i;
                    _B2E(bcount, 1) = j;
                    if (_faces[ind]._title == "Curve1") _B2E(bcount, 2) = 0;
                    else if (_faces[ind]._title == "Curve3") _B2E(bcount, 2) = 1;
                    else if (_faces[ind]._title == "Curve5") _B2E(bcount, 2) = 2;
                    else _B2E(bcount, 2) = 3;
                    bcount++;
                } else{
                    // Stores these to update I2E later
                    if (elemL[ind] == -1){
                        // Face is new and hasn't been added
                        elemL[ind] = i;
                        faceL[ind] = j;
                    } else{
                        // elemL and faceL have element info but not elemR and faceR
                        elemR[ind] = i;
                        faceR[ind] = j;
                        if (elemL[ind] > elemR[ind]){
                            std::swap(elemL[ind], elemR[ind]);
                            std::swap(faceL[ind], faceR[ind]);
                        }
                    }
                }
            }

            double s = (length[0]+length[1]+length[2])/2;
            double area = std::sqrt(s*(s-length[0])*(s-length[1])*(s-length[2]));
            // if (!isCurvedElement)
            Eigen::Matrix2d J;
            J.col(0) = node(ind2) - node(ind1);
            J.col(1) = node(ind3) - node(ind1);
            _elems.emplace_back(std::make_unique<LinearElement>(pointID, faceID, area, J));
        }
    }
    
    _B2E.conservativeResize(bcount, Eigen::NoChange);

    // Periodic boundaries, manually, only works on this one mesh
    std::deque<int> curve2Faces, curve4Faces, curve6Faces, curve8Faces;

    for (std::size_t i = 0; i < _faces.size(); i++){
        const Face& face = _faces[i];
        if (face._title == "Curve2") curve2Faces.push_back(i);
        else if (face._title == "Curve4") curve4Faces.push_front(i);
        else if (face._title == "Curve6") curve6Faces.push_back(i);
        else if (face._title == "Curve8") curve8Faces.push_front(i);
    }

    // Matching curve2 and curve4
    for (std::size_t i = 0; i < curve2Faces.size(); i++){
        int curve2ID = curve2Faces[i];
        int curve4ID = curve4Faces[i];

        if (elemL[curve2ID] < elemL[curve4ID]){
            // Left element is on row curve2ID
            elemR[curve2ID] = elemL[curve4ID];
            elemL[curve4ID] = elemL[curve2ID];
            elemR[curve4ID] = elemR[curve2ID];

            faceR[curve2ID] = faceL[curve4ID];
            faceL[curve4ID] = faceL[curve2ID];
            faceR[curve4ID] = faceR[curve2ID];
        } else{
            // Left element is on row curve4ID
            elemR[curve4ID] = elemL[curve2ID];
            elemL[curve2ID] = elemL[curve4ID];
            elemR[curve2ID] = elemR[curve4ID];

            faceR[curve4ID] = faceL[curve2ID];
            faceL[curve2ID] = faceL[curve4ID];
            faceR[curve2ID] = faceR[curve4ID]; 
        }
    }

    // Matching curve6 and curve8
    for (std::size_t i = 0; i < curve6Faces.size(); i++){
        int curve6ID = curve6Faces[i];
        int curve8ID = curve8Faces[i];

        if (elemL[curve6ID] < elemL[curve8ID]){
            // Left element is on row curve6ID
            elemR[curve6ID] = elemL[curve8ID];
            elemL[curve8ID] = elemL[curve6ID];
            elemR[curve8ID] = elemR[curve6ID];

            faceR[curve6ID] = faceL[curve8ID];
            faceL[curve8ID] = faceL[curve6ID];
            faceR[curve8ID] = faceR[curve6ID];
        } else{
            // Left element is on row curve8ID
            elemR[curve8ID] = elemL[curve6ID];
            elemL[curve6ID] = elemL[curve8ID];
            elemR[curve6ID] = elemR[curve8ID];

            faceR[curve8ID] = faceL[curve6ID];
            faceL[curve6ID] = faceL[curve8ID];
            faceR[curve6ID] = faceR[curve8ID];
        }
    }

    // Update I2E
    Eigen::Index icount = 0;
    _I2E.resize(_faces.size(), Eigen::NoChange);
    for (std::size_t i = 0; i < _faces.size(); i++){
        if (elemL[i] == -1) continue;
        _I2E(icount, 0) = elemL[i];
        _I2E(icount, 1) = faceL[i];
        _I2E(icount, 2) = elemR[i];
        _I2E(icount, 3) = faceR[i];        
    }
    _I2E.resize(icount, Eigen::NoChange);

    // Update normal vectors on each edge, always pointing from L to R
    _In.resize(icount, Eigen::NoChange);
    for (Eigen::Index i = 0; i < icount; i++){
        // Consturct a unit normal vector
        const Element& elem = *_elems[_I2E(i,0)];
        std::size_t localFaceID = _I2E(i,1);
        const Face& face = _faces[elem._faceID[localFaceID]];
        Eigen::Vector2d edge = node(face._pointID[1]) - node(face._pointID[0]);
        Eigen::Vector2d normal = Eigen::Vector2d{-edge[1], edge[0]};
        normal.normalize();

        // Find the remaining point on the "left" element and direct normal away from it
        const Eigen::Vector2d& p1 = node(elem._pointID[localFaceID]); // Point not on this edge
        const Eigen::Vector2d& p2 = node(elem._pointID[(localFaceID+1)%3]); // One of the points on this edge
        if ((p1-p2).dot(normal) > 0) normal *= -1;
        _In.row(i) = normal.transpose();
    }

    _Bn.resize(bcount, Eigen::NoChange);
    for (Eigen::Index i = 0; i < bcount; i++){
        // Consturct a unit normal vector
        const Element& elem = *_elems[_B2E(i,0)];
        std::size_t localFaceID = _B2E(i,1);
        const Face& face = _faces[elem._faceID[localFaceID]];
        Eigen::Vector2d edge = node(face._pointID[1]) - node(face._pointID[0]);
        Eigen::Vector2d normal = Eigen::Vector2d{-edge[1], edge[0]};
        normal.normalize();

        // Find the remaining point on the "left" element and direct normal away from it
        const Eigen::Vector2d& p1 = node(elem._pointID[localFaceID]); // Point not on this edge
        const Eigen::Vector2d& p2 = node(elem._pointID[(localFaceID+1)%3]); // One of the points on this edge
        if ((p1-p2).dot(normal) > 0) normal *= -1;
        _Bn.row(i) = normal.transpose();
    }
}

double TriangularMesh::length(std::size_t elemID, std::size_t localFaceID) const noexcept{
    auto& elem = *_elems[elemID];
    std::size_t globalFaceID = elem._pointID[localFaceID];
    return length(globalFaceID);
}

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