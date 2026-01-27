#include "TriangularMesh.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>

bool TriangularMesh::Face::operator==(const Face& other) const noexcept{
    return (_pointID[0] == other._pointID[0] && _pointID[1] == other._pointID[1]) || (_pointID[0] == other._pointID[1] && _pointID[1] == other._pointID[0]);
}

TriangularMesh::TriangularMesh(const std::string& fileName){
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
        std::size_t nf = std::stoi(v[1]);
        std::string title = v[2];
        for (std::size_t j = 0; j < nBFace; j++){
            splitNextLine();
            std::size_t ind1 = std::stoi(v[0]) - 1;
            std::size_t ind2 = std::stoi(v[1]) - 1;
            _faces.emplace_back(ind1, ind2, nf, title);
        }
    }

    // Fill in the elements and interior nodes
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
            _elems.emplace_back(ind1, ind2, ind3, ord, basis);
            Element& elem = _elems.back();

            // Look if its faces have already been added to faces
            for (std::size_t j = 0; j < 3; j++){
                std::size_t fInd1, fInd2;
                if (j == 0){ fInd1 = ind2; fInd2 = ind3; }
                else if (j == 1){ fInd1 = ind3; fInd2 = ind1; }
                else{ fInd1 = ind1; fInd2 = ind2; } 
                Face iface(fInd1, fInd2);

                // Assign face ID to elem and elem ID to face
                std::size_t faceID = std::find(_faces.begin(), _faces.end(), iface) - _faces.begin();
                if (faceID == _faces.size()) _faces.push_back(iface);
                elem._faceID[j] = faceID;
                if (_faces[faceID]._elemID[0] == -1) _faces[faceID]._elemID[0] = _elems.size()-1;
                else _faces[faceID]._elemID[1] = _elems.size()-1;
            }
        }
    }

    // Periodic boundaries
    splitNextLine();
    std::size_t nPG = std::stoi(v[0]);
    for (std::size_t i = 0; i < nPG; i++){
        splitNextLine();
        std::size_t nPGNode = std::stoi(v[0]);
        if (nPGNode >= 2){
            splitNextLine();
            std::size_t ind1 = std::stoi(v[0])-1;
            std::size_t ind2 = std::stoi(v[1])-1;
            Face left(ind1, ind2);

            for (std::size_t j = 1; j < nPGNode; j++){
                splitNextLine();
                ind1 = std::stoi(v[0])-1;
                ind2 = std::stoi(v[1])-1;
                Face right(ind1, ind2);

                // Locate the periodic faces
                auto it1 = std::find(_faces.begin(), _faces.end(), left);
                auto it2 = std::find(_faces.begin(), _faces.end(), right);
                it1->_periodicFaceID = it2 - _faces.cbegin();
                it1->_periodicElemID = it2->_elemID[0];

                it2->_periodicFaceID = it1 - _faces.cbegin();
                it2->_periodicElemID = it1->_elemID[0];

                left = right;
            }
        }
    }
}

Eigen::Vector2d TriangularMesh::vect(std::size_t faceID) const noexcept{
    const Face& face = _faces[faceID];
    const Eigen::Vector2d& node1 = _nodes[face._pointID[0]];
    const Eigen::Vector2d& node2 = _nodes[face._pointID[1]];
    return node2 - node1;
}

double TriangularMesh::length(std::size_t faceID) const noexcept{ return vect(faceID).norm(); }

double TriangularMesh::area(std::size_t elemID) const noexcept{
    const Element& elem = _elems[elemID];
    double a = length(elem._faceID[0]);
    double b = length(elem._faceID[1]);
    double c = length(elem._faceID[2]);
    double s = (a+b+c)/2;
    return std::sqrt(s*(s-a)*(s-b)*(s-c));
}

Eigen::Vector2d TriangularMesh::normal(std::size_t elemID, std::size_t localFaceID) const noexcept{
    // Outward normal vector to edge #localEdgeID (0, 1, 2) of element #elemID
    const Element& elem = _elems[elemID];
    Eigen::Vector2d edge = vect(elem._faceID[localFaceID]);
    Eigen::Vector2d normal{-edge[1], edge[0]};
    normal.normalize();

    const Eigen::Vector2d& p1 = _nodes[elem._pointID[(localFaceID+1)%3]]; // One of the points on this edge
    const Eigen::Vector2d& p2 = _nodes[elem._pointID[localFaceID]]; // Point not on this edge
    if ((p2-p1).dot(normal) < 0) return normal;
    return -normal;
}

void TriangularMesh::writeGri(const std::string& fileName) const noexcept{
    std::ofstream of;
    of << std::fixed << std::setprecision(15);
    std::size_t ind = fileName.find(".gri");
    std::string fileBase = (ind == std::string::npos) ? fileName : fileName.substr(0, ind);
    
    // Stored all boundary names
    std::vector<std::string> boundaries;
    bool stop = false;
    while (!stop){
        for (auto face: _faces){
            if (!face.isBoundaryFace()){
                stop = true;
                break;
            }
            else{
                auto it = std::find(boundaries.cbegin(), boundaries.cend(), face._title);
                if (it == boundaries.cend()) boundaries.push_back(face._title);
            }
        }
    }

    // I2E
    of.open(fileBase + "I2E.txt");
    for (std::size_t i = 0; i < _faces.size(); i++){
        const Face& face = _faces[i];
        if (face.isBoundaryFace() && face._periodicFaceID == -1) continue;

        if (!face.isBoundaryFace()){ // interior face
            // Elem numbers
            std::size_t elemLID = face._elemID[0];
            std::size_t elemRID = face._elemID[1];

            // Local face numbers
            const Element& elemL = _elems[elemLID];
            std::size_t faceL = std::find(elemL._faceID.cbegin(), elemL._faceID.cend(), i) - elemL._faceID.cbegin();
            const Element& elemR = _elems[elemRID];
            std::size_t faceR = std::find(elemR._faceID.cbegin(), elemR._faceID.cend(), i) - elemR._faceID.cbegin();

            of << elemLID+1 << " " << faceL+1 << " " << elemRID+1 << " " << faceR+1 << "\n";
        } else{
            // Elem numbers
            std::size_t elemID = face._elemID[0];
            std::size_t pElemID = face._periodicElemID;

            // Local face numbers
            const Element& elem = _elems[elemID];
            std::size_t faceL = std::find(elem._faceID.cbegin(), elem._faceID.cend(), i) - elem._faceID.cbegin();
            const Element& pElem = _elems[pElemID]; // does not contain this face but its periodic counterpart
            std::size_t faceR = std::find(pElem._faceID.cbegin(), pElem._faceID.cend(), face._periodicFaceID) - pElem._faceID.cbegin();

            if (elemID < pElemID) of << elemID+1 << " " << faceL+1 << " " << pElemID+1 << " " << faceR+1 << "\n";
            else of << pElemID+1 << " " << faceR+1 << " " << elemID+1 << " " << faceL+1 << "\n";
        }
    }
    of.close();

    // Periodic edges
    of.open(fileBase + "periodicEdges.txt");
    std::vector<int> periodicEdgeID;
    periodicEdgeID.reserve(_faces.size());
    stop = false;
    int i = 1;
    while (!stop){
        for (auto face: _faces){
            if (face.isBoundaryFace()){
                if (face._periodicFaceID == -1) periodicEdgeID.push_back(0);
                else{
                    periodicEdgeID.push_back(i);
                    i++;
                }
            } else{
                stop = true;
                break;
            }
        }
    }
    for (std::size_t i = 0; i < periodicEdgeID.size(); i++){
        if (periodicEdgeID[i] == 0) continue;
        const Face& face = _faces[i];
        of << periodicEdgeID[face._periodicFaceID] << "\n";
    }
    of.close();

    // Connectivity matrix B2E
    of.open(fileBase + "B2E.txt");
    for (auto face: _faces){
        if (!face.isBoundaryFace() || face._periodicFaceID != -1) continue;

        std::size_t elemID = face._elemID[0]; // Elem number
        // Local face number
        const Element& elem = _elems[elemID];
        std::size_t localFaceID = std::find(elem._faceID.cbegin(), elem._faceID.cend(), elemID) - elem._faceID.cbegin();
        // Group number
        std::size_t bGroup = std::find(boundaries.cbegin(), boundaries.cend(), face._title) - boundaries.cbegin();
        of << elemID+1 << " " << localFaceID+1 << " " << bGroup+1 << "\n";
    }
    of.close();

    // Normal vectors In
    of.open(fileBase + "In.txt");
    for (std::size_t i = 0; i < _faces.size(); i++){
        const Face& face = _faces[i];
        if (face.isBoundaryFace() && face._periodicFaceID == -1) continue;
        std::size_t elemID = face._elemID[0];
        const Element& elem = _elems[elemID];
        std::size_t localFaceID = std::find(elem._faceID.cbegin(), elem._faceID.cend(), i) - elem._faceID.cbegin();
        of << normal(elemID, localFaceID).transpose() << "\n";
    }
    of.close();

    // Normal vectors Bn
    of.open(fileBase + "Bn.txt");
    for (std::size_t i = 0; i < _faces.size(); i++){
        const Face& face = _faces[i];
        if (!face.isBoundaryFace() || face._periodicFaceID != -1) continue;
        std::size_t elemID = face._elemID[0];
        const Element& elem = _elems[elemID];
        std::size_t localFaceID = std::find(elem._faceID.cbegin(), elem._faceID.cend(), i) - elem._faceID.cbegin();
        of << normal(elemID, localFaceID).transpose() << "\n";
    }
    of.close();

    // Areas
    of.open(fileBase + "Area.txt");
    for (std::size_t i = 0; i < _elems.size(); i++) of << area(i) << "\n";
    of.close();
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