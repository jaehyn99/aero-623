#include "mesh/TriangularMesh.h"
#include <algorithm>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>

TriangularMesh::Face::Face(const TriangularMesh& mesh, std::size_t pointID1, std::size_t pointID2, std::size_t nf, const std::string& title):
    _pointID{int(pointID1), int(pointID2)},
    _length((mesh._nodes[pointID2] - mesh._nodes[pointID1]).norm()),
    _nf(nf),
    _title(title)
{}

bool TriangularMesh::Face::operator==(const Face& other) const noexcept{
    return (_pointID[0] == other._pointID[0] && _pointID[1] == other._pointID[1]) || (_pointID[0] == other._pointID[1] && _pointID[1] == other._pointID[0]);
}

TriangularMesh::Element::Element(TriangularMesh& mesh, std::size_t pointID1, std::size_t pointID2, std::size_t pointID3, std::size_t ord, const std::string& basis):
    _pointID{int(pointID1), int(pointID2), int(pointID3)},
    _order(ord),
    _basis(basis),
    _centroid((mesh.node(pointID1) + mesh.node(pointID2) + mesh.node(pointID3))/3)
{
    // Look if its faces have already been added to mesh._faces
    for (std::size_t j = 0; j < 3; j++){
        Face iface(mesh, _pointID[(j+1)%3], _pointID[(j+2)%3]);

        // Assign face ID to elem and elem ID to face
        std::size_t faceID = std::find(mesh._faces.cbegin(), mesh._faces.cend(), iface) - mesh._faces.cbegin();
        if (faceID == mesh.numFaces()) mesh._faces.push_back(iface);
        _faceID[j] = faceID;

        Face& face = mesh.face(faceID);
        if (face._elemID[0] == -1) face._elemID[0] = mesh.numElems(); // No ElemID has been updated
        else face._elemID[1] = mesh.numElems();
    }

    // Calculates area
    double a = mesh.length(_faceID[0]);
    double b = mesh.length(_faceID[1]);
    double c = mesh.length(_faceID[2]);
    double s = (a+b+c)/2;
    _area = std::sqrt(s*(s-a)*(s-b)*(s-c));

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
            _faces.emplace_back(*this, ind1, ind2, nf, title);
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
            _elems.emplace_back(*this, ind1, ind2, ind3, ord, basis);
        }
    }

    splitNextLine();
    if (v.size() > 0){
        // Works for a GRI file that has periodic boundary info
        std::size_t nPG = std::stoi(v[0]);
        for (std::size_t i = 0; i < nPG; i++){
            splitNextLine();
            std::size_t nPGNode = std::stoi(v[0]);
            if (nPGNode >= 2){
                splitNextLine();
                std::size_t ind1 = std::stoi(v[0])-1;
                std::size_t ind2 = std::stoi(v[1])-1;
                Face left(*this, ind1, ind2);

                for (std::size_t j = 1; j < nPGNode; j++){
                    splitNextLine();
                    ind1 = std::stoi(v[0])-1;
                    ind2 = std::stoi(v[1])-1;
                    Face right(*this, ind1, ind2);

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
    } else{
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

            curve2Face._periodicFaceID = curve4ID;
            curve4Face._periodicFaceID = curve2ID;
            curve2Face._periodicElemID = curve4Face._elemID[0];
            curve4Face._periodicElemID = curve2Face._elemID[0];
        }

        // Matching curve6 and curve8
        for (std::size_t i = 0; i < curve6Faces.size(); i++){
            Face& curve6Face = curve6Faces[i].get();
            Face& curve8Face = curve8Faces[i].get();
            std::size_t curve6ID = std::find(_faces.cbegin(), _faces.cend(), curve6Face) - _faces.cbegin();
            std::size_t curve8ID = std::find(_faces.cbegin(), _faces.cend(), curve8Face) - _faces.cbegin();

            curve6Face._periodicFaceID = curve8ID;
            curve8Face._periodicFaceID = curve6ID;
            curve6Face._periodicElemID = curve8Face._elemID[0];
            curve8Face._periodicElemID = curve6Face._elemID[0];
        }
    }

    // Update normal vectors on each edge, always pointing from L to R
    for (Face& face: _faces){
        // Consturct a unit normal vector
        Eigen::Vector2d edge = node(face._pointID[1]) - node(face._pointID[0]);
        face._normal = Eigen::Vector2d{-edge[1], edge[0]};
        face._normal.normalize();

        // Find the remaining point on the "left" element and direct normal away from it
        int elemID = face._elemID[0];
        const Element& elem = _elems[elemID];
        std::size_t localFaceID = 0;
        for (std::size_t i = 0; i < 3; i++){
            const Face& iface = _faces[elem._faceID[i]];
            if (face == iface){
                localFaceID = i;
                break;
            }
        }
        const Eigen::Vector2d& p1 = node(elem._pointID[localFaceID]); // Point not on this edge
        const Eigen::Vector2d& p2 = node(elem._pointID[(localFaceID+1)%3]); // One of the points on this edge
        if ((p1-p2).dot(face._normal) > 0) face._normal *= -1;
    }

    fillMatrices();
}

Eigen::Vector2d TriangularMesh::normal(std::size_t elemID, std::size_t localFaceID) const noexcept{
    const Element& elem = _elems[elemID];
    int faceID = elem._faceID[localFaceID];
    const Face& face = _faces[faceID];

    Eigen::Vector2d n = face._normal;
    if (face._elemID[1] == int(elemID)) return -n; // on element R, revert the normal vector
    return n; // on element L or on a boundary, periodic or not
}

void TriangularMesh::fillMatrices() noexcept{
    
    // Stored all boundary names
    std::vector<std::string> boundaries;
    while (true){
        for (auto face: _faces){
            if (!face.isBoundaryFace()) break;
            else{
                auto it = std::find(boundaries.cbegin(), boundaries.cend(), face._title);
                if (it == boundaries.cend()) boundaries.push_back(face._title);
            }
        }
    }

    Eigen::Index icount, bcount;
    _I2E.resize(numFaces(), Eigen::NoChange);
    _B2E.resize(numFaces(), Eigen::NoChange);
    _In.resize(numFaces(), Eigen::NoChange);
    _Bn.resize(numFaces(), Eigen::NoChange);

    for (std::size_t i = 0; i < _faces.size(); i++){
        const Face& face = _faces[i];
        if (face.isBoundaryFace() && !face.isPeriodicFace()){ // boundary face
            std::size_t elemID = face._elemID[0]; // Elem number
            // Local face number
            const Element& elem = _elems[elemID];
            std::size_t localFaceID = std::find(elem._faceID.cbegin(), elem._faceID.cend(), i) - elem._faceID.cbegin();
            // Group number
            std::size_t bGroup = std::find(boundaries.cbegin(), boundaries.cend(), face._title) - boundaries.cbegin();
            _B2E(i,0) = elemID;
            _B2E(i,1) = localFaceID;
            _B2E(i,2) = bGroup;
            _Bn.row(i) = normal(elemID, localFaceID).transpose();
            bcount++;
        } else{ // interior or periodic face
            // Elem numbers
            std::size_t elemLID = face._elemID[0];
            std::size_t elemRID = face.isPeriodicFace() ? face._periodicElemID : face._elemID[1];
            if (elemRID > elemLID) std::swap(elemLID, elemRID);

            // Local face numbers
            const Element& elemL = _elems[elemLID];
            std::size_t faceL = std::find(elemL._faceID.cbegin(), elemL._faceID.cend(), i) - elemL._faceID.cbegin();
            const Element& elemR = _elems[elemRID];
            std::size_t faceR = std::find(elemR._faceID.cbegin(), elemR._faceID.cend(), i) - elemR._faceID.cbegin();

            _I2E(i,0) = elemLID;
            _I2E(i,1) = faceL;
            _I2E(i,2) = elemRID;
            _I2E(i,3) = faceR;
            _In.row(i) = normal(elemLID, faceL).transpose();
            icount++;
        }
    }

    _I2E.conservativeResize(icount, Eigen::NoChange);
    _B2E.conservativeResize(bcount, Eigen::NoChange);
    _In.conservativeResize(icount, Eigen::NoChange);
    _Bn.conservativeResize(bcount, Eigen::NoChange);    
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