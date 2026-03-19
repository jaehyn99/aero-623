#include "TriangularMesh.h"
#include "CubicSpline.h"
#include "CurvedElement.h"
#include "CurvedFace.h"
#include "GaussLegendre1D.h"
#include "GaussLegendre2D.h"
#include "Lagrange1DBasisFunctions.h"
#include "Lagrange2DBasisFunctions.h"
#include "LinearElement.h"
#include "LinearFace.h"

#include <algorithm>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>

TriangularMesh::TriangularMesh(TriangularMesh&&) = default;
TriangularMesh::~TriangularMesh() = default;
TriangularMesh& TriangularMesh::operator=(TriangularMesh&&) = default;

TriangularMesh::TriangularMesh(const std::string& fileName, std::size_t p, std::size_t q, std::size_t r):
    _ref(p, r)
{
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
    std::vector<double> Xl, Yl, Xu, Yu;
    std::size_t nBGroup = std::stoi(v[0]); // number of boundary groups
    for (std::size_t i = 0; i < nBGroup; i++){
        splitNextLine();
        std::size_t nBFace = std::stoi(v[0]);
        std::string title = v[2];
        for (std::size_t j = 0; j < nBFace; j++){
            splitNextLine();
            std::size_t ind1 = std::stoi(v[0]) - 1;
            std::size_t ind2 = std::stoi(v[1]) - 1;
            Eigen::Vector2i pointID{ind1, ind2};
            double length = (_nodes[ind1] - _nodes[ind2]).lpNorm<2>();
            if (q == 1 || (title != "Curve1" && title != "Curve5"))
                _faces.emplace_back(std::make_unique<LinearFace>(pointID, length, title));
            else{
                _faces.emplace_back(std::make_unique<CurvedFace>(pointID, length, q, title));
                if (title == "Curve1"){
                    if (Xu.empty()){ // no points have been added into curve1
                        Xu.push_back(_nodes[ind1].x());
                        Yu.push_back(_nodes[ind1].y());
                    }
                    Xu.push_back(_nodes[ind2].x());
                    Yu.push_back(_nodes[ind2].y());                   
                }
                else{
                    if (Xl.empty()){ // no points have been added into curve5
                        Xl.push_back(_nodes[ind1].x());
                        Yl.push_back(_nodes[ind1].y());
                    }
                    Xl.push_back(_nodes[ind2].x());
                    Yl.push_back(_nodes[ind2].y());                     
                }
            }
        }
    }

    // Treatment of curve edges - do a cubic spline interpolation
    if (q > 1){
        Eigen::VectorXd Xl2 = Eigen::Map<Eigen::VectorXd>(Xl.data(), Xl.size());
        Eigen::VectorXd Yl2 = Eigen::Map<Eigen::VectorXd>(Yl.data(), Yl.size());
        Eigen::VectorXd Xu2 = Eigen::Map<Eigen::VectorXd>(Xu.data(), Xu.size());
        Eigen::VectorXd Yu2 = Eigen::Map<Eigen::VectorXd>(Yu.data(), Yu.size());
        _lower = std::make_unique<CubicSpline>(Xl2, Yl2);
        _upper = std::make_unique<CubicSpline>(Xu2, Yu2);
        GaussLegendre1D<1> GL(r);

        int iL = 0, iU = 0;
        auto& Sl = _lower->S();
        auto& Su = _upper->S();
        for (const auto& face: _faces){
            if (face->title() != "Curve1" && face->title() != "Curve5") continue;
            CurvedFace* cFace = dynamic_cast<CurvedFace*>(face.get());
            double s0 = face->title() == "Curve1" ? Su[iU] : Sl[iL];
            double s1 = face->title() == "Curve1" ? Su[iU+1] : Sl[iL+1];
            cFace->_length = s1-s0; // assign length

            // Generate the Lagrange nodes for geometry appoximation
            auto lnodes = Eigen::VectorXd::LinSpaced(q+1, s0, s1);
            cFace->_xL.resize(Eigen::NoChange, lnodes.size()-2); // ignoring the boundary nodes
            for (int i = 1; i < lnodes.size()-1; i++){
                double s = lnodes[i];
                cFace->_xL.col(i-1) = face->title() == "Curve1" ? _upper->eval(s) : _lower->eval(s);
            }

            // Use the Lagrange basis to approximate/interpolate the edge
            Lagrange1DBasisFunctions PhiLagrange(q);
            Eigen::MatrixXd Xedge(2, q+1);
            Xedge.col(0) = _nodes[cFace->pointID()[0]]; // first endpoint
            Xedge.middleCols(1, q-1) = cFace->_xL; // edge points
            Xedge.col(q) = _nodes[cFace->pointID()[1]]; // second endpoint

            // Generate the quadrature nodes
            auto qnodes = GL.getNodes(); // quadrature nodes on the domain (0, 1)
            int NQ = qnodes.size(); // number of quadrature nodes
            // Resize these matrices so that the relevant info at each quadrature node can be stored
            cFace->_xq.resize(Eigen::NoChange, NQ);
            cFace->_n.resize(Eigen::NoChange, NQ);
            cFace->_detJ.resize(NQ);
            for (int i = 0; i < NQ; i++){
                double sigma = qnodes[i];
                cFace->_xq.col(i) = PhiLagrange.funcEval(sigma, Xedge);

                Eigen::Vector2d tds_dsigma = PhiLagrange.funcXEval(sigma, Xedge);
                Eigen::Vector2d N{tds_dsigma[1], -tds_dsigma[0]};
                N.normalize();
                if (face->title() == "Curve1" && N.y() > 0) N *= -1; // normals on the "upper" curve point down
                if (face->title() == "Curve5" && N.y() < 0) N *= -1; // normals on the "lower" curve point up              
                
                cFace->_n.col(i) = N;
                cFace->_detJ[i] = tds_dsigma.lpNorm<2>();
            }

            if (face->title() == "Curve1") iU++;
            else iL++;
        }
    }

    // Fill in the elements, interior nodes, and connectivity info
    int Np = (p+1)*(p+2)/2;
    _elems.reserve(nElemTot);
    while (nElemTot > 0){
        splitNextLine();
        std::size_t nElem = std::stoi(v[0]);
        // std::size_t ord = std::stoi(v[1]);
        // std::string basis = v[2];
        nElemTot -= nElem;
        for (int i = 0; i < int(nElem); i++){
            splitNextLine();
            std::size_t ind1 = std::stoi(v[0]) - 1;
            std::size_t ind2 = std::stoi(v[1]) - 1;
            std::size_t ind3 = std::stoi(v[2]) - 1;
            Eigen::Vector3i pointID{int(ind1), int(ind2), int(ind3)};
            Eigen::Vector3i faceID;
            Eigen::Vector3d length;
            //Eigen::Vector2d centroid{ (_nodes[ind1].x() + _nodes[ind2].x() + _nodes[ind3].x())/3, (_nodes[ind1].y() + _nodes[ind2].y() + _nodes[ind3].y())/3 };

            bool isCurvedElement = false;
            for (int j = 0; j < 3; j++){
                // Look if its faces have already been added to the list of faces
                Eigen::Vector2i facePointID{pointID[(j+1)%3], pointID[(j+2)%3]};
                std::unique_ptr<Face> iface = std::make_unique<LinearFace>(facePointID, 0.0);
                auto it = std::find_if(_faces.cbegin(), _faces.cend(), [&iface](const auto& other){ return *iface == *other; });
                std::size_t ind = it - _faces.cbegin();
                faceID[j] = ind;

                // If not found, create a new face and add it to the list
                if (it == _faces.cend()){
                    iface->_length = (_nodes[pointID[(j+1)%3]] - _nodes[pointID[(j+2)%3]]).lpNorm<2>();
                    length[j] = iface->_length;
                    iface->_elemID[0] = i;
                    _faces.push_back(std::move(iface));
                } else{
                    length[j] = _faces[ind]->length();
                    if (_faces[ind]->_elemID[0] == -1) _faces[ind]->_elemID[0] = i;
                    else{
                        _faces[ind]->_elemID[1] = i;
                        if (_faces[ind]->_elemID[0] > i) std::swap(_faces[ind]->_elemID[0], _faces[ind]->_elemID[1]);
                    }
                }
                if (_faces[ind]->isCurvedFace()) isCurvedElement = true;
            }

            if (!isCurvedElement){
                double s = (length[0]+length[1]+length[2])/2;
                double area = std::sqrt(s*(s-length[0])*(s-length[1])*(s-length[2]));
                Eigen::Matrix2d J;
                J.col(0) = node(ind2) - node(ind1);
                J.col(1) = node(ind3) - node(ind1);
                _elems.emplace_back(std::make_unique<LinearElement>(pointID, faceID, area, J));
            } else{
                // Dealing with a curved element, compute the internal node and Jacobians
                auto cElem = std::make_unique<CurvedElement>(pointID, faceID, 0.0);
                // First, create a matrix of the 10 Lagrange nodes (q=3)
                Eigen::Matrix<double, 2, 10> Phi; // Lagrange nodes
                Phi.col(0) = _nodes[cElem->_pointID[0]];
                Phi.col(3) = _nodes[cElem->_pointID[1]];
                Phi.col(9) = _nodes[cElem->_pointID[2]];

                // Computes the internal node as the centroid of the surrounding 6 nodes
                for (int jj = 0; jj < 3; jj++){
                    const auto& face = _faces[cElem->_faceID[jj]];
                    if (face->isCurvedFace()){
                        // Gets the two internal nodes on a curved face
                        CurvedFace* cFace = dynamic_cast<CurvedFace*>(face.get());
                        cElem->_internal += cFace->_xL.col(0);
                        cElem->_internal += cFace->_xL.col(1);

                        int ind1, ind2, ind3, ind4; // indices of the nodes along this edge
                        if (jj == 0){
                            ind1 = 3;
                            ind2 = 6;
                            ind3 = 8;
                            ind4 = 9;
                        } else if (jj == 1){
                            ind1 = 0;
                            ind2 = 4;
                            ind3 = 7;
                            ind4 = 9;
                        } else{
                            ind1 = 0;
                            ind2 = 1;
                            ind3 = 2;
                            ind4 = 4;
                        }

                        if ((Phi(0,ind4)-Phi(0,ind1)) * (cFace->_xL(0,1)-cFace->_xL(0,0)) < 0){
                            std::swap(cFace->_pointID(0), cFace->_pointID(1));
                            cFace->_xL.col(0).swap(cFace->_xL.col(1));
                        }
                        Phi.col(ind2) = cFace->_xL.col(0);
                        Phi.col(ind3) = cFace->_xL.col(1);
                    } else{
                        // Gets the two internal nodes on a linear face
                        int ind1 = (jj+1)%3;
                        int ind2 = (jj+2)%3;
                        Eigen::Vector2d node1 = _nodes[cElem->_pointID[std::min(ind1, ind2)]];
                        Eigen::Vector2d node2 = _nodes[cElem->_pointID[std::max(ind1, ind2)]];
                        Eigen::Vector2d diff = node2-node1;
                        Eigen::Vector2d Lnode1 = node1 + diff/3;
                        Eigen::Vector2d Lnode2 = Lnode1 + diff/3;
                        cElem->_internal += Lnode1;
                        cElem->_internal += Lnode2;
                        
                        if (jj == 0){
                            Phi.col(6) = Lnode1;
                            Phi.col(8) = Lnode2;
                        } else if (jj == 1){
                            Phi.col(4) = Lnode1;
                            Phi.col(7) = Lnode2;
                        } else{
                            Phi.col(1) = Lnode1;
                            Phi.col(2) = Lnode2;
                        }
                    }
                }
                cElem->_internal /= 6;
                Phi.col(5) = cElem->_internal;

                // Compute the Jacobian at the quadrature nodes
                GaussLegendre2D<1> GL2(r);
                auto qnodes = GL2.getNodes();
                auto qweights = GL2.getWeights();
                int NQ = qweights.size();
                cElem->_area = 0.0;
                cElem->_J.resize(NQ, Eigen::Matrix2d::Zero());
                cElem->_detJ.resize(NQ);

                Lagrange2DBasisFunctions phis(q);
                for (int ii = 0; ii < NQ; ii++){
                    int Nq = (q+1)*(q+2)/2;
                    Eigen::MatrixX2d dphi(Nq, 2);
                    dphi.col(0) = phis.evalPhiX(qnodes[2*ii], qnodes[2*ii+1]);
                    dphi.col(1) = phis.evalPhiY(qnodes[2*ii], qnodes[2*ii+1]);
                    for (int jj = 0; jj < Nq; jj++) cElem->_J[ii] += Phi.col(jj)*dphi.row(jj);
                    cElem->_detJ[ii] = std::abs(cElem->_J[ii].determinant());
                    cElem->_area += cElem->_detJ[ii]*qweights[ii];
                }

                Eigen::MatrixXd M(Np, Np);
                const auto& intPhi = _ref.intPhi();
                Eigen::VectorXd intW = _ref.intW().array() * cElem->_detJ.array();
                for (int i = 0; i < Np; i++) {
                    for (int j = i; j < Np; j++) {
                        Eigen::RowVectorXd phiij = intPhi.row(i).array() * intPhi.row(j).array();
                        M(i,j) = phiij*intW;
                        M(j,i) = M(i,j);
                    }
                }
                cElem->_MLLT = M.llt();
                _elems.emplace_back(std::move(cElem));
            }
        }
    }

    // Periodic boundaries, manually, only works on this one mesh
    std::deque<int> curve2Faces, curve4Faces, curve6Faces, curve8Faces;
    for (std::size_t i = 0; i < _faces.size(); i++){
        const Face& face = *_faces[i];
        if (face._title == "Curve2") curve2Faces.push_back(i);
        else if (face._title == "Curve4") curve4Faces.push_front(i);
        else if (face._title == "Curve6") curve6Faces.push_back(i);
        else if (face._title == "Curve8") curve8Faces.push_front(i);
    }

    // Matching curve2 and curve4
    for (std::size_t i = 0; i < curve2Faces.size(); i++){
        int curve2ID = curve2Faces[i];
        int curve4ID = curve4Faces[i];
        int elemL = _faces[curve2ID]->_elemID[0];
        int elemR = _faces[curve4ID]->_elemID[0];
        if (elemL > elemR) std::swap(elemL, elemR);
        _faces[curve2ID]->_elemID[0] = elemL;
        _faces[curve2ID]->_elemID[1] = elemR;
        _faces[curve4ID]->_elemID[0] = elemL;
        _faces[curve4ID]->_elemID[1] = elemR;        
    }

    // Matching curve6 and curve8
    for (std::size_t i = 0; i < curve6Faces.size(); i++){
        int curve6ID = curve6Faces[i];
        int curve8ID = curve8Faces[i];
        int elemL = _faces[curve6ID]->_elemID[0];
        int elemR = _faces[curve8ID]->_elemID[0];
        if (elemL > elemR) std::swap(elemL, elemR);
        _faces[curve6ID]->_elemID[0] = elemL;
        _faces[curve6ID]->_elemID[1] = elemR;
        _faces[curve8ID]->_elemID[0] = elemL;
        _faces[curve8ID]->_elemID[1] = elemR;   
    }

    // Update normal vectors on each edge, always pointing from L to R
    for (int i = 0; i < int(numFaces()); i++){
        // Construct a unit normal vector
        if (LinearFace* lFace = dynamic_cast<LinearFace*>(_faces[i].get())){
            const Element& elem = *_elems[lFace->_elemID[0]];

            std::size_t localFaceID;
            if (elem._pointID[0] == i) localFaceID = 0;
            else if (elem._pointID[1] == i) localFaceID = 1;
            else localFaceID = 2;

            Eigen::Vector2d edge = node(lFace->_pointID[1]) - node(lFace->_pointID[0]);
            Eigen::Vector2d normal = Eigen::Vector2d{-edge[1], edge[0]};
            normal.normalize();

            // Find the remaining point on the "left" element and direct normal away from it
            const Eigen::Vector2d& p1 = node(elem._pointID[localFaceID]); // Point not on this edge
            const Eigen::Vector2d& p2 = node(elem._pointID[(localFaceID+1)%3]); // One of the points on this edge
            if ((p1-p2).dot(normal) > 0) normal *= -1;
            lFace->_n = normal;
        }
    }
}

double TriangularMesh::length(std::size_t faceID) const noexcept{ return _faces[faceID]->_length; }

double TriangularMesh::length(std::size_t elemID, std::size_t localFaceID) const noexcept{
    auto& elem = *_elems[elemID];
    std::size_t globalFaceID = elem._pointID[localFaceID];
    return length(globalFaceID);
}

double TriangularMesh::area(std::size_t elemID) const noexcept{ return _elems[elemID]->_area; }

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