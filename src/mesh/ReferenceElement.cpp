#include "ReferenceElement.h"
#include "GaussLegendre1D.h"
#include "GaussLegendre2D.h"
#include "Lagrange2DBasisFunctions.h"
#include <iostream>

ReferenceElement::ReferenceElement(int p, int r):
    _p(p),
    _r(r)
{
    Lagrange2DBasisFunctions PhiLagrange(_p);
    _xiL = PhiLagrange.nodes();
    Eigen::Index Np = (_p+1)*(_p+2)/2;
    
    // Computes basis functions and their derivative at internal quadrature points
    std::cout << "Computes basis functions and their derivative at internal quadrature points" << std::endl;
    GaussLegendre2D<> GL2(_r);
    _intW = GL2.getWeights();
    Eigen::Index Nq = _intW.size();
    auto intXiAlt = GL2.getNodes();
    _intXi.resize(Eigen::NoChange, Nq);
    for (int i = 0; i < Nq; i++){
        _intXi(0,i) = intXiAlt(2*i);
        _intXi(1,i) = intXiAlt(2*i+1);
    }

    _intPhi.resize(Np, Nq);
    _intPhiXi.resize(Np, Nq);
    _intPhiEta.resize(Np, Nq);

    for (int i = 0; i < Nq; i++){
        double xi = _intXi(0, i);
        double eta = _intXi(1, i);
        _intPhi.col(i) = PhiLagrange.evalPhi(xi, eta);
        _intPhiXi.col(i) = PhiLagrange.evalPhiX(xi, eta);
        _intPhiEta.col(i) = PhiLagrange.evalPhiY(xi, eta);
    }

    // Computes basis functions and their derivative at edge quadrature points
    std::cout << "Computes basis functions and their derivative at edge quadrature points" << std::endl;
    GaussLegendre1D<> GL1(_r);
    _edgeW = GL1.getWeights(); // weights from 0 to 1;
    Nq = _edgeW.size();
    auto edgeXiAlt = GL1.getNodes();
    
    for (int edge = 0; edge < 3; edge++){
        _edgeXi[edge].resize(Eigen::NoChange, Nq);
        _edgePhi[edge].resize(Np, Nq);
        _edgePhiXi[edge].resize(Np, Nq);
        _edgePhiEta[edge].resize(Np, Nq);
    }

    Eigen::Vector2d p1{0,0}, p2{1,0}, p3{0,1};
    for (int i = 0; i < Nq; i++){
        _edgeXi[0].col(i) = p2 + (p3-p2)*edgeXiAlt[i]; // edge 0 = hypotnuse
        _edgeXi[1].col(i) = p3 + (p1-p3)*edgeXiAlt[i]; // edge 1 = vertical edge
        _edgeXi[2].col(i) = p1 + (p2-p1)*edgeXiAlt[i]; // edge 2 = horizontal edge
    }

    for (int edge = 0; edge < 3; edge++){
        for (int i = 0; i < Nq; i++){
            double xi = _edgeXi[edge](0,i);
            double eta = _edgeXi[edge](1,i);
            _edgePhi[edge].col(i) = PhiLagrange.evalPhi(xi, eta);
            _edgePhiXi[edge].col(i) = PhiLagrange.evalPhiX(xi, eta);
            _edgePhiEta[edge].col(i) = PhiLagrange.evalPhiY(xi, eta);
        }
    }

    // Computes reference mass matrix using quadrature
    std::cout << "Computes reference mass matrix using quadrature" << std::endl;
    Eigen::MatrixXd M(Np, Np);
    for (int ii = 0; ii < Np; ii++) {
        for (int jj = ii; jj < Np; jj++) {
            Eigen::RowVectorXd phiij = _intPhi.row(ii).array() * _intPhi.row(jj).array();
            M(ii,jj) = phiij*_intW;
            M(jj,ii) = M(ii,jj);
        }
    }
    _MLLT = M.llt();
    if (_MLLT.info() != Eigen::Success) throw std::runtime_error("ERROR: Unable to factorize the mass matrix.");
}