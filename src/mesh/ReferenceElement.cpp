#include "ReferenceElement.h"
#include "GaussLegendre2D.h"
#include "LagrangeBasisFunctions.h"

ReferenceElement::ReferenceElement(int p, int r):
    _p(p),
    _r(r)
{
    LagrangeBasisFunctions PhiLagrange(_p);
    _xiL = PhiLagrange.getLagrangeNodes();
    Eigen::Index Np = (_p+1)*(_p+2)/2;
    
    GaussLegendre2D<> GL(_r);
    _wQ = GL.getWeights();
    Eigen::Index Nq = _wQ.size();
    auto xiQ_alt = GL.getNodes();
    _xiQ.resize(Eigen::NoChange, Nq);
    _xiQ.row(0) = xiQ_alt(Eigen::seq(0, 2*Nq, 2));
    _xiQ.row(1) = xiQ_alt(Eigen::seq(1, 2*Nq, 2));

    _phi.resize(Np, Nq);
    _phiXi.resize(Np, Nq);
    _phiEta.resize(Np, Nq);
    _M.resize(Np, Np);

    // Computes monomial coefficients
    for (int i = 0; i < Nq; i++){
        double xi = _xiQ(0, i);
        double eta = _xiQ(1, i);
        _phi.col(i) = PhiLagrange.evalPhi(xi, eta);
        _phiXi.col(i) = PhiLagrange.evalPhiX(xi, eta);
        _phiEta.col(i) = PhiLagrange.evalPhiY(xi, eta);
    }

    // Computes reference mass matrix using quadrature
    for (int ii = 0; ii < Np; ii++) {
        for (int jj = ii; jj < Np; jj++) {
            Eigen::RowVectorXd phiij = _phi.row(ii).array() * _phi.row(jj).array();
            _M(ii,jj) = phiij*_wQ;
            _M(jj,ii) = _M(ii,jj);
        }
    }
}