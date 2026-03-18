#include "DGIntegration.h"
#include "ShapeFunctions.hpp"
#include "HLLEFlux.h"
#include "RoeFlux.h"

#include <iostream>
#include <stdexcept>

// questions: where is gamma stored, are the basis and gradient already stored?, 
// should flux be replaced with Roe or HLLE flux?

static Eigen::VectorXd evalBasis(const Eigen::MatrixXd& Phi,
                                  double xi1, double xi2,
                                  int p)
{
    int Np = Phi.rows();
    Eigen::VectorXd m(Np);
    int idx = 0;
    for (int deg = 0; deg <= p; ++deg)
        for (int b = 0; b <= deg; ++b)
            m(idx++) = std::pow(xi1, deg-b) * std::pow(xi2, b);
 
    return Phi * m;   // (Np x Np) * Np = Np
}


static void evalBasisGrad(const Eigen::MatrixXd& PhiXi,
                           const Eigen::MatrixXd& PhiEta,
                           double xi1, double xi2, int p,
                           Eigen::VectorXd& dPhiXi,
                           Eigen::VectorXd& dPhiEta)
{
    int Np   = PhiXi.cols(); 
    int NpDv = PhiXi.rows(); 
 
    Eigen::VectorXd m(Np);
    int idx = 0;
    for (int deg = 0; deg <= p; ++deg)
        for (int b = 0; b <= deg; ++b)
            m(idx++) = std::pow(xi1, deg-b) * std::pow(xi2, b);
 
    dPhiXi  = PhiXi * m;   // (NpDv x Np) * Np
    dPhiEta = PhiEta * m;
}

static Eigen::Vector4d reconstructState(const Eigen::MatrixXd& modesElem,
                                         const Eigen::VectorXd& phi)
{
    return modesElem * phi;
}

static Eigen::Matrix<double,4,2> physicalFlux(const Eigen::Vector4d& U, double gamma)
{
    double gm1 = gamma - 1.0;
    double rho = U(0);
    double u = U(1)/rho;
    double v = U(2)/rho;
    double rhoE = U(3);
    double p = gm1*(rhoE - 0.5*rho*(u*u + v*v));
 
    Eigen::Matrix<double,4,2> F;

    F(0,0) = rho*u;
    F(1,0) = rho*u*u + p;
    F(2,0) = rho*u*v;
    F(3,0) = (rhoE + p)*u;

    F(0,1) = rho*v;
    F(1,1) = rho*u*v;
    F(2,1) = rho*v*v + p;
    F(3,1) = (rhoE + p)*v;
    return F;
}

// TAKE THIS OUT AND CALL GL-1D directly?
static void gaussLegendre1D(int nQE,
                             Eigen::VectorXd& pts,
                             Eigen::VectorXd& wts)
{
    pts.resize(nQE);
    wts.resize(nQE);
    switch (nQE) {
        case 1:
            pts << 0.5;
            wts << 1.0;
            break;
        case 2:
            pts << 0.5 - 0.5/std::sqrt(3.0),
                   0.5 + 0.5/std::sqrt(3.0);
            wts << 0.5, 0.5;
            break;
        case 3:
            pts << 0.5 - 0.5*std::sqrt(0.6),
                   0.5,
                   0.5 + 0.5*std::sqrt(0.6);
            wts << 5.0/18.0, 8.0/18.0, 5.0/18.0;
            break;
        default:
            throw std::runtime_error("DGIntegration: unsupported nQE (add more Gauss points)");
    }
}

static Eigen::Vector2d edgeRefPoint(int localEdge, double t)
{ // 0 < t < 1
    switch (localEdge) {
        case 0: return {t,0.0};
        case 1: return {1.0-t, t};
        case 2: return {0.0, t};
        default: throw std::runtime_error("DGIntegration: invalid local edge index");
    }
}

// take out and call stored values?
static double elemJacobian(const Eigen::Vector2d& v0,
                            const Eigen::Vector2d& v1,
                            const Eigen::Vector2d& v2)
{
    return (v1(0)-v0(0))*(v2(1)-v0(1)) - (v1(1)-v0(1))*(v2(0)-v0(0));
}

// take out and call stored values?
static Eigen::Matrix2d elemJacInvT(const Eigen::Vector2d& v0,
                                    const Eigen::Vector2d& v1,
                                    const Eigen::Vector2d& v2)
{
    double dx10 = v1(0)-v0(0), dy10 = v1(1)-v0(1);
    double dx20 = v2(0)-v0(0), dy20 = v2(1)-v0(1);
    double det  = dx10*dy20 - dy10*dx20;
    Eigen::Matrix2d Jinv;
    Jinv << dy20/det, -dy10/det,
            -dx20/det,  dx10/det;
    return Jinv.transpose();
}

void DGIntegration::getElemIntegral(int p_,
                                     double gamma,
                                     const TriangularMesh& mesh,
                                     const Eigen::MatrixXd& modes,
                                     Eigen::MatrixXd& residual)
{
    ShapeFunctions sf(p_);
    femQuadrature fq;
 
    int Np = (p_+1)*(p_+2)/2;
    int nEq = 4;                
    int qOrd = 2*p_;         
    if (qOrd < 1) qOrd = 1;
 
    Eigen::MatrixXd Phi = sf.getShapeFuncCoeffs   (p_);
    Eigen::MatrixXd PhiXi = sf.getShapeFuncXiCoeffs (p_); 
    Eigen::MatrixXd PhiEta = sf.getShapeFuncEtaCoeffs(p_); 
 
    Eigen::MatrixXd xiQ  = fq.getQuadXi(qOrd);    
    Eigen::VectorXd wQ   = fq.getQuadW (qOrd); 
    int nQV = (int)wQ.size();
 
    int nElems = (int)mesh.numElems();
 
    // Loop over elements
    for (int e = 0; e < nElems; ++e)
    {
        const auto& elem = mesh.elem(e);
        Eigen::Vector2d v0 = mesh.node(elem._pointID[0]);
        Eigen::Vector2d v1 = mesh.node(elem._pointID[1]);
        Eigen::Vector2d v2 = mesh.node(elem._pointID[2]);
 
        double detJ  = elemJacobian(v0, v1, v2);
        Eigen::Matrix2d JinvT = elemJacInvT(v0, v1, v2);
 
        Eigen::MatrixXd modesE = modes.block(0, e*Np, nEq, Np); // IS THIS THE CORRECT STORAGE METHOD

        Eigen::MatrixXd Re = Eigen::MatrixXd::Zero(nEq, Np);
        // loop over quad points
        for (int q = 0; q < nQV; ++q)
        {
            double xi1 = xiQ(q, 0);
            double xi2 = xiQ(q, 1); 
            Eigen::VectorXd phi     = evalBasis    (Phi,    xi1, xi2, p_);
            Eigen::VectorXd dPhiXi, dPhiEta;
            evalBasisGrad(PhiXi, PhiEta, xi1, xi2, p_, dPhiXi, dPhiEta);
            Eigen::MatrixXd gradPhi(2, Np);
            for (int i = 0; i < Np; ++i) {
                Eigen::Vector2d g_ref{dPhiXi(i), dPhiEta(i)};
                gradPhi.col(i) = JinvT * g_ref;
            }
            
            Eigen::Vector4d U = reconstructState(modesE, phi);
            Eigen::Matrix<double,4,2> F = physicalFlux(U, gamma); // DOES THIS NEED TO BE ROE/HLLE INSTEAD?
            Re += wQ(q) * detJ * (F * gradPhi);
        }
 
        residual.block(0, e*Np, nEq, Np) += Re;
    }
}


void DGIntegration::getLineIntegral(int p_,
                                     double gamma,
                                     const TriangularMesh& mesh,
                                     const std::vector<std::shared_ptr<BoundaryCondition>>& bc,
                                     const FVFlux& flux,
                                     const Eigen::MatrixXd& modes,
                                     Eigen::MatrixXd& residual)
{
    ShapeFunctions sf;
    int Np  = (p_+1)*(p_+2)/2;
    int nEq = 4;
    int nQE = p_ + 1;
    Eigen::VectorXd gPts, gWts;
    gaussLegendre1D(nQE, gPts, gWts);
 
    Eigen::MatrixXd Phi = sf.getShapeFuncCoeffs(p_); 
 
    int nElems = (int)mesh.numElems();
    int nFaces = (int)mesh.numFaces();
 
    // Loop over edges
    for (int f = 0; f < nFaces; ++f)
    {
        const auto& face = mesh.face(f);
        double edgeLen = face._length;
        Eigen::Vector2d nHat = face._normal;   // MAKE SURE THIS IS UNIT OUTWARD NORMAL FOR ALL CASES
 
        int eL = face._elemID[0];   //"left"  element
        int eR = face._elemID[1];   //"right" element
 
        const auto& elemL = mesh.elem(eL);
        int localEdgeL = -1;
        for (int k = 0; k < 3; ++k)
            if (elemL._faceID[k] == f) { localEdgeL = k; break; }
 
        int localEdgeR = -1;
        int eR_actual  = eR;
        bool isBoundary = face.isBoundaryFace() && !face.isPeriodicFace();
        bool isPeriodic = face.isPeriodicFace();
 
        if (!isBoundary) {
            // Interior or periodic
            int eR_for_geom = isPeriodic ? face._periodicElemID : eR;
            eR_actual = eR_for_geom;
            const auto& elemR = mesh.elem(eR_actual);
            int faceID_R = isPeriodic ? face._periodicFaceID : f;
            for (int k = 0; k < 3; ++k)
                if (elemR._faceID[k] == faceID_R) { localEdgeR = k; break; }
        }
 
        Eigen::MatrixXd modesL = modes.block(0, eL*Np, nEq, Np);
        Eigen::MatrixXd modesR;
        if (!isBoundary)
            modesR = modes.block(0, eR_actual*Np, nEq, Np);
 
        Eigen::MatrixXd ReL = Eigen::MatrixXd::Zero(nEq, Np);
        Eigen::MatrixXd ReR = Eigen::MatrixXd::Zero(nEq, Np);
 
        // Loop over quadrature
        for (int q = 0; q < nQE; ++q)
        {
            double t = gPts(q);
 
            Eigen::Vector2d xiL_ref = edgeRefPoint(localEdgeL, t);
            Eigen::VectorXd phiL = evalBasis(Phi, xiL_ref(0), xiL_ref(1), p_);
 
            Eigen::Vector4d UL = reconstructState(modesL, phiL);
 
            Eigen::Vector4d numFlux;
 
            if (isBoundary) {
                int bcGroup = (int)face._nf - 1;   // 0-based?
                numFlux = bc[bcGroup]->computeFlux(UL, nHat);
 
            } else {
                //interior/ periodic edge
                Eigen::Vector2d xiR_ref = edgeRefPoint(localEdgeR, 1.0 - t);
                Eigen::VectorXd phiR = evalBasis(Phi, xiR_ref(0), xiR_ref(1), p_);
 
                Eigen::Vector4d UR = reconstructState(modesR, phiR);
 
                // Riemann solver -- check to make sure this is right
                numFlux = flux.computeFlux(UL, UR, nHat);
                ReR -= gWts(q) * edgeLen * numFlux * phiR.transpose();
            }
            ReL += gWts(q) * edgeLen * numFlux * phiL.transpose();
        }

        residual.block(0, eL*Np, nEq, Np) -= ReL;
        if (!isBoundary)
            residual.block(0, eR_actual*Np, nEq, Np) += ReR;
    }
}

void DGIntegration::computeResidual(int p_,
                                     double gamma,
                                     const TriangularMesh& mesh,
                                     const std::vector<std::shared_ptr<BoundaryCondition>>& bc,
                                     const FVFlux& flux,
                                     const Eigen::MatrixXd& modes,
                                     Eigen::MatrixXd& residual)
{
    residual.setZero();
    getElemIntegral(p_, gamma, mesh, modes, residual);
    getLineIntegral(p_, gamma, mesh, bc, flux, modes, residual);
}