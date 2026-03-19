#include "FEAdvection.h"
#include "Element.h"
#include "Face.h"
#include "InletBC.h"
#include "InletOutletBC.h"
#include "Lagrange2DBasisFunctions.h"
#include "FVFlux.h"
#include "StateMesh.h"
#include "TriangularMesh.h"

#include <Eigen/LU>
#include <iostream>

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

Eigen::MatrixXd FEAdvection::computeResidual(const StateMesh& u) const{
    auto mesh = u.mesh();
    Eigen::MatrixXd residual = Eigen::MatrixXd::Zero(u.stateCount(), u.cellCount()*u.Np());
    if (_flux){
        std::vector<std::string> bcNames; // to keep track of boundary conditions
        const ReferenceElement& ref = mesh->reference();

        std::size_t p = u.p();
        std::size_t Np = u.Np();
        Lagrange2DBasisFunctions PhiBasis(p);

        const auto& intXi = ref.intXi();
        const auto& intW = ref.intW();
        // const auto& intPhi = ref.intPhi();
        const auto& intPhiXi = ref.intPhiXi();
        const auto& intPhiEta = ref.intPhiEta();
        std::size_t intNq = intW.size();

        const auto& edgeW = ref.edgeW();
        std::size_t edgeNq = edgeW.size();
        
        // std::cout << "Entering element loop" << std::endl;
        for (std::size_t k = 0; k < mesh->numElems(); k++){
            std::cout << k << std::endl;
            Eigen::MatrixXd cell = u.cell(k); // block matrix with basis function weight
            const Element& elem = mesh->elem(k);
            // std::cout << "Entering basis function loop" << std::endl;
            for (std::size_t i = 0; i < unsigned(Np); i++){
                // Area integral over the entire element
                // std::cout << "Computing area integral" << std::endl;
                for (std::size_t nq = 0; nq < intNq; nq++){
                    double xi = intXi(0, nq);
                    double eta = intXi(1, nq);
                    Eigen::Vector4d u = PhiBasis.funcEval(xi, eta, cell); // state values at the quadrature point
                    Eigen::Matrix<double,4,2> F = physicalFlux(u, _flux->gamma()); // get fluxes from states
                    
                    double detJ = elem.detJacobian(nq);
                    Eigen::Vector2d gRef{intPhiXi(i,nq), intPhiEta(i,nq)}; // 2-by-1
                    Eigen::Matrix2d JT = elem.jacobian(nq).transpose(); // 2-by-2
                    Eigen::Vector2d gPhi = JT.lu().solve(gRef); // 2-by-1
                    residual.col(k*Np+i) += intW[nq] * detJ * F*gPhi; // (4-by-2)*(2-by-1) = (4-by-1)
                }

                // Line integral over each of the edges
                // std::cout << "Computing line integrals" << std::endl;
                for (std::size_t edge = 0; edge < 3; edge++){
                    // std::cout << "Line #" << edge << std::endl;
                    const auto& edgeXi = ref.edgeXi(edge);
                    const auto& edgePhi = ref.edgePhi(edge);
                    // const auto& edgePhiXi = ref.edgePhiXi(edge);
                    // const auto& edgePhiEta = ref.edgePhiEta(edge);
                    double factor = (edge == 1) ? std::sqrt(2) : 1.0; // edge 0 is the hypotnuse and needs an extra weight factor

                    auto faceID = elem.faceID(edge);
                    const auto& face = mesh->face(faceID);
                    // if (face.isPeriodicFace()) std::cout << "This is a periodic face" << std::endl;
                    // else std::cout << "This is not a periodic face" << std::endl;

                    if (face.isBoundaryFace()){
                        // std::cout << "\tThis face is a boundary face." << std::endl;
                        // Boundary edge, use boundary condition to compute flux
                        auto it = std::find(bcNames.cbegin(), bcNames.cend(), face.title());
                        std::size_t boundaryID;
                        if (it == bcNames.cend()){
                            bcNames.push_back(face.title());
                            boundaryID = bcNames.size() - 1;
                        } else boundaryID = it - bcNames.cbegin();

                        auto bc = u.bc(boundaryID);
                        // Transient cases, these lines do nothing if running in steady-state
                        if (auto inlet = dynamic_cast<InletBC*>(bc.get())){
                            double y0 = mesh->node(face.pointID(0))[1];
                            double y1 = mesh->node(face.pointID(1))[1];
                            inlet->setTransientRho((y1+y0)/2);
                        } else if (auto inlet = dynamic_cast<InletOutletBC*>(bc.get())){
                            double y0 = mesh->node(face.pointID(0))[1];
                            double y1 = mesh->node(face.pointID(1))[1];
                            inlet->setTransientRho((y1+y0)/2);                       
                        }

                        for (std::size_t nq = 0; nq < edgeNq; nq++){
                            double xi = edgeXi(0, nq);
                            double eta = edgeXi(1, nq);
                            double phi = edgePhi(i, nq);
                            Eigen::Vector4d u = PhiBasis.funcEval(xi, eta, cell); // state values at the quadrature point
                            Eigen::Vector2d normal = face.normal(nq);
                            Eigen::Vector4d F = bc->computeFlux(u, normal); // get numerical flux from boundary condition                          
                            residual.col(k*Np+i) -= edgeW[nq] * face.detJ(nq) * factor * phi*F;
                        } 
                    } else{

                        // Internal edge, use numerical flux function to compute flux
                        std::size_t kn;
                        Eigen::Vector2d normal = face.normal(); // linear face so normal does not vary along the edge
                        // std::cout << face.elemID().transpose() << std::endl;
                        if (face.elemID(0) == int(k)) kn = face.elemID(1);
                        else{
                            kn = face.elemID(0);
                            normal *= -1; // Element k is the R element, has to flip normal
                        }

                        // Find the local edge index of this edge on element kn
                        std::size_t edgeN;
                        if (mesh->elem(kn).faceID(0) == faceID) edgeN = 0;
                        else if (mesh->elem(kn).faceID(1) == faceID) edgeN = 1;
                        else edgeN = 2;

                        Eigen::MatrixXd cellN = u.cell(kn);
                        Eigen::Matrix2Xd edgeXiN = ref.edgeXi(edgeN);
                        const Element& elemN = mesh->elem(kn);
                        // std::cout << "\tLocating nodes on elemL" << std::endl;
                        Eigen::Vector2d x0 = mesh->node(elem.pointID(0)); // std::cout << x0.transpose() << std::endl;
                        Eigen::Vector2d x1 = mesh->node(elem.pointID(1)); // std::cout << x1.transpose() << std::endl;
                        Eigen::Vector2d x2 = mesh->node(elem.pointID(2)); // std::cout << x2.transpose() << std::endl;
                        // std::cout << "\tLocating nodes on elemR" << std::endl;
                        Eigen::Vector2d x0n = mesh->node(elemN.pointID(0)); // std::cout << x0n.transpose() << std::endl;
                        Eigen::Vector2d x1n = mesh->node(elemN.pointID(1)); // std::cout << x1n.transpose() << std::endl;
                        Eigen::Vector2d x2n = mesh->node(elemN.pointID(2)); // std::cout << x2n.transpose() << std::endl;
                        // std::cout << "Face " << edge << " on element " << k << " is face " << edgeN << " on element " << kn << std::endl;

                        // std::cout << "\tThe actual integration" << std::endl;
                        for (std::size_t nq = 0; nq < edgeNq; nq++){
                            // std::cout << "\t Quadrature " << nq << std::endl;
                            double xi = edgeXi(0, nq);
                            double eta = edgeXi(1, nq);
                            Eigen::Vector4d uL = PhiBasis.funcEval(xi, eta, cell); // left-state values at the quadrature point

                            double xiN = edgeXiN(0, edgeNq-1-nq);
                            double etaN = edgeXiN(1, edgeNq-1-nq);
                            Eigen::Vector4d uR = PhiBasis.funcEval(xiN, etaN, cellN); // right-state values at the quadrature point

                            // std::cout << "This point on referenceL is " << xi << "\t" << eta << std::endl;
                            // std::cout << "This point on elemL is " << (x0 + (x1-x0)*xi + (x2-x0)*eta).transpose() << std::endl;
                            // std::cout << "This point on referenceL is " << xiN << "\t" << etaN << std::endl;
                            // std::cout << "This point on elemR is " << (x0n + (x1n-x0n)*xiN + (x2n-x0n)*etaN).transpose() << std::endl;

                            double phi = edgePhi(i, nq);
                            Eigen::Vector4d F = _flux->computeFlux(uL, uR, normal); // get numerical flux from boundary condition                   
                            residual.col(k*Np+i) -= edgeW[nq] * face.detJ(nq) * factor * phi*F;
                        } 
                        // std::cout << std::endl;
                    }
                }
            }
        }
    }
    return residual;
}