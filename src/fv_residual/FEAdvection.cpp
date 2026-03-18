#include "FEAdvection.h"
#include "Element.h"
#include "Face.h"
#include "InletBC.h"
#include "InletOutletBC.h"
#include "LagrangeBasisFunctions.h"
#include "FVFlux.h"
#include "StateMesh.h"
#include "TriangularMesh.h"

#include <Eigen/LU>

Eigen::MatrixXd FEAdvection::computeResidual(const StateMesh& u) const{
    auto mesh = u.mesh();
    Eigen::MatrixXd residual = Eigen::MatrixXd::Zero(u.stateCount(), u.cellCount()*u.Np());
    if (_flux){
        std::vector<std::string> bcNames; // to keep track of boundary conditions
        const ReferenceElement& ref = mesh->reference();

        Eigen::Index p = u.p();
        Eigen::Index Np = u.Np();
        LagrangeBasisFunctions PhiBasis(p);

        const auto& intXi = ref.intXi();
        const auto& intW = ref.intW();
        const auto& intPhi = ref.intPhi();
        const auto& intPhiXi = ref.intPhiXi();
        const auto& intPhiEta = ref.intPhiEta();
        Eigen::Index intNq = intW.size();

        const auto& edgeW = ref.edgeW();
        int edgeNq = edgeW.size();
        
        for (std::size_t k = 0; k < mesh->numElems(); k++){
            Eigen::MatrixXd cell = u.cell(k); // block matrix with basis function weight
            const Element& elem = mesh->elem(k);
            for (std::size_t i = 0; i < Np; i++){
                // Area integral over the entire element
                for (std::size_t nq = 0; nq < intNq; nq++){
                    double xi = intXi(0, nq);
                    double eta = intXi(1, nq);
                    Eigen::Vector4d u = PhiBasis.funcEval(xi, eta, cell); // state values at the quadrature point
                    Eigen::Matrix<double,4,2> F; // get flux from state
                    
                    double detJ = elem.detJacobian(nq);

                    Eigen::Vector2d gRef{intPhiXi(i,nq), intPhiEta(i,nq)}; // 2-by-1
                    Eigen::Matrix2d JT = elem.jacobian(nq).transpose(); // 2-by-2
                    Eigen::Vector2d gPhi = JT.lu().solve(gRef); // 2-by-1
                    residual.col(k*Np+i) += intW[nq] * detJ * F*gPhi; // (4-by-2)*(2-by-1) = (4-by-1)
                }

                // Line integral over each of the edges
                for (std::size_t edge = 0; edge < 3; edge++){
                    const auto& edgeXi = ref.edgeXi(edge);
                    const auto& edgePhi = ref.edgePhi(edge);
                    const auto& edgePhiXi = ref.edgePhiXi(edge);
                    const auto& edgePhiEta = ref.edgePhiEta(edge);
                    double factor = (edge == 1) ? std::sqrt(2) : 1.0; // edge 0 is the hypotnuse and needs an extra weight factor

                    auto faceID = elem.pointID()[edge];
                    const auto& face = mesh->face(faceID);
                    double J = face.length()*factor; // face length is the Jacobian determinant

                    if (face.isBoundaryFace()){
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
                            double y0 = mesh->node(face.pointID()[0])[1];
                            double y1 = mesh->node(face.pointID()[1])[1];
                            inlet->setTransientRho((y1+y0)/2);
                        } else if (auto inlet = dynamic_cast<InletOutletBC*>(bc.get())){
                            double y0 = mesh->node(face.pointID()[0])[1];
                            double y1 = mesh->node(face.pointID()[1])[1];
                            inlet->setTransientRho((y1+y0)/2);                       
                        }

                        for (std::size_t nq = 0; nq < edgeNq; nq++){
                            double xi = edgeXi(0, nq);
                            double eta = edgeXi(1, nq);
                            double phi = edgePhi(i, nq);
                            Eigen::Vector4d u = PhiBasis.funcEval(xi, eta, cell); // state values at the quadrature point
                            Eigen::Vector2d normal = face.normal(nq);
                            Eigen::Vector4d F = bc->computeFlux(u, normal); // get numerical flux from boundary condition                          
                            residual.col(k*Np+i) -= edgeW[nq] * J * phi*F;
                        } 
                    } else{
                        // Internal edge, use numerical flux function to compute flux
                        std::size_t kn;
                        Eigen::Vector2d normal = face.normal(); // linear face so normal does not vary along the edge
                        if (face.elemID()[0] == k) kn = face.elemID()[1];
                        else{
                            kn = face.elemID()[0];
                            normal *= -1; // Element k is the R element, has to flip normal
                        }

                        // Find the local edge index of this edge on element kn
                        std::size_t edgeN;
                        if (mesh->elem(kn).faceID()[0] == faceID) edgeN = 0;
                        else if (mesh->elem(kn).faceID()[1] == faceID) edgeN = 1;
                        else edgeN = 2;

                        Eigen::MatrixXd cellN = u.cell(kn);
                        Eigen::Matrix2Xd edgeXiN = ref.edgeXi(edgeN);

                        for (std::size_t nq = 0; nq < edgeNq; nq++){
                            double xi = edgeXi(0, nq);
                            double eta = edgeXi(1, nq);
                            Eigen::Vector4d uL = PhiBasis.funcEval(xi, eta, cell); // left-state values at the quadrature point

                            double xiN = edgeXiN(0, edgeNq-1-nq);
                            double etaN = edgeXiN(1, edgeNq-1-nq);
                            Eigen::Vector4d uR = PhiBasis.funcEval(xiN, etaN, cellN); // right-state values at the quadrature point
                            double phi = edgePhi(i, nq);
                            Eigen::Vector4d F = _flux->computeFlux(uL, uR, normal); // get numerical flux from boundary condition                          
                            residual.col(k*Np+i) -= edgeW[nq] * J * phi*F;
                        } 
                    }
                }
            }
        }
    }
    return residual;
}