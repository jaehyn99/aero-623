#include "FVAdvectionSecondOrder.h"
#include "InletBC.h"
#include "InletOutletBC.h"
#include "FVFlux.h"
#include "StateMesh.h"
#include "TriangularMesh.h"

#include <iostream>
#include <fstream>

namespace{
    using GradMatrix = std::vector<Eigen::Matrix<double,4,2>>;
    using State = Eigen::Vector4d;
    using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;
    using GradCell = Eigen::Matrix<double,4,2>;
}

static inline Eigen::Vector4d reconstructToFace(
    const Eigen::Vector4d& Uc,
    const GradCell& G,
    const Eigen::Vector2d& dx)
{
    // G.col(0) is dU/dx, G.col(1) is dU/dy
    return Uc + G.col(0)*dx(0) + G.col(1)*dx(1);
}


Eigen::MatrixXd FVAdvectionSecondOrder::computeResidual(const StateMesh& u) const
{
    auto& grad = u.getGradient();
    if (grad.size() != static_cast<std::size_t>(u.cellCount()))
        throw std::runtime_error("Gradient size mismatch with cellCount()");

    const int nElem = int(u.cellCount());
    // Unpack the gradient into gradX and gradY for easier access
    StateMatrix gradX = StateMatrix::Zero(4, nElem);
    StateMatrix gradY = StateMatrix::Zero(4, nElem);
    for (int e = 0; e < nElem; ++e) {
        gradX.col(e) = grad[e].col(0);
        gradY.col(e) = grad[e].col(1);
    }
    // Initialize the mesh from the state
    auto mesh = u.mesh();
    Eigen::MatrixXd residual = Eigen::MatrixXd::Zero(4, u.cellCount());
    if (_flux){
        std::vector<std::string> bcNames;
        std::vector<bool> visited(mesh->numFaces(), false);

        for (std::size_t i = 0; i < mesh->numFaces(); i++){
            if (visited[i]) continue;
            const auto& face = mesh->face(i);
            std::size_t e = face._elemID[0];
            Eigen::Vector2d normal = face._normal; // points outward, from elemL to elemR

            // face midpoint
            Eigen::Vector2d p0 = mesh->node(face._pointID[0]);
            Eigen::Vector2d p1 = mesh->node(face._pointID[1]);
            Eigen::Vector2d xf = 0.5*(p0 + p1);

            // left cell geometry
            Eigen::Vector2d xcL = mesh->centroid(e);
            Eigen::Vector2d dxL = xf - xcL;

            // reconstructed left state at face
            Eigen::Vector4d ULhat = reconstructToFace(u.cell(e), grad[e], dxL);

            if (face.isBoundaryFace() && !face.isPeriodicFace()){
                auto it = std::find(bcNames.cbegin(), bcNames.cend(), face._title);
                std::size_t boundaryID;
                if (it == bcNames.cend()){
                    bcNames.push_back(face._title);
                    boundaryID = bcNames.size() - 1;
                } else boundaryID = it - bcNames.cbegin();

                auto bc = u.bc(boundaryID);
                // Transient cases, these lines do nothing if running in steady-state
                if (auto inlet = dynamic_cast<InletBC*>(bc.get())){
                    double y0 = mesh->node(face._pointID[0])[1];
                    double y1 = mesh->node(face._pointID[1])[1];
                    inlet->setTransientRho((y1+y0)/2);
                } else if (auto inlet = dynamic_cast<InletOutletBC*>(bc.get())){
                    double y0 = mesh->node(face._pointID[0])[1];
                    double y1 = mesh->node(face._pointID[1])[1];
                    inlet->setTransientRho((y1+y0)/2);                       
                }
                
                Eigen::Vector4d flux = bc->computeFlux(ULhat, normal);
                if (flux.array().isNaN().any()){
                    std::cout << "UL = " << u.cell(e).transpose() << std::endl;
                    std::cout << "gradL = [" << grad[e].col(0).transpose() << "]" << std::endl;
                    std::cout << "ULhat = " << ULhat.transpose() << std::endl;
                    if (face.isPeriodicFace()) std::cout << "Boundary edge connecting on cell " << e << " have NaN flux." << std::endl << std::endl;
                }
                residual.col(e) -= flux * mesh->length(i) / mesh->area(e);
            } else{
                std::size_t ne; // neighbor element index
                if (face.isPeriodicFace()){
                    visited[i] = true;
                    visited[face._periodicFaceID] = true;
                    ne = face._periodicElemID;
                } else ne = face._elemID[1];

                // right cell geometry
                Eigen::Vector2d xcR = mesh->centroid(ne);
                Eigen::Vector2d dxR;
                if (face.isPeriodicFace()){
                    const auto& pFace = mesh->face(face._periodicFaceID);
                    p0 = mesh->node(pFace._pointID[0]);
                    p1 = mesh->node(pFace._pointID[1]);
                    dxR = 0.5*(p0 + p1) - xcR;
                }
                else dxR = xf - xcR;

                // reconstructed right state at face
                Eigen::Vector4d URhat = reconstructToFace(u.cell(ne), grad[ne], dxR);
                Eigen::Vector4d flux = _flux->computeFlux(ULhat, URhat, normal);
                if (flux.array().isNaN().any()){
                    std::cout << "UL = " << u.cell(e).transpose() << ", UR = " << u.cell(ne).transpose() << std::endl;
                    std::cout << "gradL = [" << grad[e].col(0).transpose() << "], [" << grad[e].col(1).transpose() << "]" << std::endl;
                    std::cout << "gradR = [" << grad[ne].col(0).transpose() << "], [" << grad[ne].col(1).transpose() << "]" << std::endl;
                    std::cout << "ULhat = " << ULhat.transpose() << ", URhat = " << URhat.transpose() << std::endl;
                    if (face.isPeriodicFace()) std::cout << "Edge connecting periodic cells " << e << " and " << ne << " have NaN flux." << std::endl << std::endl;
                    else std::cout << "Edge connecting interior cells " << e << " and " << ne << " have NaN flux." << std::endl << std::endl;
                }
                residual.col(e) -= flux * face._length / mesh->area(e);
                residual.col(ne) += flux * face._length / mesh->area(ne);
            }
        }
    }

    if (residual.array().isNaN().any()) throw std::runtime_error("ERROR: Cannot compute residual");
    return residual;
}