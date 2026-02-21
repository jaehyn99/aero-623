#include "FVAdvectionSecondOrder.h"
#include "BoundaryCondition.h"
#include "FVFlux.h"
#include "StateMesh.h"
#include "TriangularMesh.h"

#include "InletBC.h"
#include "InviscidWallBC.h"
#include "OutletBC.h"

#include <iostream>
#include <fstream>
#include <set>

using GradMatrix = std::vector<Eigen::Matrix<double,4,2>>;
using State = Eigen::Vector4d;
using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;
using GradCell = Eigen::Matrix<double,4,2>;

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
    auto gradients = u.getGradients();
    if (gradients.size() != static_cast<std::size_t>(u.cellCount()))
    throw std::runtime_error("gradients size mismatch with cellCount()");

    const int nElem = int(u.cellCount());
    // Unpack the gradients into gradX and gradY for easier access
    StateMatrix gradX = StateMatrix::Zero(4, nElem);
    StateMatrix gradY = StateMatrix::Zero(4, nElem);
    for (int e = 0; e < nElem; ++e) {
        gradX.col(e) = gradients[e].col(0);
        gradY.col(e) = gradients[e].col(1);
    }
    // Initialize the mesh from the state
    auto mesh = u.mesh();

    Eigen::MatrixXd residual = Eigen::MatrixXd::Zero(u.stateCount(), u.cellCount());

    if (_flux){
        std::vector<std::string> bcNames;

        // Loops through the boundary edges and add their names
        for (const auto& face: mesh->getFaces())
        {
            if (!face.isBoundaryFace()) break;
            if (face.isPeriodicFace()) continue;
            auto it = std::find(bcNames.cbegin(), bcNames.cend(), face._title);
            if (it == bcNames.cend()) bcNames.push_back(face._title);
        }

        // std::ofstream file("result.txt");

        // Loops over all elements and calculate flux contributions
        for (Eigen::Index e = 0; e < u.cellCount(); e++)
        {
            const auto& elem = mesh->elem(e);
            // file << "On element " << e << ". Area = " << mesh->area(e) << ".\n";
            
            // Loops over the three edges
            for (Eigen::Index i = 0; i < 3; i++)
            {
                const auto& face = mesh->face(elem._faceID[i]);
                // file << "\tLocal face ID " << i << ", which is face " << elem._faceID[i] << ". Length = " << mesh->length(elem._faceID[i]) << ".\n";
                Eigen::Vector2d normal = mesh->normal(e,i);
                // file << "\tThe normal vector is " << normal.transpose() << ".\n";

                double edge_length = mesh->length(elem._faceID[i]);

                Eigen::Vector4d flux;

                // face midpoint
                Eigen::Vector2d p0 = mesh->node(face._pointID[0]);
                Eigen::Vector2d p1 = mesh->node(face._pointID[1]);
                Eigen::Vector2d xf = 0.5*(p0 + p1);

                // left cell geometry
                Eigen::Vector2d xcL = mesh->centroid(e);
                Eigen::Vector2d dxL = xf - xcL;

                // reconstructed left state at face
                Eigen::Vector4d ULhat = reconstructToFace(u.cell(e), gradients[e], dxL);

                if (face.isBoundaryFace() && !face.isPeriodicFace())
                {
                    std::size_t boundaryID = std::find(bcNames.cbegin(), bcNames.cend(), face._title) - bcNames.cbegin();
                    // file << "\tThis is a non-periodic boundary face. It is on boundary " << face._title << ".";
                    // if (dynamic_cast<InletBC*>(u.bc(boundaryID).get())) file << " This is an inlet BC.\n";
                    // else if (dynamic_cast<InviscidWallBC*>(u.bc(boundaryID).get())) file << " This is an inviscid wall BC.\n";
                    // else file << " This is an outlet BC.\n";
                    flux = u.bc(boundaryID)->computeFlux(ULhat, normal);
                    residual.col(e) -= flux * face._length / mesh->area(e);
                    if (flux.array().isNaN().any()){
                        std::cout << "NaN flux on a non-periodic boundary " << face._title << " of cell " << e << std::endl;
                        std::cout << "Cell: " << u.cell(e).transpose() << std::endl;
                        std::cout << "Normal: " << normal.transpose() << std::endl;
                        std::cout << "Flux: " << flux.transpose() << std::endl; 
                        throw std::runtime_error("ERROR: Unable to evaluate flux");
                    }
                } 
                else
                {
                    Eigen::Index ne; // ne = neighbor element index
                    if (face.isPeriodicFace())
                    {
                        ne = face._periodicElemID; 
                        // file << "\tThis is a periodic boundary face on " << face._title << ".";
                        // file << "\tIt is periodic with element " << ne << " on " << mesh->face(face._periodicFaceID)._title << ".\n";
                    } 
                    else
                    {
                        ne = face._elemID[0] == e ? face._elemID[1] : face._elemID[0];
                        // file << "\tThis is an interior face. The neighboring cell is cell " << ne << ".\n";
                    }
                    
                    if (e < ne)
                    {
                        // right cell geometry
                        Eigen::Vector2d xcR = mesh->centroid(ne);
                        Eigen::Vector2d dxR = xf - xcR;

                        // reconstructed right state at face
                        Eigen::Vector4d URhat = reconstructToFace(u.cell(ne), gradients[ne], dxR);

                        flux = _flux->computeFlux(ULhat, URhat, normal);
                        residual.col(e) -= flux * face._length / mesh->area(e);
                        residual.col(ne) += flux * face._length / mesh->area(ne);
                        if (flux.array().isNaN().any()){
                            std::cout << "NaN flux on an interior boundary of cell " << e << std::endl;
                            std::cout << "Cell: " << u.cell(e).transpose() << std::endl;
                            std::cout << "Neighbor: " << u.cell(ne).transpose() << std::endl;
                            std::cout << "Normal: " << normal.transpose() << std::endl;
                            std::cout << "Flux: " << flux.transpose() << std::endl; 
                            throw std::runtime_error("ERROR: Unable to evaluate flux");
                        }
                    }
                }
                // file << "\tThe flux contribution of this face is " << (flux*face._length).transpose() << ".\n\n";
            }
        }
        // file.close();
    }
    return residual;
}