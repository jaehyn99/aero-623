#include "FVAdvectionFirstOrder.h"
#include "BoundaryCondition.h"
#include "FVFlux.h"
#include "StateMesh.h"
#include "TriangularMesh.h"

#include <iostream>

Eigen::MatrixXd FVAdvectionFirstOrder::computeResidual(const StateMesh& u) const{
    auto mesh = u.mesh();
    Eigen::MatrixXd residual = Eigen::MatrixXd::Zero(4, u.cellCount());
    if (_flux){
        std::vector<std::string> bcNames;
        // Loops through the boundary edges and add their names
        for (const auto& face: mesh->getFaces()){
            if (!face.isBoundaryFace()) break;
            if (face.isPeriodicFace()) continue;
            auto it = std::find(bcNames.cbegin(), bcNames.cend(), face._title);
            if (it == bcNames.cend()) bcNames.push_back(face._title);
        }

        // Loops over all elements and calculate flux contributions
        for (Eigen::Index e = 0; e < u.cellCount(); e++){
            const auto& elem = mesh->elem(e);
            // Loops over the three edges
            for (Eigen::Index i = 0; i < 3; i++){
                const auto& face = mesh->face(elem._faceID[i]);
                Eigen::Vector2d normal = mesh->normal(e,i);

                Eigen::Vector4d flux;
                if (face.isBoundaryFace()){
                    if (face.isPeriodicFace()){
                        Eigen::Index ne = face._periodicElemID; // ne = neighbor element index
                        flux = _flux->computeFlux(u.cell(e), u.cell(ne), normal);
                        if (flux.array().isNaN().any()){
                            std::cout << "NaN flux on a periodic boundary of cell " << e << std::endl;
                            std::cout << "Cell: " << u.cell(e).transpose() << std::endl;
                            std::cout << "Neighbor: " << u.cell(ne).transpose() << std::endl;
                            std::cout << "Normal: " << normal.transpose() << std::endl;
                            std::cout << "Flux: " << flux.transpose() << std::endl; 
                            throw std::runtime_error("ERROR: Unable to evaluate flux");
                        }
                    } else{
                        std::size_t boundaryID = std::find(bcNames.cbegin(), bcNames.cend(), face._title) - bcNames.cbegin();
                        flux = u.bc(boundaryID)->computeFlux(u.cell(e), normal);
                        if (flux.array().isNaN().any()){
                            std::cout << "NaN flux on a non-periodic boundary " << face._title << " of cell " << e << std::endl;
                            std::cout << "Cell: " << u.cell(e).transpose() << std::endl;
                            std::cout << "Normal: " << normal.transpose() << std::endl;
                            std::cout << "Flux: " << flux.transpose() << std::endl; 
                            throw std::runtime_error("ERROR: Unable to evaluate flux");
                        }
                    }
                } else{
                    Eigen::Index ne = face._elemID[0] == e ? face._elemID[1] : face._elemID[0];
                    if (ne < e) normal *= -1; // invert the normal since it's default to point from L to R
                    flux = _flux->computeFlux(u.cell(e), u.cell(ne), normal);
                    if (flux.array().isNaN().any()){
                        std::cout << "NaN flux on an interior boundary of cell " << e << std::endl;
                        std::cout << "Cell: " << u.cell(e).transpose() << std::endl;
                        std::cout << "Neighbor: " << u.cell(ne).transpose() << std::endl;
                        std::cout << "Normal: " << normal.transpose() << std::endl;
                        std::cout << "Flux: " << flux.transpose() << std::endl; 
                        throw std::runtime_error("ERROR: Unable to evaluate flux");
                    }
                }
                residual.col(e) -= flux * face._length;
            }
            residual.col(e) /= mesh->area(e);
        }
    }
    return residual;
}