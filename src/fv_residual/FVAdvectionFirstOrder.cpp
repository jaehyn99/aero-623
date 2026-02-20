#include "FVAdvectionFirstOrder.h"
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

        // std::ofstream file("result.txt");
        // Loops over all elements and calculate flux contributions
        for (Eigen::Index e = 0; e < u.cellCount(); e++){
            const auto& elem = mesh->elem(e);
            // file << "On element " << e << ". Area = " << mesh->area(e) << ".\n";
            // Loops over the three edges
            for (Eigen::Index i = 0; i < 3; i++){
                const auto& face = mesh->face(elem._faceID[i]);
                // file << "\tLocal face ID " << i << ", which is face " << elem._faceID[i] << ". Length = " << mesh->length(elem._faceID[i]) << ".\n";
                Eigen::Vector2d normal = mesh->normal(e,i);
                // file << "\tThe normal vector is " << normal.transpose() << ".\n";

                Eigen::Vector4d flux;
                if (face.isBoundaryFace() && !face.isPeriodicFace()){
                    std::size_t boundaryID = std::find(bcNames.cbegin(), bcNames.cend(), face._title) - bcNames.cbegin();
                    // file << "\tThis is a non-periodic boundary face. It is on boundary " << face._title << ".";
                    // if (dynamic_cast<InletBC*>(u.bc(boundaryID).get())) file << " This is an inlet BC.\n";
                    // else if (dynamic_cast<InviscidWallBC*>(u.bc(boundaryID).get())) file << " This is an inviscid wall BC.\n";
                    // else file << " This is an outlet BC.\n";
                    flux = u.bc(boundaryID)->computeFlux(u.cell(e), normal);
                    residual.col(e) -= flux * face._length / mesh->area(e);
                    if (flux.array().isNaN().any()){
                        std::cout << "NaN flux on a non-periodic boundary " << face._title << " of cell " << e << std::endl;
                        std::cout << "Cell: " << u.cell(e).transpose() << std::endl;
                        std::cout << "Normal: " << normal.transpose() << std::endl;
                        std::cout << "Flux: " << flux.transpose() << std::endl; 
                        throw std::runtime_error("ERROR: Unable to evaluate flux");
                    }
                } else{
                    Eigen::Index ne; // ne = neighbor element index
                    if (face.isPeriodicFace()){
                        ne = face._periodicElemID; 
                        // file << "\tThis is a periodic boundary face on " << face._title << ".";
                        // file << "\tIt is periodic with element " << ne << " on " << mesh->face(face._periodicFaceID)._title << ".\n";
                    } else{
                        ne = face._elemID[0] == e ? face._elemID[1] : face._elemID[0];
                        // file << "\tThis is an interior face. The neighboring cell is cell " << ne << ".\n";
                    }
                    
                    if (e < ne){
                        flux = _flux->computeFlux(u.cell(e), u.cell(ne), normal);
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