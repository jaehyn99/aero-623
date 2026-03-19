#include "LocalTimeStepper.h"
#include "Element.h"
#include "Face.h"
#include "BoundaryCondition.h"
#include "FVFlux.h"
#include "Lagrange2DBasisFunctions.h"
#include "StateMesh.h"
#include "TriangularMesh.h"
#include <iostream>

static double speedOfSound(const Eigen::Vector4d& U, double gamma){
    double rho = U[0];
    double u = U[1]/rho;
    double v = U[2]/rho;
    double E = U[3]/rho;
    double p = (gamma-1)*(E - 0.5*(u*u + v*v));
    return std::sqrt(gamma*p);
}

Eigen::ArrayXd LocalTimeStepper::dt(const StateMesh& u) const noexcept{
    auto mesh = u.mesh();
    int Np = u.Np();
    std::vector<std::string> bcNames;

    // Computes the wave speed on each edge
    Eigen::ArrayXd sFace(mesh->numFaces());
    for (std::size_t i = 0; i < mesh->numFaces(); i++){
        const auto& face = mesh->face(i);
        Eigen::Vector2d normal = face.normal(); // direction doesn't matter here
        Eigen::Index e = face.elemID(0);

        Eigen::Vector4d ue; // state vector at cell e
        if (Np == 1) ue = u.cell(e);
        else{
            Lagrange2DBasisFunctions Phi(u.p());
            ue = Phi.funcEval(1.0/3, 1.0/3, Eigen::MatrixXd(u.cell(e))); // computes state at the cell centroid
        }

        double cL = speedOfSound(ue, _gamma);
        if (face.isBoundaryFace()){
            std::size_t boundaryID;
            auto it = std::find(bcNames.cbegin(), bcNames.cend(), face.title());
            if (it == bcNames.cend()){
                // New boundary, not yet registered in bcNames;
                bcNames.push_back(face.title());
                boundaryID = bcNames.size()-1;
            } else boundaryID = it - bcNames.cbegin();
            Eigen::Vector4d Ub = u.bc(boundaryID)->computeBoundaryState(ue, normal); // last entry is speed of sound, not energy
            sFace[i] = _flux->computeWaveSpeed(ue, Ub, normal, cL, Ub[3]);
        } else{
            Eigen::Index ne = face.elemID(1);
            Eigen::Vector4d une; // state vector at cell e
            if (Np == 1) une = u.cell(ne);
            else{
                Lagrange2DBasisFunctions Phi(u.p());
                une = Phi.funcEval(1.0/3, 1.0/3, Eigen::MatrixXd(u.cell(ne))); // computes state at the cell centroid
            }
            double cR = speedOfSound(une, _gamma);
            sFace[i] = _flux->computeWaveSpeed(ue, une, normal, cL, cR);
        }
    }

    // Computes the wave speed on each element
    Eigen::ArrayXd sElem(mesh->numElems()*Np);
    for (std::size_t i = 0; i < mesh->numElems(); i++){
        const auto& elem = mesh->elem(i);
        double P = 0;
        for (auto faceID: elem.faceID()) P += sFace[faceID]*mesh->length(faceID);
        sElem.segment(i*Np, Np).fill(2*elem.area()*_minCFL/P);
    }
    return sElem;
}