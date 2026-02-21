#include "LocalTimeStepper.h"
#include "BoundaryCondition.h"
#include "FVFlux.h"
#include "StateMesh.h"
#include "TriangularMesh.h"

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
    std::vector<std::string> bcNames;

    // Computes the wave speed on each edge
    Eigen::ArrayXd sFace(mesh->numFaces());
    for (std::size_t i = 0; i < mesh->numFaces(); i++){
        const auto& face = mesh->face(i);
        Eigen::Vector2d normal = face._normal; // direction doesn't matter here
        Eigen::Index e = face._elemID[0];
        double cL = speedOfSound(u.cell(e), _gamma);
        if (face.isBoundaryFace()){
            if (face.isPeriodicFace()){
                Eigen::Index ne = face._periodicElemID; // ne = neighbor element index
                double cR = speedOfSound(u.cell(ne), _gamma);
                sFace[i] = _flux->computeWaveSpeed(u.cell(e), u.cell(ne), normal, cL, cR);
            } else{
                std::size_t boundaryID;
                auto it = std::find(bcNames.cbegin(), bcNames.cend(), face._title);
                if (it == bcNames.cend()){
                    // New boundary, not yet registered in bcNames;
                    bcNames.push_back(face._title);
                    boundaryID = bcNames.size()-1;
                } else boundaryID = it - bcNames.cbegin();
                Eigen::Vector4d Ub = u.bc(boundaryID)->computeBoundaryState(u.cell(e), normal); // last entry is speed of sound, not energy
                sFace[i] = _flux->computeWaveSpeed(u.cell(e), Ub, normal, cL, Ub[3]);
            }
        } else{
            Eigen::Index ne = face._elemID[1];
            double cR = speedOfSound(u.cell(ne), _gamma);
            sFace[i] = _flux->computeWaveSpeed(u.cell(e), u.cell(ne), normal, cL, cR);
        }
    }

    // Computes the wave speed on each element
    Eigen::ArrayXd sElem(mesh->numElems());
    for (std::size_t i = 0; i < mesh->numElems(); i++){
        const auto& elem = mesh->elem(i);
        double P = 0;
        for (auto faceID: elem._faceID) P += sFace[faceID]*mesh->length(faceID);
        sElem[i] = 2*elem._area*_minCFL/P;
    }
    return sElem;
}