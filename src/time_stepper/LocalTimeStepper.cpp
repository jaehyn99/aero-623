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

    // Computes the wave speed on each edge
    Eigen::ArrayXd sElem(mesh->numElems());
    auto& I2E = mesh->I2E();
    auto& B2E = mesh->B2E();
    auto& In = mesh->In();
    auto& Bn = mesh->Bn();

    // Loops over interior edges
    for (Eigen::Index i = 0; i < I2E.rows(); i++){
        Eigen::Vector2d normal = In.row(i);
        Eigen::Index e = I2E(i,0);
        double cL = speedOfSound(u.cell(e), _gamma);
        Eigen::Index ne = I2E(i,2);
        double cR = speedOfSound(u.cell(ne), _gamma);
        double sFace = _flux->computeWaveSpeed(u.cell(e), u.cell(ne), normal, cL, cR);

        double length = mesh->length(e, I2E(i,1));
        sElem[e] += sFace * length; // sElem first stores the sum of s*length
        sElem[ne] += sFace * length;
    }

    for (Eigen::Index i = 0; i < B2E.rows(); i++){
        Eigen::Vector2d normal = Bn.row(i);
        Eigen::Index e = B2E(i,0);
        double cL = speedOfSound(u.cell(e), _gamma);

        std::size_t boundaryID = B2E(i,2);
        Eigen::Vector4d Ub = u.bc(boundaryID)->computeBoundaryState(u.cell(e), normal); // last entry is speed of sound, not energy
        double sFace = _flux->computeWaveSpeed(u.cell(e), Ub, normal, cL, Ub[3]);

        double length = mesh->length(e, I2E(i,1));
        sElem[e] += sFace * length; // sElem first stores the sum of s*length
    }

    // Computes the wave speed on each element
    for (std::size_t i = 0; i < mesh->numElems(); i++){ sElem[i] = 2*mesh->area(i)*_minCFL/sElem[i]; }
    return sElem;
}