#include <iostream>
#include <fstream>

#include "Constants.h"
#include "FVAdvectionFirstOrder.h"
#include "FVSteadySolver.h"
#include "FVUnsteadySolver.h"
#include "HLLEFlux.h"
#include "InletBC.h"
#include "InletOutletBC.h"
#include "InviscidWallBC.h"
#include "LocalTimeStepper.h"
#include "OutletBC.h"
#include "RoeFlux.h"
#include "SSP_RK2.h"
#include "SSP_RK3.h"
#include "StateMesh.h"
#include "TriangularMesh.h"


int main() {
    std::shared_ptr<TriangularMesh> mesh = std::make_shared<TriangularMesh>("projects/Project-1/mesh_refined_2394.gri");

    // Inlet conditions
    double gamma = 1.4;
    double rho0 = 1;
    double a0 = 1;
    double p0 = rho0*a0*a0/gamma;
    double alpha = 50*mconst::pi/180;
    double pout = 0.7*p0;
    double M = 0.1;

    // Boundary conditions
    std::shared_ptr<BoundaryCondition> inlet = std::make_shared<InletOutletBC>(rho0, a0, alpha, pout, gamma);
    std::shared_ptr<BoundaryCondition> wall = std::make_shared<InviscidWallBC>(gamma);
    std::shared_ptr<BoundaryCondition> outlet = std::make_shared<OutletBC>(pout, gamma);
    std::vector<std::shared_ptr<BoundaryCondition>> bc{wall, inlet, wall, outlet};

    // Initialize the state mesh
    StateMesh states(mesh, bc);
    states.state(0).fill(rho0);
    states.state(1).fill(rho0*M*a0);
    states.state(2).fill(0);//(rho0*M*a0*std::sin(alpha));
    states.state(3).fill(p0/(gamma-1) + 0.5*rho0*M*M*a0*a0);

    // Solver
    std::shared_ptr<FVFlux> flux = std::make_shared<HLLEFlux>(gamma);
    std::shared_ptr<FVResidual> residual = std::make_shared<FVAdvectionFirstOrder>(flux);
    std::shared_ptr<TimeIntegrator> integrator = std::make_shared<SSP_RK2>();
    std::shared_ptr<TimeStepper> stepper = std::make_shared<LocalTimeStepper>(1, gamma, flux);
    FVSteadySolver ssSolver(residual, integrator, stepper);

    ssSolver.solve(states);

    auto results = ssSolver.getResult().back();
    auto l1norm = ssSolver.getNorm();

    std::ofstream file("projects/Project-2/coarse_mesh_steady_first_order.txt");
    for (Eigen::Index e = 0; e < states.cellCount(); e++) file << results.col(e).transpose() << "\n";
    file.close();

    file.open("projects/Project-2/coarse_mesh_steady_first_order_norm.txt");
    for (auto norm: l1norm) file << norm << "\n";
    file.close();    

    // std::ofstream f("last_iteration.txt");
    // double maxP = 0;
    // for (Eigen::Index e = 0; e < states.cellCount(); e++){
    //     f << "Cell " << e << ":\n";
    //     f << "\tStates: " << states.cell(e).transpose() << ".\n";
    //     double rhoE = states(3, e);
    //     double KE = 0.5*(states(1,e)*states(1,e) + states(2,e)*states(2,e)) / states(0,e);
    //     double P = (gamma-1) * (rhoE - KE);
    //     if (P > maxP) maxP = P;
    //     f << "\tPressure: " << P << ".\n\n";
    // }

    return 0;
}