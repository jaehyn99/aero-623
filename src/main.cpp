#include <iostream>

#include "Constants.h"
#include "FVAdvectionFirstOrder.h"
#include "FVSteadySolver.h"
#include "FVUnsteadySolver.h"
#include "HLLEFlux.h"
#include "InletBC.h"
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

    // Boundary conditions
    std::shared_ptr<BoundaryCondition> inlet = std::make_shared<InletBC>(rho0, a0, alpha, gamma);
    std::shared_ptr<BoundaryCondition> wall = std::make_shared<InviscidWallBC>(gamma);
    std::shared_ptr<BoundaryCondition> outlet = std::make_shared<OutletBC>(pout, gamma);
    std::vector<std::shared_ptr<BoundaryCondition>> bc{wall, inlet, wall, outlet};
    
    // Initialize the state mesh
    StateMesh states(mesh, bc);
    states.state(0).fill(rho0);
    states.state(1).fill(a0);
    states.state(2).fill(0);
    states.state(3).fill(p0/(gamma-1) + 0.5*rho0*a0*a0);

    // Solver
    std::shared_ptr<FVFlux> flux = std::make_shared<RoeFlux>(gamma);
    std::shared_ptr<FVResidual> residual = std::make_shared<FVAdvectionFirstOrder>(flux);
    std::shared_ptr<TimeIntegrator> integrator = std::make_shared<SSP_RK2>();
    std::shared_ptr<TimeStepper> stepper = std::make_shared<LocalTimeStepper>(0.1, gamma, flux);
    FVSteadySolver ssSolver(residual, integrator, stepper);
    ssSolver.solve(states);

    return 0;
}