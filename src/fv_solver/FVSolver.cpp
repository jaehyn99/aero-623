#include "FVSolver.h"
#include "FVResidual.h"
#include "LocalTimeStepper.h"
#include "StateMesh.h"
#include "TimeIntegrator.h"
#include "TimeStepper.h"

FVSolver::FVSolver(std::shared_ptr<FVResidual> residual, std::shared_ptr<TimeIntegrator> integrator, std::shared_ptr<TimeStepper> stepper):
    _residual(residual),
    _integrator(integrator),
    _stepper(stepper)
{
    assert(_residual && _integrator && _stepper);
}

FVSolver::FVSolver(FVSolver&&) = default;
FVSolver& FVSolver::operator=(FVSolver&&) = default; 
FVSolver::~FVSolver() = default;