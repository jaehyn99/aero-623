#include "Solver.h"
#include "Residual.h"
#include "LocalTimeStepper.h"
#include "StateMesh.h"
#include "TimeIntegrator.h"
#include "TimeStepper.h"

Solver::Solver(std::shared_ptr<Residual> residual, std::shared_ptr<TimeIntegrator> integrator, std::shared_ptr<TimeStepper> stepper):
    _residual(residual),
    _integrator(integrator),
    _stepper(stepper)
{
    assert(_residual && _integrator && _stepper);
}

Solver::Solver(Solver&&) = default;
Solver& Solver::operator=(Solver&&) = default; 
Solver::~Solver() = default;