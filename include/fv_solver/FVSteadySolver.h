#ifndef FV_STEADY_SOLVER_H
#define FV_STEADY_SOLVER_H

#include "FVSolver.h"
class FVSteadySolver: public FVSolver{
    public:
    FVSteadySolver(std::shared_ptr<FVResidual>, std::shared_ptr<TimeIntegrator>, std::shared_ptr<TimeStepper>, double tol=1e-5);
    void solve(StateMesh&) const override;

    protected:
    double _tol;
};

#endif