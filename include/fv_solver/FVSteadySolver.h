#ifndef FV_STEADY_SOLVER_H
#define FV_STEADY_SOLVER_H

#include "Solver.h"
class FVSteadySolver: public Solver{
    public:
    FVSteadySolver(std::shared_ptr<Residual>, std::shared_ptr<TimeIntegrator>, std::shared_ptr<TimeStepper>, double tol=1e-5);
    void solve(StateMesh&) const override;

    protected:
    double _tol;
};

#endif