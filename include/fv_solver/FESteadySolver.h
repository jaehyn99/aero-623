#ifndef FE_STEADY_SOLVER_H
#define FE_STEADY_SOLVER_H

#include "Solver.h"
class FESteadySolver: public Solver{
    public:
    FESteadySolver(std::shared_ptr<Residual>, std::shared_ptr<TimeIntegrator>, std::shared_ptr<TimeStepper>, double tol=1e-5);
    void solve(StateMesh&) const override;

    protected:
    double _tol;
};

#endif